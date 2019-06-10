# ############################################################################
# 
# This software is made freely available in accordance with the simplifed BSD
# license:
# 
# Copyright (c) <2017>, <Stephen McMahon>
# All rights reserved
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contacts: Stephen McMahon,	stephen.mcmahon@qub.ac.uk
# 
# ############################################################################
#
# References
#
# 1. McMahon et al, Scientific Reports, 2016, 6, 33290
# 2. McMahon et al, Scientific Reports, 2017, 7, 10790
#
# ############################################################################

import numpy as np
from scipy.integrate import odeint
from scipy.stats import poisson

from .cellcharacteriser import CellCharacteriser
from .trackcharacteriser import TrackCharacteriser

class singleCell:
	# Fixed DNA damage rate
	DSBPerMBPPerGy = 35.0/6100.0
	
	defaultDNAParams = {'sigma': 0.04187, 'NHEJFidelity': 0.98537,'MMEJFidelity': 0.4393, 
						'fastRepair': 2.081,'fastFociDelay': 8.079, 'slowRepair': 0.2604, 
						'slowFociDelay': 0.4053,'verySlowRepair': 0.008462, 
						'complexFrac': 0.4337, 'pointMutationRate': 0.04605,   'failFrac': 0.7364,
						'rNuc': 4.282 }
	defaultSurvParams= {'apoptoticRate': 0.01117, 'mitoticRate': 0.0141, 'baseRate': 0.000739}
	defaultCell = {'dna':6100.0,'chromosomes':46,'repair':0,'G1Arrest':1,'phase':0,'gene':0}

	# Utilities to cache response parameters, shared across all objects
	cellChars = CellCharacteriser(defaultDNAParams)
	trackChars = TrackCharacteriser(defaultDNAParams)
	
	# Model initialisers - start from defaults, overwrite with user-provided.
	def __init__(self, newCell={}, dnaParams={}, survParams={}):
		self.cellLine = None
		self.dnaParams = None
		self.survParams = None

		self.updateParameters(dnaParams,survParams)
		self.updateCellLine(newCell)

	def updateParameters(self,dnaParams={}, survParams={}):
		if self.dnaParams is None: self.dnaParams = singleCell.defaultDNAParams.copy()
		if self.survParams is None: self.survParams = singleCell.defaultSurvParams.copy()

		self.dnaParams.update(dnaParams)
		self.survParams.update(survParams)
		self.refreshCharacteristics()

	def updateCellLine(self,newCell={}):
		if self.cellLine is None: self.cellLine = singleCell.defaultCell.copy()
		self.cellLine.update(newCell)

		# Unroll cell-specific parameters for ease of access below
		self.dna          = self.cellLine['dna'] 
		self.repairDefect = int(self.cellLine['repair'])
		self.phase        = self.cellLine['phase']
		self.geneSize     = self.cellLine['gene']

		if self.phase<-1 or self.phase>3:
			print('Unsupported phase: ' + str(self.phase) +'. Single phases are 0 to 3 (G1-M). \
				   Asynchronous is -1.\nFalling back to G1.')
			self.phase=0
			self.cellLine['phase']=0

		self.refreshCharacteristics()

	# Rebuild characteristics, for new line or parameters. 
	def refreshCharacteristics(self):
		# Catch error if cell is not initialised
		if self.cellLine is None or self.dnaParams is None: return

		# Calculate average DSB induced per keV
		# 6.2415 is 1 Gy in keV per cubic micron
		DSBPerGray = self.DSBPerMBPPerGy*6100 # rNuc corresponds to human genome of 3.05 GBP
		self.DSBPerkeV = DSBPerGray/(6.2415*4.0/3.0*np.pi*pow(self.dnaParams['rNuc'],3))

		# Don't calculate for asynchronous cells
		if self.phase==-1:
			self.pathwayProbs = []
			self.misrepairParams = []
			self.repairFidelities = []
			self.repairRates = []
			return

		# Update characteriser objects
		singleCell.cellChars.updateParameters(self.dnaParams)
		singleCell.trackChars.updateParameters(self.dnaParams)

		# Update repair rates and fidelities
		self.pathwayProbs, self.misrepairParams = singleCell.cellChars.getCellCharacteristics(self.cellLine)
		self.repairRates = [self.dnaParams['fastRepair'], self.dnaParams['slowRepair'],
							self.dnaParams['verySlowRepair']]

		self.repairFidelities = [self.dnaParams['NHEJFidelity'], self.dnaParams['NHEJFidelity'],
								 self.dnaParams['MMEJFidelity']]

		self.eta = self.misrepairParams[0]

		repairDefect=self.repairDefect
		if self.phase>0:
			if (repairDefect // 2 ) % 2 == 0:
				self.repairFidelities[1] = 1 # Available HR is treated as accurate
			else:
				self.repairFidelities[1] = self.dnaParams['MMEJFidelity'] # HR-defective cell in S+
		if repairDefect%2==1:
			self.repairFidelities[0]=self.dnaParams['MMEJFidelity'] # NHEJ defective
			if self.phase==0:
				self.repairFidelities[1]=self.dnaParams['MMEJFidelity'] # Impacts slow repair in G1

	def printCharacteristics(self):
		print("Repair parameters")
		singleCell.cellChars.printCellCharacteristics(self.cellLine)

	# Kinetic model of rates of DSB turnover, and formation of misrepaired breaks
	# Handled by ODEint to calculate total rates
	# y is: 0-2: Fast breaks, slow breaks, very slow breaks;
	#		3-5: intra-track overlap for each break type
	#       6: Number of misrepaired breaks to date
	def breakKinetics(self,y,t,DSBRate):
		# If there are no breaks, nothing changes
		breaks = sum(y[0:3])
		if breaks==0 and DSBRate==0: return [0]*7

		# dB: change in number of breaks of each type
		# List of [fast, slow, MMEJ] rates
		dB = [ self.pathwayProbs[n]*DSBRate - self.repairRates[n]*y[n] for n in range(3)]

		# dH: High-LET clustering effect for each repair pathway
		# Increases by damage, decreases with repair
		# Stored as total effect, which is divided by breaks to get per-break for misrepair
		if self.trackEta>0:
			newOverlap = DSBRate*trackEta
			averageRepair = sum(self.repairRates[n]*y[n] for n in range(3))/breaks
			dH = [newOverlap*self.pathwayProbs[n] - (averageRepair+self.repairRates[n])*y[3+n] 
			      for n in range(3)]
		else:
			dH = [0,0,0]

		# Calculate binary misrepair, and sum contribution from each process
		if breaks>0:
			totalTrackOverlap = sum( y[3:6] )
			errorRate = self.eta*breaks + totalTrackOverlap/breaks
			correctedErrorRate = errorRate + pow(errorRate,2) # Correct for higher-order misrepair
			errorProbability = correctedErrorRate/(1+correctedErrorRate)
			dMisrepair = errorProbability*sum(self.repairRates[n]*y[n] for n in range(3))
		else:
			dMisrepair = 0

		# Return: Change in fast breaks, slow breaks & MMEJ breaks;
		# 		  change in intra-track overlap for each break type
		#		  change in misrepaired breaks
		return dB+dH+[dMisrepair]

	# Model response to a specific set of irradiations
	def modelIrradiation(self,exptCond):
		if self.phase==-1:
			print('Requesting DNA repair model in asynchronous cells. \
				   Currently unsupported - consider manual mixing.\
				   Returning G1 cell output as fallback.')
			self.phase=0

		# Get track characteristics for this irradiation
		DSBPerUm,trackParams = singleCell.trackChars.getTrackCharacteristics(exptCond,self.cellLine)
		self.trackEta = trackParams[0]
		self.eta = self.misrepairParams[0]


		# Initialise cell damage array:
		# Current numbers of breaks, clustering parameter for each pathway, and total misrepaired
		cellDamage = [0,0,0,0,0,0,0] 
		t0 = 0
		self.totalDSBs = 0

		# Get list of exposures. If it's a single value, treat as one instantaneous exposure
		exposures = exptCond['dose']
		if type(exposures) is not list: exposures = [[exposures,0]]
		# Iterate over effects of each exposure. If t=0, just increase DSB count and track overlaps
		for dose, time in exposures:
			DSBs = dose * self.dna * self.DSBPerMBPPerGy
			if self.phase==1: DSBs*=1.5
			if self.phase>=2: DSBs*=2
			self.totalDSBs+=1.0*DSBs

			if time==0:
				# For instant exposure, just increment tracks and overlap
				for n in range(3):
					cellDamage[n]  +=DSBs*self.pathwayProbs[n]
					cellDamage[n+3]+=DSBs*self.trackEta*self.pathwayProbs[n]
			else:
				# Numerical integration of breaks from t0 to t0+time. 
				# If exposure time is long, add in additional guide point at 4*slow repair time
				DSBRate = DSBs/time
				t = [t0,t0+time]
				slowRepairTime = 1.0/self.repairRates[1]
				if time>5.0*slowRepairTime: t.insert(1, t0+4.0*slowRepairTime)
				finalDamage = odeint(self.breakKinetics, cellDamage, t, args=(DSBRate,))

				# Update state for next simulation
				cellDamage = finalDamage[-1]
				t0 = t[-1]

		return cellDamage

	# Analytic expression for misrepair probability for a given break count, eta
	def misrepairProbability(self, initialBreaks, finalBreaks, eta):
		repairedBreaks = initialBreaks-finalBreaks
		if repairedBreaks < 1E-10: return 0

		atanNum = np.sqrt(3)*eta*repairedBreaks
		atanDen = 2+eta*(2*initialBreaks*finalBreaks*eta + initialBreaks+finalBreaks) 
		return 1 - 2 * np.arctan(atanNum/atanDen) / (repairedBreaks*np.sqrt(3)*eta)

	# Analytic model of repair kinetics. Can return multiple timepoints
	# Breaks repair into multiple, smaller steps.
	# Necessary to handle heterogeneous repair and misrepair rates
	def analyticMisrepair(self, timePoints, cellDamage):
		returnConditions = []
		t = 0
		breaks = sum(cellDamage[0:3])

		# Sort timepoints, moving negative to end
		timePoints.sort()
		if timePoints[0]<0: timePoints.append(timePoints.pop(0))

		# Calculate condition at each measurement time
		for repairTime in timePoints: 
			if breaks>0:
				if repairTime>=0:
					# Update count and intra-track overlaps for each break type
					timeStep = repairTime-t
					newCellDmg = [cellDamage[n]*np.exp(-self.repairRates[n]*timeStep) 
									for n in range(3)]
					avgRepair = sum(self.repairRates[n]*cellDamage[n] for n in range(3)) / breaks
					newCellDmg += [cellDamage[n+3]
					               * np.exp(-(self.repairRates[n]+avgRepair)*timeStep)
					      		   for n in range(3)] # Calculation of intra-track overlap
				else:
					# Time<0 requests full repair, so no breaks remain
					newCellDmg = [0,0,0,0,0,0]

				# Calculate binary misrepair rate
				newBreaks = sum(newCellDmg[0:3])
				trackEta = sum(cellDamage[3:6])/(breaks*breaks)  # Average track effect per break
				totalEta = self.eta+trackEta
				misrepRate = self.misrepairProbability(breaks, newBreaks, totalEta)
				newMisrepair = misrepRate  * (breaks-newBreaks)
				newCellDmg.append(cellDamage[6] + newMisrepair)

				# Update values
				cellDamage = newCellDmg
				breaks = newBreaks
			t=repairTime

			# Calculate total misrepaired DSBs, including intrinsic fidelity and binary misrepair
			currBreaks = sum(cellDamage[0:3])
			totalRepaired = self.totalDSBs - currBreaks
			misrepairedDSB = 0
			if totalRepaired>1E-12:
				binaryMisrepairProb = cellDamage[6]/totalRepaired
				# Calculate misrepair by pathway, for intrinsic misrepair and binary misrep
				for n in range(3):
					if self.repairFidelities[n]<1:
						pathRepaired = self.totalDSBs*self.pathwayProbs[n] - cellDamage[n]
						misrepairedDSB+= pathRepaired*(1-self.repairFidelities[n])
						misrepairedDSB+= pathRepaired*binaryMisrepairProb*self.repairFidelities[n]

			# Calculate rates of aberrations and deletions from misrepair
			lethalAberrations = 0.5*misrepairedDSB*self.misrepairParams[1]
			visibleAberrations = 0.5*misrepairedDSB*self.misrepairParams[2]
			deletions = misrepairedDSB * self.misrepairParams[3]
			pointMutations = (self.dnaParams['pointMutationRate'] * (self.totalDSBs-misrepairedDSB)
				              * self.geneSize/self.dna)
			totalMut = deletions+pointMutations

			# Append these results to return conditions
			returnConditions.append([currBreaks, misrepairedDSB, lethalAberrations,
				 					 visibleAberrations, totalMut, self.totalDSBs])

		return returnConditions

	# Model DNA kinetics incorporating over time
	def modelExposure(self, exptCond, trackFoci=False):
		cellDamage = self.modelIrradiation(exptCond)

		if trackFoci: return cellDamage	# Return damage status if we just want foci
		if 'time' not in exptCond: 
			exptCond['time']=0
		repairTimes = exptCond['time']
		if type(repairTimes) is not list: repairTimes = [repairTimes]

		returnConditions = self.analyticMisrepair(repairTimes,cellDamage)
		if len(returnConditions)==1: return returnConditions[0]
		return returnConditions

	# Get foci count, using multi-exponential kinetics to account for both
	# DSB repair and foci clearance
	def getFociCount(self,exptCond):
		# 'Repair' is physical stage, 'delay' is foci clearance
		fastRepair = self.dnaParams['fastRepair']
		slowRepair = self.dnaParams['slowRepair']
		mmejRepair = self.dnaParams['verySlowRepair']
		fastDelay = self.dnaParams['fastFociDelay']
		slowDelay = self.dnaParams['slowFociDelay']
		finalTime = exptCond['time']

		# Calculate breaks at end of exposure. Fraction remaining, then scale by initial number.
		initialBreaks = self.modelExposure(exptCond, trackFoci=True)[0:3]

		fastFrac = (1.0/(fastRepair-fastDelay) 
				    *( fastRepair*np.exp(-fastDelay*finalTime)
				   	  -fastDelay*np.exp(-fastRepair*finalTime) ) )
		fastFoci = initialBreaks[0]*fastFrac

		slowFrac = (1.0/(slowRepair-slowDelay)
					*( slowRepair*np.exp(-slowDelay*finalTime)
		              -slowDelay*np.exp(-slowRepair*finalTime) ) )
		slowFoci = initialBreaks[1]*slowFrac

		mmejFoci = initialBreaks[2]*np.exp(-mmejRepair*finalTime)

		return fastFoci+slowFoci+mmejFoci

	# Calculate survival for a specific single condition
	def singleSurvival(self,exptCond,phaseDist=[0.66,0,0.34,0]):
		# Async survival is a special case. Simulate clones in different phases.
		if self.phase==-1:
			tempLine = self.cellLine.copy()
			totalSurvival = 0
			for phase in range(0,4):
				if phaseDist[phase]>0:
					tempLine['phase'] = phase
					tempCell = singleCell(tempLine,self.dnaParams,self.survParams)
					totalSurvival += tempCell.singleSurvival(exptCond)*phaseDist[phase]
			return totalSurvival

		if self.phase==0:	#G1
			# If cell has functional G1 arrest, calculate number of unrepaired 
			# breaks when released into proliferation
			tempExp = exptCond.copy()
			tempExp['time']=[exptCond['time'],-1]
			damageInfo = self.modelExposure(tempExp)
			if self.cellLine['G1Arrest']==1:
				breaksAtProliferation = damageInfo[0][0]
				apoptoticSurvival = np.exp(-breaksAtProliferation*self.survParams['apoptoticRate'])
			else:
				apoptoticSurvival = 1.0

			lethalAberrations = damageInfo[1][2]
			aberrationSurvival = np.exp(-lethalAberrations)

			totalBreaks = damageInfo[0][-1]
			baseDeath = np.exp(-totalBreaks*self.survParams['baseRate'])
			return apoptoticSurvival*aberrationSurvival*baseDeath

		if self.phase==2:	#G2
			# Model number of breaks present when cell reaches mitosis.
			# This is determined either by plating delay, or G2 damage arrest (8 hours/20 DSBs)
			tempExp = exptCond.copy()
			tempExp['time']=[max(8,exptCond['time']),-1]
			damageInfo = self.modelExposure(tempExp)
			mitoticBreaks = min(damageInfo[0][0],20)
			mitoticSurvival = np.exp(-mitoticBreaks*self.survParams['mitoticRate'])

			# Model number of lethal aberrations formed during repair
			lethalAberrations = damageInfo[1][2]
			aberrationSurvival = np.exp(-lethalAberrations)

			return aberrationSurvival*mitoticSurvival

		if self.phase==3:	#Mitosis
			# Mitotic death driven by initial DSBs
			tempExp = exptCond.copy()
			tempExp['time']=0
			mitoticBreaks = self.modelExposure(tempExp)[0]
			mitoticSurvival = np.exp(-mitoticBreaks*self.survParams['mitoticRate'])

			# Aberration-induced cell death in G1
			tempLine = self.cellLine.copy()
			tempLine['phase'] = 0
			tempCell = singleCell(tempLine,self.dnaParams,self.survParams)
			G1Survival = tempCell.survival(exptCond)
			return mitoticSurvival*G1Survival
		print("Could not identify cell phase!")
		return -1

	# Wrapper method for calculating survival. Accounts for Poisson-distributed track hits, if 
	# we're modelling ions. Otherwise models a single X-ray exposure.
	def survival(self,exptCond):
		if exptCond['dose']==0: return 1

		if 'time' not in exptCond: exptCond['time']=0 # Assume no plating delay if not specified

		# If a particle isn't defined, or is 0, return photon data
		if 'LET' not in exptCond or 'particle' not in exptCond or exptCond['particle']==0:
			return self.singleSurvival(exptCond)	

		# Calculate damage and dose per track
		# rNuc_cell is scaled by ratio of cell genetic content
		DSBPerUm = self.DSBPerkeV * exptCond['LET']
		rNuc_cell = pow(self.dna/6100.0, 1.0/3.0) * self.dnaParams['rNuc']
		DSBPerTrack = (4.0/3.0) * rNuc_cell * DSBPerUm
		dosePerTrack = DSBPerTrack / (self.dna*self.DSBPerMBPPerGy)
		tracksPerNucleus = exptCond['dose'] / dosePerTrack

		# If DSB yield per track is low, can just treat as sparsely ionising
		if DSBPerTrack<0.5: return self.singleSurvival(exptCond)

		# If not, sum over 3 standard deviations of track counts, plus unirradiated cells
		tempExp  = exptCond.copy()
		tempExp['dose'] = 0
		survTotal = self.singleSurvival(tempExp)*np.exp(-tracksPerNucleus)

		# Tracks are Poisson distributed, so deviation is sqrt(n)
		minTracks = max(1,int(tracksPerNucleus-3*np.sqrt(tracksPerNucleus)-3))
		maxTracks = int(tracksPerNucleus+3*np.sqrt(tracksPerNucleus)+3)
		trackRange = list(range(minTracks,maxTracks+1))
		trackWeights = poisson.pmf(trackRange,tracksPerNucleus)
		for n,weight in zip(trackRange,trackWeights):
			tempExp['dose']=n*dosePerTrack
			survTotal = survTotal + self.singleSurvival(tempExp)*weight

		return survTotal
