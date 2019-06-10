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
import scipy.special

class CellCharacteriser:
	# Buffer to store different parameter sets, more efficient
	storedCellParameters = {}

	# Initialise object with passed in model parameters, and cache for stored results
	def __init__(self, dnaParams):
		self.updateParameters(dnaParams)

	def updateParameters(self, dnaParams):
		self.sigma       = dnaParams['sigma']
		self.complexFrac = dnaParams['complexFrac']
		self.failFrac    = dnaParams['failFrac']

	# Functions to calculate different structural properties of DNA repair.
	# See McMahon et al, Sci Rep 2016 for details.

	# Theta is the integral of recombination probability for a single random DSB, 
	# in a cell with radius = Radius.
	# Can be defined in a sub-volume of radius r, or whole cell (2*Radius, assumed if r not given).
	def theta(self,Radius=1.0,r=float('nan')):
		if(np.isnan(r)) : r=2*Radius

		termOne   = ( 8*np.sqrt(np.pi*2)*pow(Radius,3)*
					  self.sigma*scipy.special.erf(r/(np.sqrt(2)*self.sigma)) )

		termTwo   = ( -np.exp(-r*r/(2*self.sigma*self.sigma)) 
					  * ( pow(r,4)+4*pow(r,2)*(self.sigma*self.sigma-3*Radius*Radius)
			            + 16*r*pow(Radius,3)
			            + 8*self.sigma*self.sigma*(self.sigma*self.sigma-3*pow(Radius,2)) ) )

		termThree = (8*pow(self.sigma,4)-24*Radius*Radius*self.sigma*self.sigma)

		scaling   = np.pi*(self.sigma*self.sigma)/(4*pow(Radius,3))
		thetaVal = (termOne+termTwo+termThree)*scaling
		return thetaVal

	# Calculate 'eta' value, characterising probability of misrepair with a single DSB
	def calculateMisrepairCoefficient(self,radius=1.0):
		# Calculate simple average interaction rate
		rawEta = 2*self.theta(radius) / (4.0/3.0*np.pi*pow(radius,3))

		# Distribution is highly skewed, so apply skew correction to mean
		skewBase = 0.7567
		skewRate = 5.3924
		skewCorrection = skewBase+(1-skewBase)*(1-np.exp(-skewRate*self.sigma/radius))

		return rawEta/skewCorrection

	# Calculate rates of different structural changes.
	def calculateStructuralRates(self,dna,chromosomes,phase):
		# Calculate relative chromosome radius
		chromosomeRadius = pow(chromosomes,-1.0/3.0)

		# Calculate rate of inter-chromosome misrejoining
		intraChromRate = self.theta(chromosomeRadius) / self.theta()
		interChromRate = 1-intraChromRate

		# Calculate large deletion (>3MBP) rate. 
		largeDelSize   = 3.0
		largeDelRadius = pow(largeDelSize/dna,1.0/3.0)
		delFraction    = 1-self.theta(chromosomeRadius,largeDelRadius)/self.theta(chromosomeRadius)

		largeDelRate   = intraChromRate*delFraction

		# Calculate rate of inter-arm misrepair
		# Calculate length of average chromosome based on genome
		# Need additional for duplicated DNA in S/G2 phases
		chromosomeLengthMBP = dna/chromosomes
		if phase==1: chromosomeLengthMBP = chromosomeLengthMBP/1.5
		if phase==2: chromosomeLengthMBP = chromosomeLengthMBP/2.0
		# Max distance is corresponding spatial distance. 
		maxDistance = pow(chromosomeLengthMBP/2/dna,1.0/3.0)

		lastBP = 0
		interArmRate = 0

		# Sub-divide chromosome and iterate over various distances from centromere.-
		# At each distance, test probability of misrepair spanning centromere
		for n in range(1,101):
			distToCentromere = maxDistance*(n/100.0)
			basePairs  = 2*pow(distToCentromere,3)*dna
			chromosomeFrac = (basePairs-lastBP)/chromosomeLengthMBP
			rejoinAtLeastDist = 1 - (self.theta(chromosomeRadius,distToCentromere)
									/self.theta(chromosomeRadius) )
			interArmRate = interArmRate + rejoinAtLeastDist*chromosomeFrac
			lastBP = basePairs
			# Drop out once rate is low enough
			if rejoinAtLeastDist < 0.000001: break

		return [interChromRate,largeDelRate,interArmRate]

	# Calculate mutation probability for a misrepair event in a given gene
	def calculateMutationRate(self, dna, chromosomes, geneSize, 
							  interChromRate,largeDelRate,interArmRate):
		if geneSize<=0: return 0 # If gene size is zero, no mutations simulated

		# Partial deletions are any which originate within the gene, and don't kill the cell
		partialDelRate = (geneSize/dna)*(1-0.5*(interChromRate+largeDelRate))

		# For full deletions, only consider those sizes which are non-lethal and in-chromosome
		chromosomeRadius = pow(chromosomes,-1.0/3.0)
		inChromRate = self.theta(chromosomeRadius)
		largeDelSize = 3.0
		maxBreak = max(largeDelSize-geneSize,0)
	
		# Integrate over break sizes for full deletion, calculating rate of at least that size
		fullDelRate = 0
		steps    = 50
		for breakRange in np.arange(0,maxBreak,maxBreak/steps):
			minDeletion     = breakRange+geneSize
			minRadius       = pow(minDeletion/dna,1.0/3.0)
			atLeastGeneRate = ( (1-interChromRate)
							    * (1-self.theta(chromosomeRadius,minRadius)/inChromRate) )
			nonLethalRate   = atLeastGeneRate - largeDelRate # Exclude lethally large deletions
			fullDelRate     = fullDelRate + 0.5*nonLethalRate*maxBreak/steps # Only asymmetric dels

		# Correct for direction to centromere and normalise to breaks per MBP, then to whole genome
		fullDelRate = (fullDelRate*0.5/maxBreak)/dna

		return 0.5*(fullDelRate + partialDelRate) # Final scaling by 0.5 as each involves two DSBs

	# Calculate visible and lethal aberration rates per misrepaired break
	def calculateAberrations(self,phase,interChromRate,largeDelRate,interArmRate):
		# Split by phase, due to differences in visibility and deletion
		# In mitosis, no repair until G1 phase, so treat as G1
		if phase == 0 or phase==3:
			# Only asymmetric interchromosome events (dicentrics) are visible, also lethal
			aberrationRateDic = interChromRate*0.5
			visibleDic        = aberrationRateDic

			# Only large asymmetric intrachromosome events (deletions) are visible & lethal
			aberrationRateDel = largeDelRate*0.5
			visibileDel       = aberrationRateDel
		if phase == 2:
			# In G2, all inter-chromosome events are visible. Only asymmetric are lethal.
			aberrationRateDic = interChromRate*0.5
			visibleDic        = interChromRate

			# Inter-arm are always visible, and inter-chromosome/inter-arm are lethal. 
			# Intra-arm are visible if large, but non-lethal
			# Intra-arm, intra-chromatid, asymmetric are visible
			aberrationRateDel = interArmRate
			visibileDel       = interArmRate  + largeDelRate*(1-interArmRate/(1-interChromRate))
		if phase == 1:
			# For S, we take an equal mixing of G1 and G2. Crude, but an ok approximation.
			aberrationRateDic = interChromRate*0.5
			visibleDic        = interChromRate*0.75

			aberrationRateDel = largeDelRate*0.25 + interArmRate*0.5
			visibleDel = (largeDelRate*0.25
						  +0.5*(interArmRate + largeDelRate*(1-interArmRate/(1-interChromRate)) ) )

		# Sum rates and return
		totalLethal = aberrationRateDic + aberrationRateDel
		totalVisible = visibleDic + visibileDel
		aberrRates  = [ [aberrationRateDic, aberrationRateDel, totalLethal] ,
					    [visibleDic, visibileDel, totalVisible] ]

		return aberrRates

	# Build a set of cell-line specific misrepair parameters
	def characteriseMisrepair(self, dna, chromosomes, phase, geneSize):
		# Scale dna size to phase, assume mid-S
		if phase==1 : dna = dna*1.5
		if phase==2 : dna = dna*2

		# Calculate misrepair coefficient
		misRepCoefficient = self.calculateMisrepairCoefficient()
		# Calculate rates of structural DNA changes - inter-Chromosome, large deletion, inte-arm
		structuralChanges = self.calculateStructuralRates(dna,chromosomes,phase)
		# Calculate rates of lethal and visible aberrations
		aberrationRates   = self.calculateAberrations(phase,*structuralChanges)

		# Calculate rates of misrepair deletion and point mutations
		mutationRate     = self.calculateMutationRate(dna, chromosomes, geneSize, 
													  *structuralChanges)

		# Final return - misrepair coefficient, visible aberration, lethal aberration, mutation
		return [misRepCoefficient, aberrationRates[0][2], aberrationRates[1][2], mutationRate]

	# Calculate fraction of repair by each process
	def characteriseRepair(self, repairDefect, phase):
		# Basic behaviour for normal cell
		fastRepair     = (1-self.complexFrac)
		slowRepair     = self.complexFrac
		verySlowRepair = 0

		# NHEJ defective - fast component fails, in all phases
		if repairDefect % 2 == 1:
			verySlowRepair = verySlowRepair + fastRepair*self.failFrac
			fastRepair     = fastRepair*(1-self.failFrac)

		# HR defective - if we're in S or G2, slow component fails. Otherwise no effect
		if (repairDefect // 2) % 2 == 1 and (phase==1 or phase==2):
			verySlowRepair = verySlowRepair + slowRepair*self.failFrac
			slowRepair     = slowRepair * (1-self.failFrac)

		return [fastRepair,slowRepair,verySlowRepair]

	# Get cell parameters. Cache results of calculations so we don't repeat them unnecessarily
	def getCellCharacteristics(self,cellLine):
		# Build tuple of parameters to enable cell type lookup
		parameters = (self.sigma,self.complexFrac,self.failFrac,cellLine['dna'],
			          cellLine['chromosomes'], cellLine['repair'],cellLine['phase'],
			          cellLine['gene'])

		# Lookup in cache dictionary
		if parameters in self.storedCellParameters:
			return CellCharacteriser.storedCellParameters[parameters]
		else:
			# If we can't find it, calculate parameters and log for reuse
			repairProbs = self.characteriseRepair(cellLine['repair'],cellLine['phase'])
			misrepairRates = self.characteriseMisrepair(cellLine['dna'], cellLine['chromosomes'], 
														cellLine['phase'], cellLine['gene'] )
			cellParams = [repairProbs,misrepairRates]

			CellCharacteriser.storedCellParameters[parameters]=cellParams
			return cellParams

	def resetCache(self):
		CellCharacteriser.storedCellParameters.clear()
		print('Cell caches reset')

	# Print full cell line characteristics
	def printCellCharacteristics(self,cellLine):
		print('Model parameters:')
		print(('Sigma:\t',self.sigma,'\tComplexFrac:\t',self.complexFrac,
			   '\tFailFrac:\t',self.failFrac))
		print('Cell line:')
		print(cellLine)

		if cellLine['phase']==-1:
			print ("Model parameters not defined for asynchronous cells."
				   "Consider querying individual phases")
			return

		print('\nRepair probabilities (fast/slow/very slow):')
		print(self.characteriseRepair(cellLine['repair'],cellLine['phase']))

		dna = cellLine['dna']; chromosomes = cellLine['chromosomes'];
		phase = cellLine['phase']; gene = cellLine['gene']

		if phase==1: dna = dna*1.5
		if phase==2: dna = dna*2.0

		# Calculate misrepair coefficient
		print('\nMisrepair coefficient (eta):')
		print(self.calculateMisrepairCoefficient())
		print('Rate of structural changes (inter-chrom, large deletion, inter-arm)')
		structuralChanges = self.calculateStructuralRates(dna,chromosomes,phase)
		print(structuralChanges)
		print ('Aberration rates (lethal dicentrics, lethal deletions, total lethal),'
			  '(visible dicentrics, visible deletions, total)')
		print(self.calculateAberrations(phase,*structuralChanges))
		print('Mutation rate in gene of size ',gene,' MBP:')
		print(self.calculateMutationRate(dna, chromosomes, gene, *structuralChanges))