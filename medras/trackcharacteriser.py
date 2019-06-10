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
import os

from . import trackstore
from .cellcharacteriser import CellCharacteriser

class TrackCharacteriser:
	# Buffer to store different parameter sets, more efficient
	storedTracks={}

	# Static cell characteristics to calculate aberration rates
	cellChars = CellCharacteriser({'sigma':0, 'complexFrac':0, 'failFrac':0})

	# Initialise object with passed in model parameters, and cache for stored results
	def __init__(self, params):
		self.updateParameters(params)
		if TrackCharacteriser.storedTracks=={}:
			self.updateTracks()

	def updateParameters(self, params):
		self.sigma = params['sigma']
		self.rNuc  = params['rNuc']
		return

	# Locate track data and import
	def updateTracks(self):
		filePath = os.path.realpath(__file__)
		baseDir,thisFile = os.path.split(filePath)
		parentDir,thisFolder = os.path.split(baseDir)
		dataPath = parentDir+"/Data/TrackData/"
		ions = ["Proton", "Helium", "Carbon", "Nitrogen"]
		for ion in ions:
			trackstore.updateDB(dataPath+ion+"/")

	# For a shell of radius r, randomly placed in a volume radius R, what fraction
	# is within the larger volume, on average?
	def fractionInside(self, r, R):
		if r> 2*R : return 0

		return 1-(r/R)*(3-pow(r/R,2)/4)/4

	# Calculate geometric interaction probability for a break distribution
	def calculateInteraction(self, breaksAtDistance, s, rChrom):
		# For a given distance from track core, contribution is product of number of breaks, 
		# interaction rate, and probability of being in nucleus. Sum over all data. 
		interactionRates = [np.exp(-(x[0]*x[0]) / (2*s*s)) for x in breaksAtDistance]
		inNucRates = [self.fractionInside(x[0], self.rNuc) for x in breaksAtDistance]
		totalOverlap = sum([interRate*inNuc*x[1] 
			               for interRate, inNuc, x 
			               in zip(interactionRates, inNucRates, breaksAtDistance)])

		if totalOverlap==0:
			print('Warning: Calculating overlap for empty data!')
			return [0,0,0,0,0]
		
		# Calculate in-chromosome rates, similarly to above
		inChromRates   = [self.fractionInside(x[0],rChrom) for x in breaksAtDistance]
		chromOverlap   = [interRate*inChrom*x[1] 
						  for interRate, inChrom, x 
						  in zip(interactionRates, inChromRates, breaksAtDistance)]
		intraChromFrac = sum(chromOverlap) / totalOverlap

		# Find large deletion distance by interpolating in break data, then get rate
		largeDelRadius = self.rNuc * pow(3.0/self.dna,1.0/3.0)
		upperBin       = next(i for i,v in enumerate(breaksAtDistance) if v[0]>largeDelRadius)
		fractionIntoBin= ( (largeDelRadius-breaksAtDistance[upperBin-1][0])
						   /(breaksAtDistance[upperBin][0]-breaksAtDistance[upperBin-1][0]) )
		smallDelOverlap= sum(chromOverlap[0:upperBin])+chromOverlap[upperBin]*fractionIntoBin
		smallDelRate   = smallDelOverlap/sum(chromOverlap)

		# Find inter-arm rate - max separation is chromosome volume
		centromereSpan = [min(1,pow(r,3)/(0.5*pow(rChrom,3))) for r,e in breaksAtDistance]
		inChromSpan    = [span*inChrom*interRate*x[1] 
						  for span,inChrom,interRate,x 
						  in zip(centromereSpan,inChromRates,interactionRates,breaksAtDistance)]
		interArmTotal  = sum(inChromSpan)
		interArmRate   = interArmTotal / sum(chromOverlap)

		# No mutation calculation yet!
		mutationRate = 0

		# Return total overlap, inter-chromosome, deletion, inter-arm, and mutation rates
		trackRates = [2*totalOverlap, (1-intraChromFrac), (1-smallDelRate)*intraChromFrac,
				      interArmRate, mutationRate]
		return trackRates

	# Calculate intra-track effects for a given ion track.
	# Return DSB per micron, and various aberration probabilities
	def characteriseTrack(self):
		if self.particle==0: return 0,[0,0,0,0,0]

		energyAtDistance = trackstore.getTrack(self.particle,self.LET)

		# Calculate rate of DSBs based on human nucleus data
		# 0.03825 is scaling for 1 Gray deposited in a sphere of radius 1 micron
		DSBPerGyHuman = 35.0
		DSBPerkeV = (0.03825*DSBPerGyHuman)/pow(self.rNuc,3)

		breaksAtDistance = [[x[0],x[1]*DSBPerkeV] for x in energyAtDistance]
		DSBPerUm         = self.LET * DSBPerkeV

		s = self.rNuc * self.sigma
		rChrom = self.rNuc/pow(self.chromosomes,1.0/3.0)

		# Calculate track-specific interaction rates, and resulting aberrations
		trackRates = self.calculateInteraction(breaksAtDistance,s,rChrom)
		aberrRates = TrackCharacteriser.cellChars.calculateAberrations(self.phase,*trackRates[1:4])

		trackRates.append(aberrRates)

		return DSBPerUm, trackRates

	# Wrapper to cache ion track results, for efficiency
	def getTrackCharacteristics(self, exptCond, cellLine):
		self.dna = cellLine['dna']
		self.chromosomes = cellLine['chromosomes']
		self.phase = cellLine['phase']
		if 'LET' in exptCond and 'particle' in exptCond:
			self.LET = exptCond['LET'] 
			self.particle = exptCond['particle']
		else:
			self.LET = 0
			self.particle=0

		params = (self.sigma, self.rNuc, self.dna, self.chromosomes, self.phase,
			      self.LET, self.particle)

		if params in TrackCharacteriser.storedTracks:
			return TrackCharacteriser.storedTracks[params]
		else:
			trackParams = self.characteriseTrack()
			TrackCharacteriser.storedTracks[params]=trackParams
			return trackParams

	def reset(self):
		storedTracks.clear()
		print('Track cache reset')
