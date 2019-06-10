# ############################################################################
# 
# This software is made freely available in accordance with the simplifed BSD
# license:
# 
# Copyright (c) <2018>, <Stephen McMahon>
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
#
# Here, we seek to provide a basic test suite to determine if model is stable
# Runs basic fitting, gets output, and calculates parameter values and 
# correlations. 
#
# ############################################################################

import numpy as np
import scipy
import datetime

import context
from medras import medrascell

# Plot dose-response curves for cell cycle phases
def cellSurvivalCurveTest():
	doses = np.arange(0,10,0.5)
	phaseNames = ['G1','G2','M','Async']
	phaseVals  = [  0,   2,  3,      -1]

	print('\nDose response curves by cell cycle phase')
	print('Dose\t','\t'.join(map(str,doses)))
	for phase,phaseVal in zip(phaseNames,phaseVals):
		print(phase,end='\t')
		cellType = {'dna':6100, 'phase':phaseVal, 'chromosomes':46, 
					'repair':0, 'G1Arrest':1, 'gene':0}
		theCell = medrascell.singleCell(cellType)

		for d in doses:
			exptCond = {'dose':d,'time':0}
			print("%.5f" % theCell.survival(exptCond), end='\t')
		print()

# Plot dose-response curves for different LETs
def LETSurvivalCurveTest():
	doses = np.arange(0,10,0.5)
	LETs = [0.1,1,10,20,50,100,200,500]

	print('\nDose resposne curves by LET')
	print('Dose\t','\t'.join(map(str,doses)))

	cellType = {'dna':6100, 'phase':0, 'chromosomes':46, 
				'repair':0, 'G1Arrest':1, 'gene':0}
	theCell = medrascell.singleCell(cellType)
	for L in LETs:
		print(str(L)+"keV/um",end='\t')
		for d in doses:
			exptCond = {'dose':d,'time':0,'LET':L,'particle':6}
			print("%.5f" % theCell.survival(exptCond),end='\t')
		print()

# Plot kinetics of other endpoints
def repairKineticsTest():
	doses = [1, 2, 4, 6, 8]
	repairTimes = [0, 0.5,1,  2, 4, 8,  24, -1]
	print('\nImpact of repair time on mechanistic endpoints')
	endPoints = [0,1,2,3]
	endPointNames = ['Unrepaired DSBs','Misrepaired DSBs','Lethal Aberrations','Visible Aberrations']
	cellType = {'dna':6100, 'phase':0, 'chromosomes':46, 
				'repair':0, 'G1Arrest':1, 'gene':0}	
	theCell = medrascell.singleCell(cellType)
	for e in endPoints:
		print('Endpoint:\t',endPointNames[e])
		print('Time (h)\t','\t'.join([str(d)+' Gy' for d in doses]))
		for t in repairTimes:
			print(t,'\t', end=' ')
			for d in doses:
				exptCond = {'dose':d,'time':t}
				print("%.3f" % theCell.modelExposure(exptCond)[e], end='\t')
			print()
		print()

if __name__ == "__main__":
	cellSurvivalCurveTest()
	LETSurvivalCurveTest()
	repairKineticsTest()
