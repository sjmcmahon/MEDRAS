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

import os

trackDB = []

# Read in track-specific energy distribution for file
def readTrackFile(filename):
	try:
		with open(filename,'r') as fileData:
			# First value - are we reading DSBs (0), or energy in keV (1)?
			scaling=fileData.readline()
			while scaling[0][0]=='#' or len(scaling)==0:
				scaling=fileData.readline()
			scaling = float(scaling.rstrip(',').rstrip())
			if scaling>0:
				if scaling==1:
					doScaling=True
				else:
					print('Scaling header should be either 0: DSBs in file; or 1: Energy in keV')
					print('Read value was : ',scaling)
			else:
				doScaling=False

			# Particle type - Z. 0 is photon.
			particle = fileData.readline()
			while particle[0][0]=='#' or len(particle)==0:
				particle=fileData.readline()
			particle = particle.rstrip().rstrip(',')
			particle=int(particle)

			# LET/ DSB per micron value
			LET = fileData.readline()
			while LET[0][0]=='#' or len(LET)==0:
				LET=fileData.readline()
			LET = LET.rstrip().rstrip(',')
			LET = float(LET)

			# Read in: (r, energy) pairs
			energyData = []
			for row in fileData:
				if row[0][0]=='#' or len(row)==0:
					continue
				try:
					val = [float(x) for x in row.rstrip().split('\t')]
				except ValueError as e:
					val = [float(x) for x in row.rstrip().split(',')]
				energyData.append(val)

			retList = energyData

	# If we can't get the file, log and abort
	except IOError as e:
		print('Failed to read input file:\t\"',filename,'\"')
		print('Error: ',str(e))
		print('Not updating track database.')
		doScaling = False
		particle  = 0
		retList = [[0,0]]
		LET=0

	return doScaling,particle,LET,retList

# Insert track into database
def insertTrackFile(filename):
	doScaling,particle,LET,trackData = readTrackFile(filename)

	if LET==0 or particle==0: return

	# Append data to store for appropriate particle type
	try:
		subList = [x[0] for x in trackDB].index(particle)
		trackDB[subList][3].append([LET,trackData])
	# If we fail to find matching particle, initialise new particle sublist
	# Format is [particle,minLET,maxLET,[[List of LET, trackData]]]
	except ValueError: trackDB.append([particle,LET,LET,[[LET,trackData]]])

# Sort lists, and identify min and max LET
def tidyDB():
	trackDB.sort(key=lambda l: l[0])
	for n in range(len(trackDB)):
		trackDB[n][3].sort()
		trackDB[n][1]=trackDB[n][3][0][0]
		trackDB[n][2]=trackDB[n][3][-1][0]

# Read all energy/DSB files from a folder and import
def updateDB(folderName,stub="EnergyRange_",quiet=True):
	fileList = [f for f in os.listdir(folderName) if os.path.isfile(folderName+f)]
	fileList = [f for f in fileList if f.startswith(stub)]

	for f in fileList:
		insertTrackFile(folderName+f)
	tidyDB()

	if not quiet:
		for particle in trackDB:
			print(particle[0:3])
			for row in particle[3]:
				print(row)

# Get track data for a particular particle/LET combo
def getTrack(particle,LET):
	# If the particle is 0, we can just bail out as photons are expected to be sparsely ionising
	if particle==0: return []

	# If not, try and find particle data
	try:
		particleIndex = [x[0] for x in trackDB].index(particle)
		trackListing = trackDB[particleIndex]

		# If we have data for exactly for one energy, just linearly scale. Saves edge cases below.
		if len(trackListing[3])==1:
			return [ [x[0],(LET/trackListing[1])*x[1] ] for x in trackListing[3][0][1]]

		# If less than min LET, scale min LET data
		if LET < trackListing[1]:
			return [ [x[0],(LET/trackListing[1])*x[1] ] for x in trackListing[3][0][1]]
		# If greater than max LET, scale max LET data
		if LET > trackListing[2]:
			return [ [x[0],(LET/trackListing[2])*x[1] ] for x in trackListing[3][-1][1]]

		# If in between, find nearest matches and interpolate.
		upperIndex = next(i for i,v in enumerate(trackListing[3]) if v[0]>=LET)
		upperLET = trackListing[3][upperIndex][0]
		lowerLET = trackListing[3][upperIndex-1][0]
		scaleFx = (LET-lowerLET)/(upperLET-lowerLET)
		finalList = [[a[0],scaleFx*a[1]+(1-scaleFx)*b[1]] 
					  for a,b 
					  in zip(trackListing[3][upperIndex][1],trackListing[3][upperIndex-1][1])]
		return finalList

	except ValueError:
		print('Failed to find any track data for requested particle:\t',particle)
		return []

def dbList():
	for dataset in trackDB:
		particle,minLET,maxLET = dataset[0:3]
		print(particle)
		for row in dataset[3]:
			print(row[0])
