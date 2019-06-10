# Add medras directory to path.
import sys
import os

filePath = os.path.realpath(__file__)
baseDir,thisFile = os.path.split(filePath)
parentDir,thisFolder = os.path.split(baseDir)
sys.path.insert(0, parentDir)
