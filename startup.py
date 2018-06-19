# The purpose of this file is to load in all of the packages and 
# variables that we are always using to complete this project.
# By loading in this file at the start of each script, we reduce
# the need to copy/paste and continually update each individual
# script.
#
# Data directory: This is the ONE line that will be custom to 
# each individual person. Update it do reflect where your data is stored.
datadir='/Users/jkamenetzky3/data/'
#
# Import packages
# General packages 
import numpy as np
import matplotlib.pyplot as plt
# Fits and photometry
from astropy.io import fits