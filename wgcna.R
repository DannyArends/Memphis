# Calling WGCNA from python
#
# copyright (c) 2014-2020 Danny Arends - Brockmann group - HU Berlin & Lab Williams Memphis
# last modified Mar, 2015
# first written Mar, 2015

## Install package into the R environment

#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
#biocLite("preprocessCore") # Seems to be an unlisted missing dependency
#install.packages("WGCNA") 

## Download the example dataset from: http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-Data.zip
# Clean the data using R, so that it doesn't contain any annotation (just row and column names)
# options(stringsAsFactors = FALSE);
# femData = read.csv("LiverFemale3600.csv")

# datExpr0 = as.data.frame(t(femData[, -c(1:8)]))
# names(datExpr0) = femData$substanceBXH
# rownames(datExpr0) = names(femData)[-c(1:8)]
# write.table(datExpr0,"liverFemale_cleaned.csv",sep="\t")

## Python code to call WGCNA

from numpy import *
from pandas import *

import scipy as sp                            # SciPy
import rpy2.robjects as ro                    # R Objects
import pandas.rpy.common as com               # R common functions

from rpy2.robjects.packages import importr
utils = importr("utils")

## Get pointers to some common R functions
r_library       = ro.r["library"]             # Map the library function
r_options       = ro.r["options"]             # Map the options function
r_read_csv      = ro.r["read.csv"]            # Map the read.csv function
r_dim           = ro.r["dim"]                 # Map the dim function
r_c             = ro.r["c"]                   # Map the c function
r_seq           = ro.r["seq"]                 # Map the seq function
r_table         = ro.r["table"]               # Map the table function
r_names         = ro.r["names"]               # Map the names function


print(r_library("WGCNA"))                     # Load WGCNA
print(r_options(stringsAsFactors = False))

## Get pointers to some WGCNA functions
r_enableWGCNAThreads  = ro.r["enableWGCNAThreads"]              # Map the enableWGCNAThreads function
r_pickSoftThreshold   = ro.r["pickSoftThreshold"]               # Map the pickSoftThreshold function
r_blockwiseModules    = ro.r["blockwiseModules"]                # Map the blockwiseModules function
r_labels2colors       = ro.r["labels2colors"]                   # Map the labels2colors function
r_plotDendroAndColors       = ro.r["plotDendroAndColors"]       # Map the plotDendroAndColors function


## Main code
print(r_enableWGCNAThreads())                                 # Enable multi threading

# Load data
testData = r_read_csv("liverFemale_cleaned.csv", sep="\t")    # Test data was cleaned in R (see code above)
print(r_dim(testData))                                        # Should print: [1]  135 3600

# Calculate a good soft threshold
powers = r_c(r_c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), r_seq(12, 20, 2))
sft    = r_pickSoftThreshold(testData, powerVector = powers, verbose = 5)

# Create block wise modules using WGCNA
net = r_blockwiseModules(testData, power = 6, TOMType ="unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = True, pamRespectsDendro = False, saveTOMs = True, saveTOMFileBase = "femaleMouseTOM", verbose = 3)

# How many modules and how many gene per module ?
print(r_table(net[1]))

# The iconic WCGNA plot of the modules in the hanging tree
mergedColors = r_labels2colors(net[1])
r_plotDendroAndColors(net[5][0],mergedColors, "Module colors",dendroLabels = False, hang = 0.03, addGuide = True, guideHang = 0.05)