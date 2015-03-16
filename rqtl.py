# Calling R/qtl scanone from python
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin & Lab Williams Memphis, Danny Arends
# last modified Mar, 2015
# first written Mar, 2015

from numpy import *
from pandas import *

import scipy as sp                        # SciPy
import rpy2.robjects as ro                # R Objects
import pandas.rpy.common as com           # R common functions

library = ro.r["library"]                 # Map the library function

print(library("qtl"))                     # Load R/qtl

scanone     = ro.r["scanone"]             # Map the scanone function
data        = ro.r["data"]                # Map the data function
read_cross  = ro.r["read.cross"]          # Map the read.cross function
plot        = ro.r["plot"]                # Map the plot function
c           = ro.r["c"]                   # Map the c function


## Test1, use a dataset in R/qtl
res         = data("multitrait")            # Load the multitrait dataset
data1       = ro.r["multitrait"]            # Get a pointer to it

results1 = scanone(data1)
plot(results1)

## Test2, use a dataset from HDD
data2 = read_cross(file = "data/multi.data.csv", format = "csvr", genotypes=c("AA","BB"))
results2 = scanone(data2)
plot(results2)

quit()


