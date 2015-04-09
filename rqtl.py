# Calling R/qtl scanone from python
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin & Lab Williams Memphis, Danny Arends
# last modified Mar, 2015
# first written Mar, 2015

from numpy import *
from pandas import *

import scipy as sp                          # SciPy
import rpy2.robjects as ro                  # R Objects
import pandas.rpy.common as com             # R common functions

## Get pointers to some common R functions
r_library     = ro.r["library"]             # Map the library function
r_plot        = ro.r["plot"]                # Map the plot function
r_c           = ro.r["c"]                   # Map the c function
r_class       = ro.r["class"]               # Map the class function
r_list        = ro.r["list"]                # Map the list function
r_data        = ro.r["data"]                # Map the data function
r_dim         = ro.r["dim"]                 # Map the dim function
r_rownames    = ro.r["rownames"]            # Map the rownames function

r_png         = ro.r["png"]                 # Map the dim function
r_off         = ro.r["dev.off"]                 # Map the dim function

print(r_library("qtl"))                     # Load R/qtl

## Get pointers to some R/qtl functions
scanone     = ro.r["scanone"]               # Map the scanone function
read_cross  = ro.r["read.cross"]            # Map the read.cross function

## Test1 - Use a dataset in R/qtl
res         = r_data("multitrait")          # Load the multitrait dataset
data1       = ro.r["multitrait"]            # Get a pointer to it

results1 = scanone(data1)
r_plot(results1)

## Test2 - Use a dataset from the HDD
data2 = read_cross(file = "data/multi.data.csv", format = "csvr", genotypes = r_c("AA","BB"))
results2 = scanone(data2, pheno.col="the_pheno", addcovar=)
resultsInPython = np.array(results2)


r_png("test.png")
r_plot(results2)
r_off()



## Test3 - Construct a python cross object and send it to R
pydata      = DataFrame.from_csv("data/multi.data.csv", -1)

phenotypes    = pydata[isnan(pydata[1])].iloc[:,range(2, max(pydata.columns))]
rphenotypes   = com.convert_to_r_dataframe(phenotypes)
chr1          = pydata[pydata[1] == 1].iloc[:,range(2, max(pydata.columns))]
rchr1         = com.convert_to_r_dataframe(chr1)
map1          = pydata[pydata[1] == 1].iloc[:,range(0, 2)]
rmap1         = com.convert_to_r_dataframe(map1)

rgeno  = r_list(r_list(data = rchr1,map =  rmap1))
rcross = r_list(pheno = rphenotypes,geno = rgeno)
ro.globalenv["rcross"] = rcross
ro.r('class(rcross) = c("riself","cross")')                 # Update the class to be cross and riself
rcross       = ro.r["rcross"]                               # Get a pointer to it

quit()

