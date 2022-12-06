#!/usr/bin/env python
# coding: utf-8

# UMAP plot with GPU assistance
# J. Hench IfB Basel 2020
# IDAT parsing in R with feather output
# 
# Version 4: tried to invoke GPUMAP serval times within the session (GPUMAP v. 0.1.1), but also the example code for GPUMAP fails when re-invoking with a CUDA memory exception (CUDA 10.2, NVIDIA Driver 440.33.01). Not following up on this error for the moment and sticking with one round of calculation per process.
# Version 5: running with one set of top methylated probes accepting only 1 cmdline argument for the number of topN diff. meth. probes
# Version 6: implementation of subset selection through list of Sentrix_IDs (CSV file with header "Sentrix_ID", one column)
# version 8: SD modes
# version 9: encountered invalid IDAT for the first time! Contains NaN value. Filling with zero.

# configuration data
chunksize = 5000 # how many probes to transfer to GPU RAM

# dependencies
import warnings
warnings.simplefilter(action='ignore') # turn off warnings
import sys
import datetime
import feather
import numpy
import pandas
import math
#import gpumap / umap: this is being imported later depending on CMD line choice; will also run on CPU-only systems when user selects CPU mode

# read command line arguments
print ("Number of arguments: ",len(sys.argv))
if len(sys.argv) != 5 and len(sys.argv) != 6:
    raise ValueError('\n\nUsage: [THISPYTHONSCRIPT] INPUT_FILE_PATH.feather TOP_N_PROBES OUTPUT_FILE_PATH.xlsx [mode: g=GPU, c=CPU=default] [optional: SELECTIONFILE.csv - name can contain SDsubset/SDall argument to toggle SD calculation mode, SDsubset is default]\n\nPlease provide input file name and number of top differentially methylated probes to be considered.\n\n')

bVals_file=sys.argv[1] # input file path
topN = int(sys.argv[2]) # top N differentially methylated probes
umap_file=sys.argv[3] # output file path
process_mode=sys.argv[4] # processor mode, "g"=GPU, "c"=CPU
sdMode=0 # SDsubset mode (default)

print("UMAP start: ", datetime.datetime.now())
bVals_all = feather.read_dataframe(bVals_file) # load beta values
if len(sys.argv) == 6:
    print("selecting subset of cases start: ", datetime.datetime.now())
    selectedCases = pandas.read_csv(sys.argv[5]) # selection of Sentrix IDs
    print("selecting subset of cases end  : ", datetime.datetime.now())
    if sys.argv[5].count("SDall")>0:
        sdMode=1 # SDall mode
    else:
        sdMode=0 # SDsubset mode
        bVals_all = bVals_all[selectedCases['Sentrix_ID']]

print("SD start: ", datetime.datetime.now())
numberCpg,numberSamples = bVals_all.shape # determine number of rows and columns
cycles = math.ceil(numberCpg/chunksize) # determine SD calculation cycles
print("SD on CPU start: ",datetime.datetime.now())
print ("chunksize: ",chunksize)
for c in range(0,cycles):
    s = c*chunksize
    e = s+chunksize-1
    if e>=numberCpg:
        e=numberCpg-1
    if c==0:
        bVals_SD = bVals_all.loc[s:e].std(axis=1) # row-wise standard deviation, first cycle, make new cuDF dataframe
    else:
        bVals_SD = bVals_SD.append(bVals_all.loc[s:e].std(axis=1),ignore_index=True) # row-wise standard deviation, append data to cuDF dataframe
print("SD on CPU end: ",datetime.datetime.now())
print("SD sort on CPU start:",datetime.datetime.now())
bVals_SD = bVals_SD.sort_values(ascending=False) # sort the standard deviation in descending order
print("SD sort on CPU end:  ",datetime.datetime.now())
print("SD subset top N proces on CPU start:",datetime.datetime.now())
bVals_subset_index = bVals_SD[0:topN].index # get index of topN probes
print("SD subset top N proces on CPU end:  ",datetime.datetime.now())
if sdMode==1: # subset of cases but SD against all
    bVals_all = bVals_all[selectedCases['Sentrix_ID']] # change selection of cases
print("Transpose subset on CPU start:",datetime.datetime.now())
bVals_subset = bVals_all.loc[bVals_subset_index].T # get a subset of topN probes and transpose dataframe
print("Transpose subset on CPU end:  ",datetime.datetime.now())
del bVals_all # no longer needed in this session
print("Replace NaN with 0 start:  ",datetime.datetime.now())
bVals_subset = bVals_subset.fillna(0) # replace NaN with 0 - single occurance in 2 years of operations
print("Replace NaN with 0 end:  ",datetime.datetime.now())
if process_mode=="g":
    print("Processor mode: GPU")
    import gpumap # GPUMAP dependency
    print("gpumap.GPUMAP().fit_transform(bVals_subset) start: ", datetime.datetime.now())
    embedding = gpumap.GPUMAP().fit_transform(bVals_subset)
    print("gpumap.GPUMAP().fit_transform(bVals_subset) end: ", datetime.datetime.now())    
else:
    print("Processor mode: CPU")
    import umap # UMAP dependency
    print("umap.UMAP().fit_transform(bVals_subset) start: ", datetime.datetime.now())
    embedding = umap.UMAP().fit_transform(bVals_subset)
    print("umap.UMAP().fit_transform(bVals_subset) end: ", datetime.datetime.now())    
print("writing result to drive in XLSX format start:" , datetime.datetime.now())
pandas.DataFrame(embedding, index=bVals_subset.index.values).to_excel(umap_file, header=["X","Y"], index=True, index_label="Sentrix_ID", sheet_name='Sheet1')  # re-link row labels, GPUMAP omits them
print("writing result to drive in XLSX format end: " , datetime.datetime.now())
print("UMAP end: ", datetime.datetime.now())

