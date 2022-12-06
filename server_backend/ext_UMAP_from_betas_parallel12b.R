#!/usr/bin/env Rscript

# J. Hench, 2019-2020
# UMAP plot for multicore systems
# 450K and/or EPIC Illumina Methylation arrays; only 450K-equivalent probes will be considered
# from pre-filtered betas
# needs to run from Rscript with one CMD line argument, this argument is the number of top differentially methylated probes to be considere by UMAP
# v4: using integers instead of floats
# v5: try maintaining integers across programm to reduce RAM occupation
# v6: run multiple sets of UMAP in parallel on SD-sorted data
# v7: write out integer data as python pickle: fails with datasets >2 GB, unusable
# v8: feather for data exchange with external GPUAMP (Python on GPU-equipped remote computer)
# v9: use feather also for beta values (one file per case)

startTime <- Sys.time()
message("program start: ",startTime)
args = commandArgs(trailingOnly=TRUE)
numberOfProbesList <-as.numeric(args) # 25000 at DKFZ
#numberOfProbesList <- as.numeric(c("25000","50000","75000")) # TESTING ONLY
if (length(numberOfProbesList)<1){
  message("usage: Rscript SCRIPTNAME.R [number of probes 1] [number of probes 2] ... [[number of probes n]")
  message("Be aware that for each set of *number of probes* an R instance will be launched requiring the full memory for the entire dataset.")
  quit()
}else{
  for(p in 1:length(numberOfProbesList)){
    message("Task ",p,": considering top ",numberOfProbesList[p]," differentially methylated probes")
  }
}


# ------ LIBRARIES --------
library(doParallel)
library(plyr)
# ------ INPUT PARAMETERS --------
numberZeroFactor <- 1000 # multiplicator for conversion to integer = number of decimals after the comma
cpuCap_read <- 0.8 #execution cap (% CPUs of all)
#betasDir <- "/imagesets/external_internal_symlinks_beta_EPIC450Kmix/" #baseDir <- "/mnt/optane01/imagesets/IFP_TCGA_mix_mini_broken/" #" # e.g., baseDir <- "/mnt/optane01/imagesets/IFP_idat/20190211/"
betasDir <- "/mnt/optane/imagesets/betaFeatherEPIC450Kmix/"  # betasDir <- "/mnt/optane/imagesets/betaFeatherEPIC450K_mini_test/"
suffix <- "_betas_filtered.feather"
outputFilePath <- "/mnt/w/mnt/1TBraid01/applications/R/conumee_ifp_basel/Conumee_01/t-SNE_data/IFP_all/"  # output file prefix for UMAP plot
#gpuMachinelogin <-"user@mt202611.uhbs.ch" #mt202611 until 20210629
gpuMachinelogin <-"jhench@meqneuropat23.uhbs.ch"
#gpuMachineGpumapPath <- "/home/user/Documents/gpumap" #mt202611 until 20210629
gpuMachineGpumapPath <- "/applications/epidip03" #mt202611 until 20210629
# gpuMachineGpumapScript <- "gpumap_pyrap_08.py" # mt202611 until 20210629
gpuMachineGpumapScript <- "scripts/gpumap_pyrap_09.py"
#gpuMachineGpumapPython <- "source ~/miniconda3/etc/profile.d/conda.sh; conda activate base; cd  ~/Documents/gpumap/; python" #mt202611 until 20210629
gpuMachineGpumapPython <- "source /applications/epidip03/bin/activate; cd  /applications/epidip03/; python" #mt202611 until 20210629
gpuMachineGpumapInputFileName <- "betaValues.feather"
# ---- NO INPUT PARAMETERS BELOW --------
allFiles <- list.files(path=betasDir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
betasFileDf <- data.frame(allFiles[grepl(suffix,allFiles)]) # collect list of IDAT pathes, linux shell version, should be re-coded in pure R to work without the command line interpreter
colnames(betasFileDf) <- c("betasPath") # rename the betasFileDf colum to address it later on
betasFileDf$betasPath <- gsub(suffix,"",betasFileDf$betasPath) # remove the suffix (see INPUT PARAMETERS) since minfi requires base pathes
betasFileDf<-mdply(betasFileDf, function(betasPath){
  a<-unlist(strsplit(betasPath,"/"))
  b<-a[length(a)]
  return(b)
})

message("checked directory for filtered beta value FEATHER files: ",Sys.time())
cl<-makeCluster(round(detectCores()*cpuCap_read),outfile="UMAP9_mini.log") # create a cluster for parallel processing
registerDoParallel(cl) # parallel processing
iterations <- nrow(betasFileDf) # number of IDATs to process
#iterations <- 100 # outcomment for testing
message("number of IDATs to process: ",iterations)
bVals <- foreach(i=1:iterations, .errorhandling = 'pass') %dopar%{
  library(feather)
  tryCatch({
    betasFeatherFileName <- paste0(betasFileDf$betasPath[i],suffix)
    as.integer(round(unlist(read_feather(path=betasFeatherFileName))*numberZeroFactor))
  }, warning = function(w) {
    message("warning occurred",betasFileDf$betasPath[i])
    system (paste0("echo ",i,"/",iterations," warning ",Sys.time()))
  }, error = function(e) {
    message("error occurred",betasFileDf$betasPath[i])
    system (paste0("echo ",i,"/",iterations," error ",Sys.time()))
  }, finally = {
    message("cleanup",betasFileDf$betasPath[i])
    system (paste0("echo ",i,"/",iterations," cleanup ",Sys.time()))
  })
}
stopCluster(cl)
message("process feather to beta values (completed): ",Sys.time())
message("transformig beta values into dataframe: ",Sys.time())
bVals<-as.data.frame(bVals)
colnames(bVals)<-betasFileDf$V1 
rn <-readRDS("/mnt/optane/imagesets/betaFeatherEPIC450Kmix/cgIndex.rds") # read cg names from index 
rownames(bVals)<-rn
message("transformig beta values into dataframe (completed): ",Sys.time())
message("writing dataframe in feather format: ",Sys.time())
library(feather)
write_feather(bVals,path = "/mnt/optane/imagesets/betaFeatherEPIC450Kmix_all.feather")
message("writing dataframe in feather format (completed): ",Sys.time())
message("copying feather file to GPU computer: ",Sys.time())
rep <- system(paste0("scp ","/mnt/optane/imagesets/betaFeatherEPIC450Kmix_all.feather"," ",gpuMachinelogin,":",gpuMachineGpumapPath,"/"),intern=TRUE) # copy feather files to remote GPU computer
message(rep)
message("copying feather file to GPU computer (completed): ",Sys.time())

for(p in 1:length(numberOfProbesList)){
  message("Calculating GPUMAP for ",p,": top ",numberOfProbesList[p]," differentially methylated probes (remote GPU machine):", Sys.time())
  # example command line (bash): ssh user@mt202611.dyn.uhbs.ch "source ~/miniconda3/etc/profile.d/conda.sh; conda activate base; cd  ~/Documents/gpumap/; python gpumap_feather01.py bVals_SD_top_25000.feather embedding_top_25000.csv"; make sure the python / conda environment is activated for this ssh session
  rep <- system(paste0("ssh ",gpuMachinelogin," '",gpuMachineGpumapPython," ",gpuMachineGpumapScript," betaFeatherEPIC450Kmix_all.feather ",numberOfProbesList[p], " gpumap_",numberOfProbesList[p],".xlsx g'"),intern=TRUE)
  message(rep)
  message("Calculating GPUMAP for ",p,": top ",numberOfProbesList[p]," differentially methylated probes (done):", Sys.time())
  message("Copying result back from GPU machine:", Sys.time())
  rep <- system(paste0("scp ",gpuMachinelogin,":",gpuMachineGpumapPath,"/gpumap_",numberOfProbesList[p],".xlsx ",outputFilePath),intern=TRUE) # copy feather files to remote GPU computer
  message("Copying result back from GPU machine (done):", Sys.time())
}

message("program end: ",Sys.time())
message("program start was at: ",startTime)
