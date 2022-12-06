#!/usr/bin/env Rscript

# J. Hench, 2019-2020
# UMAP plot for multicore systems
# 450K and/or EPIC Illumina Methylation arrays; only 450K-equivalent probes will be considered
# from pre-filtered betas
# needs to run from Rscript with one CMD line argument, this argument is the number of top differentially methylated probes to be considere by UMAP
# v4: using integers instead of floats
# v5: try maintaining integers across programm to reduce RAM occupation
# v6: run multiple sets of UMAP in parallel on SD-sorted data

startTime <- Sys.time()
message("program start: ",startTime)
args = commandArgs(trailingOnly=TRUE)
numberOfProbesList <-as.numeric(args) # 25000 at DKFZ
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
cpuCap_read <- 0.5 #execution cap (% CPUs of all)
betasDir <- "/imagesets/external_internal_symlinks_beta_EPIC450Kmix/" #baseDir <- "/mnt/optane01/imagesets/IFP_TCGA_mix_mini_broken/" #" # e.g., baseDir <- "/mnt/optane01/imagesets/IFP_idat/20190211/"
suffix <- "_betas_filtered.rds"
outputFilePath <- "/mnt/w/mnt/1TBraid01/applications/R/conumee_ifp_basel/Conumee_01/t-SNE_data/IFP_all/UMAP6_all_bVals_top_"  # output file prefix for UMAP plot
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

message("checked directory for filtered beta value RDS files: ",Sys.time())
cl<-makeCluster(round(detectCores()*cpuCap_read),outfile="UMAP5_mini.log") # create a cluster for parallel processing
registerDoParallel(cl) # parallel processing
iterations <- nrow(betasFileDf) # number of IDATs to process
#iterations <- 100 # outcomment for testing
message("number of IDATs to process: ",iterations)
bVals <- foreach(i=1:iterations, .errorhandling = 'remove') %dopar%{
  tryCatch({
    betasRdsFileName <- paste0(betasFileDf$betasPath[i],suffix)
    c(i,as.integer(round(readRDS(betasRdsFileName)*100))) # use a memory-saving way to store beta values, first element is index and can be used later to restore sentrixID
    
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

message("process IDAT to beta values (completed): ",Sys.time())
message("cbind beta values (do.call): ",Sys.time())
bVals <- do.call(cbind,bVals)
message("cbind beta values (completed): ",Sys.time())
message("transformig beta values into dataframe: ",Sys.time())
cn <- betasFileDf$V1[bVals[1,]]
rn <- rownames(readRDS(paste0(betasFileDf$betasPath[1],suffix))) # read cg names from first RDS file in series
bVals<- as.data.frame(bVals[2:nrow(bVals),]) # the first row contains only indices for sentrix IDs
colnames(bVals) <- cn
rownames(bVals) <- rn
message("transformig beta values into dataframe (completed): ",Sys.time())

for(p in 1:length(numberOfProbesList)){
  message("Task ",p,": selecting top ",numberOfProbesList[p]," differentially methylated probes and calculating UMAP plot:", Sys.time())
}
cl<-makeCluster(length(numberOfProbesList),outfile="UMAP6_mini.log") # create a cluster for parallel processing, one process per UMAP plot
registerDoParallel(cl) # parallel processing

bVals <- foreach(numberOfProbes=numberOfProbesList, .errorhandling = 'remove') %dopar%{ #calculate UMAP plot in parallel
  message("calculate values for UMAP plot for top ",numberOfProbes, " probes (single core): ",Sys.time())
  library(umap)
  library(writexl) # to dump some dataframes etc. into generally readable files
  custom.config = umap.defaults
  bVals_top2.umap <- umap(t(bVals[names(sort(apply(bVals,1,sd), decreasing = TRUE)[1:numberOfProbes]),]), config = custom.config, method ="naive")
  bVals_top2.umap.annotated <- data.frame(rownames(bVals_top2.umap[["layout"]]),bVals_top2.umap[["layout"]])
  colnames(bVals_top2.umap.annotated)<-c("Sentrix_ID","X","Y")
  message("calculate values for UMAP plot for top ",numberOfProbes, " probes (completed): ",Sys.time())
  UMAPplotlyPlotFile <- paste0(outputFilePath,as.character(numberOfProbes)) # create path the same UMAP values
  write_xlsx(bVals_top2.umap.annotated, path = paste0(UMAPplotlyPlotFile,".xlsx"), col_names = TRUE) #save UMAP_bVals_top.data for manual lookups and as XLSX file
  save(bVals_top2.umap.annotated, file=paste0(UMAPplotlyPlotFile,".rda")) #save UMAP_bVals_top.data for manual lookups and as RDA file
}
stopCluster(cl)
message("program start was at: ",startTime)



