# generate copy number plots from Infinium 27K, 450K, and EPIC IDAT Methylome sets
# code modified from samples kindly provided by Martin Sill / Anne Schöler / David Capper 2018-5-25
# Modifications for Institut für Medizinische Genetik und Pathologie, USB, Basel, by Jürgen Hench, 2019
# including extension to 27K, 450K, and EPIC datasets
# v2: skipping recalculation of preexisting plots by a switch. 

# LIBRARIES for single core-----------------------------
library(doParallel)
# END LIBRARIES -----------------------------
# FUNCTIONS -----------------------------
getSentrix <- function(fullPath){
  a<-unlist(strsplit(fullPath,"/"))
  b<-a[length(a)]
  return(b)
}
# END FUNCTIONS ---------------------------
# INPUT PARAMETERS ---------------------------------------------------------------------------------------------------
#searchDir <- "/imagesets/IDAT_27K_450K_EPIC_symlinks/"
searchDir <- "/run/user/1001/gvfs/sftp:host=s1665.rootserver.io,user=jhench/home/jhench/Applications/shiny-server/IDAT_uploader01/uploads/" #"/mnt/w/home/jhench/patholabor/NGS_Diagnostik/Infinium_Bead_Chip_Methylation_Data/20191129/"
#searchDir <- "/mnt/w/mnt/1TBraid01/imagesets/GDC_TCGA_downloads/tcga_all_idats/"
#searchDir <- "/mnt/w/home/jhench/mac/Documents/sync/lab_journal/2019/data201903/mixed_27K_450K_EPIC_testset/"
outputDir <- "/imagesets/external_CNVplotsConumee/" # "/mnt/w/mnt/8TBraid01/GEO_cnvplots/" 
minIdatSize <- 7500000 # a typical 450K size is 8095252 bytes
saveCNVRds <- FALSE
yPlotRange <- 1.2
oneCopy <- 0.5

# INPUT PARAMETERS: if this is used as a module, these could be set as defaults
conumeeWorkingDirectory <- "/mnt/w/mnt/1TBraid01/applications/R/conumee_ifp_basel/Conumee_01" # define working directory (e.g. for temporary files)
conumeeRefEpicRdaPath <- "/mnt/w/mnt/1TBraid01/applications/R/conumee_ifp_basel/Conumee_01/conumee_examples/CNanalysis5_conumee_REF.2017-02-10.RData" # load refEPIC.data, provided by Martin Sill
conumeeAnnoEpicRdaPath <- "/mnt/w/mnt/1TBraid01/applications/R/conumee_ifp_basel/Conumee_01/conumee_examples/IlluminaArrayDBconumee_annotation_EPIC_B4.2017-06-07.RData" # load annoEPIC and annoEPICxy, provided by Martin Sill
conumeeEpicManifestRdaPath <- "/mnt/w/mnt/1TBraid01/applications/R/conumee_ifp_basel/Conumee_01/epicManifest.rda" #pre-loaded EPIC manifest CSV as dataframe below line 7
customAnnotationsCsv <- "/mnt/w/mnt/1TBraid01/applications/R/conumee_ifp_basel/Conumee_01/annotationGenes_27K.csv"
ref450Kfile <- "/mnt/w/mnt/1TBraid01/applications/R/conumee_ifp_basel/Conumee_01/refIDAT/ref450K"
cpuPlot <- 5
replot_all <- FALSE # re-plot all data irrespective of preexisting PDF files (necessary if one needs to update annotation)
# NO FURTHER MODIFYABLE / INPUT PARAMETERS BELOW ---------------------------------------------------------------------
grnSuffix <- "_Grn.idat" # collect all IDAT file locations within search directory based on the green files
idatGrn <- list.files(path=searchDir, pattern = paste0("*",grnSuffix), recursive = TRUE)
idatPathes <- ""
idatPathes <- data.frame(idatGrn)
idatPathes$idatPath <- with(idatPathes, paste0(searchDir,gsub(grnSuffix,"",idatPathes$idatGrn))) # cut away with suffix from IDAT path, leaving only Sentrix ID and plate position
message(">>> collected all files to plot ",Sys.time())
firstrun<-TRUE # will be set FALSE after initialzation for each core (or once in a normal loop)
message(">>> launching ",cpuPlot," parallel tasks for plotting ",Sys.time())
cl<-makeCluster(cpuPlot,outfile="CNV_log.txt") # select number of cores
registerDoParallel(cl) # parallel processing
#for (idatCounter in 1:length(idatPathes$idatPath)){
#for (idatCounter in 2:4){
foreach(idatCounter=1:nrow(idatPathes),.errorhandling = 'remove') %dopar%{
  library(conumee)
  library(devtools)
  library(DNAcopy)
  library(GenoGAM)
  library(GenomicRanges)
  library(GWASTools)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylation27kmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b3.hg19)
  library(minfi)
  library(minfiData)
  if (firstrun==TRUE){   # during first run (per core, load reference data to harmonize annotations to reference IDAT)
    message(">>> loading reference IDAT to create custom annotation ",Sys.time())
    setwd(conumeeWorkingDirectory)
    data(centromeres.hg19) # load centromere data (for primitive plot)
    load(conumeeRefEpicRdaPath) # load Epic control set and annotation, provided by Martin Sill
    load(conumeeAnnoEpicRdaPath)
    message(">>> EPIC control set and annotation loaded ",Sys.time())
    targetAnnotations <- read.csv(customAnnotationsCsv,header=TRUE) # create annotations from csv
    detail<- GRanges(targetAnnotations$targetChromosome, ranges=IRanges(targetAnnotations$targetCoordinateStart,targetAnnotations$targetCoordinateEnd), strand=targetAnnotations$targetStrand, name=targetAnnotations$targetGene) #taken from ncbi
    genome(detail) <- "hg19"
    myAnno <- CNV.create_anno(array_type = "450k", chrXY = TRUE, detail_regions=detail) # create new annotation object; will be modified by smallest possible array, i.e. 27K default bin_minprobes = 15, min_bin_size = 50000, max_bin_size = 5e+06
    message(">>> custom annotation created from CSV ",Sys.time())
    load(conumeeEpicManifestRdaPath) # load preprocessed EPIC manifest file
    message(">>> EPIC manifest RDA loaded ",Sys.time())
    myRawIntensitySet <- read.metharray(ref450Kfile, force=TRUE) # load a 27K IDAT to adjust annotations to this level: read in raw intensity signals
    myMethylSet <- convertArray(preprocessIllumina(myRawIntensitySet),"IlluminaHumanMethylationEPIC") # Create Methyl Set
    message(">>> 450K reference: generated methylation set from raw intensity values ",Sys.time())
    myCustomMethylSet <- mapToGenome(myMethylSet) # find overlaps between custom annotation and Mset
    myAnno@probes <- subsetByOverlaps(myAnno@probes, granges(myCustomMethylSet)) # myAnno will be used for all plots
    message(">>> 450K reference: determined overlaps between custom annotation and Mset ",Sys.time())
    firstrun<-FALSE 
  } # END of firstrun code
  epicPath <- idatPathes$idatPath[idatCounter] # start processing sample
  sampleID <- getSentrix(epicPath)
  outputPath <- outputDir # the shared path prefix for all output files, if they shall be stored in the same place as the input file. Example for epicPath: "/mnt/mydrive/myepicfolder/mysamplefolder/202135260076_R06C01"
  message (paste0("### processing ",sampleID," in ",epicPath," ",Sys.time()))
  if (file.info(paste0(epicPath,grnSuffix))$size>minIdatSize){
    
    if (file.exists(paste0(outputPath,sampleID,"_CNV_IFPBasel_annotations.pdf"))){
      message("PDF plot exists for ",sampleID,", skipping.")
    }else{
    
      myRawIntensitySet <- read.metharray(epicPath, force=TRUE)   # read in raw intensity signals
      message(">>> raw intensity signals loaded ",Sys.time())
      myMethylSet <- convertArray(preprocessIllumina(myRawIntensitySet),"IlluminaHumanMethylationEPIC")   # Create Methyl Set
      message(">>> generated methylation set from raw intensity values ",Sys.time())
      myCustomMethylSet <- mapToGenome(myMethylSet)   # find overlaps between custom annotation and Mset
      message(">>> determined overlaps between custom annotation and Mset ",Sys.time())
      myCNV.data <- CNV.load(myMethylSet) # load CNV data
      mysampleCust.fit <- CNV.fit(myCNV.data[sampleID],refEPIC.data, myAnno)   # find CNVs (custom annotation)
      mysampleCust.CNVs <- CNV.segment(CNV.detail(CNV.bin(mysampleCust.fit)))
      message(">>> determined CNVs ",Sys.time())
      if (saveCNVRds==TRUE){
        saveRDS(myProbesCg,file=paste0(outputPath,sampleID,"_CNV_myProbesCg.rds"))
        message(">>> saved myProbes as RDS")
        saveRDS(mysampleCust.CNVs,file=paste0(outputPath,sampleID,"_mysampleCust_CNVs.rds"))
        message(">>> saved mysampleCust.CNVs as RDS")
      }
      message(">>> plotting classical conumee annotations")
      pdf(paste0(outputPath,sampleID,"_CNV_IFPBasel_annotations.pdf"), height = 8, width = 17)
      CNV.genomeplot(mysampleCust.CNVs)
      dev.off()
      
      message(">>> done with sample ",Sys.time())
    }
  }else{
    message(">>> sample is not a 450K or EPIC idat (size too small) ",Sys.time())
  }
} # end of main for loop that goes through IDAT sets
stopCluster(cl)
message("done plotting ",Sys.time())
