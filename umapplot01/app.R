# plotly plotter to create interactive tSNE/UMAP plots with shiny and CNV plots
# 2019, J. Hench, Institut für Medizinische Genetik und Pathologie Basel
# 4.2b fixed bug with global variables; created config file; switched to generic version
rm(list=ls())
gc()

versionString <- "4.2"

# DEPENDENCIES -------------------------------------------------------
library(dplyr)
library(ggplot2)
library(gplots)
library(gridExtra)
library(RSpectra)
library(reshape)
library(RefFreeEWAS)
library(RSVGTipsDevice) # for plotting with tooltips
library(SemiPar) # for plotting with tooltips
library(shiny)
library(readxl) # to read XLSX input, e.g. from google sheets
library(stringr) # to parse sheet data
library(plyr) # to merge dataframes and keep their order
library(plotly) # according to help from M. Sill, DKFZ, 2018'
library(compiler)
library(shinyBS) # Additional Bootstrap Controls
library(htmlwidgets) # detect window dimensions

# FUNCTIONS -----------------------------------------------------------------------------------------------------------
getSampleID <- function(totalPath){ #determine sample ID from EPIC path (typically the enclosing folder name)
  b<-unlist(strsplit(totalPath,"/"))
  c<-b[length(b)-1] #2nd last position is enclosing folder name
  return(c) 
}

getSentrixID <- function(totalPath){ #determine Sentrix ID from EPIC path (typically the enclosing folder name)
  b<-unlist(strsplit(totalPath,"/"))
  c<-b[length(b)] #last position is enclosing folder name
  return(c) 
}

merge.with.order <- function(x,y, ..., sort = T, keep_order){ # taken from https://www.r-statistics.com/2012/01/merging-two-data-frame-objects-while-preserving-the-rows-order/ # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
  add.id.column.to.data <- function(DATA){
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  order.by.id...and.remove.it <- function(DATA){   # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5))); gets in a data.frame with the "id..." column.  Orders by it and returns it
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]
  }
  if(!missing(keep_order)){   # tmp <- function(x) x==1; 1	# why we must check what to do if it is missing or not...  # tmp()
    if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
    if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
    warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")     # if you didn't get "return" by now - issue a warning.
  }else{
    return(merge(x=x,y=y,..., sort = sort))
  }
}
# NO FURTHER FUNCTIONS BELOW  ------------------------------------------------------------------------------------------

# INPUT PARAMETERS ----------------------------------------------------------------------------------------------------
numberOfProbes <- 25000
configData<-read.csv("config.csv")
sampleSheetKey<<-configData$sampleSheetKey[1]
cloudsFilePathPrefix<<-configData$cloudsFilePathPrefix[1]

sampleSheetURL <<- paste0("https://docs.google.com/spreadsheets/d/",sampleSheetKey,"/export?format=xlsx")
targetSheet <<- c("All_IfP_IDAT_Annotation")
cloudsFile <<- paste0(cloudsFilePathPrefix,as.character(numberOfProbes),".xlsx")
cnvUrlSuffix <<- "_CNV_IFPBasel_annotations.pdf"
# NO FURTHER MODIFYABLE / INPUT PARAMETERS BELOW ----------------------------------------------------------------------
# collect metadata for the Sentrix IDs (global)
system(paste0("wget -O LandBsamlpeSheet.xlsx ",sampleSheetURL)) # download from google drive (read access required; local installation of wget required)
annotation <<- read_excel(path="LandBsamlpeSheet.xlsx",sheet=targetSheet[1],skip=1,col_names=TRUE)
system("rm LandBsamlpeSheet.xlsx") # local installation of rm required to delete file

# SHINY APP -----------------------
ui <- fluidPage(
  titlePanel("",windowTitle="All Methylomes (UMAP) and CNV plots"),
  sidebarLayout(
    sidebarPanel(
      htmlOutput("versionView"),
      htmlOutput("UMAPtext"),
      selectInput("numberOfProbes", label = "Number of probes", choices = c("25000","50000","75000"), selected = "25000", multiple = FALSE, selectize = TRUE),
      textInput("searchbox", label = "Search annotation for:", value = "ID, Sentrix ID, diagnosis, etc."),
      actionButton("start","Start search"),
      htmlOutput("hr1"),
      sliderInput("plotYslider",label = "plot height [px]", min = 300, max=2000, value=800),
      htmlOutput("pdfViewSmall"),
      textInput("annoKey", label = "License key:", value = sampleSheetKey),
      bsTooltip(id = "numberOfProbes", title = "25000 (default), 50000, 75000", placement = "bottom", trigger = "hover"),
      bsTooltip(id = "searchbox", title = "search for ID, Sentrix ID, diagnosis, etc.", placement = "bottom", trigger = "hover"),
      bsTooltip(id = "start", title = "initiate search (can take a while)", placement = "bottom", trigger = "hover"),
      bsTooltip(id = "plotYslider", title = "adjust UMAP plot size to screen height", placement = "bottom", trigger = "hover"),
      bsTooltip(id = "pdfViewSmall", title = "A larger version of the plot (PDF file) can be accessed through the CNV tab next to UMAP Plot", placement = "bottom", trigger = "hover")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("UMAP Plot",plotlyOutput("plot",width = "100%")),
        tabPanel("CNV Profile", htmlOutput("pdfView"))
      )
    )
  )
)
server <- function(input, output, session) {
  searchBoxString<-reactiveValues(data="FooBarDoesNotExist") # see https://shiny.rstudio.com/articles/action-buttons.html
  observeEvent(input$start,{
    searchBoxString$data<-input$searchbox
  })
  output$UMAPtext <- renderPrint({
    showModal(modalDialog("Loading application and data, please wait.", footer=NULL))
    HTML(paste0("<b>UMAP plot</b> - top differentially methylated probes with search filter<hr>"))
  })
  output$plot <- renderPlotly({
    numberOfProbes <- input$numberOfProbes
    if (sampleSheetKey!=input$annoKey){ # only reload annotation if it has changed
      sampleSheetKey <<- input$annoKey
      sampleSheetURL <<- paste0("https://docs.google.com/spreadsheets/d/",sampleSheetKey,"/export?format=xlsx")
      system(paste0("wget -O LandBsamlpeSheet.xlsx ",sampleSheetURL)) # download from google drive (read access required; local installation of wget required)
      annotation <<- read_excel(path="LandBsamlpeSheet.xlsx",sheet=targetSheet[1],skip=1,col_names=TRUE)
      system("rm LandBsamlpeSheet.xlsx") # local installation of rm required to delete file
    }
    cloudsFile <<- paste0(cloudsFilePathPrefix,as.character(numberOfProbes),".xlsx")
    cloudsData <<- read_excel(path=cloudsFile,sheet=1,col_names=TRUE)
    cloudsColNames <<- colnames(cloudsData)
    cloudsDataRounded <<-data.frame(cloudsData$Sentrix_ID,round(cloudsData$X,digits = 3),round(cloudsData$Y,digits = 3)) # to search for points in plotly
    colnames(cloudsDataRounded) <<- cloudsColNames
    cloudsDataAnnotated <<- merge.with.order(cloudsData,annotation,by='Sentrix_ID', all.x = TRUE, sort=FALSE , keep_order = 1)
    cn <<- colnames(cloudsDataAnnotated)    #remove evil characters from dataframe column names
    cn <<- str_replace_all(cn, "\n", "")
    cn <<- str_replace_all(cn, " ", "")
    colnames(cloudsDataAnnotated) <<- cn
    cloudsDataAnnotated$ClassifierDiagnosis <<- with(cloudsDataAnnotated, gsub('(.{1,30})(\\s|$)', '\\1<br>', cloudsDataAnnotated$ClassifierDiagnosis)) # insert line breaks after 30 characters for legibility (https://stackoverflow.com/questions/2351744/insert-line-breaks-in-long-string-word-wrap) #does not work, SVG does not recognize such line breaks
    cloudsDataAnnotated$ClassifierDiagnosis[is.na(cloudsDataAnnotated$ClassifierDiagnosis)] <<- "-"
    cloudsDataAnnotated$MethylationClass[is.na(cloudsDataAnnotated$MethylationClass)] <<- "-"
    cloudsDataAnnotated$plotlyString <<- with(cloudsDataAnnotated,paste0(cloudsDataAnnotated$Sentrix_ID,":<br>",cloudsDataAnnotated$ClassifierDiagnosis))
    arrowlist <-list() #generate annotation for plot
    enableJIT(3)
    for (caseline in 1:length(cloudsDataAnnotated$MethylationClass)){
      # arrows
      if (grepl(searchBoxString$data,paste(toString(cloudsDataAnnotated$plotlyString[caseline]),toString(cloudsDataAnnotated$MethylationClass[caseline])))==TRUE){ # selectively annotate one (or more) MCs with arrows
        al <- list(
          x = cloudsDataAnnotated$X[caseline],
          y = cloudsDataAnnotated$Y[caseline],
          text = toString(cloudsDataAnnotated$Sample_Name[caseline]),
          xref = "x",
          yref = "y",
          font = list(family = "Courier New", size = 10),
          showarrow = TRUE,
          arrowhead = 4,
          ax = 10,
          ay = -20,
          arrowsize = 0.75
        )
        arrowlist<- c(arrowlist, list(al)) #append to list
      }
    }
    enableJIT(0)
    # plot as interactive HTML with Plotly     # dots only
    plot_ly(x=cloudsDataAnnotated$X,
            
            y=cloudsDataAnnotated$Y,
            name=cloudsDataAnnotated$MethylationClass,
            text=cloudsDataAnnotated$plotlyString,
            type="scatter",
            source="myPlot",
            height=input$plotYslider,
            marker = list(size = 4)) %>% layout(showlegend = TRUE, annotations = arrowlist)
  })  
  output$pdfView <- renderPrint({
    d <- event_data("plotly_click", source="myPlot")
    if (is.null(d)) HTML("Copy number profiles appear here; click on case dots in plot above.<hr>")
    else HTML(paste0('Copy number profile<br><iframe style="width:100%; height:500px" src="', data.frame(filter(cloudsDataRounded,X==round(as.numeric(d[3]),digits = 3) & Y==round(as.numeric(d[4]),digits = 3)))$Sentrix_ID[1] ,cnvUrlSuffix,'#view=FitV"></iframe>')) #files need to be readable from a directory called "www" within the root directory of this shiny app. It can be a symlink to another directory.
  })
  output$pdfViewSmall <- renderPrint({
    d <- event_data("plotly_click", source="myPlot")
    if (is.null(d)) HTML("<hr>Copy number profiles appear here; click on case dots in UMAP plot")
    else HTML(paste0('<hr>Copy number profile<br><iframe style="width:100%;',input$plotYslide,'px" src="', data.frame(filter(cloudsDataRounded,X==round(as.numeric(d[3]),digits = 3) & Y==round(as.numeric(d[4]),digits = 3)))$Sentrix_ID[1] ,cnvUrlSuffix,'#view=FitV"></iframe>')) #files need to be readable from a directory called "www" within the root directory of this shiny app. It can be a symlink to another directory.
  })
  output$versionView <- renderPrint({
    removeModal()
    window_height <- as.numeric(JS('window.innerHeight'))
    window_width <- as.numeric(JS('window.innerWidth'))
    #HTML(paste0("<b>Methylome & CNV Viewer</b><br>by J. Hench, IfP Basel, <i>version ",versionString,"</i><hr>window width/height: ",toString(window_width),",",toString(window_height),"<hr>"))
    HTML(paste0("<b>Methylome & CNV Viewer</b><br>by J. Hench, IfP Basel, <i>version ",versionString,"</i><hr>"))
  })
  output$hr1 <- renderPrint({
    HTML("<hr>")
  })
}

shinyApp(ui, server, options = list(display.mode = "showcase"))



