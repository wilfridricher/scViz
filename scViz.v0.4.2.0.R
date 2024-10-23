rm(list=ls())

# Packages needed to run the scViz:
#install.packages("shiny",dep=TRUE)
#install.packages("shinyjs",dep=TRUE)
#install.packages("shinyWidgets",dep=TRUE)
#install.packages("colourpicker",dep=TRUE)
#install.packages("ggplot2",dep=TRUE)
#install.packages("Seurat",dep=TRUE)
#install.packages("clustree",dep=TRUE)
#install.packages("BiocManager")
#BiocManager::install("schex")
	
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(colourpicker)
library(ggplot2)
library(Seurat)
library(clustree)
library(schex)

# Define UI for data upload app ----
ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  # App title ----
  titlePanel("scViz (v0.4.2.0)"),
  h4("RShiny application to visualize single cell data"),
  h6("wilfrid.richer@curie.fr"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
		h4("Options:"),
		tabsetPanel(id = "optnavbar",
			tabPanel("Global",

	  			p(""),
	  			
	  			uiOutput("assayObs"),

	  			# Input: Select gene ----
	  			radioButtons("autocompletion", "Insert Gene(s)",
                   choices = c("Autocompletion mode" = TRUE,
                               "Exact term mode" = FALSE),
                   selected = TRUE),
      		uiOutput("genes"),
      		uiOutput("itemObs")
			),

			tabPanel("Feature",

	  			p(""),

      		# Input: Select output format ----
      		radioButtons("feature.format", "Feature format",
                   choices = c("t-SNE"	= "tsne",
                               "UMAP"	= "umap",
                               "PCA"	= "pca"),
                   selected = "umap"),

      		# Input: Specification of range within an interval ----
          #slider version
#      			sliderInput("x.range", "x axis range:",
#                  min = -50, max = 50,
#                  value = c(-10,10)),
          numericRangeInput("x.range", label = "x axis range:", value=c(-10,10),separator=" to "),
#slider version
#      			sliderInput("y.range", "y axis range:",
#                  min = -50, max = 50,
#                  value = c(-10,10)),
          numericRangeInput("y.range", label = "y axis range:", value=c(-10,10),separator=" to "),

     			# Horizontal line ----
      		tags$hr(),
	  			
      		# Input: Specification of point size ----
      		sliderInput("point.size", "size of points:",
                  min = 0.5, max = 3,
                  value = 1)

	  	),
	  		
			tabPanel("Violin",

      		# Input: Checkbox for parameters ----
      		checkboxInput("violindots", "Show cells for Violin Plots", FALSE),	  
	       					
				  # Horizontal line ----
      		tags$hr(),
	  			
	  			# Input: Select clusters ----
	  			textInput("clusters", "Selection of clusters",value = NULL)

			),

      tabPanel("Signature",
               
          # Input: Select number of genes as Ctrl for AddModuleScore() ----
          numericInput('nb.ctrlFeatures', 'Number of control features to check expression levels', 100, min = 2, max = 1000),
          
          # Input: Select number of bins for AddModuleScore() ----
          numericInput('nb.bins', 'Number of bins of aggregate expression levels', 24, min = 5, max = 50)
      )         
		)
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Show a plot of the generated distribution
     tabsetPanel(id = "mainnavbar",
      tabPanel("", icon = icon("home", lib = "glyphicon"),
	  		p(""),
      		   fileInput("file", label = "R file"), 
      		   actionButton(inputId="btnLoad","Load"),
      		   p("Be patient, it can take time..."),
      		   verbatimTextOutput("out")), 
      tabPanel("ClusterPlot", 
	  		p(""),
      		   fluidRow(
    				column(4,uiOutput("selection2")),
                    column(4,
                    checkboxInput("legend", "Include legend", TRUE),
                    checkboxInput("label", "Include labels", FALSE))
    		   ),
    		   a(id = "toggleAdvanced", "Show/hide advanced options", href = "#"),
    		   uiOutput("observData1"),
    		   plotOutput("cca"),
	  		    fluidRow(
	  		      column(6,checkboxInput("percbp", "Using Percentage", TRUE))
	  		   ),
    		   plotOutput("barplotcca"),
    		   downloadButton('downloadCca', 'Download Plot'), 
    		   downloadButton('downloadBarplotcca', 'Download Barplot')
      ), 
      tabPanel("ClustreePlot", 
	  		p(""),
	  		fluidRow(
	  		    uiOutput("resObs")
	  		   ),
	  		   plotOutput("clustree.plot"),
    		   downloadButton('downloadClustree', 'Download Clustree Plot') 
      ), 
      tabPanel("DensityPlot", 
        p(""),
            fluidRow(
                    #column(4, checkboxInput("legend_dp", "Include legend", TRUE)),
                    column(8, checkboxInput("rmFltCluster_dp", "Include advanced options from ClusterPlot", FALSE)),
                    column(4, selectInput("legend_dp", "Legend Location", 
        								  choices=c("right","left","bottom","top","none")))
             ),
             plotOutput("density.plot"),
               downloadButton('downloadDensity', 'Download Plot')),
      tabPanel("FeaturePlot", 
	  		p(""),
      		   fluidRow(
			# Selection of format of the feature plot output
               		radioButtons("format_fp", "Feature options",
                   		choices = c("by default"				  = "default",
                               		"order by expression"	= "express",
                   		            "randomly"            = "randomly",
                               		"hexagonal"				    = "schex"),
                   		selected = "default", inline=T)),
#                  	column(4, checkboxInput("legend_fp", "Include legend", TRUE)),
                  	column(8, checkboxInput("rmFltCluster_fp", "Include advanced options from ClusterPlot", FALSE)),
                  	column(4, selectInput("legend_fp", "Legend Location", 
        								  choices=c("right","left","bottom","top","none"))),
      		   fluidRow(                
      				tags$hr(),
					h5("Options for 'default' and 'order by expression' feature plot"), 
      				column(1,
      					colourInput("fp_col1", "", "gray",showColour="both"),
      					colourInput("fp_col2", "", "red",showColour="both")),
      		   		column(6,sliderInput("cutoff.scale", "Cutoff values:",
                  				min = -10, max = 20,
                  				value = c(0,15), step = 0.1))),
      		   fluidRow(                
      				tags$hr(),
      				h5("Options for 'hexagonal' feature plot"), 
                  	radioButtons("schex_action", "Calculation of bins:",
                   		choices = c("proportion of observations"= "prop_0",
                               		"means of observations"		= "mean",
                               		"medians of observations"	= "median"),
                   		selected = "mean", inline=T),
                  	radioButtons("schex_palette", "Palette of colors:",
                   		choices = c("magma"		= "magma",
                               		"inferno"	= "inferno",
                               		"plasma"	= "plasma",
                               		"viridis"	= "viridis",
                               		"cividis"	= "cividis"),
                   		selected = "magma", inline=T)                  		
               ),
      		   plotOutput("feature"),
      		   downloadButton('downloadFeature', 'Download Plot')),
      tabPanel("ViolinPlot", 
	  		p(""),
      		   fluidRow(
                  	column(4,
                  	checkboxInput("legend_vp", "Include legend", TRUE))
                ),
      		   plotOutput("violin"),
      		   downloadButton('downloadViolinPlot', 'Download Plot')
      ),        
      tabPanel("Heatmap", 
	  		p(""),
      		   fluidRow(
      				column(1, colourInput("hm_col1", "", "blue",showColour="both")),
      				column(1, colourInput("hm_col2", "", "white",showColour="both")),
      				column(1, colourInput("hm_col3", "", "red",showColour="both")),
      				column(8,
                  	checkboxInput("legend_hm", "Include legend", TRUE),
                  	checkboxInput("topgenes", "Selection of top5 highest mean ratio\n(cautious: help to open up a idea)", FALSE))
                ),
      		   plotOutput("heatmap.plot"),
      		   downloadButton('downloadHeatmap', 'Heatmap')
      )        
	 )
	)

  )
)

# Define server logic to read selected file ----
server <- function(input, output, session) {

  options(shiny.maxRequestSize=10*1024^3)

###
#  shinyjs process

  observe({
    shinyjs::hide(selector = "#mainnavbar li a[data-value=ClusterPlot]")
    shinyjs::hide(selector = "#mainnavbar li a[data-value=ClustreePlot]")
    shinyjs::hide(selector = "#mainnavbar li a[data-value=DensityPlot]")
    shinyjs::hide(selector = "#mainnavbar li a[data-value=ViolinPlot]")
    shinyjs::hide(selector = "#mainnavbar li a[data-value=FeaturePlot]")
    shinyjs::hide(selector = "#mainnavbar li a[data-value=Heatmap]")
  })
  
  observe({
	shinyjs::toggleState("btnLoad", !is.null(input$file))
  })
  
  observeEvent(input$btnLoad, {
    shinyjs::toggle(selector = "#mainnavbar li a[data-value=ClusterPlot]")
    shinyjs::toggle(selector = "#mainnavbar li a[data-value=ClustreePlot]")
    shinyjs::toggle(selector = "#mainnavbar li a[data-value=DensityPlot]")
    shinyjs::toggle(selector = "#mainnavbar li a[data-value=ViolinPlot]")
    shinyjs::toggle(selector = "#mainnavbar li a[data-value=FeaturePlot]")
    shinyjs::toggle(selector = "#mainnavbar li a[data-value=Heatmap]")
  })
   
###
#  load(file1$datapath, envir=environment())

  load_Rdata <- function(){
    if(is.null(input$file)) return(NULL)
    inFile <- isolate({ input$file })

    # Create new environment to load data into
	if(length(grep(".rdata$",inFile$datapath,ignore.case=TRUE))>0){	
    	tmp.SingleCells.Normdata <- load(inFile$datapath, envir=.GlobalEnv)
		SingleCells.Normdata <<- get(tmp.SingleCells.Normdata)
	}else{
    	SingleCells.Normdata <<- readRDS(inFile$datapath)
	}
	# Verify the version of Seurat object and transform in Seurat object version 3
	if(as.numeric(unlist(SingleCells.Normdata@version)[1])<4){
    	SingleCells.Normdata <<- UpdateSeuratObject(SingleCells.Normdata)
    }

	  if(!is.null(dim(SingleCells.Normdata@assays$integrated))){
		DefaultAssay(object = SingleCells.Normdata) <<- input$select2
		obsDataMatrix<<-SingleCells.Normdata@assays$integrated
	  }else{
		obsDataMatrix<<-SingleCells.Normdata@assays$RNA
	  }
    output$out <- renderPrint({ cat(paste("Data are loaded!\nData are composed by:\n",
      										nrow(SingleCells.Normdata@meta.data), "cells\n",
      										paste(names(table(SingleCells.Normdata@meta.data[,"orig.ident"])),collapse=",")
      							))  })
  }

  observeEvent(input$btnLoad,{
    load_Rdata()
  	updateSelectizeInput(session, "gene", server = TRUE, choices = rownames(SingleCells.Normdata@assays$RNA)) 
  	updateSelectizeInput(session,"select1", server = TRUE, 
  						 choices=names(which(sapply(colnames(SingleCells.Normdata@meta.data),
													function(x){length(table(SingleCells.Normdata@meta.data[,x]))})<=50)))
  	updateSelectizeInput(session,"select2", server = TRUE, 
  	                     choices=names(SingleCells.Normdata@assays))
  	updateSelectizeInput(session,"select4", server = TRUE, 
  	                     choices=names(table(gsub("[0-9.].+","",colnames(SingleCells.Normdata@meta.data)[grep("res",colnames(SingleCells.Normdata@meta.data))]))))
  	
  })

###  	

  output$genes <- renderUI({
	if(input$autocompletion==FALSE){
		textInput("gene", "Gene(s)", placeholder ='Type a gene symbol')
	}else{
		if(exists("SingleCells.Normdata")){
			selectizeInput("gene", label = "Gene(s)", choices = rownames(SingleCells.Normdata@assays$RNA), options = list(
			         placeholder = 'Type a gene symbol', maxOptions = 20),multiple=TRUE, selected=NULL)
		} else {
			selectizeInput("gene", label = "Gene(s)", choices = NULL, options = list(
			         placeholder = 'Type a gene symbol', maxOptions = 20),multiple=TRUE, selected=NULL)		
		}
	}
  })
  
  output$itemObs <- renderUI({
  	if(exists("SingleCells.Normdata")){
		selectizeInput("select1", label = "Information to observe",
			choices=names(which(sapply(colnames(SingleCells.Normdata@meta.data),
									   function(x){length(table(SingleCells.Normdata@meta.data[,x]))})<=50)),
			options = list(placeholder = 'None'),multiple=FALSE, selected=NULL)
	}else{
		selectizeInput("select1", label = "Information to observe",
			choices=NULL,
			options = list(placeholder = 'None'),multiple=FALSE, selected=NULL)
	}
 })

  output$assayObs <- renderUI({
    if(exists("SingleCells.Normdata")){
      selectizeInput("select2", label = "Selection of assay",
                     choices=names(SingleCells.Normdata@assays),
                     options = list(placeholder = 'None'),multiple=FALSE, selected="RNA")
    }else{
      selectizeInput("select2", label = "Selection of assay",
                     choices="RNA",
                     options = list(placeholder = 'None'),multiple=FALSE, selected="RNA")
    }
  })

  output$resObs <-  renderUI({
    if(exists("SingleCells.Normdata")){
      selectizeInput("select4", label = "Selection of assay",
                     choices=names(table(gsub("[0-9.].+","",colnames(SingleCells.Normdata@meta.data)[grep("res",colnames(SingleCells.Normdata@meta.data))]))),
                     options = list(placeholder = 'None'),multiple=FALSE, selected=NULL)
    }else{
      selectizeInput("select4", label = "Selection of assay",
                     choices=NULL,
                     options = list(placeholder = 'None'),multiple=FALSE, selected=NULL)
    }
  })
  
# ViolinPlot Panel

  violinPlot <- reactive({

	if(input$violindots==TRUE){
		point.size=1
	}else{
		point.size=0
	}
	
	selected.res=input$select1
	DefaultAssay(object = SingleCells.Normdata) <<- input$select2
	refDataMatrix=SingleCells.Normdata@assays[[input$select2]]
	slct.assay = input$select2

	if(input$clusters==""){
		select.clusters=NULL
	}else{
		select.clusters=unlist(strsplit(input$clusters,",|;|\\.|\ |\\-"))
		Idents(SingleCells.Normdata) <- SingleCells.Normdata@meta.data[,selected.res]		
	}
			
	validate(
		need(!is.null(input$gene),"Please, give a gene")
	)

  tmp.gene=unlist(strsplit(input$gene,",|;|\ "))
  	
	if(length(tmp.gene)<=1){
		gene=rownames(refDataMatrix)[which(tmp.gene==rownames(refDataMatrix))]
		if(length(gene)==0){
		  validate(need(length(gene)!=0,"No data found for this gene, give another gene"))
		}
		item=SingleCells.Normdata
	}else{
		ntmp.gene<-rownames(refDataMatrix)[unlist(lapply(tmp.gene,function(x){which(x==rownames(refDataMatrix))}))]
		if(length(ntmp.gene)<=1){
		  gene<-ntmp.gene
		  item=SingleCells.Normdata
		  if(length(gene)==0){
		    validate(need(length(gene)!=0,"No data found for this list of genes, give another list of genes"))
		  }
		}else{
		  prefix = paste(input$select2,"_", sep="")
		  modulScore <- AddModuleScore(object 	= SingleCells.Normdata,
		                               features	= list(ntmp.gene),
		                               ctrl		  = input$nb.ctrlFeatures,
		                               assay		= slct.assay,
		                               nbin     = input$nb.bins,
		                               name		  = paste(prefix,"Signature.",sep=""))
		  gene=paste(prefix,"Signature.1",sep="")
		  item=modulScore
		}
	}
	
	validate(
		need(length(gene)!=0,"Please, give a gene")
	)

	vp<-VlnPlot(object = item,
  			features  = gene,
  			pt.size = point.size,
  			idents  = select.clusters,
			group.by = selected.res)
  	
	if(length(grep("Signature",gene))==1){
	  vp <- vp + ggtitle(paste(input$select2,"_Signature (",length(ntmp.gene)," genes)", sep=""))
	}else{
	  vp <- vp + ggtitle(paste(input$select2,gene, sep="_"))
	}
	
  	if(input$legend_vp==TRUE){
		print(vp)
	}else{
		print(vp + NoLegend())
	}

  })

  output$violin <- renderPlot({
  	print(violinPlot())
  })

  output$downloadViolinPlot <- downloadHandler(
    filename = function() { paste("violinplot", format(Sys.time(), format="%H:%M:%S.%d%h%y"), "png", sep=".") },
    content = function(file) {
      ggsave(file,violinPlot(), device = "png",dpi = 600, width = 5, height = 4)
    }
  )
  
# FeaturePlot Panel

  feature <- reactive({
	
	validate(
		need(!is.null(input$gene),"Please, give a gene")
	)
	
  DefaultAssay(object = SingleCells.Normdata) <<- input$select2
  refDataMatrix = SingleCells.Normdata@assays[[input$select2]]
  slct.assay = input$select2
    
  if(input$rmFltCluster_fp==TRUE){
  	tmp.SingleCells.Normdata=subset(SingleCells.Normdata,
				cells=unlist(sapply(input$rmCluster1,function(x){which(SingleCells.Normdata@meta.data[,input$select1]==x)})))
	}else{
		tmp.SingleCells.Normdata=SingleCells.Normdata
	}

##
  tmp.gene=unlist(strsplit(input$gene,",|;|\ "))
  
	if(length(tmp.gene)<=1){
	  gene=rownames(refDataMatrix)[which(tmp.gene==rownames(refDataMatrix))]
		if(length(gene)==0){
		  validate(need(length(gene)!=0,"No data found for this gene, give another gene"))
		}
		item=tmp.SingleCells.Normdata
	}else{
		ntmp.gene<-rownames(refDataMatrix)[unlist(lapply(tmp.gene,function(x){which(x==rownames(refDataMatrix))}))]
		if(length(ntmp.gene)<=1){
			gene<-ntmp.gene
			item=SingleCells.Normdata
			if(length(gene)==0){
			  validate(need(length(gene)!=0,"No data found for this list of genes, give another list of genes"))
			}
		}else{
		  prefix = paste(input$select2,"_", sep="")
		  modulScore <- AddModuleScore(object 	= tmp.SingleCells.Normdata,
		                               features	= list(ntmp.gene),
		                               ctrl		  = input$nb.ctrlFeatures,
		                               assay		= slct.assay,
		                               nbin     = input$nb.bins,
		                               name		= paste(prefix,"Signature.",sep=""))
		  gene=paste(prefix,"Signature.1",sep="")
		  item=modulScore		
		}
	}
	
	validate(
		need(length(gene)!=0,"Please, give a gene")
	)
	
	if(input$format_fp=="schex"){
		item<-make_hexbin(item,dimension_reduction = input$feature.format, nbins=80)
		fp<-schex::plot_hexbin_feature(item, feature = gsub("rna_","",gene), type="data", mod=slct.assay, action=input$schex_action, xlab="UMAP_1", ylab="UMAP_2") + 
			scale_fill_viridis_c(option=input$schex_palette)
	}else{
		if(input$format_fp=="default"){ local.format.fp=FALSE ; slct.cells= NULL }
    if(input$format_fp=="express"){ local.format.fp=TRUE ; slct.cells= NULL }
    if(input$format_fp=="randomly"){ local.format.fp=FALSE ; slct.cells= sample(Cells(item)) }

		fp<-FeaturePlot(object = item,
						features    = gene,
						cells       = slct.cells,
						cols        = c(input$fp_col1, input$fp_col2), 
	    			reduction   = input$feature.format,
	    			min.cutoff  = input$cutoff.scale[1], max.cutoff=input$cutoff.scale[2], 
	    			order       = local.format.fp,
	    			pt.size     = input$point.size
	    			)
	}
	
	if(length(grep("Signature",gene))==1){
		fp <- fp + ggtitle(paste(input$select2,"_Signature (",length(ntmp.gene)," genes)", sep=""))
	}else{
	  fp <- fp + ggtitle(paste(input$select2,gene, sep="_"))
	}

	print(fp + 
		  coord_cartesian(xlim = c(input$x.range[1], input$x.range[2]), ylim = c(input$y.range[1], input$y.range[2]), expand = FALSE)) +
		  theme(legend.position = input$legend_fp, legend.justification='center')
	
  })

  output$feature <- renderPlot({
  	print(feature())
  })

  output$downloadFeature <- downloadHandler(
    filename = function() { paste("feature", format(Sys.time(), format="%H:%M:%S.%d%h%y"), "png", sep=".") },
    content = function(file) {
      ggsave(file,feature(), device = "png",dpi = 600, width = 5, height = 4)
    }
  )
  
  # DensityPlot Panel
  
  density.plot <- reactive({
    
    if(input$rmFltCluster_dp==TRUE){
      tmp.SingleCells.Normdata=subset(SingleCells.Normdata,
                                      cells=unlist(sapply(input$rmCluster1,function(x){which(SingleCells.Normdata@meta.data[,input$select1]==x)})))
    }else{
      tmp.SingleCells.Normdata=SingleCells.Normdata
    }
    
    tmp.data<-as.data.frame(Embeddings(object = tmp.SingleCells.Normdata, 
                                       reduction = input$feature.format))
    if(input$feature.format=="umap") tmpDensity.ggplot<-ggplot(tmp.data, aes(x = UMAP_1, y = UMAP_2))
    if(input$feature.format=="tsne") tmpDensity.ggplot<-ggplot(tmp.data, aes(x = tSNE_1, y = tSNE_2))
    if(input$feature.format=="pca") tmpDensity.ggplot<-ggplot(tmp.data, aes(x = PC_1, y = PC_2))
    
    tmpDensity.ggplot<-tmpDensity.ggplot + geom_point(colour="#00000000") +
      stat_density_2d(aes(fill = stat(level)), geom = "polygon", bins=25) +
      scale_fill_gradientn(colors = c("#4169E100","royalblue", "darkolivegreen3","goldenrod1","red")) +
      xlim(c(input$x.range[1], input$x.range[2])) + ylim(c(input$y.range[1], input$y.range[2])) + 
      theme_classic() + ggtitle("Density") +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))

    print(tmpDensity.ggplot + theme(legend.position = input$legend_dp, legend.justification='center'))
    
  })
  
  output$density.plot <- renderPlot({
    print(density.plot())
  })
  
  output$downloadDensity <- downloadHandler(
    filename = function() { paste("plot", format(Sys.time(), format="%H:%M:%S.%d%h%y"), "png", sep=".") },
    content = function(file) {
      ggsave(file,density.plot(), device = "png",dpi = 600, width = 5, height = 4)
    }
  )
  
  # ClustreePlot Panel

  clustree.plot <- reactive({

	res.motif=input$select4
	print(clustree(SingleCells.Normdata@meta.data, prefix=paste(res.motif,".",sep="")) + NoLegend())
	
  })

  output$clustree.plot <- renderPlot({
  	print(clustree.plot())
  })

  output$downloadClustree <- downloadHandler(
    filename = function() { paste("plot", format(Sys.time(), format="%H:%M:%S.%d%h%y"), "png", sep=".") },
    content = function(file) {
      ggsave(file,clustree.plot(), device = "png",dpi = 600, width = 5, height = 4)
    }
  )

# Heatmap Panel

  heatmap.plot <- reactive({

    DefaultAssay(object = SingleCells.Normdata) <<- input$select2
    Idents(object = SingleCells.Normdata) <- SingleCells.Normdata@meta.data[,input$select1]

  	if(input$topgenes==TRUE){
		tmp.colnames=levels(SingleCells.Normdata@active.ident)
		mean.matrix=matrix(NA,ncol=length(levels(SingleCells.Normdata@active.ident)),nrow=nrow(obsDataMatrix))
		colnames(mean.matrix)=tmp.colnames;rownames(mean.matrix)=rownames(obsDataMatrix)
		for(i in colnames(mean.matrix)){
			mean.matrix[,which(colnames(mean.matrix)==i)]=apply(obsDataMatrix[,which(SingleCells.Normdata@active.ident==i)],1,mean)	
		}

		#filter the 500 expressed genes with the higher value
		mean.matrix=mean.matrix[names(sort(apply(mean.matrix,1,max),decreasing=TRUE))[1:500],]

		ratio.matrix=matrix(NA,ncol=length(levels(SingleCells.Normdata@active.ident)),nrow=nrow(mean.matrix))
		colnames(ratio.matrix)=tmp.colnames;rownames(ratio.matrix)=rownames(mean.matrix)
		for(i in colnames(ratio.matrix)){
			ratio.matrix[,which(colnames(ratio.matrix)==i)]=mean.matrix[,which(tmp.colnames==i)]/apply(mean.matrix[,which(tmp.colnames!=i)],1,mean)
		}	
		
		rm(tmp.gene)
		for(i in colnames(ratio.matrix)){
			if(!exists("tmp.gene")){
				tmp.gene=names(head(sort(ratio.matrix[,which(colnames(ratio.matrix)==i)],decreasing=TRUE),5))
			}else{
				tmp.gene=c(tmp.gene,names(head(sort(ratio.matrix[,which(colnames(ratio.matrix)==i)],decreasing=TRUE),5)))
			}
		}	
	}else{
		validate(
			need(!is.null(input$gene),"Please, give a gene")
		)	
  		tmp.gene=unlist(strsplit(input$gene,",|;|\ "))
	}

	visualgenes=ceiling(sum(table(SingleCells.Normdata@active.ident))/500)
	
	for(ident in levels(SingleCells.Normdata@active.ident)){
		if(ident==levels(SingleCells.Normdata@active.ident)[1]){
			new.sample=sample(which(SingleCells.Normdata@active.ident==ident),
								ceiling(length(which(SingleCells.Normdata@active.ident==ident))/visualgenes))
			count.cells=length(new.sample)
			tmp.sample=new.sample
		}else{
			new.sample=sample(which(SingleCells.Normdata@active.ident==ident),
								ceiling(length(which(SingleCells.Normdata@active.ident==ident))/visualgenes))
			count.cells=c(count.cells,length(new.sample))
			tmp.sample=c(tmp.sample,new.sample)
		}
	}

  	if(input$legend_hm==TRUE){
		print(DoHeatmap(object = SingleCells.Normdata, features = tmp.gene, cells=names(SingleCells.Normdata@active.ident)[tmp.sample], size=0, angle=90) + 
				scale_fill_gradientn(colors = c(input$hm_col1, input$hm_col2, input$hm_col3)))
	}else{
		print(DoHeatmap(object = SingleCells.Normdata, features = tmp.gene, cells=names(SingleCells.Normdata@active.ident)[tmp.sample], angle=90) + 
				NoLegend() + scale_fill_gradientn(colors = c(input$hm_col1, input$hm_col2, input$hm_col3)))
	}
	
  })

  output$heatmap.plot <- renderPlot({
  	print(heatmap.plot())
  })

  output$downloadHeatmap <- downloadHandler(
    filename = function() { paste("plot", format(Sys.time(), format="%H:%M:%S.%d%h%y"), "png", sep=".") },
    content = function(file) {
      ggsave(file,heatmap.plot(), device = "png",dpi = 600, width = 6, height = 6)
    }
  )

# ClusterPlot Panel

  output$selection2 <- renderUI({
	selectInput("select3","Information to observe in the barplot",
		choices=names(which(sapply(colnames(SingleCells.Normdata@meta.data),
									function(x){length(table(SingleCells.Normdata@meta.data[,x]))})<=50)))
  })

  output$observData1 <- renderUI({
	
	c.list=names(table(SingleCells.Normdata@meta.data[,input$select1]))
	if(!anyNA(as.numeric(c.list))){
		c.list=as.character(sort(as.numeric(c.list)))
	}
        shinyjs::hidden(
          div(id = "advanced",
          checkboxGroupInput("rmCluster1", "cluster(s) to select",
                       c.list, 
                       selected = c.list,
                       inline = TRUE)
          )
        )
  })
  
  shinyjs::onclick("toggleAdvanced",
                     shinyjs::toggle(id = "advanced", anim = TRUE))    
  				
  cca <- reactive({
	
	#Condition to remove some clusters
	if(length(input$rmCluster1)==length(names(table(SingleCells.Normdata@meta.data[,input$select1])))){
		ccap<-DimPlot(SingleCells.Normdata, 
					  reduction = input$feature.format, group.by=input$select1, pt.size = input$point.size,
					  label = input$label, label.size = 5) 
	}else{
		ccap<-DimPlot(subset(SingleCells.Normdata,
					  cells=unlist(sapply(input$rmCluster1,function(x){which(SingleCells.Normdata@meta.data[,input$select1]==x)}))), 
					  reduction = input$feature.format, group.by=input$select1, pt.size = input$point.size,
					  label = input$label, label.size = 5) 
	}
	
	if(input$legend==TRUE){
		print(ccap + 
			coord_cartesian(xlim = c(input$x.range[1], input$x.range[2]), ylim = c(input$y.range[1], input$y.range[2]), expand = FALSE))
	}else{
		print(ccap + NoLegend() + 
			coord_cartesian(xlim = c(input$x.range[1], input$x.range[2]), ylim = c(input$y.range[1], input$y.range[2]), expand = FALSE))
	}
	
  })

  output$cca <- renderPlot({
  	print(cca())
  })

  output$downloadCca <- downloadHandler(
    filename = function() { paste("plot", format(Sys.time(), format="%H:%M:%S.%d%h%y"), "png", sep=".") },
    content = function(file) {
      ggsave(file,cca(), device = "png",dpi = 600, width = 5, height = 4)
    }
  )

  barplotcca <- reactive({

	info1=input$select1
	info2=input$select3

	validate(
		need(info1!=info2,"The two items need to be differents")
	)
	
	if(exists("memory.count")) rm("memory.count")
	for(obs1 in input$rmCluster1){
		for(obs2 in names(table(SingleCells.Normdata@meta.data[,info2]))){
	
			tmp.row=c(obs1,obs2,length(which(
					(SingleCells.Normdata@meta.data[,info1]==obs1)&
      	     		(SingleCells.Normdata@meta.data[,info2]==obs2))))
    		if(exists("memory.count")){
    			memory.count=rbind(memory.count,tmp.row)
    		}else{
    			memory.count=tmp.row
    		}
    	}
	}
	colnames(memory.count)=c(info1,info2,"sum")
	rownames(memory.count)=1:nrow(memory.count)

		new.count=as.matrix(memory.count)
		for(tmp.rescluster in names(table(new.count[,info1]))){
			new.count[which(new.count[,info1]==tmp.rescluster),"sum"]=
			(as.numeric(as.character(new.count[which(new.count[,info1]==tmp.rescluster),"sum"]))/
			sum(as.numeric(as.character(new.count[which(new.count[,info1]==tmp.rescluster),"sum"]))))*100
		}

		new.count=as.data.frame(new.count)
		if(!is.na(as.numeric(levels(factor(new.count[,info1])))[1])){
			new.count[,info1]<-factor(
				x=new.count[,info1],
				levels= sort(as.numeric(levels(factor(new.count[,info1])))))
		}
		
		if(input$percbp==TRUE){
		  ggplot(data=as.data.frame(new.count), aes(x=get(info1), y=as.numeric(as.character(sum)), fill=get(info2))) +
		    geom_bar(stat="identity") + theme(legend.position="bottom") + labs(y = "Percentage", x = info1) +
		    guides(fill=guide_legend(title=info2)) + theme_bw()
		}else{
		  ggplot(data=as.data.frame(memory.count), aes(x=get(info1), y=as.numeric(as.character(sum)), fill=get(info2))) +
		    geom_bar(stat="identity") + theme(legend.position="bottom") + labs(y = "Number of cells", x = info1) +
		    guides(fill=guide_legend(title=info2)) + theme_bw()
		}

  })
  
  output$barplotcca <- renderPlot({
  	print(barplotcca())
  })

  output$downloadBarplotcca <- downloadHandler(
    filename = function() { paste("barplot", format(Sys.time(), format="%H:%M:%S.%d%h%y"), "png", sep=".") },
    content = function(file) {
      ggsave(file,barplotcca(), device = "png",dpi = 600, width = 5, height = 4)
    }
  )
    
}

# Create Shiny app ----
shinyApp(ui, server)

