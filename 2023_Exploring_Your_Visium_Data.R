# install library commands ------------------------------------------------------

## These are the libraries you will need to perform the example analysis. 

install.packages("remotes")
remotes::install_github('JEFworks-Lab/STdeconvolve')
remotes::install_version("Seurat", version = "4.3.0")
install.packages("hdf5r")
install.packages("tidyverse")

# load libraries ---------------------------------------------------------

## In this section we load the required library files. 

library(tidyverse)
library(Seurat)
library(STdeconvolve)

# download input files ---------------------------------------------------

## In this section we download the required input files for the analysis from the 10x Genomics Analysis Guides GitHub repository.

download.file("https://raw.githubusercontent.com/10XGenomics/analysis_guides/main/2023_Exploring_Your_Visium_Data_input_files/VisiumFFPE_Mouse_Brain_Transgenic_Age_17p9_Rep_1.h5", "VisiumFFPE_Mouse_Brain_Transgenic_Age_17p9_Rep_1.h5")
download.file("https://raw.githubusercontent.com/10XGenomics/analysis_guides/main/2023_Exploring_Your_Visium_Data_input_files/astro_markers.csv","astro_markers.csv")
download.file("https://raw.githubusercontent.com/10XGenomics/analysis_guides/main/2023_Exploring_Your_Visium_Data_input_files/spatial_cord_subset_17p9_rep1.csv", "spatial_cord_subset_17p9_rep1.csv")
download.file("https://raw.githubusercontent.com/10XGenomics/analysis_guides/main/2023_Exploring_Your_Visium_Data_input_files/optlDA.17p9_rep1_astrogenes.rds","optlDA.17p9_rep1_astrogenes.rds")
list.files(path = ".")

# define functions ---------------------------------------------------------------

## In this section we define functions used in the analysis. 

## This function returns a tibble containing highly expressed genes and their log2 fold change. 
## Input - gene expression matrix from the optimal lda model (Gexp), topic number from the optimal lda model (topic), gene expression cut-off for considering genes highly expressed (exp_value)

marker_gene_list<-function(topic,exp_value,Gexp){
  ## highly expressed in cell-type of interest
  highgexp <- names(which(Gexp[topic,] > exp_value))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(Gexp[topic,highgexp]/colMeans(Gexp[-topic,highgexp])), decreasing=TRUE)
  return(tibble(Gene=names(log2fc),log2fc=log2fc))
}

## This function generates a spatial plot colored by topic proportion.
## Input - a table consisting of the spot proportions for a topic of interest (plot_data), a suffix defining the image type (suffix1), output directory (dir), and the topic of interest number (i)   

plot_spatial<-function(plot_data=plot_data,suffix1='_prop.jpg',dir=dir,i){
  p1<-ggplot(plot_data, aes(x, y,fill = prop)) +
    geom_point(shape=21,size=4) + 
    guides(size="none")+
    labs(title=str_c("topic ",i))+
    scale_fill_viridis_c() +
    theme_void()+
    coord_equal()
  
  ggsave(plot = p1, paste(dir,"topic_",i,suffix1, sep=''), 
         height=5, width=5, units='in', dpi=300)
 
}

## This is a plotting function that generates a QC plot for fitted models.
## Input -  the model object (ldas), starting K (starting_k) and ending K (ending_k) value. 

run_me_QC<-function(ldas,starting_k=2,ending_k=22,dir){
  if(!dir.exists(dir)){dir.create(dir)}
  alpha<-map(.x = (starting_k-1):(ending_k-1),~ldas$model[[.x]]@alpha)%>%unlist
  
  plot_df<-tibble(K=starting_k:ending_k,alpha=alpha,perplexities=ldas$perplexities,rare=ldas$numRare)
  
  p1<-ggplot(data = plot_df) +
    geom_line(mapping = aes(x = K,y = perplexities), color="red3",size=2) +
    geom_point(mapping = aes(x = K,y = perplexities),shape=21, color="black", fill=ifelse(alpha > 1, "white", "red3"), size=6)+
    theme_linedraw(base_size = 16,base_rect_size =2,base_line_size = 2)+ylab("perplexity")+
    ylim(min(plot_df$perplexities)-10,10+max(plot_df$perplexities))
  p2<-ggplot(data = plot_df) +
    geom_point(mapping = aes(x = K,y = rare),shape=21, color="black", fill="blue", size=4)+ 
    geom_line(mapping = aes(x = K,y = rare), color="blue",size=2)+
    theme_linedraw(base_size = 16,base_rect_size =2,base_line_size = 2)+ylab("cell−types with mean proportion < 5%")
  
  p3<-ggplot(data = plot_df) +
    geom_line(mapping = aes(x = K,y = alpha), color="darkgreen",size=2) +ylim(c(0,1))+
    geom_point(mapping = aes(x = K,y = alpha),shape=21, color="black", fill=ifelse(alpha > 1, "white", "darkgreen"), size=6)+
    theme_linedraw(base_size = 16,base_rect_size =2,base_line_size = 2)+ylab("alpha")
 
  ggsave(plot = p1+p2+p3,filename = str_c(dir,"merged_QC_plot.jpg"),height=5, width=12, units='in', dpi=300)
  
}

## This is a wrapper function that generates the spatial plots for each topic for the optimal model and exports the log2 fold change for highly expressed genes for each topic for the optimal model.
## Input - the optimal model number (opt), lda object returned from the fitLDA function (ldas), and an output directory (dir) 

run_me_results<-function(opt,
                         dir,ldas ){
  optLDA <- optimalModel(models = ldas, opt = opt)
  results <- getBetaTheta(optLDA,
                          perc.filt = 0.05,
                          betaScale = 1000)
  deconProp <- results$theta
  deconGexp <- results$beta
  
  if(!dir.exists(dir)){dir.create(dir)}
  
  for(i in 1:dim(deconProp)[2]){
    plot_data<-merge(pos,deconProp[,i],by = 0)
    names(plot_data)<-c("barcode","x","y","prop")
    
    plot_spatial(plot_data=plot_data,suffix1='_prop.jpg',dir=dir,i=i)
    
  }
  
  marker_gene_output<-map(.x = 1:dim(deconGexp)[1],
                          ~marker_gene_list(topic = .x,exp_value = 2,Gexp = deconGexp))
  names(marker_gene_output)<-str_c("topic_genes_exp2.",1:dim(deconGexp)[1],".csv")
  
  map2(.x = names(marker_gene_output),
       .y = marker_gene_output,
       ~write_csv(x = .y,file = paste(dir,.x)))
  
 
}

# load and preprocess data ---------------------------------------------------------------

## In this section we load and preprocess the data.
## note: check R_input_files folder for input data

# loading input data

counts<-Read10X_h5(filename = "VisiumFFPE_Mouse_Brain_Transgenic_Age_17p9_Rep_1.h5")
spatial_barcodes<-read_csv("spatial_cord_subset_17p9_rep1.csv")

# subset the input data to focus on the region of interest

counts_subset <- counts[,colnames(counts)%in%spatial_barcodes$barcode]
pos<-as.data.frame(spatial_barcodes)
rownames(pos)<-pos[,1]
pos<-pos[,5:6]
names(pos)<-c("x","y")

## filter count matrix to remove low quality spots and poorly expressed genes

counts_subset_clean <- cleanCounts(counts = counts_subset,
                                   min.lib.size = 100,
                                   min.reads = 1,
                                   min.detected = 1,
                                   verbose = TRUE)
## selecting genes for the model
## first, the overdispersed genes are determined 

odGenes <- getOverdispersedGenes(as.matrix(counts_subset_clean),
                                 gam.k=5,
                                 alpha=0.05,
                                 plot=FALSE,
                                 use.unadjusted.pvals=FALSE,
                                 do.par=TRUE,
                                 max.adjusted.variance=1e3,
                                 min.adjusted.variance=1e-3,
                                 verbose=FALSE, details=TRUE)


genes <- odGenes$ods
length(genes)

## second, we load a list of canonical astrocyte markers and then merge them with the list of overdispersed genes

astro<-read_csv(file = "astro_markers.csv")
astro_overlap_fit<-rownames(counts_subset_clean)%in%astro$`Astrocyte Markers`
astro_overlap<-rownames(counts_subset_clean)[astro_overlap_fit]
gene_astro<-c(genes,astro_overlap)%>%unique()

## the merged overdispersed gene and canonical astrocyte marker list is used to generate the corpus for the model
corpus<-preprocess(t(as.matrix(counts_subset_clean)),
                   selected.genes = gene_astro,plot=FALSE,
                   min.reads = 1, 
                   min.lib.size = 100, 
                   min.detected = 1,
                   ODgenes = FALSE, 
                   verbose = TRUE)

# fit LDA model for a range of topics (K values) --------

## In this section we fit a series of LDA models for different K values

## model fitting for a range of K values

ldas <- fitLDA(corpus$corpus, Ks = seq(2, 22, by = 1),
               perc.rare.thresh = 0.05,
               plot=FALSE,
               ncores=2,
               verbose=TRUE)

## exporting QC metrics for the models

run_me_QC(ldas,starting_k=2,ending_k=22,dir="output_18/")

# the following # commented commands allow you to save the "ldas" variable (line 186 above) 
# to an RDS file
# more information about RDS files: 
# https://rstudio-education.github.io/hopr/dataio.html#saving-r-files
# remove the # (comment) to run the command:

# saveRDS(object = ldas,file = "optlDA.17p9_rep1_astrogenes.rds")

# here is the command to load our RDS file (downloaded on line 26 above) 
# remove the # (comment) to run the command:

# ldas<-readRDS(file = "optlDA.17p9_rep1_astrogenes.rds")

## exporting spatial plots and topics for the optimal model
## This function assume a gene expresison cut-off of 2.

run_me_results(opt=18,dir = "output_18/",ldas=ldas)

