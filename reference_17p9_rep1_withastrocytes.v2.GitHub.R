# install libraries ------------------------------------------------------

install.packages("remotes")
remotes::install_github('JEFworks-Lab/STdeconvolve')
remotes::install_version("Seurat", version = "4.3.0")
install.packages("hdf5r")
install.packages("tidyverse")

# load libraries ---------------------------------------------------------

library(tidyverse)
library(Seurat)
library(STdeconvolve)

# download input files ---------------------------------------------------

download.file("https://raw.githubusercontent.com/10XGenomics/analysis_guides/main/input_files/VisiumFFPE_Mouse_Brain_Transgenic_Age_17p9_Rep_1.h5", "VisiumFFPE_Mouse_Brain_Transgenic_Age_17p9_Rep_1.h5")
download.file("https://raw.githubusercontent.com/10XGenomics/analysis_guides/main/input_files/astro_markers.csv","astro_markers.csv")
download.file("https://raw.githubusercontent.com/10XGenomics/analysis_guides/main/input_files/spatial_cord_subset_17p9_rep1.csv", "spatial_cord_subset_17p9_rep1.csv")
download.file("https://raw.githubusercontent.com/10XGenomics/analysis_guides/main/input_files/optlDA.17p9_rep1_astrogenes.rds","optlDA.17p9_rep1_astrogenes.rds")
list.files(path = ".")

# functions ---------------------------------------------------------------
marker_gene_list<-function(topic=1,exp_value=2,Gexp){
  ## highly expressed in cell-type of interest
  highgexp <- names(which(Gexp[topic,] > exp_value))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(Gexp[topic,highgexp]/colMeans(Gexp[-topic,highgexp])), decreasing=TRUE)
  return(tibble(Gene=names(log2fc),log2fc=log2fc))
}

plot_spatial<-function(plot_data=plot_data,suffix1='_prop.eps',suffix2='_prop_lim1.eps',y_max=1,dir=dir,i){
  p1<-ggplot(plot_data, aes(x, y,fill = prop)) +
    geom_point(shape=21,size=4) + 
    guides(size="none")+
    labs(title=str_c("topic ",i))+
    scale_fill_viridis_c() +
    theme_void()+
    coord_equal()
  
  p2<-ggplot(plot_data, aes(x, y,fill = prop)) +
    geom_point(shape=21,size=4) + 
    guides(size="none")+
    labs(title=str_c("topic ",i))+
    scale_fill_viridis_c(limits=c(0,y_max)) +
    theme_void()+
    coord_equal()
  ggsave(plot = p1, paste(dir,"topic_",i,suffix1, sep=''), 
         height=5, width=5, units='in', dpi=300)
  ggsave(plot= p2, paste(dir,"topic_",i,suffix2, sep=''), 
         height=5, width=5, units='in', dpi=300)
}

run_me_results<-function(opt,
                         dir ){
  optLDA <- optimalModel(models = ldas, opt = opt)
  results <- getBetaTheta(optLDA,
                          perc.filt = 0.05,
                          betaScale = 1000)
  deconProp <- results$theta
  deconGexp <- results$beta
  
  dir.create(dir)
  
  for(i in 1:dim(deconProp)[2]){
    plot_data<-merge(pos,deconProp[,i],by = 0)
    names(plot_data)<-c("barcode","x","y","prop")
    
    plot_spatial(plot_data=plot_data,suffix1='_prop.eps',suffix2='_prop.jpg',y_max=1,dir=dir,i=i)
    
  }
  
  marker_gene_output<-map(.x = 1:dim(deconGexp)[1],
                          ~marker_gene_list(topic = .x,exp_value = 2,Gexp = deconGexp))
  names(marker_gene_output)<-str_c("topic_genes_exp2.",1:dim(deconGexp)[1],".csv")
  
  map2(.x = names(marker_gene_output),
       .y = marker_gene_output,
       ~write_csv(x = .y,file = paste(dir,.x)))
  
  alpha<-map(.x = 1:21,~ldas$model[[.x]]@alpha)%>%unlist
  
  plot_df<-tibble(K=2:22,alpha=alpha,perplexities=ldas$perplexities,rare=ldas$numRare)
  plot_df
  
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
  
  print(p1+p2+p3)
  
  ggsave(plot = p1+p2+p3,filename =str_c(dir,"merged_QC_plot.eps") ,height=5, width=12, units='in', dpi=300)
  ggsave(plot = p1+p2+p3,filename = str_c(dir,"merged_QC_plot.jpg"),height=5, width=12, units='in', dpi=300)
  
}

# load data and preprocess data ---------------------------------------------------------------

#note: check R_input_files folder for input data
counts<-Read10X_h5(filename = "VisiumFFPE_Mouse_Brain_Transgenic_Age_17p9_Rep_1.h5")
counts

spatial_barcodes<-read_csv("spatial_cord_subset_17p9_rep1.csv")

colnames(counts)%in%spatial_barcodes$barcode%>%summary

counts_subset <- counts[,colnames(counts)%in%spatial_barcodes$barcode]

pos<-as.data.frame(spatial_barcodes)
rownames(pos)<-pos[,1]
pos<-pos[,5:6]
names(pos)<-c("x","y")
head(pos)


counts_subset_clean <- cleanCounts(counts = counts_subset,
                                   min.lib.size = 100,
                                   min.reads = 1,
                                   min.detected = 1,
                                   verbose = TRUE)
## feature select for genes
##I am going to use the over dispersed genes and an astrocyte subset.
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
astro<-read_csv(file = "astro_markers.csv")
astro_overlap_fit<-rownames(counts_subset_clean)%in%astro$`Astrocyte Markers`
astro_overlap<-rownames(counts_subset_clean)[astro_overlap_fit]
gene_astro<-c(genes,astro_overlap)%>%unique()

corpus<-preprocess(t(as.matrix(counts_subset_clean)),
                   selected.genes = gene_astro,plot=FALSE,
                   min.reads = 1, 
                   min.lib.size = 100, 
                   min.detected = 1,
                   ODgenes = FALSE, 
                   verbose = TRUE)

ldas <- fitLDA(corpus$corpus, Ks = seq(2, 22, by = 1),
               perc.rare.thresh = 0.05,
               plot=FALSE,
               ncores=2,
               verbose=TRUE)


saveRDS(object = ldas,file = "optlDA.17p9_rep1_astrogenes.rds")
#ldas<-readRDS(file = "optlDA.17p9_rep1_astrogenes.rds")

run_me_results(opt=18,dir = "output_18_Tran_17p9_rep1_astrogenes_astro_newplots/")
#run_me_results(opt=12,dir = "output_12_Tran_17p9_rep1_astrogenes_astro/")

