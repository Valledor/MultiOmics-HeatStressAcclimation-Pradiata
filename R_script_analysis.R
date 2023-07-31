#### Import dataset and loading processomics ----
#### Load independent datasets
transcripts <- read.csv("transcripts.csv",row.names = 1)
proteins <- read.csv("proteins.csv",row.names = 1)
metabolites <- read.csv("metabolites.csv",row.names = 1)
physiology <- read.csv("physiology.csv",row.names = 1)

#### Create object with all datasets for working with pRocessomics
original.dataset <- list(transcripts=transcripts, proteins=proteins, metabolites=metabolites, physiology=physiology)
class(original.dataset) <- "pRoDS"

#### Load pRocessomics (available at: https://github.com/Valledor/pRocessomics). Version 0.1.14 - 19.01.2023
library(pRocessomics)

#### Preprocess the multiomic dataset ----

#### Preprocess of the dataset using the wizard
preprocess_wizard(original.dataset)
  # Missing value imputation: Yes
  # Imputation threshold: 0.34 (1/3)
  # Imputation algorithm: KNN
  # Number of neighbours: 3
  # Apply consistency pre-filter: Yes
  # Consistency threshold: 0.45 (4/9)
  # Dataset balancing: Average Intensity
  # Omic levels preprocessed: transcripts, proteins, metabolites

#### Feature selection. The employ of the whole dataset led to unclear analyses (code not shown).
#### This step was performed to select the more variable variables, excluding the noisiest, 
#### over previously generated "original.dataset_preprocessed" dataset.
featureselection_wizard()
  # Method: Variation Coefficient
  # Omic datasets: transcripts, proteins, metabolites
  # Threshold: 85 (0.85)
  # Keep variables: Below
featureselection_wizard()
  # Dataset: "original.dataset_preprocessed_VarCoef_featsel"
  # Method: Stats
  # Omic datasets: transcripts, proteins, metabolites
  # Threshold: 0.3
  # Employ parametric analysis: FALSE
  # Statistical to be used: p-value
  # Keep variables: Below
original.dataset_featured <- original.dataset_preprocessed_VarCoef_featsel_Stat_featsel
rm(original.dataset_preprocessed_VarCoef_featsel_Stat_featsel,original.dataset_preprocessed_VarCoef_featsel)

#### Univariate analyses ####

#### Calculation of Average, SD, ANOVA, FDR, Post-Hoc
univariate_wizard()
  # Dataset: "original.dataset_featured"
  # Paramecitry test: Yes
  # Transformation: Box-Cox-Lambda
  # Test to apply: Parametric (ANOVA)
  # Perform post-hoc analysis: Yes
  # FDR calculation: Yes
  # Save plots and XLS: No

export_table(original.dataset_featured_univariate) #Export supplementary tables S2-S5. Descriptions were added later in Microsoft Excel.

#### Multivariate: K-means and PCA ----
####PCA transcriptome dataset (Figure 1A-B; Table S6)
pca_analysis_wizard()
  # Dataset: "original.dataset_featured"
  # Omic levels: transcripts
  # Add annotations: No
  # Save PCA Results: Yes #Generation of Table S6
  # Plots: scree, score, loadings
pca_scoreplot <- pca_plot(original.dataset_featured_pca,1,1,2,plottype = "Score") #Generation of Figure 1A
export_plot(pca_scoreplot, filename = "pca_scoreplot.pdf")

#Plot some supplementary figures
pca_topscoring_1 <- pca_plot(original.dataset_featured_pca,treatment = 1,compX = 1,plottype = "Topscoring",fortopscoring =c(1,25,"abs"))
export_plot(pca_topscoring_1, filename = "pca_topscoring_1.pdf")
pca_topscoring_2 <- pca_plot(original.dataset_featured_pca,treatment = 1,compX = 1,plottype = "Topscoring",fortopscoring =c(2,25,"abs"))
export_plot(pca_topscoring_2, filename = "pca_topscoring_2.pdf")

#k-means transcriptome dataset (Figure 1C; Table S7)
kmeans_wizard()
  # Dataset: "original.dataset_featured"
  # Treatment for splitting data: 1
  # Scaling method: scaling and not centering
  # Omic levels: transcripts
  # Max clusters: 18
  # Min clusters: 15

export_plot(original.dataset_featured_kmeans, filename = "kmeans.pdf") #Export plot of Figure 1B

export_table(original.dataset_featured_kmeans) #Export supplementary table S7. Descriptions were added later in Microsoft Excel.

#### Determination of the functional enrichment of the different clusters
library(gprofiler2) #Required for the analysis
library(ggplot2) #Required for plotting
library(clusterProfiler)
library(enrichplot)
library(DOSE)

PinusRadiataToken <- "gp__GBCz_Eh2n_bII" #Token of the employed GMT file
clustno <- 16 #We select 16 clusters
aux.clustno <- paste(clustno, "Clusters", sep=" ")
aux.clust.index <- which(names(original.dataset_featured_kmeans$kmeans_list) == aux.clustno)
clustering <- original.dataset_featured_kmeans$kmeans_list[[aux.clust.index]]
b <- ncol(clustering)
kme <- clustering[,c(b-1,b)] #kme are the two last columns, which correspond to groups and ID

#### Functional enrichment analysis by cluster
usthre <- 0.01
cluster_funct_enrichment <- list()
for (i in 1:max(kme[,1])){ 
  cluster <- kme[kme[,1]==i,2]
  cluster_funct_enrichment[[i]] <- gprofiler2::gost(query = cluster, 
                                    organism = PinusRadiataToken, 
                                    user_threshold = usthre, 
                                    significant = T)
}
names(cluster_funct_enrichment) <- paste("Cluster", 1:i)

cluster_funct_enrichment[sapply(cluster_funct_enrichment, is.null)] <- NULL

#Initial data parsing
cluster_funct_enrichment_for_plot <- lapply(cluster_funct_enrichment, function(x) x$result[,c("query", "source", "term_id",
                                                              "term_name", "p_value", "query_size", 
                                                              "intersection_size", "term_size", 
                                                              "effective_domain_size")])
#Adding GeneRatio and BgRatio
for (i in 1:length(cluster_funct_enrichment_for_plot)){
  cluster_funct_enrichment_for_plot[[i]] <- cbind(cluster_funct_enrichment_for_plot[[i]], 
                           "GeneRatio"= paste0(cluster_funct_enrichment_for_plot[[i]]$intersection_size,  "/", cluster_funct_enrichment_for_plot[[i]]$query_size),
                           "BgRatio" = paste0(cluster_funct_enrichment_for_plot[[i]]$term_size, "/", cluster_funct_enrichment_for_plot[[i]]$effective_domain_size))
  #Adding new colnames
  names(cluster_funct_enrichment_for_plot[[i]]) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                             "query_size", "Count", "term_size", "effective_domain_size", 
                             "GeneRatio", "BgRatio")
  #rownames
  row.names(cluster_funct_enrichment_for_plot[[i]]) = cluster_funct_enrichment_for_plot[[i]]$ID
  #Clusters
  cluster_funct_enrichment_for_plot[[i]]$Cluster <- do.call(c, strsplit(names(cluster_funct_enrichment_for_plot)[[i]], split = " "))[2]
}

cluster_funct_enrichment_for_plot <- do.call(rbind, cluster_funct_enrichment_for_plot) #Return from list to df


#### Multi-Plot

### All clusters in one plot: personal favourite option
# define as compareClusterResult object
###ERROR: el codigo siguiente no funciona si no hay NAs y debe ser omitido:
#cluster_funct_enrichment_for_plot_no_NA <- cluster_funct_enrichment_for_plot[-which(cluster_funct_enrichment_for_plot$Description == "NA"),] 
cluster_funct_enrichment_for_plot_no_NA <- cluster_funct_enrichment_for_plot
cluster_funct_enrichment_plot_data= new("compareClusterResult", compareClusterResult = cluster_funct_enrichment_for_plot_no_NA)

View(cluster_funct_enrichment_plot_data)
#multi-cluster plot
pdf("FEkme16fixeds5.pdf",15,23)
enrichplot::dotplot(cluster_funct_enrichment_plot_data, showCategory = 5, includeAll = F, font.size = 11) + 
  ggplot2::scale_color_gradientn(colours = grDevices::colorRampPalette(rev(c("darkgreen", "yellow","darkred"))) (20))
dev.off()

####Multivariate BlockPLS-DA####
#Models were previously set up (code not shown)
list.keepX <- list(transcripts=c(60,60),proteins=c(50,50),metabolites=c(30,30),physiology=c(5,5))
original.dataset_featured_bda <- bda_analysis(original.dataset_featured,initialrow = 1,initialcolumn = 2,treatment1col = 1, keepX = list.keepX,omiclevel="all")
bda_distance <- bda_plot(original.dataset_featured_bda_keepX2,plottype = "Distance") #Figure 4A
bda_cim <- bda_plot(original.dataset_featured_bda_keepX2,plottype = "Cim") #Figure 4B
export_plot(bda_distance,filename = "bda_dist.pdf")
export_plot(bda_cim, filename = "bda_cim.pdf")
export_table(original.dataset_featured_bda) #Table S8
original.dataset_featured_bda$bda$prop_expl_var #Check of explained variance
network_heatstress<-mixOmics::network(original.dataset_featured_bda$bdanetwork, blocks = c(1,2,3,4))
export_table(original.dataset_featured_bda,filename = "network.xlsx") #Export network in xlsx
igraph::write.graph(network_heatstress$gR,file="network.gml",format="gml") #Export network in gml


####Volcano plots transcriptomics dataset ----
library(plotly)
data_volcanos <- list("C-T1" = original.dataset_featured$transcripts[c(1:3,4:6), 2:ncol(original.dataset_featured$transcripts)],
                      "C-T3" = original.dataset_featured$transcripts[c(1:3,7:9), 2:ncol(original.dataset_featured$transcripts)],
                      "T1-T3" = original.dataset_featured$transcripts[c(4:6,7:9), 2:ncol(original.dataset_featured$transcripts)])

#Removal of empty columns and qualitative variables
data_volcanos_clean <-lapply(data_volcanos, function(x) t(RemoveEmptyColumns(x,1,1)))
data_volcanos_clean <- lapply(data_volcanos_clean, function(x) RemoveQualVars(x,3)$"QuantitativeVars")

#Fold change and statistics
G1 <- lapply(data_volcanos_clean , function(x) apply(x[,1:3], 1, mean))
G2 <- lapply(data_volcanos_clean , function(x) apply(x[,4:6], 1, mean))
Foldchange <- mapply("/", G2, G1)
Log2foldChange <- lapply(Foldchange, log2)
pvalues <- lapply(data_volcanos_clean, function(x) apply(x, 1, ttest, grp1 = c(1:3), grp2 = c(4:6)))
data_volcanos_stats <- mapply(cbind, Log2foldChange, pvalues)
data_volcanos_stats <- lapply(data_volcanos_stats, as.data.frame) #required to maintain numeric vectors in the next steps
data_volcanos_stats <- mapply(cbind, data_volcanos_stats, lapply(data_volcanos_stats, updownregulated2) , SIMPLIFY = FALSE) #Up or down regulated
data_volcanos_stats <- lapply(data_volcanos_stats, function(x) cbind(rownames(x),x)) #add accesions to tables
data_volcanos_stats <- lapply(data_volcanos_stats, "colnames<-", c("accession","log2FoldChange","pvalue","diffexpressed"))
rm(data_volcanos,data_volcanos_clean,G1,G2,Foldchange,Log2foldChange,pvalues) #time to clean the environment

#Export dataset
writexl::write_xlsx(data_volcanos_stats,"data_volcanos_stats.xlsx")

####Plotting volcanos and exporting tables and figures
#Palette
pal.new <-list("C-T1" = c("#BF2549","black","#fae650") ,
               "C-T3" = c("#BF2549","black","#0082c8"),
               "T1-T3" = c("#fae650","black","#0082c8"))

#Plotly 
VolcanoListPlots <- mapply(function(results, n, pal.new){
  afig <- plot_ly(type = 'scatter', mode = 'markers', data = results, colors = pal.new, text =~results[,1]) 
  afig <- afig %>% add_markers(x = ~log2FoldChange, y = ~-log10(pvalue), color = ~diffexpressed, opacity = 0.7, showlegend = T)
  afig <- afig %>% plotly::layout(annotations = list(x = 0.05, y = 1, text = paste(n,"contrast", sep=" "), showarrow = F, xref='paper', yref='paper'))
},
data_volcanos_stats, names(data_volcanos_stats), pal.new, SIMPLIFY = FALSE)

#Export plots
lapply(seq_along(VolcanoListPlots), function(i) orca(VolcanoListPlots[[i]],paste(names(VolcanoListPlots)[[i]],".pdf",sep="")))

#Get unified plot
Volcano_plot_unified <- subplot(style(VolcanoListPlots[[1]],showlegend=T),
               style(VolcanoListPlots[[2]],showlegend=T),
               style(VolcanoListPlots[[3]],showlegend=T),
               nrows=1, shareY = T, shareX=T, 
               margin = 0.05)
Volcano_plot_unified <- unify_legend(Volcano_plot_unified)
orca(Volcano_plot_unified, "Volcano_plot_unified.pdf")


####Auxiliary functions----
#Functions required for cleaning the data
RemoveEmptyColumns <- function(matriz, initialrow, initialcolumn){
  matrizfiltrada <-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  initial<-ncol(matrizfiltrada)
  vectorseleccion <-which(colSums(matrizfiltrada,na.rm=TRUE)==0)
  end<-length(vectorseleccion)
  cat(paste("\nOut of ",initial,"initial variables, ",end, "have been found to be empty and dropped from the analysis"))
  if(length(vectorseleccion)==0){
    return(matriz)}
  else {matrizfiltrada <- matrizfiltrada[,-vectorseleccion]}
  if(initialrow==1){
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    return(matrizfiltrada)
  } else if(initialrow>=2) {
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    matrizfiltrada <- rbind(matriz[1:initialrow-1,],matrizfiltrada)
    
    return(matrizfiltrada)
  } else {
    return(cat("\nError. Please revise your data input."))
  }
} #Removes empty columns
RemoveQualVars <- function(matriz,n){
  Amat <- apply(matriz[,1:n], 1, sum)
  Bmat <- apply(matriz[,(n+1):ncol(matriz)], 1, sum)
  if(any(Amat == 0)|(any(Bmat == 0))){
    QLindex <- c(which(Amat == 0), which(Bmat == 0))
    
    cat(paste("\n\nThere are", length(QLindex), 
              "qualitative variables in your dataset. \nThese will be discarted for Volcano plots.","\nFirst treament has", 
              length(which(Amat == 0)), 
              "and second treatment has", 
              length(which(Bmat == 0)),".", sep=" "))
  } else stop
  
  QuanVars <- matriz[-QLindex,]
  QualVars <- cbind(names(QLindex), 
                    c(rep(c("G1"), length(which(Amat == 0))), 
                      rep(c("G2"), length(which(Bmat == 0)))))
  colnames(QualVars) <- c("IDENTIFIER","Treatment")
  finalObj <- list("QualitativeVars"=QualVars, "QuantitativeVars"=QuanVars)
  return(finalObj)
} #Removes qualitative variables
#Statistical functions
ttest <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  results = t.test(x, y)
  results$p.value
}
updownregulated <- function(x){
  #x[1]=log2foldchange, x[2]=pvalue
  if(x[1] > 2 & x[2] < 0.05) {y <- "UP"}
  else if(x[1] < -2 & x[2] < 0.05) {y <- "DOWN"}
  else {y <- "NO"}
  return(y)
}
updownregulated2 <- function(x){
  y <- apply(x, 1, updownregulated)
  return(as.data.frame(y))
}
#Plotly
unify_legend <- function(x) {
  x <- plotly::plotly_build(x)
  lgroup <- lapply(x$x$data, "[[", "legendgroup")
  #change NULL to "" or " " or "NA 
  lgroup <- lapply(lgroup, function(x){
    if(is.null(x)) x <- " "
    else x <- x # I don't understand why I need this line
  })
  pltype <- sapply(x$x$data, "[[", "type")
  showit <- sapply(x$x$data, "[[", "showlegend")
  showit <- !duplicated(unlist(lgroup))
  
  for (i in seq(x$x$data)) { # 
    x$x$data[[i]]$showlegend <- showit[i] 
  }
  x
} #Unify plotly legends
