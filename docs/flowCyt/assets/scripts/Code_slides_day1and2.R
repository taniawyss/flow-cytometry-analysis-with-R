#######################################
# This code comes from the slides     #
#######################################


## Load packages


# clear environment
rm(list = ls())

set.seed(123)

# Basic structures for flow cytometry data
library(flowCore)
# Provides S4 data structures and basic functions to deal with flow cytometry data.


# Visualization tools for flow cytometry data.
library(ggcyto)
library(flowViz)

# Per-channel variance stabilization from a collection of flow cytometry samples by Bertlett test for homogeneity of variances. The approach is applicable to microarrays data as well
library(flowVS)

# Automatic and interactive quality control for flow cytometry data
library(flowAI)
# The package is able to perform an automatic or interactive quality control on FCS data acquired using flow cytometry instruments. By evaluating three different properties: 1) flow rate, 2) signal acquisition, 3) dynamic range, the quality control enables the detection and removal of anomalies.

# Peak-based selection of high quality cytometry data
library(PeacoQC)
# This is a package that includes pre-processing and quality control functions that can remove margin events, compensate and transform the data and that will use PeacoQCSignalStability for quality control. This last function will first detect peaks in each channel of the flowframe. It will remove anomalies based on the IsolationTree function and the MAD outlier detection method. This package can be used for both flow- and mass cytometry data.

# Automated quality control based on changepoint analysis
library(flowClean)

# Cytometry dATa anALYSis Tools
library(CATALYST)
# Mass cytometry (CyTOF) uses heavy metal isotopes rather than fluorescent tags as reporters to label antibodies, thereby substantially decreasing spectral overlap and allowing for examination of over 50 parameters at the single cell level. While spectral overlap is significantly less pronounced in CyTOF than flow cytometry, spillover due to detection sensitivity, isotopic impurities, and oxide formation can impede data interpretability. We designed CATALYST (Cytometry dATa anALYSis Tools) to provide a pipeline for preprocessing of cytometry data, including i) normalization using bead standards, ii) single-cell deconvolution, and iii) bead-based compensation.

# S4 Classes for Single Cell Data
library(SingleCellExperiment)

# Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction
library(uwot)

# Differential testing
library(diffcyt)



#############################################
# Starting to work with flow cytometry data #
#############################################

?read.FCS
## Read a FCS file into a flowFrame
FCS_file <- flowCore::read.FCS(filename = "course_datasets/FR_FCM_Z4KT/T_cells_REU270_alive_T cells.fcs",
                     transformation = FALSE,
                     truncate_max_range = FALSE)

## summarize a flowFrame
FCS_file
summary(FCS_file)


## Access data elements in a flowFrame

# Matrix of expression
FCS_file@exprs
exprs(FCS_file) # alternative

# column names in the expression matrix
colnames(FCS_file)


## Access data elements in a flowFrame

# metadata
pData(FCS_file@parameters)
pData(parameters(FCS_file))

## Replace channel names by antigen names

# copy the metadata to a data frame
panel <- pData(FCS_file@parameters)

# copy the names to a new column
pData(FCS_file@parameters)$channel <- panel$name

# replace the names by antigens
colnames(FCS_file)[!(is.na(panel$desc))] <- panel$desc[!is.na(panel$desc)]

# check
head(exprs(FCS_file))[,10:15]

## Read a list of FCS files into a flowSet


fcs_data <- flowCore::read.flowSet(path="course_datasets/FR_FCM_Z4KT/", 
                         pattern="*.fcs", 
                         transformation = FALSE, 
                         truncate_max_range = FALSE) 
# Access 1 of the flowFrame
fcs_data@frames$`T_cells_REU267_alive_T cells.fcs`

## Methods applied to a flowSet
# list sample names
sampleNames(fcs_data)

# change sample names
sampleNames(fcs_data) <- c("REU267","REU268","REU269","REU270",
                           "REU271_12_july","REU271_13_april",
                           "REU271_14_april","REU271_7_apr",
                           "REU271_9_april","REU271","REU272_12_july",
                           "REU272_13_april","REU272_14_april",
                           "REU272_7_apr","REU272_9_apri","REU272")

## Access phenotypic data of samples
pData(fcs_data)
class(pData(fcs_data)) # "data.frame"


## Add a new column to the phenotypic data
pData(fcs_data)$gender <- c(rep("male",8), rep("female",8))

gender<-c("male", "female", "female")

# check
pData(fcs_data)
fcs_data@phenoData@data

## Manipulating a flowSet

# extract a flowFrame from a flowSet using the [[ operator
fcs_data[[1]]

# create a new flowSet object by subsetting with the [ operator
fcs_data[1:5]

# subset a flowSet based on a condition
fcs_data_males <- fcs_data[pData(fcs_data)$gender=="male"]
fcs_data_females <- subset(fcs_data,pData(fcs_data)$gender=="female") # alternative

# split a flowSet baed on a condition
fcs_data_split <- flowCore::split(fcs_data, pData(fcs_data)$gender)

# check
names(fcs_data_split)

# combine several flowSet objects (or flowSets and flowFrames)
fcs_data_combined <- rbind2(fcs_data_split$female, fcs_data_split$male)

# check
pData(fcs_data_combined)

## Vizualizing a single flowFrame within a flowSet
?autoplot
# scatterplot
ggcyto::autoplot(object = fcs_data[[5]], x="FSC-H", y="FJComp-BUV496-A", bins = 2^7)

# univariate density plot
ggcyto::autoplot(object = fcs_data[[5]], x="FSC-H")


## In-line transformation
# univariate density plot
ggcyto::autoplot(fcs_data[[5]], x="FJComp-BUV496-A") # original scale
ggcyto::autoplot(fcs_data[[5]], x="FJComp-BUV496-A") + scale_x_flowjo_fasinh() # flowJo inverse hyperbolic sine
ggcyto::autoplot(fcs_data[[5]], x="FJComp-BUV496-A") + scale_x_flowjo_fasinh(m=1/10) # example with a cofactor of 100

# Customize using ggplot2 functions:
ggcyto::autoplot(fcs_data[[5]], x="FJComp-BUV496-A") + ggtitle("original scale")

# bivariate density plot
ggcyto::autoplot(object=fcs_data[[5]], x="FSC-H", y="FJComp-BUV496-A", bins = 2^7) + 
  scale_y_flowjo_fasinh()
ggcyto::autoplot(object=fcs_data[[5]], x="FSC-H", y="FJComp-BUV496-A", bins = 2^7) + 
  scale_y_flowjo_fasinh(m=1/10) # example with a cofactor of 100


## Visualizing a flowSet

ggcyto::autoplot(object = fcs_data[10:15],
         x="FSC-H",
         y="SSC-H",
         bins=2^7) 



#####################
# Transformation    #
#####################


## Constructing a data frame for the panel


# retrieve the list of channels and corresponding antigens
fcs_colname <- colnames(fcs_data) 
antigen <- pData(parameters(fcs_data[[1]]))$desc

# setting marker classes
marker_class <- rep("none", ncol(fcs_data[[1]])) # setting all to "none"
marker_class[c(8:31,33:36,38)] <- "type" # type markers
marker_class[32] <- "state" # state markers
marker_class <- factor(marker_class,levels=c("type","state","none"))

# put everything together in a data frame
panel <- data.frame(fcs_colname, antigen, marker_class, row.names = NULL) 


## Downsample the data


# define function
Downsampling_FlowSet <- function(x, samplesize , replace=TRUE, prob=NULL){
  if(missing(samplesize))
    samplesize <- min(flowCore::fsApply(x,nrow))
  flowCore::fsApply(x, function(ff){
    i <- sample(nrow(ff), size = samplesize, replace=replace, prob)
    ff[i,]
  })
}

# Create a downsampled flowSet
fcs_data_small <- Downsampling_FlowSet(x=fcs_data, 
                                       samplesize = 2000)

## Arcsimh transformation with flowVS

# select markers to be transformed
markerstotransf <- panel$fcs_colname[panel$marker_class!="none"]

# Estimate cofactors based on the downsampled data
# This takes very long to run
# Please load the resulting vector instead 
if(FALSE) {
cofactors <- estParamFlowVS(fcs_data_small, 
                             channels=markerstotransf)
 
 save(cofactors, file="course_datasets/FR_FCM_Z4KT/cofactors.RData") } else {
load("course_datasets/FR_FCM_Z4KT/cofactors.RData") }


# check the cofactors
cofactordata <- data.frame(markerstotransf, cofactors)
cofactordata


## Arcsinh transformation with with flowVS and the estimated cofactors

# transform
fcs_transform <- transFlowVS(fcs_data, 
                             channels = markerstotransf,cofactordata$cofactors)


## Check transformation with FlowViz: 


# Density plots
densityplot(~`FJComp-BUV496-A`, fcs_data) #density plot before transformation, you can replace `FJComp-BUV496-A` by . to view all markers.

densityplot(~`FJComp-BUV496-A`, fcs_transform) # density plot after transformation


##  Alternative: Arcsinh transformation with fixed cofactors

# set cofactor vector
cofactor <- 3000
l <- length(markerstotransf)             
cofactors <- rep(cofactor, l)

# transform
fcs_transform <- transFlowVS(fcs_data,
                             channels = markerstotransf,
                             cofactors)

# save to a file for downstream analysis
save(fcs_transform,file = "course_datasets/FR_FCM_Z4KT/fcs_transform.RData")


# Check transformation with density plots
densityplot(~`FJComp-BUV496-A`, fcs_data) # before transformation
densityplot(~`FJComp-BUV496-A`, fcs_transform) # density plot after transformation


##  Logicle transformation with flowCore


# estimate parameters and transform

# create an empty list for the transformed flowFrames
fcs_list <- list() 

# iterate over flowFrames
for(i in 1:16){
  
  # select the flowFrame
  ff <- fcs_data[[i]] 
  
  # estimate the parameters
  algcl <- estimateLogicle(ff, 
                           channels = markerstotransf, 
                           m=6,
                           t = 4E6) 
  
  # transform and add to list of flowFrames
  fcs_list[[i]] <- transform(ff, algcl) 
  
  
}



# recreate the transformed flowSet
names(fcs_list) <- sampleNames(fcs_data)
fcs_transform <- as(fcs_list, "flowSet")


# save
save(fcs_transform, file = "course_datasets/FR_FCM_Z4KT/fcs_transform_Logicle.RData")

# Density plots
flowViz::densityplot(~`FJComp-BUV496-A`, fcs_data) # density plot before transformation
flowViz::densityplot(~`FJComp-BUV496-A`, fcs_transform) # density plot after transformation


#########################################
# Automated Quality Control             #
#########################################

# Load the arcsinh transformed flowSet
load("course_datasets/FR_FCM_Z4KT/fcs_transform.RData")

# run flowAI
# This step takes very long to run
# Please load the resulting flowSet instead 
if(FALSE) {
  fcs_QC <- flow_auto_qc(fcs_transform, folder_results = "course_datasets/FR_FCM_Z4KT/flowAI_results/") 
 save(fcs_QC,file = "course_datasets/FR_FCM_Z4KT/fcs_QC.RData") } else {
load("course_datasets/FR_FCM_Z4KT/fcs_QC.RData") }


# If pre-gated, skip FR step
# Do not run
# fcs_QC <- flow_auto_qc(fcs_transform, remove_from = "FS_FM")  

## Automated Quality control with peacoQC

# run peacoQC first on one file to optimize parameters
peacoqc_res <- PeacoQC(ff= fcs_transform[[1]], 
                       channels=markerstotransf,
                       determine_good_cells = "all", 
                       save_fcs = FALSE, 
                       plot=TRUE, 
                       output_directory = "course_datasets/FR_FCM_Z4KT/PeacoQCresults",
                       IT_limit = 0.65, 
                       MAD=8)

# after choosing the right parameters, apply to all samples
for(i in 1:16){
  peacoqc_res <- PeacoQC(fcs_transform[[i]], 
                         markerstotransf, 
                         determine_good_cells = "all",
                         IT_limit=0.55, 
                         MAD=5, 
                         save_fcs = TRUE, 
                         plot=TRUE, 
                         output_directory = "course_datasets/FR_FCM_Z4KT/PeacoQCresults")
} 

# construct new flowSet from the cleaned fcs files
fcs_QC <- read.flowSet(path= "course_datasets/FR_FCM_Z4KT/PeacoQCresults/PeacoQC_results/", 
                       transformation=FALSE,
                       truncate_max_range = FALSE)


## Automated Quality control with flowClean
# This takes very long time
# Please don't run
# # create an empty list for the transformed flowFrames
# fcs_list <- list() 
# 
# # iterate over flowFrames
# for(i in 1:16){
#   
#   # select the flowFrame
#   ff <- fcs_transform[[i]] 
#   
#   # transform and add to list of flowFrames
#   fcs_list[[i]]  <- clean(fF = fcs_transform[[i]] ,
#                           vectMarkers = match(markerstotransf,colnames(fcs_transform))) 
#
#   
# }
# 
# # recreate the flowSet with an added parameter "GoodVsBad"
# names(fcs_list) <- sampleNames(fcs_transform)
# fcs_QC <- as(fcs_list, "flowSet")


###############################
# Dimensionality Reduction    #
###############################


## Example dataset from the CATALYST package

data(PBMC_fs, PBMC_panel, PBMC_md)

# check
PBMC_fs
View(PBMC_md)
View(PBMC_panel)


## Create the sce object

sce_PBMC <- prepData(PBMC_fs,
                md=PBMC_md, 
                panel= PBMC_panel,
                transform = TRUE,
                FACS = FALSE,
                features = NULL)


# check name of the matrix with transformed values
assayNames(sce_PBMC)


## UMAP with the uwot R package

# Extract the expression matrix and transpose
exprs_PBMC <- assay(sce_PBMC, "exprs")
exprs_PBMC <- t(exprs_PBMC)

# subset to markers you want to use for clustering
marker_type <- PBMC_panel$antigen[PBMC_panel$marker_class=="type"]
exprs_PBMC <- exprs_PBMC[,c(marker_type)]

# compute the UMAP
set.seed(1234)
umap_PBMC <- umap(exprs_PBMC)

# add UMAP coordinates to sce object
reducedDim(sce_PBMC, "UMAP") <- umap_PBMC

# Plot
plotDR(sce_PBMC, 
       color_by="sample_id")

plotDR(sce_PBMC,
       dr = "UMAP",
       assay = "exprs",
       color_by = "CD3",
       facet_by =  "sample_id")

# try different values for parameters

set.seed(1234)
umap_PBMC <- umap(exprs_PBMC, 
                  n_neighbors=5)
reducedDim(sce_PBMC, "UMAP") <- umap_PBMC
plotDR(sce_PBMC, 
       color_by="sample_id")


set.seed(1234)
umap_PBMC <- umap(exprs_PBMC, 
                  min_dist = 0.5)
reducedDim(sce_PBMC, "UMAP") <- umap_PBMC
plotDR(sce_PBMC, 
       color_by="sample_id")




# plot expression of given markers


# redo UMAP with default parameter values
set.seed(1234)
umap_PBMC <- umap(exprs_PBMC)
reducedDim(sce_PBMC, "UMAP") <- umap_PBMC


# CD45
plotDR(sce_PBMC, 
       "UMAP",
       color_by="CD45")

# CD3
plotDR(sce_PBMC, 
       "UMAP",
       color_by="CD3")


## UMAP with the CATALYST package

# run dimensionality reduction
set.seed(1234)
sce_PBMC <- runDR(sce_PBMC,
             assay="exprs",
             dr="UMAP",
             cells = NULL,
             features="type",
             n_neighbors = 10,
             min_dist = 0.05)

# plot UMAP with expression of CD3 by sample
plotDR(sce_PBMC,
       dr="UMAP",
       assay="exprs",
       color_by = "CD3",
       facet_by = "sample_id")


#############################
# Clustering and Annotation #
#############################


# Umsupervised clustering with CATALYST
sce_PBMC <- cluster(sce_PBMC, 
               features="type", 
               xdim = 10,
               ydim = 10,
               maxK=20, 
               seed=5024)

# names of clustering schemes
names(cluster_codes(sce_PBMC))


# plot UMAP with SOM clusters
plotDR(sce_PBMC, "UMAP", 
       color_by="som100")

# plot UMAP with 20 metaclusters
plotDR(sce_PBMC, "UMAP", 
       color_by="meta20")


# Heatmap of the median expression per marker and metacluster
plotExprHeatmap(sce_PBMC, 
                features = "type",
                by = "cluster_id", 
                k = "meta20", 
                scale = "first", 
                q = 0.01, 
                perc = TRUE, 
                row_clust = FALSE, 
                col_dend = TRUE)

# Ridge plots of expression per marker and metacluster
plotClusterExprs(sce_PBMC, 
                 k = "meta20", 
                 features = "type")


## Manual cluster renaming / merging

# Create a 2 column data.frame containing old_cluster and new_cluster
# data(merging_table) # an example of this table is also available with the CATALYST package
merging_table <- data.frame(old_cluster = 1:20, 
                            new_cluster = c("B-cells IgM+","surface-","NK cells",
                                            "CD8 T-cells","B-cells IgM-","monocytes",   
                                            "monocytes","CD8 T-cells","CD8 T-cells",
                                            "monocytes","monocytes","CD4 T-cells",
                                            "DC","CD8 T-cells","CD4 T-cells","DC",
                                            "CD4 T-cells","CD4 T-cells","CD4 T-cells",
                                            "CD4 T-cells"))


# merge / rename clusters
sce_PBMC <- mergeClusters(sce_PBMC, 
                          k = "meta20", 
                          table = merging_table, 
                          id = "final_annotation")

# plot UMAP with final annotation
plotDR(sce_PBMC, "UMAP", 
       color_by="final_annotation")


# heatmap of the median expression per markers and metacluster
plotExprHeatmap(sce_PBMC, 
                features = "type",
                by = "cluster_id", 
                k = "final_annotation", 
                scale = "first", 
                q = 0.01, 
                perc = TRUE, 
                row_clust = FALSE,
                col_dend = TRUE)

#save(sce_PBMC,file =  "course_datasets/FR_FCM_Z4KT/DA_example_sce_PBMC.RData")


###########################################
# Differential testing with diffcyt       #
###########################################

# example dataset
load("course_datasets/FR_FCM_Z4KT/DA_example_sce_PBMC.RData")

plotDR(sce_PBMC, color_by="final_annotation", facet_by="condition")

# check experimental information
ei(sce_PBMC)

# set design matrix
design <- createDesignMatrix(ei(sce_PBMC), cols_design = "condition")
design

# Set the contrast matrix
contrast <- createContrast(c(0,1))
contrast


## Differential abundance (DA) analysis

# Plot relative population abundances (stacked bar plot)
plotAbundances(sce_PBMC,
               k="final_annotation",
               by="sample_id",
               group_by="condition")

# Plot relative population abundances (box plot)
plotAbundances(sce_PBMC,
               k="final_annotation",
               by="cluster_id",
               group_by="condition",
               shape_by="patient_id")

# test for differences in abundance between conditions
res_DA <- diffcyt(sce_PBMC,
                  clustering_to_use = "final_annotation",
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  design = design,
                  contrast = contrast)

# extract results
tbl_DA <- rowData(res_DA$res)
tbl_DA


# plot results
plotDiffHeatmap(sce_PBMC,
                tbl_DA,
                fdr = 0.05,
                lfc = 1,
                top_n = 20,
                all=TRUE,
                normalize = TRUE,
                col_anno = "condition")


## Differential state (DS) analysis

# test
res_DS <- diffcyt(sce_PBMC,
                  clustering_to_use = "final_annotation",
                  analysis_type = "DS",
                  method_DS = "diffcyt-DS-limma",
                  design = design,
                  contrast = contrast)

# extract results table
tbl_DS <- rowData(res_DS$res)
tbl_DS


# Plot results
plotDiffHeatmap(sce_PBMC, 
                tbl_DS,
                fdr = 0.05,
                sort_by = "lfc",
                col_anno = "condition")

