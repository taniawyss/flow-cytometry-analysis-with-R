#############################
# Let’s practice – 4        #
# Dimensionality reduction  #
# Solutions                 #
#############################


# clear the environment
rm(list = ls())
gc()

# load libraries
library(flowCore)
library(CATALYST)


# 1) load the flowSet object from previous exercise of flowAI output
load("course_datasets/FR_FCM_Z3WR/fcs_clean.RData")
fcs_clean

# 2) Downsample to 2'000 cells 

# source downsampling function
source("course_datasets/function_for_downsampling_flowSets.R")

# downsample
set.seed(1234)
fcs_small <- Downsampling_flowSet(fcs_clean,samplesize = 2000)

# 3) Create a sce object from the downsampled flowSet
# We need the panel data frame
panel <- read.csv("course_datasets/FR_FCM_Z3WR/panel_with_marker_classes.csv")

# We also need a metadata dataframe
md <- pData(fcs_small)

# create sce
sce <- CATALYST::prepData(fcs_small, 
                          md = md,
                          md_cols = list(file="name", id = "name", factors = "time_point"),
                          panel = panel,
                          panel_cols = list(channel = "channels", antigen = "antigen", class = "marker_class"),
                          transform = FALSE,
                          FACS=TRUE,
                          features = panel$channels[panel$marker_class!="none"])
# Overview of the object:
sce
# change the assay name to "exprs" because it was already transformed
assayNames(sce) <- "exprs"
head(assays(sce)$exprs[,1:10])

# Phenotypic data
colData(sce)
# Channel parameters
rowData(sce)

# save
save(sce,file = "course_datasets/FR_FCM_Z3WR/sce.RData")

# 4) UMAP with default parameters (n_neighbors=15 and min_dist = 0.01)
set.seed(1601)
?runDR
sce <- runDR(sce, 
             assay = "exprs", 
             dr = "UMAP", 
             pca = 10,
             features = "type", 
             cells = NULL)

# plot
plotDR(sce,
       dr =  "UMAP",
       assay = "exprs", 
       color_by="CD3", 
       facet_by = "time_point" ) + ggtitle("CD3, n=15, dist=0.01")

reducedDims(sce)
# List of length 1
# names(1): UMAP 

# save
save(sce,file = "course_datasets/FR_FCM_Z3WR/sce_UMAP.RData" )


# 4) UMAP with n_neighbors=15 (defaults) and min_dist = 0.5
set.seed(1601)
sce <- runDR(sce, 
             assay = "exprs", 
             dr = "UMAP", 
             features = "type", 
             cells = NULL,
             min_dist = 0.5)
# the previous UMAP coordinates are overwritten, we still have only 1 
# reducedDims element:
reducedDims(sce)

# plot
plotDR(sce,
       dr =  "UMAP",
       assay = "exprs", 
       color_by="CD3", 
       facet_by = "time_point" ) + ggtitle("CD3, n=15, dist=0.5")


# 4) UMAP with n_neighbors=5 and min_dist = 0.01 (default)
set.seed(1601)
sce <- runDR(sce, 
             assay = "exprs", 
             dr = "UMAP", 
             features = "type", 
             cells = NULL,
             n_neighbors = 5)

# plot
plotDR(sce,
       dr =  "UMAP",
       assay = "exprs", 
       color_by="CD3", 
       facet_by = "time_point" ) + ggtitle("CD3, n=5, dist=0.01")


# 5) UMAP with n_neighbors=2 and min_dist = 0.5
set.seed(1601)
sce <- runDR(sce, 
             assay = "exprs", 
             dr = "UMAP", 
             features = "type", 
             cells = NULL,
             n_neighbors = 2,
             min_dist = 0.5)


# plot
plotDR(sce,
       dr =  "UMAP",
       assay = "exprs", 
       color_by="CD3", 
       facet_by = "time_point" ) + ggtitle("CD3, n=2, dist=0.5")

# Try tSNE! It is slower than UMAP
if(TRUE) {  # change this to FALSE after running it so you can just quickly
  # load the object with TSNE afterwards and save time
set.seed(1601)
sce <- runDR(sce, 
             assay = "exprs", 
             dr = "TSNE", 
             features = "type", 
             cells = NULL)
reducedDims(sce)
# List of length 2
# names(2): UMAP TSNE
save(sce,file = "course_datasets/FR_FCM_Z3WR/sce_TSNE.RData" )
} else {
  load("course_datasets/FR_FCM_Z3WR/sce_TSNE.RData")
  reducedDims(sce)
  # List of length 2
  # names(2): UMAP TSNE
}
# plot
plotDR(sce,
       dr =  "TSNE",
       color_by="sample_id",
       facet_by = "time_point") + ggtitle("TSNE")

# Try PCA, using all features and a max number of components
table(rowData(sce)$marker_class)

sce<-runDR(sce, 
           dr="PCA",
           features = NULL,
           ncomponents = 15)
reducedDims(sce)
# List of length 3
# names(3): UMAP TSNE PCA


plotDR(sce, dr="PCA", color_by = "time_point")
plotDR(sce, dr="PCA", color_by = "CD3", facet_by = "time_point")

# ElbowPlot of the top 15 components, could be used as input for 
# the runDR function with dr = "UMAP" and pca = 10.
pcs <- SingleCellExperiment::reducedDim(x = sce, type = "PCA")
variance <- apply(X = pcs, MARGIN = 2, FUN = var)

plot(x = 1:15, 
     y = variance, 
     xlab = "Principal components", 
     ylab = "Variance")




