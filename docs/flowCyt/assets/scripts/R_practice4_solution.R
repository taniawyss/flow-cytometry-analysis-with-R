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


# 1) load the flowSet object from previous exercise
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
sce <- runDR(sce, 
             assay = "exprs", 
             dr = "UMAP", 
             features = "type", 
             cells = NULL)

# plot
plotDR(sce,
       dr =  "UMAP",
       assay = "exprs", 
       color_by="CD3", 
       facet_by = "time_point" )

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
# the previous UMAP coordinates are overwritten, we still have only 1 reducedDims element:
reducedDims(sce)

# plot
plotDR(sce,
       dr =  "UMAP",
       assay = "exprs", 
       color_by="CD3", 
       facet_by = "time_point" )


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
       facet_by = "time_point" )


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
       facet_by = "time_point" )

# Try tSNE! It is slower than UMAP
set.seed(1601)
sce <- runDR(sce, 
             assay = "exprs", 
             dr = "TSNE", 
             features = "type", 
             cells = NULL)
reducedDims(sce)
# List of length 2
# names(2): UMAP TSNE

# plot
plotDR(sce,
       dr =  "TSNE",
       color_by="sample_id",
       facet_by = "time_point")






