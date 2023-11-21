
In this section, you will find the R code that we will use during the course. We will explain the code and output during correction of the exercises.

<!-- This is commented text -->

## Dimensionality reduction, clustering and differential testing - packages

[CATALYST](https://bioconductor.org/packages/release/bioc/html/CATALYST.html) provides functions for preprocessing of cytometry data such as FACS, CyTOF, and IMC, as well as functions for dimensional reduction, clustering and methods for differential composition and expression analysis. The CATALYST package provides a function to first cluster data with [FlowSOM](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html) clustering and then apply [ConsensusClusterPlus](https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html) metaclustering.

CATALYST requires the data be contained within an object of class [SingleCellExperiment](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) of Bioconductor.

An R implementation of the Uniform Manifold Approximation and Projection (UMAP) method for dimensionality reduction is available in the CRAN package [uwot](https://cran.r-project.org/web/packages/uwot/index.html). 

The [diffcyt](https://bioconductor.org/packages/release/bioc/html/diffcyt.html) package provides     statistical methods for differential discovery analyses in high-dimensional cytometry data (including flow cytometry and mass cytometry).

## Let's practice - 4

n this exercise we will continue with the clean flowSet from the last exercise. We will use the CATALYST package to create a SingleCellExperiment (sce) object,  perform dimensionality reduction (UMAP) and use the UMAP to plot the expression of markers.

Create a new script in which you will:

1) Load the clean flowSet from last exercise («fcs_clean.Rdata»).

2) Downsample the flowSet to 2’000 cells per flowFrame (source the file «function_for_downsampling_flowSets.R»)

3) Create a sce object from the downsampled flowSet.

4) Create a UMAP with default parameters, based on the expression of the «type» markers. Show the expression of CD3 by time point.

5) Check the effect of changing parameters «min_dist» and «n_neighbors» from the default values.


??? done "Answer"
	```r
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
                panel_cols = list(channel = "channels", antigen = "antigen", class =    "marker_class"),
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

    ```

## Let's practice - 5

In this exercise we will apply the FlowSom method for unsupervised clustering of cells, followed by ConsensusClusterPlus metaclustering. We then check the expression of markers by metacluster. Finally, we will rename / merge the metaclusters to annotate major cell populations.

Create a new script in which you will:

1) Load the sce object with UMAP from the previous exercise ("course_datasets/FR_FCM_Z3WR/sce_UMAP.RData").

2) Apply FlowSOM clustering + ConsensusClusterPlus metaclustering.

3) Plot a UMAP showing the location of metaclusters; marker expression heatmap and ridge plots. Use 8 metaclusters.

4) Rename / merge metaclusters as major cell populations according to the expression of markers.

5) Plot a UMAP showing the major cell populations.

??? done "Answer"
	```r
    # load libraries
    library(flowCore)
    library(CATALYST)

    # 1) Load the sce object with UMAP from previous exercise
    load("course_datasets/FR_FCM_Z3WR/sce_UMAP.RData")

    # 2) Apply FlowSOM clustering + ConsensusClusterPlus metaclustering

    set.seed(1234)
    sce <- cluster(sce, 
               features = "type")

    names(cluster_codes(sce))
    # [1] "som100" "meta2"  "meta3"  "meta4"  "meta5"  "meta6"  "meta7"  "meta8"  "meta9"  "meta10"
    # [11] "meta11" "meta12" "meta13" "meta14" "meta15" "meta16" "meta17" "meta18" "meta19" "meta20"

    # 3) Plot UMAP with clusters, expression heatmap and ridge plots

    # UMAP
    plotDR(sce,
       dr = "UMAP",
       color_by = "meta8" )

    # Heatmap 
    plotExprHeatmap(sce,row_clust = F, col_clust = F,
                features = "type",
                by="cluster_id",
                k="meta8")

    # Ridge plots
    plotClusterExprs(sce,
                 k="meta10")

    # 4) Rename / merge clusters
    # create a merging table
    merging_table <- data.frame(old_cluster = 1:8,
                            new_cluster = c("CD4 T cells","Other","Monocytes","Other","CD8 T cells","CD8 T cells","Other","B cells"))

    # write to file
    write.csv(merging_table,file = "course_datasets/FR_FCM_Z3WR/merging_table.csv",quote = F,row.names = F)

    # annotate clusters
    sce <- mergeClusters(sce,
                     k="meta8",
                     table = merging_table,
                     id = "Major_cell_populations")

    # 5) UMAP with major cell populations
    plotDR(sce,
       dr = "UMAP",
       color_by = "Major_cell_populations")

    # save
    save(sce,file = "course_datasets/FR_FCM_Z3WR/sce_annotated.RData" )

    ```

## Let's practice - 6

In this exercise we will test if cell populations have significantly different abundances between two time points (D14 compared to D0).

Create a new script in which you will:

1) Load the sce object from the previous exercise ("sce_annotated.RData").

2) Plot relative cell population abundances by sample and time point.

3) Set up the design and contrast matrices.

4) Test for differences in abundances between D14 and D0.

5) View table of results

??? done "Answer"
	```r
    # load libraries
    library(CATALYST)
    library(diffcyt)

    # 1) Load the sce object with annotation from previous exercise
    load("course_datasets/FR_FCM_Z3WR/sce_annotated.RData")

    # 2) Plot relative population abundances by sample, grouped by time point

    sce$time_point <- factor(sce$time_point, levels = c("D0","D7","D14"))

    # stacked bar plot
    plotAbundances(sce, 
               k = "Major_cell_populations", 
			         by = "sample_id", 
			         group_by = "time_point")

    # 3) Set up designand contrast  matrices.

    # Set design matrix
    design <- createDesignMatrix(ei(sce),
                             cols_design = c("time_point"))
                             
    # check
    design

    # Set contrast matrix (D14 vs D0)
    contrast <- createContrast(c(0,0,1))

    # 4) Compute differential abundance and show top differentially abundant cell populations

    # compute DA
    res_DA <- diffcyt(sce, 
                  clustering_to_use = "Major_cell_populations",
				analysis_type = "DA", 
				method_DA = "diffcyt-DA-edgeR",
    				design = design, 
				contrast = contrast)

    # show top differentially abundant cell populations
    tbl_DA <- rowData(res_DA$res)
    tbl_DA

    topTable(res_DA, format_vals = T)

    ```

## Let's practice - 7

In this exercise we will test if markers were differentially expressed between two time points (D14 compared to D0).

Create a new script in which you will:

1) Load the sce object from the previous exercise ("sce_annotated.RData").

2) Set up the design and contrast matrices.

3) Test for differences in marker expression between D14 and D0.

4) View table of results.

??? done "Answer"
	```r
    # load libraries
    library(CATALYST)
    library(diffcyt)

    # 1) Load the sce object with annotation from previous exercise
    load("course_datasets/FR_FCM_Z3WR/sce_annotated.RData")

    # 2) Set up design and contrast matrices (D14 vs D0)

    design <- createDesignMatrix(ei(sce),
                             cols_design = c("time_point"))

    contrast <- createContrast(c(0,0,1))

    # 3) Compute differential state expression

    res_DS <- diffcyt(sce, 
                  clustering_to_use = "Major_cell_populations",
				analysis_type = "DS", 
				method_DS = "diffcyt-DS-limma",
    				design = design, 
				contrast = contrast)

    # Are there any differentially expressed markers ?

    topTable(res_DS)
    ```

**End of Day 2, good job!** :medal:

<!--
## Feedback :sparkle:
-->









