
In this section, you will find the R code that we will use during the course. We will explain the code and output during correction of the exercises.

<!-- This is commented text -->

## Normalization, gating - packages

[CytoNorm](https://github.com/saeyslab/CytoNorm) is a package that allows to perform batch correction, i.e. normalization of samples when technical replicates of a single and same sample were run on the several batches used to run the experimental samples.

[flowWorkspace](https://bioconductor.org/packages/release/bioc/html/flowWorkspace.html) provides the GatingSet class of objects as an efficient data structure to store, query and visualize gated flow data.

[CytoML](https://bioconductor.org/packages/release/bioc/html/CytoML.html) uses platform-specific implementations of the GatingML2.0 standard to exchange gated cytometry data. 

[flowClust](https://www.bioconductor.org/packages/release/bioc/html/flowClust.html) can be used for automated gating. It can help in identifying cell populations in flow cytometry data. Robust Model-based Clustering of Flow Cytometry Data (Lo et al. 2008)

[openCyto](https://bioconductor.org/packages/release/bioc/html/openCyto.html) implements a hierachical gating pipeline for flow cytometry data. 

[flowGate](https://bioconductor.org/packages/release/bioc/html/flowGate.html) provides interactive cytometry gating in R, based on a shiny app (web application using R). It is especially geared toward wet-lab cytometerists looking to take advantage of R without having a lot of experience. 

## Let's practice - 8

In this exercise, we will perform normalization using CytoNorm, and estimating quantiles and splines using a technical replicate distributed across batches. The data comes from the FlowRepository accession .

1) Create a flowSet called fcs_data of all samples within the `/course_dataset/FR_FCM_Z4KT` folder

2) Generate a panel data.frame using `colnames(fcs_data)` and antigen names extracted with `pData(parameters(fcs_data[[1]]))$desc`. Create a new column called marker_class that will contain the type of markers: all the ones that are not NA should be labeled as "type", except PD-1 which should be labeled as "state". Make sure that the antigen "Zombie UV" is labeled as "none" and not as "type". Save the panel to an Excel file using `write.xlsx2()`.

3) Transform the data: extract a vector from the panel data.frame which are the channels to be transformed, which are the ones that are not labeled with "none". Perform asinh transformation with a cofactor of 3000 for all channels to be transformed, using `transFlowVS()` from the flowVS package.

4) Split the flowSet resulting from transformation into a training flowSet containing all flowFrames from the sample "REU271", and a flowSet with the rest of the flowFrames not corresponding to sample "REU271". Use the `grep()` function on the `sampleNames` of the flowSet.

5) Perform pre-clustering with flowSOM with function `prepareFlowSOM()`, providing the flowSet with the training flowFrames, the vector of channels to transform, and `FlowSOM.params=list(xdim=10, ydim=10,
nClus=20, scale=FALSE)`.

6) Test the coefficient of variation within clusters with the `testCV()` function.

7) Import the metadata with the batch label of each sample contained in the excel file [md.xlsx](https://github.com/taniawyss/flow-cytometry-analysis-with-R/tree/master/docs/flowCyt/assets/data), using `read.xlsx2()`. Create 2 vectors using the column "batch" in the md.xlsx file. One vector contains the batch labels of the samples that correspond to sample "REU271", and another vector contains the batch labels of the other samples (i.e. not "REU271").

8) Estimate quantiles from the training flowSet using `CytoNorm.train()`. Use `FlowSOM.params = list(nCells = 6000, xdim = 10, ydim = 10, nClus = 5, scale = FALSE)`.

9) Normalize the rest of the samples using `CytoNorm.normalize()`, and using `outputDir = "course_datasets/FR_FCM_Z4KT/Normalized"`  ; Make sure this is a new folder.

10) Choosing one channel, create a ridge plot of its distribution within samples before normalization (without the training samples), and one for the normalized samples. For this, you need to create a new flowSet with the created "Norm_" fcs files within the newly created output folder. Use the `densityplot()` function for each flowSet, storing the output in 2 objects, then use the cowplot `plot_grid()` function to plot one ridge plot above the other.

??? done "Answer"
	```r
    # load libraries
    library(flowCore)
    library(CytoNorm)
    library(xlsx) # package to export and import data within Excel
    library(cowplot) # package with add-on functions for ggplot2 plotting

    # 1) generate a flowSet with the fcs files from the FR_FCM_Z4KT data
    # path to the directory (folder) with the fcs files
    fcs.dir<- file.path("course_datasets/FR_FCM_Z4KT/")

    # read fcs files into a flowSet
    fcs_data <- read.flowSet(path=fcs.dir, 
                         pattern="*.fcs", 
                         transformation = FALSE, 
                         truncate_max_range = FALSE) #fcs_data will be a FlowSet object
    # Check the object:
    fcs_data

    # retrieve the list of channels and corresponding antigens
    fcs_colname <- colnames(fcs_data) 
    antigen <- pData(parameters(fcs_data[[1]]))$desc

    # setting marker classes. In this panel, only PD-1 is a state marker:
    marker_class <- rep("none", ncol(fcs_data[[1]])) # setting all to "none"
    marker_class[which(!is.na(antigen))] <- "type" # type markers
    marker_class[which(antigen=="PD-1")] <- "state" # state markers
    # change UV marker to "none"
    marker_class[which(antigen=="Zombie UV")] <- "none" # state markers
    marker_class <- factor(marker_class,levels=c("type","state","none"))

    # put everything together in a data frame
    panel <- data.frame(fcs_colname, antigen, marker_class, row.names = NULL) 
    # check that the type/none/state column is correct:
    panel
    # Export to excel file, which can be read in with read.xlsx2()
    xlsx::write.xlsx2(panel, file="course_datasets/FR_FCM_Z4KT/panel_Z4KT.xlsx", 
                  sheetName="panel_Z4KT")

    # Import an Excel sheet as a data.frame:
    # Using sheet index
    panel<-read.xlsx2("course_datasets/FR_FCM_Z4KT/panel_Z4KT.xlsx", sheetIndex = 1)
    # using sheet name:
    panel<-read.xlsx2("course_datasets/FR_FCM_Z4KT/panel_Z4KT.xlsx", sheetName =  "panel_Z4KT")

    ## Arcsinh transformation with a fixed cofactor of 3000
    # Select markers to be transformed:
    markerstotransf <- panel$fcs_colname[panel$marker_class!="none"]

    # set cofactor vector
    cofactor <- 3000
    l <- length(markerstotransf)             
    cofactors <- rep(cofactor, l)

    # transform
    fcs_transform <- transFlowVS(fcs_data,
                             channels = markerstotransf,
                             cofactors)
    # the output is a flowSet:
    fcs_transform

    # save to a file for downstream analysis
    save(fcs_transform, file = "course_datasets/FR_FCM_Z4KT/fcs_transform.RData")

    # Used transformed data for CytoNorm:
    # Separate the files according to training and validation set:
    # sample REU271 was measured on several days:
    train_files <- fcs_transform[grep("REU271", sampleNames(fcs_transform))]
    validation_files <- fcs_transform[-c(grep("REU271", sampleNames(fcs_transform)))]

    # Pre-clustering with FlowSOM:
    ?CytoNorm::prepareFlowSOM

    fsom <- prepareFlowSOM(train_files, 
                       colsToUse = markerstotransf, 
                       transformList = NULL, 
                       FlowSOM.params = list(xdim=10,
                                             ydim=10, 
                                             nClus=20, scale=FALSE))
    fsom

    # Check coefficient of variation within cluster:
    # Function to inspect whether all control samples contain a 
    # similar percentage of cells in all FlowSOM clusters
    cvs <- CytoNorm::testCV(fsom, cluster_values = c(5,10,15,20))
    range(cvs$cvs$`20`) # 0.05758965 1.43114512

    # If the clusters are impacted by batch effects, CV values of >1.5 or 2 will
    # occur, then you can choose to put FlowSOM.params to NULL 
    # and skip clustering.

    # Import sample metadata and batch info, which should include Sample_ID
    # and batch columns:
    md <- xlsx::read.xlsx2("course_datasets/FR_FCM_Z4KT/md.xlsx", sheetIndex = 1)
    head(md)

    # Check that the order of fcs files within folder is the same as within the metadata:
    file_names <- list.files(fcs.dir, pattern = ".fcs")

    summary(md$file_name == file_names) 
    #    Mode    TRUE 
    # logical      16 

    # extract the batch labels of the training sample:
    labels_train <- md$batch[grep("REU271", md$Sample_ID)]
    labels_train
    # [1] "B" "C" "D" "E" "F" "A"

    model <- CytoNorm.train(files = train_files,
                        labels = labels_train,
                        channels = markerstotransf,
                        transformList = NULL,
                        FlowSOM.params = list(nCells = 6000, 
                                              xdim = 10,
                                              ydim = 10,
                                              nClus = 5,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)
    # View the quantile values:
    model$clusterRes$`5`

    # Normalize the rest of the files:
    label_norm <- md$batch[-c(grep("REU271", md$Sample_ID))]

    CytoNorm.normalize(model = model,
                   files = validation_files,
                   labels = label_norm,
                   transformList = NULL,
                   transformList.reverse = NULL,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "course_datasets/FR_FCM_Z4KT/Normalized",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE)

    fcs.dir<- file.path("course_datasets/FR_FCM_Z4KT/Normalized")
    fcs_norm <- read.flowSet(path=fcs.dir, 
                         pattern="*.fcs", 
                         transformation = FALSE, 
                         truncate_max_range = FALSE)

    # Compare the distribution a marker across replicates
    colnames(fcs_norm)
    # before normalization (without plotting the training samples)
    p1<-densityplot(~`FJComp-BUV496-A`, 
            fcs_transform[-c(grep("REU271", sampleNames(fcs_transform)))])
    # after normalization
    p2<-densityplot(~`FJComp-BUV496-A`, fcs_norm)

    cowplot::plot_grid(p1, p2, nrow = 2)

    ```

## Let's practice - 9

In this exercise we will do some gating using flowGate and data from the FR_FCM_Z3WR of the FlowRepository.

Create a new script in which you will:

1) Create a flowSet of all samples within the `/course_dataset/FR_FCM_Z3WR` folder

2) Perform asinh transformation with a cofactor of 3000 for all channels not labeled with "none", using `transFlowVS()` from the flowVS package. Use the csv file with the panel previously created `"/course_datasets/FR_FCM_Z3WR/panel_with_marker_classes.csv"`.

3) Convert the flowSet to a GatingSet

4) Using flowGate, create a gating hierarchy according to the scheme depicted below. Don't forget to check your gating with scatter or density plots.

5) Do necessary adjustements so that your gating hierarchy looks like the one depicted below.

6) What is the percentage of CD8+ T cells among T cells ("CD3").

Gating hierarchy:

<figure>
  <img src="../../assets/images/gates9.png" width="700"/>
</figure>



??? done "Answer"
	```r
    # load libraries
    library(flowCore)
    library(flowWorkspace)
    library(flowGate)
    library(flowVS)

    # 1) Create a flowSet from the FCS files in "course_datasets/FR_FCM_Z3WR/"

    fs <- read.flowSet(path="course_datasets/FR_FCM_Z3WR/",
                   pattern = "*.fcs",
                   transformation = FALSE,
                   truncate_max_range = FALSE)

    # 2 ) FlowVS Arcsinh transformation with fixed factors (3000). 
    # Use the csv file containing the panel and marker classes previously created 
    # ("course_datasets/FR_FCM_Z3WR/panel_with_marker_classes.csv")


    panel <- read.csv("course_datasets/FR_FCM_Z3WR/panel_with_marker_classes.csv")
    markerstotransform <- as.character(panel$channels)[!is.na(panel$antigen)]

    fs <- transFlowVS(fs,
                  channels = markerstotransform,
                  cofactors = rep(3000,length(markerstotransform)))


    # 3) Convert the flowSet to a GatingSet

    gs <- GatingSet(fs)

    # save for next exercise
    #save_gs(gs, path = "course_datasets/FR_FCM_Z3WR/gs_preprocessed")


    # 4) Using flowGate, create the following gates

    # Polygon gate ("Leukocytes")
    gs_gate_interactive(gs,
                    filterId = "Leukocytes",
                    dims = list("FSC-H", "SSC-H"))

    autoplot(gs[[1]], gate = "Leukocytes")
    plot(gs)

    # span gate ("CD3")
    gs_gate_interactive(gs,
                    filterId = "CD3",
                    dims = "BV510-A",
                    subset = "Leukocytes")

    autoplot(gs[[1]],gate = "CD3")
    plot(gs)

    # quadrant gate ("CD4 CD8")
    # BUV615-A = CD4
    # BUV805-A = CD8
    gs_gate_interactive(gs,
                    filterId = "CD4 CD8",
                    dims = list("BUV615-A","BUV805-A"),
                    subset = "CD3")

    plot(gs)
    autoplot(gs[[1]],gate = gs_pop_get_children(gs, "CD3"))

    # 5) Make the necessary adjustements so that the gate tree loks like the one depicted.

    # Check node names
    gs_get_pop_paths(gs, path = 1)

    # Rename the CD4 CD8 gates to CD4+, CD8+, DNT and DPT
    gs_pop_set_name(gs,"BUV615-A-BUV805-A+","CD8+")
    gs_pop_set_name(gs,"BUV615-A+BUV805-A-","CD4+")
    gs_pop_set_name(gs,"BUV615-A+BUV805-A+","DPT")
    gs_pop_set_name(gs,"BUV615-A-BUV805-A-","DNT")

    plot(gs)

    # 6) What is the percentage of CD8+ T cells among T cells ("CD3")
    gs_pop_get_stats(gs, nodes = "CD8+", type = "percent")
    
    ```

## Let's practice - 10

In this exercise we will repeat the gating previously done with flowGate on the data from FR_FCM_Z3WR of the FlowRepository, but this time using automated gating with openCyto.

Create a new script in which you will:

1) Repeat steps 1 to 4 from the previous exercise (loading data from fcs files and preprocessing). You can also load the preprocessed GatingSet from `/course_dataset/FR_FCM_Z3WR/gs_preprocessed/`.

2) Using the `gs_add_gating_method()` function (i.e., without a template), create a gating hierarchy according to the scheme depicted below. Don't forget to check your gating with scatter or density plots.

3) Do necessary adjustements so that your gating hierarchy looks like the one depicted below (hide the "BUV805-A+" and "BUV615-A+" nodes from the tree, and rename the CD4, CD8, DNT and DPT nodes).

4) Create boxplots showing the percentage of CD8+ T cells among T cells ("CD3") as a function of time points.

Gating hierarchy:

<figure>
  <img src="../../assets/images/gates10.png" width="700"/>
</figure>


??? done "Answer"
	```r
    # load libraries
    library(flowCore)
    library(flowWorkspace)
    library(openCyto)
    library(ggcyto)

    # 1) Repeat steps 1 to 4 from previous exercise (loading data, preprocessing)

    # load the preprocessed ungated data
    gs <- load_gs("course_datasets/FR_FCM_Z3WR/gs_preprocessed")

    # 2) Apply the gating as in the scheme

    # gate "Leukocytes"
    # use flowClust
    # set a K=3, to split in three populations and select the one with the highest "peak"
    gs_add_gating_method(gs,
                     alias = "Leukocytes",
                     parent = "root",
                     dims = "FSC-H,SSC-H",
                     gating_method = "flowClust",
                     gating_args = "K=3" )

    # check
    plot(gs)

    autoplot(gs[[1]], gate="Leukocytes")

    # gate "CD3"
    # use minDensity
    gs_add_gating_method(gs,
                     alias = "CD3",
                     parent = "Leukocytes",
                     dims = "BV510-A",
                     gating_method = "gate_mindensity",
                     pop = "+")

    # check
    plot(gs)

    autoplot(gs[[1]], gate="CD3")

    # gate "CD4 CD8"
    # use minDensity
    gs_add_gating_method(gs,
                     alias = "*",
                     parent = "CD3",
                     dims = "BUV615-A,BUV805-A",
                     gating_method = "gate_mindensity",
                     pop = "-/++/-")

    # check
    plot(gs)

    gs_pop_get_children(gs,"CD3")
    autoplot(gs[[1]], gate=gs_pop_get_children(gs,"CD3"))
    autoplot(gs[[1]], gate=gs_pop_get_children(gs,"CD3")[3:6])

    # 3) 

    # Hide the "BUV805-A+" and "BUV615-A+" nodes from the tree

    gs_pop_set_visibility(gs,"BUV805-A+", FALSE)
    gs_pop_set_visibility(gs,"BUV615-A+", FALSE)

    # check
    plot(gs)

    # 4) Rename the CD4 CD8 gates to CD4+, CD8+, DNT and DPT

    # Rename the CD4 CD8 gates to CD4+, CD8+, DNT and DPT
    gs_pop_set_name(gs,"BUV615-A-BUV805-A+","CD8+")
    gs_pop_set_name(gs,"BUV615-A+BUV805-A-","CD4+")
    gs_pop_set_name(gs,"BUV615-A+BUV805-A+","DPT")
    gs_pop_set_name(gs,"BUV615-A-BUV805-A-","DNT")

    # check
    plot(gs)

    # 5) Create a boxplot showing the percentage of CD8+ T cells among T cells ("CD3") as a function of time point

    # extract the proportions
    my_proportions <- gs_pop_get_stats(gs, nodes = "CD8+", type = "percent")

    # add the timepoints
    my_proportions$time_point <- rep(c(0,14,7), 5)

    # convert to factor (order in the plot)
    my_proportions$time_point <- factor(my_proportions$time_point, levels = c(0,7,14))

    # create the boxplot
    boxplot(percent ~ time_point, data = my_proportions)
    
    ```


**End of Day 3, good job!**
