
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

1) Create a flowSet called fcs_data of all samples within the /course_dataset/FR_FCM_Z4KT folder

2) Generate a panel data.frame using colnames(fcs_data) antigen names extracted with pData(parameters(fcs_data[[1]]))$desc. Create a new column called marker_class that will contain the type of markers: all that are not NA should be labeled as "type", except PD-1 which should be labeled as "state". Make sure that the antigen "Zombie UV" is labeled as "none" and not as "type". Save the panel to an Excel file using write.xlsx2().

3) Transform the data: extract a vector from the panel data.frame which are the channels to be transformed, which are not labeled with "none". Perform asinh transformation with a cofactor of 3000 for all channels to be transformed, using transFlowVS() from the flowVS package.

4) Split the flowSet resulting from transformation into a training flowSet containing all flowFrames from the sample "REU271", and a flowSet with the rest of the flowFrames not corresponding to sample "REU271".

5) Perform pre-clustering with flowSOM with function prepareFlowSOM(), providing the flowSet with the training flowFrames, the vector of channels to transform, and FlowSOM.params =list(xdim=10, ydim=10,
nClus=20, scale=FALSE)

6) Test the coefficient of variation within clusters with the testCV() function.

7) Import the metadata with the batch label of each sample contained in the excel file md.xlsx, using read.xlsx2(). Create 2 vectors using the column "batch" in the md.xlsx file. One vector contains the batch labels of the samples that correspond to sample "REU271", and another vector contains the batch labels of the other samples (i.e. not "REU271").

8) Estimate quantiles from the training flowSet using CytoNorm.train(). Use FlowSOM.params = list(nCells = 6000, xdim = 10, ydim = 10, nClus = 5, scale = FALSE)

9) Normalize the rest of the samples using CytoNorm.normalize(), and using outputDir = "course_datasets/FR_FCM_Z4KT/Normalized"  ; Make sure this is a new folder.


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

1) Import a xml workspace from FlowJo

2) Plot the gating hierarchy

3) Plot the gates

4) Add a gate

5) Statistics

6) Export flow data


??? done "Answer"
	```r
    # load libraries
    library(flowCore)
    
    ```

## Let's practice - 10

??? done "Answer"
	```r
    # load libraries
    library(flowCore)
    
    ```




**End of Day 3, good job!**
