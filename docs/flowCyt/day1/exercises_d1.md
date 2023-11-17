
In this section, you will find the R code that we will use during the course. We will explain the code and output during correction of the exercises.

<!-- This is commented text -->

## Starting to work with flow cytometry data in R - packages

One package that provides data structures and basic functions to deal with flow cytometry data is [flowCore](https://bioconductor.org/packages/release/bioc/html/flowCore.html). flowCore allows to import data contained within a FCS file and store it in a flowFrame object. Data combined from several FCS files are imported and stored in a flowSet object. 

Packages that provide functions for automatic quality control are [flowAI](https://bioconductor.org/packages/release/bioc/html/flowAI.html) and [PeacoQC](https://bioconductor.org/packages/release/bioc/html/PeacoQC.html). 

For visualization, the package [ggcyto](https://www.bioconductor.org/packages/release/bioc/html/ggcyto.html
) provides plotting functions as an interface to ggplot2 plots using flowFrame or flowSet objects. One example is the autoplot() function.

Another package available for visualization is [flowViz](https://bioconductor.org/packages/release/bioc/html/flowViz.html), which includes the densityplot() function.


## Let's practice - 1

In this exercise we will use a 36-color spectral flow cytometry dataset from a study performed in the context of Covid-19 research. Only a subset from 5 healthy donors will be used. For each healthy donors, there are three time points, as indicated in the FCS file names. Data was downloaded through the Flow Repository database [(FR-FCM-Z3WR)](https://flowrepository.org/id/FR-FCM-Z3WR). FCS files were pre-gated on live CD3+CD19- T cells in FlowJo.

Create a new script in which you will:

1) Import the FCS files (within course_datasets/FR_FCM_Z3WR/) into a flowSet. Do not transform or truncate the values.

2) Create a data frame with the list of channels and corresponding antigens, and view it. Hint: get the antigens from the parameters of one of the flowFrame in the set

3) Add a new column to the phenotypic data with the time point of the sample. View the phenotypic data

4) Convert the channel names in the expression matrices to the corresponding antigen names (where applicable).

5) Create a bivariate density plot showing «FSC-H» againts «HLA-DR» for all samples from day 0. Apply a flowJo inverse hyperbolic sine scale to the y axis («HLA-DR»)


??? done "Answer"
	```r
    # load libraries
    library(flowCore)
    library(ggcyto)
    
    # 1) Import the FCS files (course_datasets/FR_FCM_Z3WR/). 
    # Do not transform or truncate the values. 

    # path to the directory with the fcs files
    fcs.dir<- file.path( "course_datasets/FR_FCM_Z3WR/")

    # read fcs files into a floSwet
    fcs_data <- read.flowSet(path=fcs.dir, pattern="*.fcs", transformation = FALSE, truncate_max_range = FALSE) 

    #2) Create a data frame with the list of channels and corresponding antigens, and view it.
    #Hint: get the antigens from the parameters of one of the flowFrame in the set
    channels <- colnames(fcs_data)
    antigen <- pData(parameters(fcs_data[[1]]))$desc
    panel <- data.frame(channels = channels, antigen= antigen)

    # view the panel
    panel

    # write panel to csv file
    write.csv(panel,file = "course_datasets/FR_FCM_Z3WR/panel.csv", quote = FALSE, row.names = FALSE)

    #3) Add a new column to the phenotypic data with the time point of the sample

    # check sample names
    sampleNames(fcs_data)
    # [1] "0E1F8E_0.fcs"  "0E1F8E_14.fcs" "0E1F8E_7.fcs"  "180E1A_0.fcs"  "180E1A_14.fcs" "180E1A_7.fcs" 
    # [7] "1A9B20_0.fcs"  "1A9B20_14.fcs" "1A9B20_7.fcs"  "61BBAD_0.fcs"  "61BBAD_14.fcs" "61BBAD_7.fcs" 

    # add column with time point
    pData(fcs_data)$time_point <- rep(c("D0","D14","D7"),5)

    # View the phenotypic data
    pData(fcs_data)

    # save flowSet for next exercise
    save(fcs_data,file="course_datasets/FR_FCM_Z3WR/fcs_data.RData")

    #4) Convert the channel names in the expression matrices to the corresponding antigen names (where applicable)

    colnames(fcs_data)[!is.na(antigen)] <- antigen[!is.na(antigen)] 

    # check
    head(exprs(fcs_data[[1]])[,c(5:10)])

    # 5) Create a bivariate density plot showing "FSC-H" against "HLA-DR" for all samples from day 0. 
    # Apply a flowJO inverse hyperbolic sine scale to the y axis ("HLA-DR")

    # split by time point 
    fcs_data.split <- split(fcs_data, pData(fcs_data)$time_point)

    # create the bivariate density plot
    autoplot(fcs_data.split$D0, x="FSC-H",y="HLA-DR", bins = 64) + 
     scale_x_flowjo_biexp() + 
      scale_y_flowjo_fasinh()

    ```

## Let's practice - 2
We will use the flowSet created in the previous exercise, and transform the data using two sets of cofactors: fixed and estimated using a function from the flowVS package.

Create a new script in which you will:

1) Load the flowSet object saved at the end of the previous exercise.

2) Read the «course_datasets/FR_FCM_Z3WR/panel.csv» file into a data frame. The last column contains the marker classes («none», «type» or «state»).

3) Downsample the flowSet to 2’000 cells per flowFrame (you can find the downsampling function in the «course_datasets/function_for_downsampling_flowSets.R» file).

4) Transform the «type» and «state» markers using  both Logicle (hints: use the downsampled flowSet for parameter estimation; start with default parameters, and adjust if needed) and arcsinh transformations (fixed cofactors of 3000).

5) Compare the transformation in the first flowFrame using density plots.

??? done "Answer"
	```r
    # load the libraries
    library(flowCore)
    library(flowVS)
    library(flowViz)

    # 1) load the flowSet object from previous exercise
    load("course_datasets/FR_FCM_Z3WR/fcs_data.RData")

    # 2) Add marker_class to panel
    # load panel from previous exercise
    panel <- read.csv("course_datasets/FR_FCM_Z3WR/panel.csv")

    # Set the marker classes
    panel$marker_class <- rep("none", nrow(panel))
    panel$marker_class[c(7:10,11:15,17:18,20,21,23:29,31:36,38,41,42)] <- "state"
    panel$marker_class[c(16,19,22,37,39,40)] <- "type"
    panel$marker_class <- factor(panel$marker_class,levels = c("type","state","none"))

    # write new panel to csv file
    write.csv(panel,file = "course_datasets/FR_FCM_Z3WR/panel_with_marker_classes.csv", quote = FALSE, row.names = FALSE)

    # View the panel
    panel

    #3) Downsample the flowSet to 2'000 cells per flowFrame for parameter estimation

    # load the function for downsampling a flowset
    source("course_datasets/function_for_downsampling_flowSets.R")

    # downsample to 2000 cells
    fcs_data_small <- Downsampling_flowSet(x = fcs_data, samplesize = 2000)

    #4) Transform using Logicle and Arcsinh transformation (fixed cofactors)

    # select markers to be transformed
    markerstotransform <- panel$channels[panel$marker_class!="none"]

    # transform with Logicle
    fcs_list <- list()
    for(i in 1:length(fcs_data)){
  
      algcl <- estimateLogicle(fcs_data_small[[i]],
                           channels = markerstotransform, m=6, t=4E6)
  
      fcs_list[[i]] <- transform(fcs_data[[i]], algcl)
  
    }

    fcs_transform_logicle <- as(fcs_list, "flowSet")
    sampleNames(fcs_transform_logicle) <- sampleNames(fcs_data)
    pData(fcs_transform_logicle) <- pData(fcs_data)

    # transform with fixed cofactors
    fcs_transform_arcsinh <- transFlowVS(fcs_data, 
                                     channels = markerstotransform, 
                                     rep(3000, length(markerstotransform)))

    sampleNames(fcs_transform_arcsinh) <- sampleNames(fcs_data)
    pData(fcs_transform_arcsinh) <- pData(fcs_data)

    # 5) Density plots 
    densityplot( ~ ., fcs_transform_logicle[[1]]) # worst
    densityplot( ~ ., fcs_transform_arcsinh[[1]]) # worst
    
    # save
    save(fcs_transform_logicle, markerstotransform,file = "course_datasets/FR_FCM_Z3WR/fcs_transform_logicle.RData")

    ```

## Let's practice - 3

We will continue with the Logicle transformed flowSet created in the last exercise, and apply the flowAI quality control algorithm to remove low quality cells.

Create a new script in which you will:

1) Load the flowSet object from exercice 2 («/course_datasets/FR_FCM_Z3WR/fcs_transform_logicle.RData»).

2) Run the flowAI quality control algorithm. Set the output directory to «course_datasets/FR_FCM_Z3WR/flowAI_res».

3) Load the «Qcmini.txt» report created by flowAI and view it.

4) Check the html report for sample 1A9B20_0. What happened ?

??? done "Answer"
	```r
    # load libraries
    library(flowCore)
    library(flowAI)

    # 1) load the flowSet object from previous exercise and 
    load("course_datasets/FR_FCM_Z3WR/fcs_transform_logicle.RData")

    # 2) Run the flowAI quality control algorith. 
    # Output the results to a"course_datasets/FR_FCM_Z3WR/flowAI_res/"
    
    fcs_clean <- flow_auto_qc(fcs_transform_logicle, 
                          folder_results = "course_datasets/FR_FCM_Z3WR/flowAI_res/")

    # save clean flowSet
    save(fcs_clean, file = "course_datasets/FR_FCM_Z3WR/fcs_clean.RData")

    # 3) Load and view the report created by flowAI

    # load
    QCmini <- read.delim("course_datasets/FR_FCM_Z3WR/flowAI_res/QCmini.txt")

    # change the names of the columns
    names(QCmini) <- gsub("X..","% ", names(QCmini))

    # View
    QCmini

    ```


**End of Day 1, good job!** :sparkle:

<!--
## Feedback :sparkle:
-->









