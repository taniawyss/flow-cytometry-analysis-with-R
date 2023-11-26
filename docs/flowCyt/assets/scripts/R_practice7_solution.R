#######################################
# Let’s practice – 7                  #
# Differential abundance testing      #
# Solutions                           #
#######################################


# clear the environment
rm(list = ls())


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


## Adding a covariate such as patient id to perform paired analysis:
?diffcyt
# method_DS = c("diffcyt-DS-limma", "diffcyt-DS-LMM"),
?testDS_limma
# we can use the block_id argument which has to be a vector of patient IDs

# create a vector of patient IDs for block design:
patient_id <- ei(sce)$sample_id
patient_id <- gsub("_0.fcs", "", patient_id)
patient_id <- gsub("_14.fcs", "", patient_id)
patient_id <- gsub("_7.fcs", "", patient_id)

head(patient_id)
# [1] "0BF51C" "0BF51C" "0BF51C" "0E1F8E" "0E1F8E" "0E1F8E"

# 2) Set up design and contrast matrices

design <- createDesignMatrix(ei(sce), cols_design = c("time_point"))

contrast <- createContrast(c(0,0,1))

# 3) Compute differential state expression using a vector
# of patient IDs as block_id argument:

res_DS_paired <- diffcyt(sce, 
                  clustering_to_use = "Major_cell_populations",
                  analysis_type = "DS", 
                  method_DS = "diffcyt-DS-limma",
                  design = design, 
                  contrast = contrast,
                  block_id = patient_id)

# Are there any differentially expressed markers ?

topTable(res_DS_paired)


# Heatmap with time points re-ordered:
colData(sce)$sample_id <- factor(colData(sce)$sample_id, 
                          levels=c(ei(sce)$sample_id[grep("_0", ei(sce)$sample_id)],
                                   ei(sce)$sample_id[grep("_7", ei(sce)$sample_id)],
                                   ei(sce)$sample_id[grep("_14", ei(sce)$sample_id)]))

plotDiffHeatmap(sce, rowData(res_DS_paired$res), all=T, sort_by = "lfc", col_anno ="time_point")

