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

