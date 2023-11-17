#######################################
# Let’s practice – 6                  #
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


