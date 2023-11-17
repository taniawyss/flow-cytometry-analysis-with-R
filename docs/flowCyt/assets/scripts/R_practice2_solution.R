#########################
# Let’s practice – 2    #
# Transformation        #
# Solutions             #
#########################


# clear the environment
rm(list = ls())


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


# plot the panel
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








