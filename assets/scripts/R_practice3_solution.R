#########################
# Let’s practice – 3    #
# Quality control       #
# Solutions             #
#########################


# clear the environment
rm(list = ls())


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

# 3) Load and plot the report created by flowAI

# load
QCmini <- read.delim("course_datasets/FR_FCM_Z3WR/flowAI_res/QCmini.txt")

# change the names of the columns
names(QCmini) <- gsub("X..","% ", names(QCmini))

# plot
QCmini





