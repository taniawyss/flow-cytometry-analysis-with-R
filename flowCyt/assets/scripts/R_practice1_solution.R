#########################
# Let’s practice – 1    #
# Import data           #
# Solutions             #
#########################

# In this exercise we will use a 36-color spectral flow cytometry dataset from a 
# study performed in the context of Covid-19 research. 
# Only a subset from 4 healthy donors will be used. 
# For each healthy donor, there are three time points, as indicated in FCS file names. 
# Data was downloaded through the Flow Repository database (FR-FCM-Z3WR) 
# at https://flowrepository.org/id/FR-FCM-Z3WR. 
# FCS files were pre-gated on living CD3+CD19-T cells in FlowJo.

setwd("/export/scratch/twyss/SIB_training/flowCyt_2023/flowCyt_nov2023/data/")
# clear the environment
rm(list = ls())
gc()


# load libraries
library(flowCore)
library(ggcyto)

# 1) Import the FCS files (course_datasets/FR_FCM_Z3WR/). 
# Do not transform or truncate the values. 


# path to the directory (folder) with the fcs files
fcs.dir<- file.path("course_datasets/FR_FCM_Z3WR/")

# read fcs files into a flowSet
fcs_data <- flowCore::read.flowSet(path = fcs.dir, 
                         pattern = "*.fcs", 
                         transformation = FALSE, truncate_max_range = FALSE) 
# Explore the object:
fcs_data

summary(fcs_data[[1]])

#2) Create a data frame with the list of channels and corresponding antigens, and show it.
#Hint: get the antigens from the parameters of one of the flowFrame in the set
channels <- colnames(fcs_data)
antigen <- pData(parameters(fcs_data[[1]]))$desc
# fcs_data@frames$`0BF51C_0.fcs`@parameters$desc
panel <- data.frame(channels = channels, antigen= antigen)

fcs_data@frames$`0BF51C_0.fcs`@parameters@data

# show the panel
panel

# write panel to csv file
write.csv(panel,file = "course_datasets/FR_FCM_Z3WR/panel.csv", quote = FALSE, row.names = FALSE)

#3) Add a new column to the phenotypic data with the time point of the sample

# check sample names
sampleNames(fcs_data)
# [1] "0E1F8E_0.fcs"  "0E1F8E_14.fcs" "0E1F8E_7.fcs"  "180E1A_0.fcs"  "180E1A_14.fcs" "180E1A_7.fcs" 
# [7] "1A9B20_0.fcs"  "1A9B20_14.fcs" "1A9B20_7.fcs"  "61BBAD_0.fcs"  "61BBAD_14.fcs" "61BBAD_7.fcs" 
# [13] "61BBAD_0.fcs"  "61BBAD_14.fcs" "61BBAD_7.fcs"
# add column with time point
pData(fcs_data)$time_point <- rep(c("D0","D14","D7"),5)

# show the phenotypic data
pData(fcs_data)

# save flowSet for next exercise
save(fcs_data,file="course_datasets/FR_FCM_Z3WR/fcs_data.RData")


#4) Convert the channel names in the expression matrices to the corresponding 
# antigen names (where applicable)
colnames(fcs_data)[!is.na(antigen)] <- antigen[!is.na(antigen)] 

# check that the antigen name change was effective:
head(exprs(fcs_data[[1]])[,c(5:10)])


# 5) Create a bivariate density plot showing "FSC-H" against "HLA-DR" for all samples from day 0. 
# Apply a flowJO inverse hyperbolic sine scale to the y axis ("HLA-DR")

# split by time point 
fcs_data.split <- split(fcs_data, pData(fcs_data)$time_point)
class(fcs_data.split) # list
class(fcs_data.split$D0) # flowSet

# create the bivariate density plot
ggcyto::autoplot(fcs_data.split$D0, x="FSC-H",y="HLA-DR", bins = 64) + 
  ggcyto::scale_x_flowjo_biexp() + 
  ggcyto::scale_y_flowjo_fasinh()

# FSC-H forward scatter height
# FSC-A forward  scatter area
# SSC-H side scatter height
# SSC-A side scatter area



