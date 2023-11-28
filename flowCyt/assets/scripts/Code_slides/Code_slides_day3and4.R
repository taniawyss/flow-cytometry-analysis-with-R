##############################
# Code from slides           #
# Days 3 and 4               #
##############################

rm(list = ls())


###################################
# cytoML                          #
# Import a workspace from FlowJo  #
###################################

# First install and load flowWorkspaceData
#BiocManager::install("flowWorkspaceData")
library(flowWorkspaceData)

# We will use a xml file with workspace exported from FlowJo
# Available from the flowWorkspaceData package
# dataDir <- system.file("extdata",package="flowWorkspaceData")
# wsfile <- list.files(dataDir, pattern="manual.xml",full=TRUE)
# I copied this file to course_datasets/flowWorspaceData/manual.xml

# Install and load CytoML
#BiocManager::install("CytoML")
library(CytoML)


# import workspace from FlowJo
# ws <- open_flowjo_xml(wsfile)
ws <- open_flowjo_xml("course_datasets/flowWorspaceData/manual.xml")


# convert workspace to a GatingSet object
# we will only use the first two FCS files
# which were copied to course_datasets/flowWorspaceData
gs <- flowjo_to_gatingset(ws, name = "T-cell")


###############################
# flowWorspace                #
# basics on GatingSets        #
###############################

library(flowWorkspace)

# check names of samples in object
sampleNames(gs)

# check metadata
pData(gs)

# add metadata
pData(gs)$condition <- c("treatment","control")
pData(gs)

# Subset a GatingSet by metadata column
# If you run this command you will keep only the sencond sample
#subset(gs, subset = condition == "control")

# Retrieve a GatingHierarchical (one sample)​
gh <- gs[[1]]
gh


# plot the gating hierarchy
plot(gs)


# list nodes (cell populations)
gs_get_pop_paths(gs, path = 2)
gs_get_pop_paths(gs, path = "full")
gs_get_pop_paths(gs, path = "auto")

# Retrieve data from all nodes as a cytoset
cs <- gs_pop_get_data(gs)
class(cs)

# Convert the cytoset to a flowSet
fs <- cytoset_to_flowSet(cs)

# Check the number of cells in the flowSet
fsApply(fs,nrow)

# retrieve data associated to one node (gate)
cs <- gs_pop_get_data(gs, "CD4")
fs <- cytoset_to_flowSet(cs)
fsApply(fs,nrow)

# get membership indices with respect to a gate
gh_pop_get_indices(gs[[1]], "CD4")
table(gh_pop_get_indices(gs[[1]], "CD4"))

# Delete a gate
gs_pop_remove(gs, "DPT")
plot(gs)

# save / load a GatingSet
save_gs(gs, path = "course_datasets/flowWorspaceData/gs")
gs2 <- load_gs("course_datasets/flowWorspaceData/gs")


###############################
# flowWorspace                #
# manual gating               #
###############################

rm(list = ls())

library(flowWorkspace)
library(flowCore)
library(CATALYST)
library(ggcyto)
library(flowVS)


# start from a flowSet
fs <- read.flowSet(path="course_datasets/flowWorspaceData/", pattern = "*.fcs")

# transform
panel <- pData(parameters(fs[[1]]))
markerstotransf <- as.character(panel$name[!is.na(panel$desc)])
fs <- transFlowVS(fs,
                  channels = markerstotransf,
                  cofactors = rep(3000,length(markerstotransf)))

# convert flowSet to GatingSet
gs <- GatingSet(fs) 


# ALTERNATIVE: start from a single cell experience object
# load("course_datasets/FR_FCM_Z4KT/DA_example_sce_PBMC.RData") # load the sce
# fs <- sce2fcs(sce, split_by = "sample_id",  assay = "exprs") # convert sce to flowSet (CATALYST)
# gs_PBMC <- GatingSet(fs) # convert flowSet to GatingSet


# Create rectange gate "NotDebris"
rg1  <- rectangleGate("FSC-A"=c(60000,260000), 
                      "SSC-A"=c(1, 250000), 
                      filterId="NotDebris")

# Add the gate to the GatingSet
gs_pop_add(gs, rg1, parent = "root")

# recompute the GatingSet
recompute(gs)

# plot the gate
autoplot(gs[[1]], gate = "NotDebris")

# plot the gating hierarchy
plot(gs)


# get statistics
gs_pop_get_stats(gs, "NotDebris") # counts
gs_pop_get_stats(gs, "NotDebris",type = "percent") # proportions


## Create a polygon gate

# Define the vertices of the polygon
my_vertices <-matrix(c(1,0.6,1,2,2.3,2.2,
                       25000,65000,120000,120000,75000,25000),
                     ncol=2,nrow=6)

colnames(my_vertices) <- c("V450-A","SSC-A")

# Create polygon gate "CD3"
rg2  <- polygonGate( boundaries= my_vertices,
                     filterId="CD3")

# Add the gate to the GatingSet
gs_pop_add(gs, rg2, parent = "NotDebris")

# Recompute the GatingSet
recompute(gs)

# Check
autoplot(gs[[1]], gate = "CD3")
plot(gs)


## Create a quadrant gate

# Create quadrant gate "CD4 CD8"
rg3  <- quadGate("B710-A"= 1.5, "R780-A"= 3, filterId = "CD4 CD8")

# Add the gate to the GatingSet
gs_pop_add(gs, rg3, parent = "CD3")

# Recompute the GatingSet
recompute(gs)

# Check
gs_get_pop_paths(gs)
autoplot(gs[[1]], gate = gs_get_pop_paths(gs)[4:7])
plot(gs)



# rename node names
gs_pop_set_name(gs,"B710-A-R780-A+","CD8+")
gs_pop_set_name(gs,"B710-A+R780-A-","CD4+")
gs_pop_set_name(gs,"B710-A-R780-A-","DNT")
gs_pop_set_name(gs,"B710-A+R780-A+","DPT")
plot(gs)

gs_pop_remove(gs, "DNT") 
plot(gs)

# retrieve the flow data for a node
fs_CD8 <- gs_pop_get_data( gs, "CD8+")


###############################
# flowGate                    #
###############################

rm(list = ls())

#BiocManager::install("flowGate")
library(flowGate)
library(flowWorkspace)
library(flowCore)
library(flowVS)

# read the flow data from FCS files
# start from a flowSet
fs <- read.flowSet(path="course_datasets/flowWorspaceData/", pattern = "*.fcs")

# transform
panel <- pData(parameters(fs[[1]]))
markerstotransf <- as.character(panel$name[!is.na(panel$desc)])
fs <- transFlowVS(fs,
                  channels = markerstotransf,
                  cofactors = rep(3000,length(markerstotransf)))

# convert flowSet to GatingSet
gs <- GatingSet(fs) 


# Create a rectangular gate
# Choose a type of gate (e.g. "rectangle")
# Draw a gate and click "done"
gs_gate_interactive(gs,
                    filterId = "NotDebris",
                    dims = list("FSC-A", "SSC-A"))

# Plot the data with the new gate
autoplot(gs[[1]], gate = "NotDebris")

# Plot hierarchy
plot(gs)


# Create a polygon gate for singlets
gs_gate_interactive(gs,
                    filterId = "singlets",
                    dims = list("FSC-A","FSC-H"),
                    subset = "NotDebris")


# check
autoplot(gs[[1]], gate = "singlets")
plot(gs)


# create a 1-D (span) gate
gs_gate_interactive(gs,
                    filterId = "CD3",
                    dims = "V450-A",
                    subset = "singlets")

# check
autoplot(gs[[1]],gate="CD3")
plot(gs)


gs_gate_interactive(gs,
                    filterId = "CD4 CD8",
                    dims = c("B710-A","R780-A"),
                    subset = "CD3")

# check
my_nodes <- gs_pop_get_children(gs, "CD3")
autoplot(gs[[1]],my_nodes )

plot(gs)



###########################
# openCyto                #
###########################

rm(list = ls())

library(flowCore)
library(flowWorkspace)
library(flowWorkspaceData)
library(ggcyto)
library(CytoML)
#BiocManager::install("openCyto")
library(openCyto)

# example of gating Hierarchy
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
wsfile <- list.files(flowDataPath, pattern="manual.xml",full = TRUE)
ws <- open_flowjo_xml("course_datasets/flowWorspaceData/manual.xml")
gs_final <- flowjo_to_gatingset(ws, name= "T-cell", subset =1)

# Check complete gating hierarchy​
plot(gs_final[[1]])

# plot gates
autoplot(gs_final[[1]])



# create the gatingTemplate from a csv file
# This is a gating scheme for a `T cell` panel, which tries to identify `T cell` 
# sub-populations
#gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
# I copied this file to /course_datasets/flowWorspaceData/
gt_tcell <- gatingTemplate("course_datasets/flowWorspaceData/tcell.csv")


# Visualize the gatingTemplate as a scheme
plot(gt_tcell)

# examine the structure of the gatingTemplate
# don't need to run this
# my_gt <- read.csv("course_datasets/flowWorkspaceData/tcell.csv")
# View(my_gt)


## Run the gating pipeline

# Load the preprocess but ungated data
# The code used to preprocess this data is avalable in
# /Code_slides/Code_preprocessing_data_openCyto.R
gs <- load_gs("course_datasets/flowWorspaceData/gs_preprocessed")

# check
plot(gs)

# Apply the gatingTemplate to the ungated data
gt_gating(gt_tcell, gs)


# check after gating
plot(gs[[1]])


# hide populations we are not interested in
nodesToHide <- c("cd8+", "cd4+"
                 , "cd4-cd8-", "cd4+cd8+"
                 , "cd4+cd8-/HLA+", "cd4+cd8-/CD38+"
                 , "cd4-cd8+/HLA+", "cd4-cd8+/CD38+"
                 , "CD45_neg/CCR7_gate", "cd4+cd8-/CD45_neg"
                 , "cd4-cd8+/CCR7+", "cd4-cd8+/CD45RA+"
)

# apply gs_pop_set_visibility to all the nodes
lapply(nodesToHide, function(thisNode) gs_pop_set_visibility(gs, thisNode, FALSE))


# rename some cell populations
gs_pop_set_name(gs, "cd4+cd8-", "cd4")
gs_pop_set_name(gs, "cd4-cd8+", "cd8")

# check renaming
plot(gs[[1]])

# visualize the gates
autoplot(gs[[1]])


# gating without template
gs_add_gating_method(gs,alias = "non-activated cd4",
                     pop = "--",parent = "cd4",
                     dims = "CD38,HLA",
                     gating_method = "tailgate")


plot(gs[[1]])

