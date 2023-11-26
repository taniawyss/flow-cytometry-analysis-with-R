#######################################
# Let’s practice – 9                  #
# Manual gating with flowGate.        #
# Solutions                           #
#######################################


# clear the environment
rm(list = ls())


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
