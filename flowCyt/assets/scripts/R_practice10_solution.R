#######################################
# Let’s practice – 10                 #
# Gating with openCyto                #
# Solutions                           #
#######################################


# clear the environment
rm(list = ls())


# load libraries
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)


# 1) Repeat steps 1 to 4 from previous exercise (loading data, preprocessing)

# load the preprocessed ungated data
gs <- load_gs("course_datasets/FR_FCM_Z3WR/gs_preprocessed")

# 2) Apply the gating as in the scheme

# gate "Leukocytes"
# use flowClust
# set a K=3, to split in three populations and select the one with the highest "peak"
gs_add_gating_method(gs,
                     alias = "Leukocytes",
                     parent = "root",
                     dims = "FSC-H,SSC-H",
                     gating_method = "flowClust",
                     gating_args = "K=3" )

# check
plot(gs)

autoplot(gs[[1]], gate="Leukocytes")


# gate "CD3"
# use minDensity
gs_add_gating_method(gs,
                     alias = "CD3",
                     parent = "Leukocytes",
                     dims = "BV510-A",
                     gating_method = "gate_mindensity",
                     pop = "+")

# check
plot(gs)

autoplot(gs[[1]], gate="CD3")

# gate "CD4 CD8"
# use minDensity
gs_add_gating_method(gs,
                     alias = "*",
                     parent = "CD3",
                     dims = "BUV615-A,BUV805-A",
                     gating_method = "gate_mindensity",
                     pop = "-/++/-")

# check
plot(gs)

gs_pop_get_children(gs,"CD3")
autoplot(gs[[1]], gate=gs_pop_get_children(gs,"CD3"))
autoplot(gs[[1]], gate=gs_pop_get_children(gs,"CD3")[3:6])


# 3) 

# Hide the "BUV805-A+" and "BUV615-A+" nodes from the tree

gs_pop_set_visibility(gs,"BUV805-A+", FALSE)
gs_pop_set_visibility(gs,"BUV615-A+", FALSE)

# check
plot(gs)

# 4) Rename the CD4 CD8 gates to CD4+, CD8+, DNT and DPT

# Rename the CD4 CD8 gates to CD4+, CD8+, DNT and DPT
gs_pop_set_name(gs,"BUV615-A-BUV805-A+","CD8+")
gs_pop_set_name(gs,"BUV615-A+BUV805-A-","CD4+")
gs_pop_set_name(gs,"BUV615-A+BUV805-A+","DPT")
gs_pop_set_name(gs,"BUV615-A-BUV805-A-","DNT")

# check
plot(gs)

# 5) Create a boxplot showing the percentage of CD8+ T cells among T cells ("CD3") as a function of time point

# extract the proportions
my_proportions <- gs_pop_get_stats(gs, nodes = "CD8+", type = "percent")

# add the timepoints
my_proportions$time_point <- rep(c(0,14,7), 5)

# convert to factor (order in the plot)
my_proportions$time_point <- factor(my_proportions$time_point, levels = c(0,7,14))

# create the boxplot
boxplot(percent ~ time_point, data = my_proportions)
