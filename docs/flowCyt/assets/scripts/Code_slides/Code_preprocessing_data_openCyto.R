###############################################
# OpenCyto                                    #
# Code for preprocessing flow cytometry data  #
###############################################

rm(list = ls())

library(CytoML)
library(openCyto)
library(data.table)


# get the data for the compensation
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
wsfile <- list.files(flowDataPath, pattern="manual.xml",full = TRUE)
ws <- open_flowjo_xml("course_datasets/flowWorspaceData/manual.xml")
gs_raw <- flowjo_to_gatingset(ws, name= "T-cell", subset =1)
gh <- gs_raw[[1]]

# Load the raw data
fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
cs  <- load_cytoset_from_fcs(fcsFiles)
cf <- realize_view(cs[[1]])
gs <- GatingSet(cs)

# compensate
compMat <- gh_get_compensations(gh)
compensate(gs, compMat)

# transform
chnls <- parameters(compMat)
trans <- estimateLogicle(gs[[1]], channels = chnls)
gs <- transform(gs, trans)

# save
save_gs(gs,"course_datasets/flowWorspaceData/gs_preprocessed")
