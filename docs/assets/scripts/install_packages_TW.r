#!/usr/bin/env Rscript

# R script to install requirements for exercises ------


#################################
# Non-CRAN package installation #
#################################

install.packages("devtools")
install.packages("BiocManager")


###############
# R markdown  #
###############

install.packages("knitr")
install.packages("rmarkdown")
install.packages("pander")


#####################
# Excel files       #
#####################

install.packages("readxl")
install.packages("xlsx")
install.packages("WriteXLS")


#################
# Visualization #
#################

install.packages("ggplot2")
BiocManager::install("ggcyto")
install.packages("manipulate")
install.packages("ggrepel")
install.packages("ggpubr") # done
install.packages("RColorBrewer")
install.packages("gridExtra")
install.packages("cowplot")
install.packages("ggsignif")
BiocManager::install("plotly") # done
install.packages("ggridges")
install.packages("scales")
BiocManager::install("ComplexHeatmap")
# done

#################
# Data handling #
#################

install.packages("reshape2")
install.packages("matrixStats")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("tibble")


###############################################
# Statistical functions                       #
###############################################

install.packages("lme4")
install.packages("multcomp")
install.packages("rstatix")
install.packages("DescTools")
install.packages("statmod") # done
BiocManager::install("edgeR") # done
install.packages("MASS") # done
BiocManager::install("diffcyt")
install.packages('sfsmisc')
install.packages("rms")


#########################################
# Libraries for flow cytometry analysis #
#########################################

BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("flowDensity")
BiocManager::install("MetaCyto")
BiocManager::install("scDataviz")
BiocManager::install("flowViz")
BiocManager::install("flowVS")
BiocManager::install("flowAI")
BiocManager::install("PeacoQC")
BiocManager::install("flowClean")
BiocManager::install("CATALYST")
devtools::install_github('saeyslab/CytoNorm')
BiocManager::install("SingleCellExperiment")
install.packages("uwot")
BiocManager::install("FlowSOM")
BiocManager::install("ConsensusClusterPlus")
install.packages("Rtsne")
BiocManager::install("scater")
devtools::install_github("RGLab/scamp")
devtools::install_github("RGLab/FAUST")


############################
# Survival analysis        #
############################

install.packages("survival")
install.packages("survminer")


