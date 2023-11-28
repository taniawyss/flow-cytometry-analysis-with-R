###################
# FAUST           #
###################

rm(list = ls())

# these packages must be installed
library(knitr)
library(rmarkdown)
#install.packages("ggdendro")
library(ggdendro)
library(remotes)
library(lme4)
library(multcomp)
library(ggplot2)

# install FAUST from github
# install vignettes
#remotes::install_github("RGLab/FAUST", force = TRUE, build_vignettes = TRUE)
library(faust)


# Available Vignettes
#vignette('faustIntro')
#vignette('faustTuning')
#vignette('faustPFDA')


# An example workflow with simulated data

# load a simulated flowSet
# A simple dataset consisting of 15 samples with three markers, 
# which we generically call V1, V2, and V3.
load("course_datasets/FAUST/Simulated_fs.RData")

# Check
pData(fs)
# name
# sample001 sample001
# sample002 sample002
# sample003 sample003
# sample004 sample004
# sample005 sample005
# sample006 sample006
# sample007 sample007
# sample008 sample008
# sample009 sample009
# sample010 sample010
# sample011 sample011
# sample012 sample012
# sample013 sample013
# sample014 sample014
# sample015 sample015

pData(parameters(fs[[1]]))
# name desc range  minRange maxRange
# $P1   V1   V1    14 -4.889966       13
# $P2   V2   V2    12 -3.195790       11
# $P3   V3   V3    13 -2.797357       12
# $P4   V4   V4    13 -3.260996       12
# $P5   V5   V5    14 -3.620577       13

# convert to a GatingSet
gs <- flowWorkspace::GatingSet(fs)


# Annotation thresholds - depth score tuning

set.seed(123)
generateAnnotationThresholds(
  gatingSet           = gs,
  startingCellPop     = "root",
  projectPath         = "course_datasets/FAUST/",
  depthScoreThreshold = 0.85,
  selectionQuantile   = 0.5,
  plottingDevice      = "png",
  threadNum           = 4,
  activeChannels      = c("V1","V4","V5")
)


# Discovering phenotypes
discoverPhenotypes(
  gatingSet   = gs,
  projectPath = "course_datasets/FAUST/",
  threadNum   = 4)

# Using the faust function in R
?faust

faust(gatingSet = gs,
      startingCellPop = "root",
      depthScoreThreshold = 0.85,
      selectionQuantile   = 0.5,
      projectPath = "course_datasets/FAUST/",
      activeChannels      = c("V1","V4","V5"),
      annotationsApproved = TRUE,
      threadNum = 4)


# Examine output
count_df <- as.data.frame(readRDS("course_datasets/FAUST/faustData/faustCountMatrix.rds"))
count_df

# names of samples 
snVec <- pData(gs)$name
  
  
# produces an annotation embedding for samples analyzed by FAUST
annoEmbed <- faust::makeAnnotationEmbedding(          
  projectPath="course_datasets/FAUST/",
  sampleNameVec=snVec
)

# plot the UMAP
ggplot(annoEmbed[annoEmbed$sampleOfOrigin=="sample001",],aes(x=umapX,y=umapY,color=faustLabels))+
  geom_point()+
  theme_classic()+
  xlab("Annotation Embedding X")+
  ylab("Annotation Embedding Y")+
  facet_grid(vars(sampleOfOrigin))+
  theme(legend.position="bottom")




