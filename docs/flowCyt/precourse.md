
## R and RStudio

### Previous knowledge / Competencies

We expect participants to have previous knowledge in:

* R beginner level (Rstudio, install a library, data frame manipulation, import data from csv file). An introduction to R is available [here](https://taniawyss.github.io/flow-cytometry-analysis-with-R/introR/material/).

### Technical

This course will be streamed, you are thus required to have your own computer with an internet connection, and with latest the version of [R](https://cran.r-project.org/)
and the free version of [RStudio](https://www.rstudio.com/products/rstudio/download/) installed. Admin rights may be needed to install the necessary packages.

The packages we will need are hosted on [CRAN](https://cran.r-project.org/web/packages/available_packages_by_name.html), [Bioconductor](https://bioconductor.org/) and [Github](https://github.com/). You can install the necessary packages using:

```r
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
install.packages("ggpubr")
install.packages("RColorBrewer")
install.packages("gridExtra")
install.packages("cowplot")
install.packages("ggsignif")
BiocManager::install("plotly")
install.packages("ggridges")
install.packages("scales")
BiocManager::install("ComplexHeatmap")
install.packages("circlize")
install.packages("cowplot")

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
install.packages("statmod")
BiocManager::install("edgeR")
install.packages("MASS")
BiocManager::install("diffcyt")
install.packages('sfsmisc')
install.packages("rms")


#########################################
# Libraries for flow cytometry analysis #
#########################################

BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("flowWorkspaceData")
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

############################
# Gating                   #
############################

BiocManager::install("flowClust")
BiocManager::install("CytoML")
BiocManager::install("openCyto")

```

After installation, packages can be loaded using:

```r

#################################
# Non-CRAN package installation #
#################################

# install packages
library(devtools) # install.packages("devtools")
library(BiocManager) # install.packages("BiocManager")


###############
# R markdown  #
###############

library(knitr) # install.packages("knitr")
library(rmarkdown) # install.packages("rmarkdown")
library(pander) # install.packages("pander")


#####################
# Excel files       #
#####################

library(readxl) # install.packages("readxl")
library(xlsx) # install.packages("xlsx")
library(WriteXLS) # install.packages("WriteXLS")


#################
# Visualization #
#################

library(ggplot2) # install.packages("ggplot2")
library(ggcyto) # BiocManager::install("ggcyto")
library(manipulate) # install.packages("manipulate")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(RColorBrewer) # install.packages("RColorBrewer")
library(gridExtra) # install.packages("gridExtra")
library(cowplot) # install.packages("cowplot")
library(ggsignif) # install.packages("ggsignif")
library(plotly) # BiocManager::install("plotly")
library(ggridges) # install.packages("ggridges")
library(scales) # install.packages("scales")
library(ComplexHeatmap) # BiocManager::install("ComplexHeatmap")
library(circlize) # install.packages("circlize")
library(cowplot) # install.packages("cowplot")


#################
# Data handling #
#################

library(reshape2) # install.packages("reshape2")
library(matrixStats) # install.packages("matrixStats")
library(tidyverse) # install.packages("tidyverse")
library(dplyr) # install.packages("dplyr")
library(tibble) # install.packages("tibble")


###############################################
# Statistical functions                       #
###############################################

library(lme4) # install.packages("lme4")
library(multcomp) # install.packages("multcomp")
library(rstatix) # install.packages("rstatix")
library(DescTools) # install.packages("DescTools")
library(statmod) # install.packages("statmod")
library(edgeR) # BiocManager::install("edgeR")
library(MASS) # install.packages("MASS")
library(diffcyt) # BiocManager::install("diffcyt")
library(sfsmisc) # install.packages('sfsmisc')
library(rms) # install.packages("rms")

#########################################
# Libraries for flow cytometry analysis #
#########################################

library(flowCore) # BiocManager::install("flowCore")
library(flowWorkspace) # BiocManager::install("flowWorkspace")
library(flowWorkspaceData) # BiocManager::install("flowWorkspaceData")
library(flowDensity) # BiocManager::install("flowDensity")
library(MetaCyto) # BiocManager::install("MetaCyto")
library(scDataviz) # BiocManager::install("scDataviz")
library(flowViz) # BiocManager::install("flowViz")
library(flowVS) # BiocManager::install("flowVS")
library(flowAI) # BiocManager::install("flowAI")
library(PeacoQC) # BiocManager::install("PeacoQC")
library("flowClean") # BiocManager::install("flowClean")
library(CATALYST) # BiocManager::install("CATALYST")
library(CytoNorm) # install_github('saeyslab/CytoNorm')
library(SingleCellExperiment) # BiocManager::install("SingleCellExperiment")
library(uwot)  # install.packages("uwot")
library(FlowSOM) # BiocManager::install("FlowSOM")
library(ConsensusClusterPlus) # BiocManager::install("ConsensusClusterPlus")
library(Rtsne) # install.packages("Rtsne")
library(scater) # BiocManager::install("scater")
library(scamp) # devtools::install_github("RGLab/scamp")
library(faust) # devtools::install_github("RGLab/FAUST")


############################
# Survival analysis        #
############################

library(survival) # install.packages("survival")
library(survminer) # install.packages("survminer")

############################
# Gating                   #
############################

library(flowClust) # BiocManager::install("flowClust")
library(CytoML) # BiocManager::install("CytoML")
library(openCyto) # BiocManager::install("openCyto")


```



