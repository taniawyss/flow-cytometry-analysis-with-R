
In this section, you will find the R code that we will use during the course. We will explain the code and output during correction of the exercises.

<!-- This is commented text -->

## Example of workflow

Here we present an R markdown file that allows to perform a full flow cytometry data analysis workflow, including clustering, dimensional reduction and trajectory analysis.

The source and inspiration of this workflow was published by [Melsen et al, 2020](https://journals.aai.org/jimmunol/article/205/3/864/60873/A-Comprehensive-Workflow-for-Applying-Single-Cell). In the paper, Melsen et al also used some tools that were not implemented in R. In our workflow, all steps are done within R.

We need some packages that we have not used before, so they need to be installed.

```r

#########################################
# Libraries for flow cytometry analysis #
#########################################

BiocManager::install("flowStats")
BiocManager::install("destiny")
BiocManager::install("slingshot")
devtools::install_github("JinmiaoChenLab/cytofkit2", dependencies=TRUE)

```

The data are available for download at the bottom of the [material](https://taniawyss.github.io/flow-cytometry-analysis-with-R/flowCyt/material/#day-5) page. 

[Download solution Rmd file](../../assets/scripts/Workflow_clustering_pseudotime.Rmd){: .md-button }


## Feedback










