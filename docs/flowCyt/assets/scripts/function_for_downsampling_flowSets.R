#######################################
# Function for downsampling a flowSet #
#######################################

# create a function for downsampling the data
Downsampling_flowSet <- function(x, samplesize, replace=TRUE, prob=NULL) {
  
  
  if(missing(samplesize)) samplesize <- min(flowCore::fsApply(x,nrow))
  
  flowCore::fsApply(x, function(ff) {
    
    i <- sample(nrow(ff), size = samplesize, replace=replace, prob)
    ff[i,]
  })    
  
}