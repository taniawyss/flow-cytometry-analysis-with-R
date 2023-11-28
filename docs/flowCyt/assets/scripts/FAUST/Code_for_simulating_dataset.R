#####################################
# Code for simulating a dataset     #
#####################################

rm(list = ls())

library(flowCore)
library(faust)

# functions (copied from the Vignette)

Posdef <- function (n, ev = runif(n, 1, 2))
{
  #function written by Ravi Varadhan, from r-help mailing list
  #Thu Feb 7, 20:02:30 CET 2008
  #Generates a positive definine matrix of dimension n
  #ev bounds the variance by 2
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

labelMeanVector <- function(meanVector) {
  baseStr <- paste0("V",seq(length(meanVector)))
  tokenStr <- rep("+",length(meanVector))
  tokenStr[which(meanVector==0)] <- "-"
  return(paste0(paste0(baseStr,tokenStr),collapse=""))
}

genClusterCentersWithLabels <- function(possibleCenterMat,
                                        nPop=5,
                                        seedVal=0)
{
  pcMat <- t(possibleCenterMat)
  allPops <- expand.grid(split(pcMat,rep(seq(nrow(pcMat)),ncol(pcMat))))
  if (nPop > nrow(allPops)) {
    print("Too many clusters relative to specMat/possibleCenterMat.")
    stop("Reduce number of clusters or increase possible mean vectors.")
  }
  clusterIndex <- sample(seq(nrow(allPops)),nPop)
  outMat <- allPops[clusterIndex,,drop=FALSE]
  outLab <- as.character(apply(outMat,1,labelMeanVector))
  outList <- list(fixedMeanMatrix=outMat,fixedLabelVector=outLab)
  return(outList)
}

getSimTrueCountMatrix <- function(truthList)
{
  uniqueNames <- c()
  for (i in seq(length(truthList))) {
    uniqueNames <- append(uniqueNames,names(truthList[[i]]))
  }
  uniqueNames <- sort(unique(uniqueNames))
  truthMat <- matrix(0,nrow=length(truthList),ncol=length(uniqueNames))
  colnames(truthMat) <- uniqueNames
  rownames(truthMat) <- names(truthList)
  for (i in seq(length(truthList))) {
    cv <- truthList[[i]]
    for (cName in colnames(truthMat)) {
      lookup <- which(names(cv) == cName)
      if (length(lookup)) {
        truthMat[i,cName] <- cv[lookup]
      }
    }
  }
  return(truthMat)
}

simSample <- function(sampleDim,
                      sampleSize,
                      transformationType,
                      mixtureType="gaussianOnly",
                      fixedMeanMatrix=NA,
                      fixedLabelVector=NA,
                      noiseDim=0,
                      probVecSample,
                      isKnockout=FALSE,
                      isSpikeIn=FALSE,
                      sRegime=0,
                      targetRanks=c(length(probVecSample)-1,length(probVecSample)),
                      hasBatchEffect=FALSE,
                      batchEffect=0
)
{
  numClusters <- nrow(fixedMeanMatrix)
  probVec <- probVecSample
  knockoutStatus <- "No knockout"
  if (isKnockout) {
    targetProbs <- sort(probVec)
    targetProb <- targetProbs[ceiling(length(targetProbs)/2)]
    targetLookup <- which(probVec == targetProb)
    modProbLookups <- setdiff(seq(length(probVec)),targetLookup)
    knockoutStatus <- fixedLabelVector[targetLookup]
    probVec[modProbLookups] <- (probVec[modProbLookups] + (targetProb/length(modProbLookups)))
    probVec[targetLookup] <- 0
  }
  if ((sRegime) && (isSpikeIn)) {
    currentSpikeMass <- probVec[targetRanks[2]]
    targetSpikeMass <- probVec[targetRanks[1]]
    massIncrement <- currentSpikeMass/(length(probVec)-1)
    probVec <- probVec + massIncrement
    probVec[targetRanks[2]] <- 0
    probVec <- probVec - (targetSpikeMass * probVec)
    probVec[targetRanks[2]] <- targetSpikeMass
    if (abs(sum(probVec) - 1) > 1e-9) {
      print(probVec)
      print(sum(probVec))
      print(abs(sum(probVec) - 1))
      stop("Error in probability reapportionment")
    }
    if (min(probVec) == probVec[targetRanks[2]]) {
      stop("Error in spiking in population")
    }
  }
  sampleSizeVec <- as.vector(t(rmultinom(1,sampleSize,probVec)))
  outData <- matrix(nrow=0,ncol=(sampleDim+noiseDim))
  outLabel <- c()
  for (clusterNumber in seq(numClusters)) {
    if (sampleSizeVec[clusterNumber] == 0) {
      next
    }
    currentMu <- as.numeric(fixedMeanMatrix[clusterNumber,,drop=TRUE])
    if (hasBatchEffect) {
      currentMu <- currentMu + batchEffect 
    }
    subjectShift <- round(rnorm(length(currentMu),mean=0,sd=(1/sqrt(2))))
    currentMu <- currentMu + subjectShift 
    currentLabel <- fixedLabelVector[clusterNumber]
    outLabel <- append(outLabel,rep(currentLabel,sampleSizeVec[clusterNumber]))
    currentSigma<- Posdef(sampleDim)
    if (mixtureType == "tPlusGauss") {
      if ((clusterNumber %% 2) == 0) {
        currentSample <- mvrnorm(sampleSizeVec[clusterNumber],mu=currentMu,Sigma=currentSigma)
      }
      else {
        currentSample <- rmvt(sampleSizeVec[clusterNumber],delta=currentMu,sigma=currentSigma,df=2)
      }
    }
    else if (mixtureType == "tOnly") {
      currentSample <- rmvt(sampleSizeVec[clusterNumber],delta=currentMu,sigma=currentSigma,df=2)
    }
    else  {
      currentSample <- mvrnorm(sampleSizeVec[clusterNumber],mu=currentMu,Sigma=currentSigma)
    }
    if (is.vector(currentSample)) {
      currentSample <- t(as.matrix(currentSample))
    }
    if (noiseDim > 0) {
      currentSigma<- Posdef(noiseDim)
      noiseSample <- mvrnorm(nrow(currentSample),mu=rep(5,noiseDim),Sigma=currentSigma)
      if (is.vector(noiseSample)) {
        noiseSample <- t(as.matrix(noiseSample))
      }
      currentSample <- cbind(currentSample,noiseSample)
    }
    outData <- rbind(outData,t(apply(currentSample,1,transformationType)))
  }
  colnames(outData) <- paste0("V",seq(ncol(outData)))
  outList <- list(sampleMatrix=outData,sampleLabels=outLabel,knockoutInfo=knockoutStatus)
  return(outList)
}

simulateExperiment <- function(
    meanVectorBoundsMatrix, 
    numSamples=100, 
    randomSeed=0,
    transformationList=list(function(x){return(x)},
                            function(x){return(x)},
                            function(x){return(x)}), 
    noiseDimension=5,
    probVecIn, 
    minSampleSize=5000,
    maxSampleSize=40000,
    tncp=10000,
    knockoutNumber=0, 
    spikeInFlag=FALSE, 
    targetRanks=c(length(probVecIn)-1,length(probVecIn)),
    useBatchFlag=FALSE, 
    batchEffectShift=0.75,
    fixedResFlag=FALSE, 
    responderStatusVec=c(NA)
)
{
  probVecForSim <- sort(probVecIn,decreasing=TRUE)
  sampleSpecs <- genClusterCentersWithLabels(
    possibleCenterMat=meanVectorBoundsMatrix,
    nPop=length(probVecForSim),
    seedVal=randomSeed
  )
  currentTransformation <- transformationList[[1]]
  currentSample <- 1
  startKnockoutNum <- numSamples - knockoutNumber + 1
  isKnockoutSample <- FALSE
  regime <- rep(0,numSamples)
  regime[sample(seq(numSamples),(numSamples/2))] <- 1
  if (fixedResFlag) {
    regime <- rep(0,numSamples)
    regime[which(responderStatusVec==1)] <- 1
  }
  nextBatchEffect <- rep((-1*batchEffectShift),ncol(sampleSpecs$fixedMeanMatrix))
  regimeList <- knockoutList <- labelsList <- flowList <- truthList <- list()
  while (currentSample <= numSamples) {
    nextSampleSize <- min(max(minSampleSize,round(rt(1,df=3,ncp=tncp))),maxSampleSize)
    nextRegime <- regime[currentSample]
    if (currentSample >= startKnockoutNum) {
      isKnockoutSample <- TRUE
    }
    else {
      isKnockoutSample <- FALSE
    }
    if ((useBatchFlag) && (((currentSample - 1) %% 10) == 0)) {
      nextBatchEffect <- nextBatchEffect + batchEffectShift
      if ((currentSample - 1) > floor(numSamples/3)) {
        currentTransformation <- transformationList[[2]]
      }
      if ((currentSample - 1) > floor((2*(numSamples/3)))) {
        currentTransformation <- transformationList[[3]]
      }
    }
    
    sampleData <- simSample(
      sampleDim=ncol(sampleSpecs$fixedMeanMatrix),
      sampleSize=nextSampleSize,
      transformationType=currentTransformation,
      fixedMeanMatrix=sampleSpecs$fixedMeanMatrix,
      fixedLabelVector=sampleSpecs$fixedLabelVector,
      noiseDim=noiseDimension,
      probVecSample=probVecForSim,
      isKnockout=isKnockoutSample,
      isSpikeIn=spikeInFlag,
      sRegime=nextRegime,
      targetRanks=targetRanks,
      hasBatchEffect=useBatchFlag,
      batchEffect=nextBatchEffect
    )
    if (currentSample < 10) {
      outName <- paste0("sample00",currentSample)
    }
    else if (currentSample < 100) {
      outName <- paste0("sample0",currentSample)
    }
    else {
      outName <- paste0("sample",currentSample)
    }
    ff <- flowCore::flowFrame(sampleData$sampleMatrix)
    flowList <- append(flowList,ff)
    names(flowList)[length(flowList)] <- outName
    labelsList <- append(labelsList,list(sampleData$sampleLabels))
    names(labelsList)[length(labelsList)] <- outName
    truthList <- append(truthList,list(table(sampleData$sampleLabels)))
    names(truthList)[length(truthList)] <- outName
    knockoutList <- append(knockoutList,list(table(sampleData$knockoutInfo)))
    names(knockoutList)[length(knockoutList)] <- outName
    regimeList <- append(regimeList,list(nextRegime))
    names(regimeList)[length(regimeList)] <- outName
    currentSample <- currentSample + 1
  }
  truthMat <- getSimTrueCountMatrix(truthList)
  
  outputList <- list(
    "flowFrameList"=flowList,
    "truthValueList"=truthList,
    "truthMat"=truthMat,
    "labelsList"=labelsList,
    "knockoutList"=knockoutList,
    "spikedPop"=sampleSpecs$fixedLabelVector[targetRanks[2]],
    "spikedInPopList"=regimeList
  )
  return(outputList)
}
dimension <- 5 
sampleSize <- 500 
ssLB <- 500 
ssUB <- 500
specMat <- matrix(0,nrow=2,ncol=dimension)
specMat[2,] <- c(8,8,8,8,8)
dataTransformations <- list(function(x){return((x))},
                            function(x){return((x))},
                            function(x){return((x))})
startProb <- 0.1425
refClusterProbVec <- c(
  startProb,rep(startProb/2,2),rep(startProb/4,4),rep(startProb/8,8),rep(startProb/16,16),rep(startProb/32,32),rep(startProb/64,64)
)
numberOfClusters <- 5
clusterProbVec <- sort(refClusterProbVec[seq(numberOfClusters)],decreasing=TRUE)
if (sum(clusterProbVec) < 1) {
  residualPV <- 1-sum(clusterProbVec)
  clusterProbVec <- (clusterProbVec + (residualPV/length(clusterProbVec)))
}


# simulated data

simmedExp <- simulateExperiment(
  meanVectorBoundsMatrix=specMat,
  numSamples=15, # nr of samples
  noiseDimension = 0, # nr of markers
  transformationList=dataTransformations,
  probVecIn=clusterProbVec,
  minSampleSize=ssLB,
  maxSampleSize=ssUB,
  randomSeed=12345,
  tncp=sampleSize
)
fs <- as(simmedExp[["flowFrameList"]], "flowSet")
save(fs,file =  "course_datasets/FAUST/Simulated_fs.RData")
