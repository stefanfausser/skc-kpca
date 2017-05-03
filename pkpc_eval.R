## Used solely in pkpc_experiments_eval_basic.R

simpleCap <- function(x)
{
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}

getNumberClusterSeq <- function(nClasses)
{
    s <- c(nClasses - 4, nClasses - 2, nClasses, nClasses + 2, nClasses + 4)
    s <- c(s, round(nClasses / 4), round(nClasses / 3), round(nClasses / 2), round(nClasses), round(nClasses * 2), round(nClasses * 3))
    
    # remove duplicated number of clusters
    s <- sort(s[!duplicated(s)])
    
    # min 2 clusters
    s <- s[s >= 2]
    
    return(s)
}

pkpc.eval.func <- function(seedSamples, seedLabels, seed, datasetname, datasetfunc, M, labelFactorSeq = 0, a_K = 2, degree = 1, offset = 0, sigma = 0, scaleFeatures = 1, normalizeKernel = 0, removeOutliers = 1, removeFeatures = 1, a_method = 5, fuzzifier = 1.25, R = 35, kFoldCrossValidation = 2, normalizeAssignmentsKFCM = 0, maxDataChunksSeq = 0, weightFactorSeq = 0)
{
    ret <- datasetfunc(removeOutliers = removeOutliers, removeFeatures = removeFeatures)
    
    N <- ret$N
    x <- ret$x
    labels <- ret$labels
    
    if(sigma > 0)
    {
        # Gaussian Kernel
        kernelfunc <- kernel.gauss
        kernelparam <- list(sigma=sigma)
    }
    else if(degree > 0)
    {
        # Polynomial Kernel (for degree = 1: Linear Kernel)
        kernelfunc <- kernel.poly
        kernelparam <- list(offset=offset, degree=degree, normalizeKernel=normalizeKernel)
    }
    
    dirname <- "labelsTmp/"
    if(!dir.exists(dirname))
    {
        dir.create(dirname)
    }
    
    savefile <- paste(dirname, datasetname, "eval", sep="")
    
    a_mu <- c(1)
    a_M <- M
    a_L <- c(0)
    
    pkpca_params <-  list(verbose=0, maxRepeats=5, maxRowsKmatOut=1000, maxDataChunks=0, allSamplesFromLast=0, maxSamplesValidate=2000, getHardClusteringsKFCM=0, normalizeAssignmentsKFCM=normalizeAssignmentsKFCM, limitSamplesByHardAssignments=0, localMinimum=0, softInitialAssignmentsApproxPseudoCentresRNG=0, hardAssignmentsApproxPseudoCentresRNG=0, excludeSamplesFromAPCInGreedySelection=0)
    
    a_pkpcparam <- rbind(pkpca_params, NULL)
    
    fileToSave <- paste("results-eval-", datasetname, "_M", M[1], "_method", a_method, "_normalizeAssignmentsKFCM", normalizeAssignmentsKFCM, "_fuzzifier", fuzzifier, "_degree", degree, "_offset", offset, "_sigma", sigma, "_normalize", scaleFeatures, "_normalizeKernel", normalizeKernel, sep="")
    # prefix results-eval/ directory
    path <- "results-eval/"
    if(!dir.exists(path))
    {
        dir.create(path)
    }
    fileToSave <- paste(path, fileToSave, sep="")
    
    epsilon <- 0.000001
    T <- 100 # number of maximum iterations per data chunk
    
    labelledUnlabelledFactorSeq <- 1
    
    printf("Create samplesMat\n")
    
    samplesMat <- saveOrRestoreSampleMat(seedSamples, R, N)
    
    printf("Randomly zero-out labels\n")
    
    ## only performed once
    createAndSaveRandomZeroLabels(seedLabels, savefile, labels, labelFactorSeq, R)
    
    printf("Start pkpc.main.test\n")
    
    # start the kernel clustering without labels
    pkpc.main.test(seed, fileToSave, savefile, a_pkpcparam, a_mu, a_M, a_K, a_L, a_method, x, kernelfunc, kernelparam, labelledUnlabelledFactorSeq, labelFactorSeq, samplesMat, kFoldCrossValidation, T, N, labels, fuzzifier, epsilon, scaleFeatures, maxDataChunksSeq, weightFactorSeq)
}

pkpc.start.eval <- function(seedSamples, seedLabels, seed, datasetname, datasetfunc, M, method, a_K, a_sigma, scaleFeatures = 1, linKernel = 1, polyKernel = 1, gaussKernel = 1, ...)
{
    if(linKernel)
    {
        # Linear kernel
        
        pkpc.eval.func(seedSamples, seedLabels, seed, datasetname = datasetname, datasetfunc = datasetfunc, M = M, a_K = a_K, degree = 1, offset = 0, sigma = 0, scaleFeatures = scaleFeatures, normalizeKernel = 0, a_method = method, ...)
    }
    
    if(polyKernel)
    {
        pkpc.eval.func(seedSamples, seedLabels, seed, datasetname = datasetname, datasetfunc = datasetfunc, M = M, a_K = a_K, degree = 2, offset = 0, sigma = 0, scaleFeatures = scaleFeatures, normalizeKernel = 1, a_method = method, ...)
        pkpc.eval.func(seedSamples, seedLabels, seed, datasetname = datasetname, datasetfunc = datasetfunc, M = M, a_K = a_K, degree = 3, offset = 0, sigma = 0, scaleFeatures = scaleFeatures, normalizeKernel = 1, a_method = method, ...)
    }
    
    if(gaussKernel)
    {
        # Gaussian Kernel
        for(sigma in a_sigma)
        {
            pkpc.eval.func(seedSamples, seedLabels, seed, datasetname = datasetname, datasetfunc = datasetfunc, M = M, a_K = a_K, degree = 0, offset = 0, sigma = sigma, scaleFeatures = scaleFeatures, normalizeKernel = 0, a_method = method, ...)
        }
    }
}

pkpc.start.eval.evaluation <- function(datasetname, a_sigma, legendLocationARI = "topleft", legendLocationARI2 = "topleft")
{
    # Linear
    
    path <- "results-eval/"
    
    file <- list.files(pattern=paste("results-eval-", datasetname, ".*method5", ".*degree1", sep=""), path=path)
    file <- paste(path, file, sep="")
    
    if(!file.exists(file))
        return(NULL)
    
    resLinear <- read.table(file)    
    
    # Polynomial kernel, degree 2    
    
    file <- list.files(pattern=paste("results-eval-", datasetname, ".*method5", ".*degree2", ".*normalizeKernel1", sep=""), path=path)
    file <- paste(path, file, sep="")
    
    if(!file.exists(file))
        return(NULL)
    
    resPolynomialNormalized2 <- read.table(file)
    
    # Polynomial kernel, degree 3
    
    file <- list.files(pattern=paste("results-eval-", datasetname, ".*method5", ".*degree3", ".*normalizeKernel1", sep=""), path=path)
    file <- paste(path, file, sep="")
    
    if(!file.exists(file))
        return(NULL)
    
    resPolynomialNormalized3 <- read.table(file)
    
    # Gaussian kernel
    
    file <- list.files(pattern=paste("results-eval-", datasetname, ".*method5", ".*degree0", ".*sigma", a_sigma[1], sep=""), path=path)
    file <- paste(path, file, sep="")
    
    if(!file.exists(file))
        return(NULL)
    
    resGaussian1 <- read.table(file)
    
    file <- list.files(pattern=paste("results-eval-", datasetname, ".*method5", ".*degree0", ".*sigma", a_sigma[2], sep=""), path=path)
    file <- paste(path, file, sep="")
    
    if(!file.exists(file))
        return(NULL)
    
    resGaussian2 <- read.table(file)
    
    file <- list.files(pattern=paste("results-eval-", datasetname, ".*method5", ".*degree0", ".*sigma", a_sigma[3], sep=""), path=path)
    file <- paste(path, file, sep="")
    
    if(!file.exists(file))
        return(NULL)
    
    resGaussian3 <- read.table(file)
    
    xLabelSeq <- c(resLinear$K, resPolynomialNormalized2$K)    
    xLabelSeq <- xLabelSeq[!duplicated(xLabelSeq)]
    
    res <- matrix(0, 6, length(xLabelSeq))
    
    res[1,1:length(resLinear$predAdjustedRandPseudoCentresMean)] <- resLinear$predAdjustedRandPseudoCentresMean
    res[2,1:length(resPolynomialNormalized2$predAdjustedRandPseudoCentresMean)] <- resPolynomialNormalized2$predAdjustedRandPseudoCentresMean
    res[3,1:length(resPolynomialNormalized3$predAdjustedRandPseudoCentresMean)] <- resPolynomialNormalized3$predAdjustedRandPseudoCentresMean
    res[4,1:length(resGaussian1$predAdjustedRandPseudoCentresMean)] <- resGaussian1$predAdjustedRandPseudoCentresMean
    res[5,1:length(resGaussian2$predAdjustedRandPseudoCentresMean)] <- resGaussian2$predAdjustedRandPseudoCentresMean
    res[6,1:length(resGaussian3$predAdjustedRandPseudoCentresMean)] <- resGaussian3$predAdjustedRandPseudoCentresMean
    
    nSuccesses <- c(resLinear$nSuccesses, resPolynomialNormalized2$nSuccesses, resPolynomialNormalized3$nSuccesses, resGaussian1$nSuccesses, resGaussian2$nSuccesses, resGaussian3$nSuccesses)
    
    cat("Min successes: ", min(nSuccesses), "\n")
    
    ylim <- c(min(res),max(res) + (max(res) - min(res)) * 0.1)
    
    ylabel <- "ARI"
    
    nVals <- 5        
    yStep <- round((max(ylim) - min(ylim)) / (nVals - 1),digits=2)
    yMax <- round(max(ylim),digits=2)
    yMin <- yMax - (nVals - 1) * yStep        
    yLabelSeq <- seq(yMax, yMin, -yStep)
    
    setEPS()
    
    postscript(paste("results-eval-", datasetname, "-", ylabel, "1.eps", sep=""))
    # margins: top, left, bottom, right. Default: 5,4,4,2 + 0.1
    par(mar=c(5, 4 + 0.5, 4, 2))
    plot(xLabelSeq, res[1,],type='o',ylim=ylim, lwd=2, lty = "solid", col='grey', main=sapply(datasetname, simpleCap), xlab="number of clusters", ylab=ylabel, pch=21, yaxt="n", xaxt="n", cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=1.7, cex=2.0)
    axis(1, at=xLabelSeq, labels=xLabelSeq, cex.axis=2.0)
    axis(2, at=yLabelSeq, labels=yLabelSeq, cex.axis=2.0)
    lines(xLabelSeq, res[4,],type='o', lwd=2, lty="longdash", col='grey', pch=23, cex=2.0)
    lines(xLabelSeq, res[5,],type='o', lwd=2, lty="solid", col='black', pch=24, cex=2.0)
    lines(xLabelSeq, res[6,],type='o', lwd=2, lty="dotted", col='black', pch=25, cex=2.0)
    legend(legendLocationARI, lwd=c(2,2,2,2), pch=c(21,23,24,25), col=c('grey', 'grey', 'black', 'black'), lty=c("solid","longdash","solid","dotted"), legend=c("linear",paste("gauss., ", a_sigma[1], sep=""),paste("gauss., ", a_sigma[2], sep=""), paste("gauss., ", a_sigma[3], sep="")), cex=1.5)
    dev.off()
    
    postscript(paste("results-eval-", datasetname, "-", ylabel, "2.eps", sep=""))
    # margins: top, left, bottom, right. Default: 5,4,4,2 + 0.1
    par(mar=c(5, 4 + 0.5, 4, 2))
    plot(xLabelSeq, res[1,],type='o',ylim=ylim, lwd=2, lty = "solid", col='grey', main=sapply(datasetname, simpleCap), xlab="number of clusters", ylab=ylabel, pch=21, yaxt="n", xaxt="n", cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=1.7, cex=2.0)
    axis(1, at=xLabelSeq, labels=xLabelSeq, cex.axis=2.0)
    axis(2, at=yLabelSeq, labels=yLabelSeq, cex.axis=2.0)
    lines(xLabelSeq, res[2,],type='o', lwd=2, lty="longdash", col='grey', pch=23, cex=2.0)
    lines(xLabelSeq, res[3,],type='o', lwd=2, lty="solid", col='black', pch=24, cex=2.0)
    legend(legendLocationARI2, lwd=c(2,2,2), pch=c(21,23,24), col=c('grey', 'grey', 'black'), lty=c("solid","longdash","solid","dotted"), legend=c("linear", "poly., norm., 2","poly., norm., 3"), cex=1.5)
    dev.off()
}
