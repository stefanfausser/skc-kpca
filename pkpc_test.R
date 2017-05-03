## Used solely in pkpc_experiments.R

simpleCap <- function(x)
{
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}

pkpc.test.unsupervised.framework <- function(param, basic = TRUE, evalWeight = FALSE, evalWeightTwoDataChunks = FALSE, evalMaxDataChunks = FALSE)
{
    # Unsupervised
    
    if(evalWeight)
    {
        maxDataChunksSeq <- param$maxDataChunksEval
        a_L <- param$a_L[length(param$a_L)]
        
        weightFactorSeq <- param$weightFactorEvalSeq
        
        a_method_kpca <- 2 # KKMEANS
        subfix <- "-evalWeight"
    }
    else if(evalWeightTwoDataChunks)
    {
        a_L <- param$a_L[length(param$a_L)]
        weightFactorSeq <- param$weightFactorEvalSeq
        a_method_kpca <- 2 # KKMEANS
        subfix <- "-evalWeightTwoDataChunks"
    }
    else if(evalMaxDataChunks)
    {
        maxDataChunksSeq <- param$maxDataChunksEvalSeq
        a_L <- param$a_L[length(param$a_L)]
        weightFactorSeq <- param$weightFactor
        a_method_kpca <- 2 # KKMEANS
        subfix <- "-evalMaxDataChunks"
    }
    else
    {
        maxDataChunksSeq <- param$maxDataChunks
        a_L <- param$a_L
        
        weightFactorSeq <- param$weightFactor
        a_method_kpca <- 0:2
        a_method <- 3:5
        subfix <- ""
    }
    
    if(!evalWeight && !evalWeightTwoDataChunks && basic)
        pkpc.test.func(param$seedSamples, param$seedLabels, param$seed, param$dataset, param$datasetfunc, "-unsupervised-2fold-basic", a_method = a_method, a_M = param$MBasicSeq, K = param$K, degree = param$degree, offset = param$offset, sigma = param$sigma, normalizeKernel = param$normalizeKernel, fuzzifier = param$fuzzifier)
    
    if(evalWeightTwoDataChunks)
    {
        pkpc.test.func(param$seedSamples, param$seedLabels, param$seed, param$dataset, param$datasetfunc, paste("-unsupervised-2fold-kpca", subfix, sep=""), a_method = a_method_kpca, a_M = param$M, K = param$K, a_L = a_L, degree = param$degree, offset = param$offset, sigma = param$sigma, normalizeKernel = param$normalizeKernel, fuzzifier = param$fuzzifier, weightFactorSeq = weightFactorSeq, maxDataChunks = 2)
        
        pkpc.test.func(param$seedSamples, param$seedLabels, param$seed, param$dataset, param$datasetfunc, paste("-unsupervised-2fold-all", subfix, sep=""), a_method = a_method_kpca, a_M = param$M, K = param$K, a_L = 0, allSamplesFromLast = 1, degree = param$degree, offset = param$offset, sigma = param$sigma, normalizeKernel = param$normalizeKernel, fuzzifier = param$fuzzifier, weightFactorSeq = weightFactorSeq)
    }
    else if(!basic)
        pkpc.test.func(param$seedSamples, param$seedLabels, param$seed, param$dataset, param$datasetfunc, paste("-unsupervised-2fold-kpca", subfix, sep=""), a_method = a_method_kpca, a_M = param$M, K = param$K, a_L = a_L, degree = param$degree, offset = param$offset, sigma = param$sigma, normalizeKernel = param$normalizeKernel, fuzzifier = param$fuzzifier, weightFactorSeq = weightFactorSeq, maxDataChunks = maxDataChunksSeq)
}

pkpc.test.semisupervised.framework <- function(param, evalLabelledUnlabelledFactor = FALSE)
{
    # Semi-Supervised
    
    if(evalLabelledUnlabelledFactor)
    {
        labelFactorSeq <- param$labelFactorEvalSeq
        labelledUnlabelledFactorSeq <- param$labelledUnlabelledFactorEvalSeq
        a_method <- 5 # KKMEANS
        subfix <- "-evalLabelledUnlabelledFactor"
    }
    else
    {
        labelFactorSeq <- param$labelFactorSeq
        labelledUnlabelledFactorSeq <- param$labelledUnlabelledFactor
        a_method <- 3:5
        subfix <- ""
    }
    
    pkpc.test.func(param$seedSamples, param$seedLabels, param$seed, param$dataset, param$datasetfunc, paste("-unsupervised-2fold-skc", subfix, sep=""), a_method = a_method, a_M = param$M_semisupervised, K = param$K, labelFactorSeq = labelFactorSeq, labelledUnlabelledFactorSeq = labelledUnlabelledFactorSeq, degree = param$degree, offset = param$offset, sigma = param$sigma, normalizeKernel = param$normalizeKernel, fuzzifier = param$fuzzifier)
}

pkpc.test.func <- function(seedSamples, seedLabels, seed, datasetname, datasetfunc, filesuffix, a_method, a_M, a_L = 0, labelFactorSeq = 0, K = 2, degree = 1, offset = 0, sigma = 0, scaleFeatures = 1, normalizeKernel = 0, removeOutliers = TRUE, removeFeatures = TRUE, fuzzifier = 1.25, labelledUnlabelledFactorSeq = 1, R = 350, kFoldCrossValidation = 2, normalizeAssignmentsKFCM = 0, getHardClusteringsKFCM = 0, allSamplesFromLast = 0, maxDataChunksSeq = 0, weightFactorSeq = 0.5, epsilon = 0.000001, T = 100, hardAssignmentsApproxPseudoCentresRNG = 1, softInitialAssignmentsApproxPseudoCentresRNG = 0, excludeSamplesFromAPCInGreedySelection = 0)
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
    
    a_mu <- c(1)
    
    a_K <- K    
    
    dirname <- "labelsTmp/"
    if(!dir.exists(dirname))
    {
        dir.create(dirname)
    }
    
    savefile <- paste(dirname, datasetname, sep="")
    
    pkpca_params <-  list(verbose=0, maxRepeats=10, maxRowsKmatOut=1000, allSamplesFromLast=allSamplesFromLast, maxSamplesValidate=2000, getHardClusteringsKFCM=getHardClusteringsKFCM, normalizeAssignmentsKFCM=normalizeAssignmentsKFCM, limitSamplesByHardAssignments=0, localMinimum=0, hardAssignmentsApproxPseudoCentresRNG=hardAssignmentsApproxPseudoCentresRNG, softInitialAssignmentsApproxPseudoCentresRNG=softInitialAssignmentsApproxPseudoCentresRNG,  excludeSamplesFromAPCInGreedySelection=excludeSamplesFromAPCInGreedySelection)
    
    a_pkpcparam <- rbind(pkpca_params, NULL)
    
    fileToSave <- paste("results-", datasetname, filesuffix, "_", date(), sep="")
    
    # write parameter file
    
    fileToSaveParameters <- paste("parameters-", datasetname, filesuffix, "_", date(), sep="")
    
    parameterList <- list(kFoldCrossValidation=kFoldCrossValidation, scaleFeatures=scaleFeatures, removeFeatures=removeFeatures, removeOutliers=removeOutliers, normalizeKernel=normalizeKernel, degree=degree, offset=offset, sigma=sigma, fuzzifier=fuzzifier, T=T, epsilon=epsilon, normalizeAssignmentsKFCM=normalizeAssignmentsKFCM, getHardClusteringsKFCM=getHardClusteringsKFCM, hardAssignmentsApproxPseudoCentresRNG=hardAssignmentsApproxPseudoCentresRNG, softInitialAssignmentsApproxPseudoCentresRNG=softInitialAssignmentsApproxPseudoCentresRNG,  excludeSamplesFromAPCInGreedySelection=excludeSamplesFromAPCInGreedySelection)
    
    write.table(parameterList, fileToSaveParameters)
    
    samplesMat <- saveOrRestoreSampleMat(seedSamples, R, N)
    
    ## only performed once
    createAndSaveRandomZeroLabels(seedLabels, savefile, labels, labelFactorSeq, R)
    
    pkpc.main.test(seed, fileToSave, savefile, a_pkpcparam, a_mu, a_M, a_K, a_L, a_method, x, kernelfunc, kernelparam, labelledUnlabelledFactorSeq, labelFactorSeq, samplesMat, kFoldCrossValidation, T, N, labels, fuzzifier, epsilon, scaleFeatures, maxDataChunksSeq, weightFactorSeq)
}

pkpc.plot.get.ylim <- function(vals, nVals = 5, nDigits = 2, spaceTopFactor = 0.1)
{
    ylim <- c(min(vals),max(vals) + (max(vals) - min(vals)) * spaceTopFactor)
    
    yStep <- round((max(ylim) - min(ylim)) / (nVals - 1),digits=nDigits)
    yMax <- round(max(ylim),digits=nDigits)
    yMin <- yMax - (nVals - 1) * yStep        
    yLabelSeq <- seq(yMax, yMin, -yStep)
    
    return(list(yLabelSeq=yLabelSeq, ylim=ylim))
}


pkpc.semisupervised.evaluation <- function(datasetname, date, legendLocationARI = "topleft", legendLocationDBI = "topright", shiftDBI=NULL, shiftARI=NULL, spaceTopFactor=0.35)
{
    if(length(date) != length(datasetname))
    {
        cat("Length of datasetname and date must match\n")
        return(NULL)
    }
    
    iSeq <- 1:length(date)
    
    valsARI <- NULL
    valsDBI <- NULL
    legendSeq <- NULL
    
    for(i in iSeq)
    {
        file <- list.files(pattern = paste("results-", datasetname[i], ".+", date[i], sep=""))
        
        if(!file.exists(file))
        {
            cat("File ", file, " does not exist\n")
            return(NULL)
        }
        
        res <- read.table(file)
        
        xVal <- res$labelledUnlabelledFactor
        xLabel <- "gamma"
        mainLabel <- "SKC(KKM)"
        
        valsARICurr <- res$predAdjustedRandPseudoCentresMean
        valsDBICurr <- res$errDaviesBouldinOrgRmsePseudoCentresMean
        
        if(!is.null(shiftARI))
            valsARICurr <- valsARICurr + shiftARI[i]
        if(!is.null(shiftDBI))
            valsDBICurr <- valsDBICurr + shiftDBI[i]
        
        valsARI <- rbind(valsARI, valsARICurr)
        valsDBI <- rbind(valsDBI, valsDBICurr)
        legendSeq <- c(legendSeq, sapply(datasetname[i],simpleCap))
    }
    
    lwdSeq <- c(2,2,2,2,2)
    pchSeq <- c(21,23,24,25,22)
    colSeq <- c('grey', 'grey', 'black', 'black','black')
    ltySeq <- c("solid","longdash","solid","dotted",'longdash')
    
    for(j in 1:2)
    {
        if(j == 1)
        {
            # ARI
            ylabel <- "ARI"            
            vals <- valsARI
            legendLocation <- legendLocationARI
        }
        else
        {
            # DBI
            ylabel <- "DBI"
            vals <- valsDBI
            legendLocation <- legendLocationDBI
        }            
        
        xLabelSeq <- xVal        
        
        ret <- pkpc.plot.get.ylim(vals, spaceTopFactor = spaceTopFactor)
        
        ylim <- ret$ylim
        yLabelSeq <- ret$yLabelSeq
        
        setEPS()            
        postscript(paste("results-semisupervised-eval-", ylabel, ".eps", sep=""))        
        for(i in iSeq)
        {
            if(i == 1)
            {
                # margins: top, left, bottom, right. Default: 5,4,4,2 + 0.1
                par(mar=c(5, 4 + 0.5, 4, 2))
                
                plot(xLabelSeq, vals[1,], type='o',ylim=ylim, lwd=lwdSeq[i], lty=ltySeq[i], col=colSeq[i], main=mainLabel, xlab=xLabel, ylab=ylabel, pch=pchSeq[i], yaxt="n", xaxt="n", cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=1.7, cex=2.0)
                axis(1, at=xLabelSeq, labels=xLabelSeq, cex.axis=2.0)
                axis(2, at=yLabelSeq, labels=yLabelSeq, cex.axis=2.0)
            }
            else
                lines(xLabelSeq, vals[i,], type='o', lwd=lwdSeq[i], lty=ltySeq[i], col=colSeq[i], pch=pchSeq[i], cex=2.0)
        }
        
        legend(legendLocation, lwd=lwdSeq[iSeq], pch=pchSeq[iSeq], col=colSeq[iSeq], lty=ltySeq[iSeq], legend=legendSeq[iSeq], cex=1.5, ncol=2)
        dev.off()
    }
}

pkpc.test.evaluation.maxDataChunks <- function(datasetname, date, legendLocation = "top")
{
    
    file <- list.files(pattern = paste("results-", datasetname, ".+", date, sep=""))
    
    if(!file.exists(file))
    {
        cat("File ", file, " does not exist\n")
        return(NULL)
    }
    
    ret <- read.table(file)
    
    valsARI <- ret$predAdjustedRandPseudoCentresMean
    valsSSE <- ret$errPseudoCentreQuantMean
    
    xLabelSeq <- ret$maxDataChunks
    
    ret <- pkpc.plot.get.ylim(valsARI, nVals = 4, nDigits = 3)
    yLabelSeq <- ret$yLabelSeq
    ylim <- ret$ylim
    
    ret2 <- pkpc.plot.get.ylim(valsSSE, nVals = 4, nDigits = 0)
    yLabelSeq2 <- ret2$yLabelSeq
    ylim2 <- ret2$ylim
    
    setEPS()            
    postscript(paste("results-", datasetname, "-maxChunks.eps", sep=""))        
    
    # margins: top, left, bottom, right. Default: 5,4,4,2 + 0.1
    par(mar=c(5, 4.5, 4, 4) + 0.1)
    
    plot(xLabelSeq, valsARI, type='o', ylim=ylim, lwd=2, lty = "solid", col='grey', main=paste(sapply(datasetname,simpleCap), ", KPC-A(KKM)", sep=""), xlab="number of data chunks", ylab="ARI", pch=21, yaxt="n", xaxt="n", cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=1.7, cex=2.0)
    axis(1, at=xLabelSeq, labels=xLabelSeq, cex.axis=2.0)
    axis(2, at=yLabelSeq, labels=yLabelSeq, cex.axis=2.0)
    
    par(new=TRUE)
    
    plot(xLabelSeq, valsSSE, type='o', ylim=ylim2, lwd=2, lty = "solid", col='black', xlab="", ylab="", pch=23, yaxt="n", xaxt="n", cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=1.7, cex=2.0)
    mtext("SSE",side=4,line=3,cex=2)
    axis(4, at=yLabelSeq2, labels=yLabelSeq2, cex.axis=2.0)
    
    legend(legendLocation, lwd=c(2,2), pch=c(21,23), col=c("gray","black"), lty=c("solid","solid"), legend=c("ARI","SSE"), cex=1.5)
    
    dev.off()
}

pkpc.test.evaluation.dist <- function(datasetname, date, namelen = 5)
{
    iSeq <- 1:length(date)
    
    vals <- NULL    
    
    for(i in iSeq)
    {
        file <- list.files(pattern = paste("results-", datasetname[i], ".+", date[i],sep=""))
        
        if(!file.exists(file))
        {
            cat("File ", file, " does not exist\n")
            return(NULL)
        }
        
        ret <- read.table(file)
        
        if(i == 1)
            Lseq <- unique(ret$L)
        
        ind <- NULL
        for(l in Lseq)
            ind <- c(ind, which(ret$L == l))
        
        ret <- ret[ind,]
        
        ind2 <- ret$method == 2
        vals <- cbind(vals, ret[ind2,]$bestDistsMean)        
    }
    
    ylim <- c(0,0.02)
    
    setEPS()            
    postscript(paste("results-dists.eps", sep=""))
    
    # margins: top, left, bottom, right. Default: 5,4,4,2 + 0.1
    par(mar=c(5, 4.5, 4, 4) + 0.1)
    
    barplot(vals, names.arg=substr(datasetname,1,namelen), ylab="error", col=c("white","gray","black"), beside=TRUE, ylim=ylim, cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=1.7, cex=1.9, main="KPC-A(KKM)")
    legend("topleft", legend = Lseq, fill = c("white","gray","black"), cex=1.5)
    
    dev.off()
}

pkpc.test.evaluation <- function(datasetname, date, semisupervised = FALSE, baselineSamples = 0)
{
    nDigitsExternal <- 3
    nDigitsInternal <- 3
    nDigitsInternalQuant <- 0
    nDigitsInternalQuantSd <- 1
    
    str <- NULL
    
    if((length(date) > 3 && semisupervised) || (length(date) > 2 && !semisupervised))
    {
        printf("Wrong length of date (%i)\n", length(date))
        return(NULL)
    }
    
    if(semisupervised && baselineSamples <= 0)
    {
        printf("Must set baselineSamples when semisupervised\n")
        return(NULL)
    }
    
    iSeq <- 1:length(date)
    
    index1Str <- "DBI"            
    index2Str <- "SSE"
    
    if(semisupervised)    
        str <- paste("Method \t\t & La. & NMI & ARI & ", index1Str, " & ", index2Str, "\\\\\n", sep="")
    else
        str <- paste("Method \t\t & L & NMI & ARI & ", index1Str, " & ", index2Str, "\\\\\n", sep="")
    
    # four indexes and three methods
    indexesBaseline <- matrix(NA, 3, 4)
    indexesBaselineSd <- matrix(NA, 3, 4)
    hasIndexesBaseline <- rep(FALSE, 3)
    indexesBaselineR <- rep(0,3)
    
    for(i in iSeq)
    {
        file <- list.files(pattern = paste("results-", datasetname, ".+", date[i],sep=""))
        
        if(!file.exists(file))
        {
            cat("File ", file, " does not exist\n")
            return(NULL)
        }
        
        res <- read.table(file)
        
        for(r in 1:nrow(res))
        {
            if(i == 1 && baselineSamples > 0)
            {
                if(res[r,]$M != baselineSamples)
                    next # for
            }
            
            if(res[r,]$method == 0)
                method <- "KPC-A(RNG)"
            else if(res[r,]$method == 1)
                method <- "KPC-A(KFCM)"
            else if(res[r,]$method == 2)
                method <- "KPC-A(KKM)"
            if(res[r,]$method == 3)
                method <- "RNG"
            else if(res[r,]$method == 4)
                method <- "KFCM"
            else if(res[r,]$method == 5)
            {
                if(i == 3)
                    method <- "SS-KKM"
                else
                    method <- "KKM"
            }
            
            doTtest <- FALSE
            haveBaseline <- FALSE
            if(res[r,]$method >= 3 && res[r,]$method <= 5)
            {
                method <- paste(method, "(\\num{", res[r,]$M ,"})", sep="")
                k <- "NA"
                
                baseLineMethod <- res[r,]$method - 2
                
                if(hasIndexesBaseline[baseLineMethod])
                {
                    doTtest <- TRUE
                    baseLineM <- baseLineMethod
                }
                else
                    haveBaseline <- TRUE
            }
            else
            {
                k <- paste(res[r,]$L, sep="")
                doTtest <- TRUE
                baseLineM <- res[r,]$method + 1
            }
            
            nDigitsInternal2 <- nDigitsInternal
            nDigitsInternal2Sd <- nDigitsInternal - 1
            
            extIndex1 <- res[r,]$predNMIPseudoCentresMean
            extIndex1Sd <- res[r,]$predNMIPseudoCentresStd
            extIndex2 <- res[r,]$predAdjustedRandPseudoCentresMean
            extIndex2Sd <- res[r,]$predAdjustedRandPseudoCentresStd
            
            index1 <- res[r,]$errDaviesBouldinOrgRmsePseudoCentresMean
            index1Sd <- res[r,]$errDaviesBouldinOrgRmsePseudoCentresStd
            
            index2 <- res[r,]$errPseudoCentreQuantMean
            index2Sd <- res[r,]$errPseudoCentreQuantStd
            nDigitsInternal2 <- nDigitsInternalQuant
            nDigitsInternal2Sd <- nDigitsInternalQuantSd
            
            strSignif1 <- rep("",4)
            strSignif2 <- rep("",4)
            
            if(doTtest)
            {
                signif <- rep(FALSE, 4)
                
                R1 <- indexesBaselineR[baseLineM]
                R2 <- res[r,]$nSuccesses
                
                ret <- myttest(c(extIndex1,indexesBaseline[baseLineM,1]),c(extIndex1Sd,indexesBaselineSd[baseLineM,1]),c(R1,R2))
                
                if(ret$p.value < 0.05 && extIndex1 > indexesBaseline[baseLineM,1])
                    signif[1] <- TRUE
                
                ret <- myttest(c(extIndex2,indexesBaseline[baseLineM,2]),c(extIndex2Sd,indexesBaselineSd[baseLineM,2]),c(R1,R2))
                
                if(ret$p.value < 0.05 && extIndex2 > indexesBaseline[baseLineM,2])
                    signif[2] <- TRUE
                
                ret <- myttest(c(index1,indexesBaseline[baseLineM,3]),c(index1Sd,indexesBaselineSd[baseLineM,3]),c(R1,R2))
                
                if(ret$p.value < 0.05 && index1 < indexesBaseline[baseLineM,3])
                    signif[3] <- TRUE
                
                ret <- myttest(c(index2,indexesBaseline[baseLineM,4]),c(index2Sd,indexesBaselineSd[baseLineM,4]),c(R1,R2))
                
                if(ret$p.value < 0.05 && index2 < indexesBaseline[baseLineM,4])
                    signif[4] <- TRUE
                
                for(j in 1:4)
                {
                    if(signif[j])
                    {
                        strSignif1[j] <- "\\bfseries "
                        strSignif2[j] <- ""
                    }
                }                
            }
            
            if(haveBaseline)
            {
                if(!hasIndexesBaseline[baseLineMethod])
                {
                    indexesBaseline[baseLineMethod,1] <- extIndex1
                    indexesBaseline[baseLineMethod,2] <- extIndex2
                    indexesBaseline[baseLineMethod,3] <- index1
                    indexesBaseline[baseLineMethod,4] <- index2
                    
                    indexesBaselineSd[baseLineMethod,1] <- extIndex1Sd
                    indexesBaselineSd[baseLineMethod,2] <- extIndex2Sd
                    indexesBaselineSd[baseLineMethod,3] <- index1Sd
                    indexesBaselineSd[baseLineMethod,4] <- index2Sd
                    
                    indexesBaselineR[baseLineMethod] <- res[r,]$nSuccesses
                    
                    hasIndexesBaseline[baseLineMethod] <- TRUE
                }
            }
            
            if(semisupervised)
                str2 <- res[r,]$labelFactor
            else
                str2 <- k
            
            str <- paste(str, method, " \t & ", str2,
                         " & ", strSignif1[1], round(extIndex1, digits=nDigitsExternal), " \\pm ", round(extIndex1Sd, nDigitsExternal), strSignif2[1],
                         " & ", strSignif1[2], round(extIndex2, digits=nDigitsExternal), " \\pm ", round(extIndex2Sd, nDigitsExternal), strSignif2[2],
                         " & ", strSignif1[3], round(index1, digits=nDigitsInternal), " \\pm ", round(index1Sd, nDigitsInternal), strSignif2[3],
                         " & ", strSignif1[4], round(index2, digits=nDigitsInternal2), " \\pm ", round(index2Sd, nDigitsInternal2), strSignif2[4],
                         "\\\\ \n", sep="")
        } # for r
    } # for i
    
    cat("Latex table entries:\n", str,"\n", sep="")    
}
