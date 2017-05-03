dyn.load("kernelclustering.so")
dyn.load("getD_1.so")
dyn.load("getD_2.so")
dyn.load("kernel_dist_prototypes3.so")
dyn.load("greedy_approximation.so")
dyn.load("contingency_matrix.so")
dyn.load("gauss_km_flexible.so")
dyn.load("poly_km_flexible.so")
dyn.load("scale_features.so")

library(doRNG)
library(doParallel)
library(R.utils)

registerDoParallel(cores=10)

# input (excerpt):
# samplesTrain: (permutated) sample indices of training set
# mu: number of parallelizations
# M: number of samples per (parallelized) clustering
# K: number of cluster
# L: number of samples for k-approximation
# method: 0 = kernel batch neural gas, 1 = kernel fuzzy c-means, 2 = kernel k-means, 3 = arbitrary choose M samples and cluster with kernel batch neural gas, 4 = arbitrary choose M samples and cluster with kernel fuzzy c-means, 5 = arbitrary choose M samples and cluster with kernel k-means
# kernelfunc: kernel function for calculating the kernel matrix
# kernelparam: kernel function parameters
# fuzzifier: parameter for kernel fuzzy c-means
# epsilon: parameter for kernel fuzzy c-means
# T: parameter for kernel batch neural gas, number of iterations

pkpc <- function(pkpcparam, samplesTrain, mu, M, K, L, method, data, kernelfunc, kernelparam, fuzzifier, epsilon, T, nClasses, reducedLabels, labelledUnlabelledFactor, labelFactor, maxDataChunks, weightFactor)
{
    if(labelFactor == 0)
    {
        supervised <- 0
        
        printf("Completely unsupervised clustering\n")
    }
    else
    {
        supervised <- 1
        printf("Semi-Supervised Clustering, (%i) classes are known, (%f) labels shall be known, labelledUnlabelledFactor (%f)\n", nClasses - 1, labelFactor, labelledUnlabelledFactor)        
    }
    
    N <- length(samplesTrain)
    
    printf("Training set size (%i)\n", N)
    
    if(supervised)
        printf("Cluster labels are detected in each step in the kernel clustering method\n")
    
    maxDataChunksTmp <- 0
    
    if(method >= 0 && method <= 2)
    {
        if(pkpcparam$allSamplesFromLast)
            maxDataChunksTmp <- 2
        else if(maxDataChunks > 0)
            maxDataChunksTmp <- maxDataChunks
    }
    else
        maxDataChunksTmp <- 1
    
    if(M > N)
        Mtmp <- N
    else
        Mtmp <- M
    
    nDataChunksReal <- floor(N / Mtmp)
    
    Nmax <- Mtmp * (nDataChunksReal - 1)
    
    if(maxDataChunksTmp > 0)
        nDataChunks <- min(nDataChunksReal, maxDataChunksTmp)
    else
        nDataChunks <- nDataChunksReal
    
    lastsampleI <- Mtmp * (nDataChunksReal - nDataChunks)
    
    printf("Index of first sample of first d-chunk (%i), index of first sample of last d-chunk (%i), number of d-chunks (%i)\n", lastsampleI + 1, Nmax + 1, nDataChunks)
    
    Nmaxtmp <- Nmax
    
    haveLp <- 0
    
    breakCondition <- 0
    
    w <- NULL
    prototypes <- NULL
    t <- 0
    lastPrototypePercentages <- NULL
    lastNumberBestSamples <- NULL
    
    bestDistsMean <- NULL
    bestDistsMin <- NULL
    
    nSamples <- NULL
    
    lambda_start <- K / 2
    lambda_end <- 0.01
    
    if(method >= 0 && method <= 2)
    {
        while(lastsampleI < N)
        {
            summed_w <- NULL
            summed_reduced_samples <- NULL
            
            wbak <- w
            for(p in 1:mu)
            {
                # is the current data chunk the last data chunk or are there no more data chunks?
                if(p > 1 && lastsampleI == Nmaxtmp)
                    break # for
                
                additionalSamplesIndex <- (lastsampleI+1):(lastsampleI+Mtmp)
                newSamples <- samplesTrain[additionalSamplesIndex]
                actualSamples <- c(prototypes, newSamples)
                w <- c(wbak, rep(1, length(additionalSamplesIndex)))
                
                ## with the following lines the K \cdot k samples from the last data chunk have no sample label
                actualLabels <- c(rep(0,length(prototypes)), reducedLabels[newSamples])
                
                nLabelledSamplesClass <- rep(0, nClasses)
                
                if(supervised)
                {
                    # some (or all) samples are labelled
                    for(m in 1:nClasses)
                        nLabelledSamplesClass[m] <- sum(reducedLabels[newSamples] == (m - 1))
                }
                else
                {
                    # all samples are unlabelled
                    nLabelledSamplesClass[1] <- length(newSamples)
                }
                
                printf("Number of samples with matching label (first value = number of unlabelled samples):\n");
                print(nLabelledSamplesClass)
                
                lastsampleI <- lastsampleI + length(additionalSamplesIndex)
                
                printf("t = (%i), p = (%i), samples to cluster = (%i)\n", t, p, length(actualSamples))
                
                nRepeats <- 0
                while(1)
                {
                    # start clustering
                    ret <- pkpc.cluster(pkpcparam, actualSamples, w, K, L, method, data, kernelfunc, kernelparam, actualLabels, nLabelledSamplesClass, supervised, fuzzifier, epsilon, lambda_start, lambda_end, T, lastNumberBestSamples, labelledUnlabelledFactor, FALSE, weightFactor)
                    
                    if(!ret$emptyCluster)
                        break
                    
                    printf("\t Got empty cluster(s). Repeat the clustering.\n")
                    
                    nRepeats <- nRepeats + 1
                    
                    if(nRepeats > pkpcparam$maxRepeats)
                    {
                        breakCondition <- 1
                        break
                    }
                }
                
                if(breakCondition)
                    break # for
                
                bestDistsMean <- c(bestDistsMean, mean(ret$bestDists))
                bestDistsMin <- c(bestDistsMin, min(ret$bestDists))
                
                if(p == 1)
                {
                    if(t > 0)
                        nSamples <- c(nSamples, length(prototypes))
                    
                    lastSamples <- actualSamples
                    lastW <- w
                    lastAlpha <- ret$alpha
                    lastF <- ret$f
                    lastNumberBestSamples <- ret$nBestSamples
                }
                
                summed_reduced_samples <- c(summed_reduced_samples, ret$reduced_samples)
                summed_w <- c(summed_w, ret$w)
                
                realmu <- p
            } # for p in 1:mu
            
            if(breakCondition)
                break # while
            
            prototypes <- summed_reduced_samples
            w <- summed_w / realmu
            
            if(haveLp)
            {
                percent <- length(intersect(prototypes,lp)) / length(prototypes)
                
                lastPrototypePercentages <- c(lastPrototypePercentages, percent)
                
                if(pkpcparam$verbose)
                    cat("Have ", percent * 100, "% prototypes of last prototypes\n")
            }
            else
                haveLp <- 1
            
            lp <- prototypes
            
            # prepare last data chunk
            if(mu > 1 && lastsampleI == Nmax)
            {
                if(pkpcparam$verbose)
                    printf("Preparing last data chunk (will be clustered serialized)\n")
                
                Nmaxtmp <- N
            }
            
            t <- t + 1
        } # for all data chunks
        
        if(breakCondition)
            return( list(emptyCluster=1) )
    }
    else
    {
        method <- method - 3
        
        # choose last data chunk
        ind <- (Nmax + 1):N
        
        actualSamples <- samplesTrain[ind]
        actualLabels <- reducedLabels[actualSamples]
        
        # each sample has the same weight
        w <- rep(1,length(actualSamples))
        
        nLabelledSamplesClass <- rep(0, nClasses)
        
        if(supervised)
        {
            # some (or all) samples are labelled
            for(m in 1:nClasses)
                nLabelledSamplesClass[m] <- sum(actualLabels == (m - 1))                
        }
        else
        {
            # all samples are unlabelled
            nLabelledSamplesClass[1] <- length(actualSamples)                
        }
        
        printf("Number of samples with matching label (first value = number of unlabelled samples):\n");
        print(nLabelledSamplesClass)
        
        nRepeats <- 0
        cat("(randomly chosen) samples to cluster =",length(actualSamples),"\n")
        
        while(1)
        {
            # start clustering
            
            ret <- pkpc.cluster(pkpcparam, actualSamples, w, K, 0, method, data, kernelfunc, kernelparam, actualLabels, nLabelledSamplesClass, supervised, fuzzifier, epsilon, lambda_start, lambda_end, T, NULL, labelledUnlabelledFactor, TRUE, 0)
            
            if(!ret$emptyCluster)
                break
            
            printf("\t Got empty cluster(s). Repeat the clustering.\n")
            
            nRepeats <- nRepeats + 1
            
            if(nRepeats > pkpcparam$maxRepeats)
            {
                breakCondition <- 1
                break
            }
        }
        
        if(breakCondition)
            return( list(emptyCluster=1) )
        
        nSamples <- NULL
        lastSamples <- actualSamples
        lastW <- w
        lastAlpha <- ret$alpha        
        lastF <- ret$f
        
        bestDistsMean <- NA
        bestDistsMin <- NA
        nSamples <- NA
    }
    
    printf("Clustering: done. Starting internal / external cluster validations\n")
    
    if(pkpcparam$verbose)
    {
        w_k <- rep(0,K)
        for(ki in 1:K)
            w_k[ki] <- sum(lastW * lastF[ki,])
        
        printf("Final number of samples in cluster:\n")
        print(w_k)
    }
    
    if(pkpcparam$verbose >= 3)
    {
        printf("Final last samples:\n")
        print(lastSamples)
        printf("Final weights:\n")
        print(lastW)
        printf("Final sample-to-cluster weights:\n")
        print(lastAlpha)
        printf("Final sample-to-cluster assignments:\n")
        print(lastF)
    }
    
    if(pkpcparam$verbose >= 2)
    {
        printf("Final weights sum: %f\n", sum(lastW))
        printf("List of last prototypes percentages (divided by 100):\n")
        print(lastPrototypePercentages)    
    }
    
    # input, given by procedures above: lastSamples, lastW, lastAlpha, lastF
    
    return(list(emptyCluster = 0, lastSamples = lastSamples, lastF = lastF, lastW = lastW, lastAlpha = lastAlpha, nSamplesMean = mean(nSamples), bestDistsMean = mean(bestDistsMean), bestDistsMin = mean(bestDistsMin)))
}

pkpc.calculateIndexes <- function(pkpcparam, samplesValidate, K, data, kernelfunc, kernelparam, labels, lastSamples, lastF, lastW, lastAlpha)
{
    printf("Validation set size (%i)\n", length(samplesValidate))
    
    errPseudoCentresQuant <- 0
    
    summed_diameterPseudoCentres <- rep(0,K)
    
    # calculate input kernel matrix
    KmatInPseudoCentres <- kernelfunc(data, kernelparam, lastSamples, lastSamples, FALSE, TRUE)
    
    ret2 <- pkpc.distToPseudoCentres1(KmatInPseudoCentres, lastF, lastW, lastAlpha)
    
    largeDAll <- NULL
    
    Nvalidate <- length(samplesValidate)
    lastsampleI <- 0
    
    while(lastsampleI < Nvalidate)
    {
        sampleIndices <- (lastsampleI+1):min((lastsampleI+pkpcparam$maxRowsKmatOut),Nvalidate)
        
        validateSampleIndices <- samplesValidate[sampleIndices]
        
        lastsampleI <- lastsampleI + length(sampleIndices)
        
        KmatDiag <- kernelfunc(data, kernelparam, validateSampleIndices, NULL, TRUE, FALSE)
        
        ## calculate distance of samples validateSampleIndices (subset of all samples) to the pseudo cluster centres
        
        # calculate output kernel matrix
        KmatOutPseudoCentres <- kernelfunc(data, kernelparam, validateSampleIndices, lastSamples, FALSE, FALSE)
        
        ret <- pkpc.distToPseudoCentres2(KmatOutPseudoCentres, KmatDiag, ret2$g, ret2$amount_points, ret2$connection_strength)
        
        errPseudoCentresQuant <- errPseudoCentresQuant + ret$err
        summed_diameterPseudoCentres <- summed_diameterPseudoCentres + ret$summed_diameter
        
        largeDAll <- cbind(largeDAll, ret$d)                
    } # while lastsampleI < Nvalidate
    
    if(pkpcparam$verbose)
        printf("(pseudo-centres) err = (%f)\n", errPseudoCentresQuant)
    
    # build hard clusters
    
    largeAAll <- matrix(0,K,Nvalidate)
    for(j in 1:Nvalidate)
        largeAAll[which.min(largeDAll[,j]),j] <- 1
    
    nSamplesInClusterPseudoCentres <- rep(0,K)
    
    for(ki in 1:K)
        nSamplesInClusterPseudoCentres[ki] <- sum(largeAAll[ki,])
    
    # check nSamplesInClusterPseudoCentres
    if(sum(nSamplesInClusterPseudoCentres) != Nvalidate)
    {
        retList <- list(emptyCluster=1)
        printf("sum(nSamplesInClusterPseudoCentres) != Nvalidate, (%i) != (%i)\n", sum(nSamplesInClusterPseudoCentres), Nvalidate)
        return(retList)
    }
    
    # empty cluster?
    if(sum(nSamplesInClusterPseudoCentres == 0) > 0)
        return(list(emptyCluster=1))    
    
    values <- .C("kernel_dist_prototypes3", din=as.numeric(largeDAll), A=as.integer(largeAAll), dout=as.numeric(matrix(0,K,K)), K=as.integer(K), N=as.integer(Nvalidate))    
    
    ret3 <- list(d = matrix(values$dout, values$K, values$K))        
    
    distancesBetweenPseudoCentres <- ret3$d
    
    # Use the average distances of the samples of one pseudo-centre to another pseudo-centre as distances between clusters (may be the best approach)
    
    errDaviesBouldinOrgPseudoCentres <- mydaviesbouldin(nSamplesInClusterPseudoCentres, summed_diameterPseudoCentres, distancesBetweenPseudoCentres, FALSE)
    
    errDaviesBouldinOrgRmsePseudoCentres <- mydaviesbouldin(nSamplesInClusterPseudoCentres, summed_diameterPseudoCentres, distancesBetweenPseudoCentres, TRUE)
    
    if(pkpcparam$verbose)
    {
        printf("Calculating the Davies–Bouldin Index with Pseudo-Centres:\n");
        printf("(pseudo-centres) Distances between the cluster prototypes:\n")
        print(distancesBetweenPseudoCentres)
        printf("(pseudo-centres) Summed diameter:\n")
        print(summed_diameterPseudoCentres)
        printf("Number of samples in cluster:\n")
        print(nSamplesInClusterPseudoCentres)
    }
    
    if(pkpcparam$verbose)
        printf("(pseudo-centres) Davies–Bouldin Index: (%f)\n", errDaviesBouldinOrgPseudoCentres)
    
    # Use all labels (and not the reduced set of labels)
    sample_labels <- labels[samplesValidate[1:Nvalidate]]
    
    ## with (approximate) pseudo-centres        
    
    contingencyMatrix <- matrix(0,max(labels),K)
    
    # calculate contingency matrix
    values <- .C("contingency_matrix", f=as.numeric(largeAAll), sample_labels=as.integer(sample_labels), K=as.integer(K), N=as.integer(Nvalidate), M=as.integer(max(labels)), c=as.integer(contingencyMatrix))
    
    contingencyMatrix <- matrix(values$c, values$M, values$K)
    
    if(sum(as.numeric(contingencyMatrix)) != Nvalidate)
    {
        printf("### Internal error, sum of entries in contingencyMatrix is not equal to N\n")
        return(list(emptyCluster=1))
    }
    
    # calculate NMI with all samples
    predNMI <- swc.calculateNMI(contingencyMatrix, 0)
    
    # calculate Purity with all samples
    predPurity <- 0
    for(ki in 1:K)
        predPurity <- predPurity + max(contingencyMatrix[,ki])
    
    predPurity <- 1 / Nvalidate * predPurity
    
    # calculate adjusted Rand Index
    predAdjustedRand <- calculateAdjustedRandIndex(contingencyMatrix)
    
    if(pkpcparam$verbose)
    {
        printf("NMI Pseudo-Centres: (%f)\n", predNMI)
        printf("Purity Pseudo-Centres: (%f)\n", predPurity)
        printf("Adjusted Rand Index Pseudo-Centres: (%f)\n", predAdjustedRand)
        
        printf("Contingency matrix:\n");
        print(contingencyMatrix)
    }
    
    # Number of samples from the last data chunk
    nsamples <- length(lastSamples)
    
    nSamplesClusterPseudoCentresMin <- min(nSamplesInClusterPseudoCentres)
    
    retList <- list(emptyCluster=0, nsamples=nsamples,
                    nSamplesClusterPseudoCentresMin=nSamplesClusterPseudoCentresMin,
                    errPseudoCentresQuant=errPseudoCentresQuant,
                    errDaviesBouldinOrgPseudoCentres=errDaviesBouldinOrgPseudoCentres,
                    errDaviesBouldinOrgRmsePseudoCentres=errDaviesBouldinOrgRmsePseudoCentres,
                    predPurityPseudoCentres=predPurity,
                    predNMIPseudoCentres=predNMI,
                    predAdjustedRandPseudoCentres=predAdjustedRand)
    
    return(retList)
}

pkpc.distToPseudoCentres1 <- function(KmatIn, f, w, alpha)
{
    values <- .C("getD_1", KmatIn=as.numeric(KmatIn), f=as.numeric(f), K=as.integer(nrow(f)), Nin=as.integer(nrow(KmatIn)), win=as.numeric(w), alpha=as.numeric(alpha), g=as.numeric(rep(0,nrow(f))), amount_points=as.numeric(rep(0,nrow(f))), connection_strength=as.numeric(matrix(0,nrow(f),nrow(KmatIn))))
    
    ret <- list(g = values$g, amount_points = values$amount_points, connection_strength = matrix(values$connection_strength,values$K,values$Nin))
    
    return(ret)
}

pkpc.distToPseudoCentres2 <- function(KmatOut, KmatDiag, g, amount_points, connection_strength)
{
    K <- nrow(connection_strength)
    
    values <- .C("getD_2", KmatOut=as.numeric(KmatOut), KmatDiag=as.numeric(KmatDiag), K=as.integer(K), Nin=as.integer(ncol(KmatOut)), Nout=as.integer(nrow(KmatOut)),
                 g=as.numeric(g), amount_points=as.numeric(amount_points), connection_strength=as.numeric(connection_strength), d=as.numeric(matrix(0,K,nrow(KmatOut))), err=as.numeric(0), summed_diameter=as.numeric(rep(0,K)))
    
    ret <- list(err = values$err, d = matrix(values$d,values$K,values$Nout), summed_diameter=values$summed_diameter)
    
    return(list(d=ret$d,err=ret$err,summed_diameter=ret$summed_diameter))
}

pkpc.initializeClusterAssignments.hard <- function(K, N, equalFac=0.5)
{
    f <- matrix(0,K,N)
    
    Nequally <- floor(N * equalFac)
    Nrest <- N - Nequally
    
    ind <- sample(1:N, Nequally)
    
    # assign 'Nequally' samples equally to the clusters (avoids empty clusters)
    
    k <- 1
    for(i in 1:Nequally)
    {
        f[k,ind[i]] <- 1
        
        k <- (k + 1) %% K
        if(k == 0)
            k <- K
    }
    
    ind2 <- setdiff(1:N, ind)
    
    # randomly assign the rest of the samples to the clusters
    for(i in ind2)
    {
        ki <- round(runif(1, min=1, max=K))
        f[ki,i] <- 1
    }
    
    return(f)
}

# input: 
# Kmat: kernel matrix for samples to cluster, size: number of samples x number of samples
# w: weight matrix for samples to cluster, size: number of cluster x number of samples
# K: number of cluster
# method: 0 = kernel batch neural gas, 1 = kernel fuzzy c-means
# fuzzifier: parameter for kernel fuzzy c-means
# epsilon: parameter for kernel fuzzy c-means
# T: maximum iterations to cluster

pkpc.cluster <- function(pkpcparam, samples, w, K, L, method, data, kernelfunc, kernelparam, labels, nLabelledSamplesClass, supervised, fuzzifier, epsilon, lambda_start, lambda_end, T, lastNumberBestSamples, labelledUnlabelledFactor, isBasicClustering, weightFactor)
{
    # calculate kernel matrix
    Kmat <- kernelfunc(data, kernelparam, samples, samples, FALSE, TRUE)
    
    N <- nrow(Kmat)
    Nnew <- N - sum(lastNumberBestSamples)    
    nPrototypesIn <- length(lastNumberBestSamples)
    
    if(pkpcparam$verbose)
    {
        printf("Number of input prototypes from last data chunk: %f\n", nPrototypesIn);
        print(lastNumberBestSamples)
        
        printf("Sample labels:\n")
        print(table(labels))
    }
    
    if(pkpcparam$verbose >= 2)
    {
        printf("w:\n");
        print(w)
        printf("Summed weights input = (%f)\n", sum(w))
        
        printf("Sample labels (vector):\n")
        print(labels)
    }
    
    if(pkpcparam$verbose >= 2)
        printf("Initialize f randomly\n")
    
    # initialize sample-to-cluster assignment matrix f
    
    f <- NULL
    
    # k-approximation need to stay together -> same initialization
    if(nPrototypesIn > 0)
    {
        # The previous prototypes shall be assigned to different clusters
        for(ki in 1:nPrototypesIn)
        {
            if(pkpcparam$softInitialAssignmentsApproxPseudoCentresRNG && method == 0)
                val <- exp(-1.0 * sample(1:K-1) / lambda_start)
            else
            {
                val <- rep(0,K)
                val[ki] <- 1
            }
            
            for(ka in 1:lastNumberBestSamples[ki])
                f <- cbind(f, val)
        }
        
        NnewTmp <- Nnew
    }
    else
        NnewTmp <- N
    
    if(method == 0)
    {
        for(i in 1:NnewTmp)
            f <- cbind(f, exp(-1.0 * sample(1:K-1) / lambda_start))
        
        # normalize f (input)
        for(i in 1:N)
            f[,i] <- (f[,i] / sum(f[,i]))
    }
    else if (method == 1)
    {
        for(i in 1:NnewTmp)
            f <- cbind(f, runif(K, min=0, max=1))
        
        # normalize f (input)
        for(i in 1:N)
            f[,i] <- (f[,i] / sum(f[,i]))
    }
    else if(method == 2)
    {
        fnew <- pkpc.initializeClusterAssignments.hard(K, NnewTmp)
        
        f <- cbind(f, fnew)
    }
    
    # for the first time clustering there are no old samples forming the pseudo cluster centres (missing K * k samples)
    
    if(nPrototypesIn > 0)
        distToPseudoCentreMethodTmp <- 1
    else
        distToPseudoCentreMethodTmp <- 0
    
    decreaseLambda <- 1
    verbose <- 0
    
    nClasses <- length(nLabelledSamplesClass)
    
    if(labelledUnlabelledFactor < 0)
        supervisedTmp <- 0
    else
        supervisedTmp <- supervised
    
    mySeed <- runif(1, min=0, max=0x7FFFFFFF)
    
    values <- .C("kernelclustering", seed=as.integer(mySeed), retval=as.integer(0), method=as.integer(method), Kmat=as.numeric(Kmat), w=as.numeric(w), alpha=as.numeric(matrix(1,K,N)), cMat=as.numeric(matrix(0,K,nClasses)), l=as.integer(labels), nLabelledSamplesClass=as.numeric(nLabelledSamplesClass), M=as.integer(nClasses), bSupervised=as.integer(supervisedTmp), winner_label=as.integer(rep(0,K)), winner_count=as.numeric(rep(0,K)), f=as.numeric(f), K=as.integer(K), mu=as.integer(nPrototypesIn), T=as.integer(T), lambda_start=as.numeric(lambda_start), lambda_end=as.numeric(lambda_end), N=as.integer(N), d=as.numeric(matrix(0,K,N)), errSoft=as.numeric(0), verbose=as.integer(verbose), decreaseLambda=as.integer(decreaseLambda), distToPseudoCentre=as.integer(distToPseudoCentreMethodTmp), nBestSamples=as.integer(lastNumberBestSamples), labelledUnlabelledFactor=as.numeric(labelledUnlabelledFactor), fuzzifier=as.numeric(fuzzifier), epsilon=as.numeric(epsilon), emptyCluster=as.integer(0), getHardClusteringsKFCM=as.integer(pkpcparam$getHardClusteringsKFCM), normalizeAssignmentsKFCM=as.integer(pkpcparam$normalizeAssignmentsKFCM), hardAssignmentsApproxPseudoCentresRNG=as.integer(pkpcparam$hardAssignmentsApproxPseudoCentresRNG))
    
    ret <- list(retval = values$retval, K = values$K, E = values$errSoft, f = matrix(values$f,K,N), d = matrix(values$d,K,N), emptyCluster = values$emptyCluster, winner_label = values$winner_label, cMat=matrix(values$cMat,K,nClasses), winner_count = values$winner_count, alpha = matrix(values$alpha,K,N))
    
    if(ret$emptyCluster)
        return(list(emptyCluster=1))
    
    if(ret$retval)
        return(list(emptyCluster=1))
    
    alpha <- ret$alpha
    
    if(pkpcparam$verbose)
    {
        printf("Cluster labels:\n")
        print(ret$winner_label)
        
        printf("Cluster winner-counts:\n")
        print(ret$winner_count)
    }
    
    if(pkpcparam$verbose >= 3)
    {
        printf("Contingency matrix:\n");
        print(ret$cMat)
    }
    
    if(pkpcparam$verbose >= 2)
        printf("Error returned by clustering method: (%f)\n", ret$E)
    
    f <- ret$f
    d <- ret$d
    
    # build hard clusters
    A <- matrix(0,K,N)
    for(j in 1:N)
        A[which.min(d[,j]),j] <- 1
    
    # empty cluster ?
    emptyCluster <- 0
    for(ki in 1:K)
    {
        if(sum(A[ki,])==0)
            emptyCluster <- 1
    }
    
    if(emptyCluster)
        return(list(emptyCluster=1))
    
    if(pkpcparam$verbose >= 2)
    {
        nSamplesInClusterCurrent <- rep(0,K)
        for(ki in 1:K)
            nSamplesInClusterCurrent[ki] <- sum(A[ki,])
        
        printf("Current number of samples in cluster (without previous samples):\n")
        print(nSamplesInClusterCurrent)
    }
    
    reduced_samples <- NULL
    nBestSamples <- rep(0,K)
    
    if(pkpcparam$allSamplesFromLast && !isBasicClustering)
    {
        # Use all samples
        
        BestSamples <- matrix(0,K,N)
        BestSamplesDist <- matrix(0,K,N) # first data chunk: L samples, second data chunk: N samples
        
        indWithoutAPC <- (sum(lastNumberBestSamples) + 1):N
        
        for(ki in 1:K)
        {
            nBestSamples[ki] <- sum(A[ki,indWithoutAPC])
            BestSamples[ki,1:nBestSamples[ki]] <- which(A[ki,indWithoutAPC] == 1)
        }
    }
    else if(!isBasicClustering)
    {
        # Nstart is the index of the first sample, starting with index zero
        
        Nstart <- 0
        if(pkpcparam$excludeSamplesFromAPCInGreedySelection)
            Nstart <- sum(lastNumberBestSamples)
        
        printf("Starting greedy sample-selection for the APCs, Nstart (%i)\n", Nstart)
        
        uniqueSampleIndices <- !duplicated(samples)
        
        values <- .C("greedy_approximation", Kmat=as.numeric(Kmat), f=as.numeric(f), A=as.integer(A), w=as.numeric(w), alpha=as.numeric(alpha), K=as.integer(K), nSamplesMax=as.integer(L), N=as.integer(N), best_samples=as.integer(matrix(-1,K,L)), uniqueSampleIndices=as.integer(uniqueSampleIndices), best_samples_dist=as.numeric(matrix(0,K,L)), number_samples=as.integer(rep(0,K)), limitSamplesByHardAssignments=as.integer(pkpcparam$limitSamplesByHardAssignments), localMinimum=as.integer(pkpcparam$localMinimum), Nstart=as.integer(Nstart), ret=as.integer(0))
        
        ret2 <- list(best_samples = matrix(values$best_samples,K,L), best_samples_dist = matrix(values$best_samples_dist,K,L), number_samples = values$number_samples, ret = values$ret)
        
        if(ret2$ret)
        {
            printf("greedy_approximation failed!\n")
            return(list(emptyCluster=1))            
        }
        
        nBestSamples <- ret2$number_samples
        BestSamples <- ret2$best_samples + 1
        BestSamplesDist <- ret2$best_samples_dist
    }
    
    bestDists <- rep(0,K)
    
    wnew <- NULL
    
    if(!isBasicClustering)
    {
        if(pkpcparam$verbose)
        {
            printf("Number best samples per cluster:\n")
            print(nBestSamples)
            print(sum(nBestSamples))
            
            printf("Distances to cluster centre:\n")
            for(ki in 1:K)
                print(BestSamplesDist[ki,1:nBestSamples[ki]])
        }
        
        for(ki in 1:K)
        {
            reduced_samples <- c(reduced_samples, samples[BestSamples[ki,1:nBestSamples[ki]]])
            bestDists[ki] <- BestSamplesDist[ki,nBestSamples[ki]] 
        }
        
        if(pkpcparam$verbose >= 3)
        {
            printf("reduced samples:\n")
            print(reduced_samples)
        }
        
        # Calculate the weights of the reduced_samples
        
        for(ki in 1:K)
        {
            nSamplesInCluster <- sum(w * f[ki,]) / nBestSamples[ki]
            
            wnew <- c(wnew, rep(nSamplesInCluster, nBestSamples[ki]))
        }
        
        # normalize the new weights
        
        summedWeights <- Nnew
        
        if(weightFactor > 0 && weightFactor <= 1)
            wnew <- wnew / sum(wnew) * summedWeights * weightFactor
        else
            wnew <- wnew / sum(wnew) * summedWeights
        
        if(pkpcparam$verbose >= 2)
        {
            printf("Total number of samples in cluster (with previous samples):\n")
            print(nSamplesInCluster)
        }
        
        if(pkpcparam$verbose >= 2)
        {
            printf("reduced_samples:\n")
            print(reduced_samples)
            printf("wnew:\n")
            print(wnew)
        }
        
        printf("(Normalized or discounted:) Summed weights output = (%f)\n", sum(wnew))
    }
    
    return(list(emptyCluster=0, reduced_samples=reduced_samples, w=wnew, alpha=alpha, f=f, nBestSamples=nBestSamples, bestDists=bestDists))
}

standardizeFeatures.get <- function(data)
{
    D <- ncol(data)
    
    offset <- rep(0,D)
    divisor <- rep(0,D)
    
    for(d in 1:D)
    {
        offset[d] <- - mean(as.numeric(data[,d]))
        divisor[d] <- sd(as.numeric(data[,d]))
        if(divisor[d] == 0)
            divisor[d] <- 1
        
        # zero mean and unit variance: (data + offset) * factor
    }
    
    return(list(offset=offset, divisor=divisor))
}

scaleFeatures.get <- function(data)
{
    D <- ncol(data)
    
    offset <- rep(0,D)
    divisor <- rep(0,D)
    
    for(d in 1:D)
    {
        offset[d] <- - min(data[,d])
        divisor[d] <- max(data[,d]) - min(data[,d])
    }
    
    return(list(offset=offset, divisor=divisor))
}

scaleFeatures.apply <- function(data, offset, divisor)
{
    values <- .C("scale_features", data=as.numeric(as.matrix(data)), N=as.integer(nrow(data)), D=as.integer(ncol(data)), offset=as.numeric(offset), divisor=as.numeric(divisor))
    
    return(matrix(values$data, values$N, values$D))
}

showDataScaling <- function(data)
{
    means <- 0
    for( i in 1:ncol(data) )
    {
        cat( "dimension",i,": min =",min(data[,i]),"max =",max(data[,i]),"mean = ",mean(data[,i]),"variance =",var(data[,i]),"\n" )
        means <- c(means,mean(data[,i]))
    }
    
    means <- means[-1]
    cat( "mean of all means =",mean(means),"\n" )
    cat( "variance of all means = ",var(means),"\n" )
}

kernel.gauss <- function(data, kernelparam, indices_row, indices_col, isDiag, isSymmetric)
{
    # RBF-Kernel
    # already given: sigma, data
    
    if(isDiag)
    {
        # the matrix diagonal is all unity, i.e. 1
        KmatDiag <- rep(1,length(indices_row))
        
        return(KmatDiag)
    }
    
    data_row <- as.matrix(data[indices_row,])
    data_col <- as.matrix(data[indices_col,])
    
    values <- .C("gauss_km_flexible",data_row=as.numeric(data_row), data_col=as.numeric(data_col), D=as.integer(ncol(data_row)), N_row=as.integer(nrow(data_row)), N_col=as.integer(nrow(data_col)), Kmat=as.numeric(matrix(0,nrow(data_row),nrow(data_col))), sigma=as.numeric(kernelparam$sigma), isSymmetric=as.integer(isSymmetric))
    
    return( matrix(values$Kmat,nrow(data_row),nrow(data_col)) )
}

kernel.poly <- function(data, kernelparam, indices_row, indices_col, isDiag, isSymmetric)
{
    # Polynomial Kernel
    
    if(isDiag)
    {
        data_row <- as.matrix(data[indices_row,])
        
        values <- .C("poly_km_flexible",data_row=as.numeric(data_row), data_col=as.numeric(0), D=as.integer(ncol(data_row)), N_row=as.integer(nrow(data_row)), N_col=0, Kmat=as.numeric(matrix(0,nrow(data_row),1)), offset=as.numeric(kernelparam$offset), degree=as.numeric(kernelparam$degree), normalizeKernel=as.integer(kernelparam$normalizeKernel), isDiag=as.integer(isDiag), isSymmetric=as.integer(0))
        
        return(as.numeric(values$Kmat,nrow(data_row)))
    }
    
    data_row <- as.matrix(data[indices_row,])
    data_col <- as.matrix(data[indices_col,])
    
    values <- .C("poly_km_flexible",data_row=as.numeric(data_row), data_col=as.numeric(data_col), D=as.integer(ncol(data_row)), N_row=as.integer(nrow(data_row)), N_col=as.integer(nrow(data_col)), Kmat=as.numeric(matrix(0,nrow(data_row),nrow(data_col))), offset=as.numeric(kernelparam$offset), degree=as.numeric(kernelparam$degree), normalizeKernel=as.integer(kernelparam$normalizeKernel), isDiag=as.integer(isDiag), isSymmetric=as.integer(isSymmetric))
    
    return( matrix(values$Kmat,nrow(data_row),nrow(data_col)) )
}

saveOrRestoreSampleMat <- function(seed, R, N)
{
    # set the seed for the random number generator for drawing the sample indices
    set.seed(seed)
    
    samplesMat <- matrix(0, R, N)
    for(r in 1:R)
        samplesMat[r,] <- sample(1:N)
    
    return(samplesMat)
}

removeRandomLabels <- function(labels, labelFactor)
{
    nonLabelFactor <- 1.0 - labelFactor
    if(nonLabelFactor > 0)
        labels[sample(1:length(labels),nonLabelFactor * length(labels))] <- 0
    
    return(labels)
}

createAndSaveRandomZeroLabels <- function(seed, savfile, labels, labelFactorSeq, R, forceOverwrite = TRUE)
{
    # set the seed for the random number generator for randomly removing the class labellings
    set.seed(seed)
    
    for(labelFactor in labelFactorSeq)
    {
        for(r in 1:R)
        {
            file <- paste(savfile, "_labels_r", r, "labelFactor", labelFactor, sep="")
            if(!file.exists(file) || forceOverwrite)
            {
                tmpLabels <- removeRandomLabels(labels, labelFactor)
                write.table(tmpLabels, file=file)
            }
        }
    }
}

restoreLabels <- function(savfile, r, labelFactor)
{
    file <- paste(savfile, "_labels_r", r, "labelFactor", labelFactor, sep="")
    labels <- as.numeric(t(read.table(file=file)))
    
    return(labels)
}

mydaviesbouldin <- function(nSamplesCluster, summedDistSamplesCluster, distClusterCluster, rmse)
{
    # for the variation of the DB index see 'An extensive comparative study of cluster validity indices'
    
    K <- length(nSamplesCluster)
    
    D <- rep(0,K)
    for(k in 1:K)
    {
        vals <- rep(-Inf,K)
        
        Sk <- summedDistSamplesCluster[k] / nSamplesCluster[k]
        
        for(l in (1:K)[-k])
        {
            Sl <- summedDistSamplesCluster[l] / nSamplesCluster[l]
            d <- distClusterCluster[k,l]  # M_{kl}
            
            if(rmse)
                vals[l] <- (sqrt(Sk) + sqrt(Sl)) / sqrt(d)
            else
                vals[l] <- (Sk + Sl) / d
        }
        
        D[k] <- max(vals[-k])
    }
    
    errDaviesBouldin <- mean(D)
    
    return(errDaviesBouldin)
}

calculateAdjustedRandIndex <- function(contingencyMatrix)
{
    # this function (calculateAdjustedRandIndex) produces
    # the same result as mclust:adjustedRandIndex,
    # i.e. function test => Passed
    
    K <- ncol(contingencyMatrix)
    M <- nrow(contingencyMatrix)
    N <- sum(as.numeric(contingencyMatrix))
    
    A <- 0
    hc <- rep(0, M)
    for(m in 1:M)
    {
        hc[m] <- sum(contingencyMatrix[m, ])
        if(hc[m] >= 2)
        {
            A <- A + choose(hc[m],2)
        }
    }
    
    B <- 0
    hk <- rep(0, K)
    for(k in 1:K)
    {
        hk[k] <- sum(contingencyMatrix[, k])
        if(hk[k] >= 2)
        {
            B <- B + choose(hk[k],2)
        }
    }
    
    Index <- 0
    for(m in 1:M)
    {
        for(k in 1:K)
        {
            if(contingencyMatrix[m,k] >= 2)
            {
                Index <- Index + choose(contingencyMatrix[m,k],2)
            }
        }
    }
    
    val <- (Index - (A * B) / choose(N,2)) / (1/2 * (A + B) - (A * B) / choose(N,2))
    
    return(val)
}

swc.calculateNMI <- function(contingencyMatrix, verbose=0)
{
    # NMI by Strehl & Ghosh
    
    K <- ncol(contingencyMatrix)
    M <- nrow(contingencyMatrix)
    N <- sum(as.numeric(contingencyMatrix))
    
    if(verbose)
    {
        printf("calculateNMI: K (%i) M (%i) N (%i)\n", K, M, N)
    }
    
    hc <- rep(0, M)
    for(m in 1:M)
    {
        hc[m] <- sum(contingencyMatrix[m, ])
    }
    
    hk <- rep(0, K)
    for(k in 1:K)
    {
        hk[k] <- sum(contingencyMatrix[, k])
    }
    
    HC <- 0
    for(m in 1:M)
    {
        if(hc[m] != 0)
            HC <- HC + hc[m] / N * log10(hc[m] / N) 
    }
    HC <- - HC
    
    HK <- 0
    for(k in 1:K)
    {
        if(hk[k] != 0)
            HK <- HK + hk[k] / N * log10(hk[k] / N)
    }
    HK <- - HK
    
    HC_K <- 0
    for(m in 1:M)
    {
        for(k in 1:K)
        {
            if(contingencyMatrix[m, k] != 0)
                HC_K <- HC_K + contingencyMatrix[m, k] / N * log10(contingencyMatrix[m, k] / hk[k])
        }
    }
    HC_K <- - HC_K
    
    if(verbose)
    {
        printf("HC (%f), HK (%f), HC_K (%f)\n", HC, HK, HC_K)
    }
    
    NMI <- (HC - HC_K) / sqrt(HC * HK)
    
    return(NMI)
}

myttest <- function(sampleMean, sampleSd, sampleSize)
{
    # Welch's t-test (unpaired t-test with unequal variances)
    
    t <- (sampleMean[1] - sampleMean[2]) / sqrt(sampleSd[1]^2 / sampleSize[1] + sampleSd[2]^2 / sampleSize[2])
    
    df <- (sampleSd[1]^2 / sampleSize[1] + sampleSd[2]^2 / sampleSize[2])^2 / (sampleSd[1]^4 / (sampleSize[1]^2 * (sampleSize[1] - 1)) + sampleSd[2]^4 / (sampleSize[2]^2 * (sampleSize[2] - 1)))
    
    p <- 2 * pt(-abs(t), df=df)
    
    return(list(t=t, df=df, p.value=p))
}

pkpc.main.test <- function(seed, fileToSave, savefileLabels, a_pkpcparam, a_mu, a_M, a_K, a_L, a_method, data, kernelfunc, kernelparamSeq, labelledUnlabelledFactorSeq, labelFactorSeq, samplesMat, kFoldCrossValidation, T, N, labels, fuzzifier, epsilon, scaleFeatures, maxDataChunksSeq, weightFactorSeq)
{
    minSamplesPerCluster <- 5
    
    if(min(a_M) / 2 < max(a_K) * minSamplesPerCluster)
    {
        cat("Increase a_M or decrease a_K\n")
        return (NULL)
    }
    
    pkpc_results <- NULL
    
    R <- nrow(samplesMat)
    
    if(!is.matrix(kernelparamSeq))
        kernelparamSeq <- t(as.matrix(kernelparamSeq))
    
    if(!is.matrix(a_pkpcparam))
        a_pkpcparam <- t(as.matrix(a_pkpcparam))
    
    for(kernelparamI in 1:nrow(kernelparamSeq))
    {
        kernelparam <- kernelparamSeq[kernelparamI,]
        
        # save data
        dataBak <- data        
        
        for(pkpcparamI in 1:nrow(a_pkpcparam))
        {
            pkpcparam <- a_pkpcparam[pkpcparamI,]
            
            for(mu in a_mu)
            {
                for(M in a_M)
                {
                    for(labelFactor in labelFactorSeq)
                    {
                        for(labelledUnlabelledFactor in labelledUnlabelledFactorSeq)
                        {                        
                            if(labelledUnlabelledFactor != labelledUnlabelledFactorSeq[1] && labelFactor == 0)
                                next
                            
                            for(K in a_K)
                            {
                                for(L in a_L)
                                {
                                    if(L < 0)
                                        L <- M * abs(L)
                                    
                                    for(maxDataChunks in maxDataChunksSeq)
                                    {
                                        for(weightFactor in weightFactorSeq)
                                        {
                                            for(method in a_method)
                                            {
                                                # parameter 'mu = parallelizations' are unused for basic clustering methods
                                                if(method >= 3 && (mu != a_mu[1]))
                                                    next
                                                
                                                # parameter 'pkpcparam' are unused for basic clustering methods
                                                if(method >= 3 && (pkpcparamI != 1))
                                                    next
                                                
                                                # parameter 'L = number of samples for approximate pseudo-centres' are unused for basic clustering methods
                                                if(method >= 3 && (L != a_L[1]))
                                                    next
                                                
                                                # parameter maxDataChunks are unused for basic cluster methods
                                                if(method >= 3 && (maxDataChunks != maxDataChunksSeq[1]))
                                                    next
                                                
                                                # parameter weightFactor are unused for basic cluster methods
                                                if(method >= 3 && (weightFactor != weightFactorSeq[1]))
                                                    next
                                                
                                                if(kFoldCrossValidation <= 0)
                                                    kFoldCrossValidation <- 1
                                                
                                                cat("Setting seed: ", seed, "\n")
                                                
                                                # set the seed for the random number generator
                                                set.seed(seed)
                                                
                                                Rnew <- R * kFoldCrossValidation
                                                
                                                nVals <- 12
                                                
                                                tmp <- foreach(r = 1:Rnew, .combine=cbind) %dorng%
                                                {
                                                    offset <- 0
                                                    vals <- rep(0, nVals)
                                                    
                                                    # s in {1,..., kFoldCrossValidation}
                                                    s <- ((r - 1) %% kFoldCrossValidation) + 1
                                                    
                                                    # r2 in {1,..., R}
                                                    r2 <- floor((r - 1) / kFoldCrossValidation) + 1
                                                    
                                                    if(kFoldCrossValidation > 1)
                                                    {
                                                        # training set has (S - 1) / S of the samples, validation set has the rest
                                                        # for S = 10, training set has 90% of the samples and validation set 10%
                                                        factorTrainingSet <- (kFoldCrossValidation - 1) / kFoldCrossValidation
                                                        factorValidationSet <- 1.0 - factorTrainingSet
                                                        nSamplesValidationSet <- ncol(samplesMat) * factorValidationSet
                                                        
                                                        offset <- (s - 1) * nSamplesValidationSet
                                                        indicesValidationSet <- (1 + offset):(nSamplesValidationSet + offset)
                                                        
                                                        samplesValidate <- samplesMat[r2,indicesValidationSet]
                                                        samplesTrain <- samplesMat[r2,-indicesValidationSet]
                                                    }
                                                    else
                                                    {
                                                        # training set = validation set
                                                        samplesTrain <- samplesMat[r2,]
                                                        
                                                        samplesValidate <- samplesTrain                                    
                                                    }
                                                    
                                                    reducedLabels <- restoreLabels(savefileLabels, r2, labelFactor)
                                                    
                                                    # restore data
                                                    data <- dataBak
                                                    
                                                    # get offset and divisor from the training set
                                                    if(scaleFeatures == 1)
                                                        ret2 <- standardizeFeatures.get(data[samplesTrain,])
                                                    else if(scaleFeatures == 2)
                                                        ret2 <- scaleFeatures.get(data[samplesTrain,])
                                                    
                                                    if(scaleFeatures)
                                                    {
                                                        # apply offset and divisor to training set and to validation set
                                                        data[samplesTrain,] <- scaleFeatures.apply(data[samplesTrain,], ret2$offset, ret2$divisor)
                                                        
                                                        data[samplesValidate,] <- scaleFeatures.apply(data[samplesValidate,], ret2$offset, ret2$divisor)
                                                    }
                                                    
                                                    # remove the samples that will not be clustered (last data chunk) (before or after rescaling the features?)
                                                    
                                                    N <- length(samplesTrain)
                                                    
                                                    if(M > N)
                                                        Mtmp <- N
                                                    else
                                                        Mtmp <- M
                                                    
                                                    Nnew <- Mtmp * floor(N / Mtmp)
                                                    samplesTrain <- samplesTrain[1:Nnew]
                                                    
                                                    startTime <- as.integer(format(Sys.time(),"%s"))
                                                    
                                                    ret1 <- pkpc(pkpcparam, samplesTrain, mu, M, K, L, method, data, kernelfunc, kernelparam, fuzzifier, epsilon, T, max(labels) + 1, reducedLabels, labelledUnlabelledFactor, labelFactor, maxDataChunks, weightFactor)
                                                    
                                                    endTime <- as.integer(format(Sys.time(),"%s"))
                                                    
                                                    if(ret1$emptyCluster)
                                                    {
                                                        vals <- rep(NA, nVals)
                                                        # if one of the values is NA then the mean of these values will also be NA so give up immediately
                                                    }
                                                    else
                                                    {
                                                        ret <- pkpc.calculateIndexes(pkpcparam, samplesValidate, K, data, kernelfunc, kernelparam, labels, ret1$lastSamples, ret1$lastF, ret1$lastW, ret1$lastAlpha)                                                                                                        
                                                        
                                                        if(ret$emptyCluster)
                                                        {
                                                            vals <- rep(NA, nVals)
                                                            # if one of the values is NA then the mean of these values will also be NA so give up immediately
                                                        }                                                                                                        
                                                        else
                                                        {
                                                            vals <- vals + c(ret$errPseudoCentresQuant, 
                                                                             ret$errDaviesBouldinOrgPseudoCentres,
                                                                             ret$errDaviesBouldinOrgRmsePseudoCentres,
                                                                             ret$predPurityPseudoCentres, ret$predNMIPseudoCentres, ret$predAdjustedRandPseudoCentres, 
                                                                             ret1$nSamplesMean, 
                                                                             ret$nsamples, ret$nSamplesClusterPseudoCentresMin, 
                                                                             ret1$bestDistsMean, ret1$bestDistsMin,
                                                                             endTime - startTime)
                                                        }
                                                    }
                                                    
                                                    vals
                                                }
                                                
                                                tmp <- as.matrix(tmp)
                                                
                                                nEmptyCluster <- sum(is.na(tmp[1,]))
                                                if(nEmptyCluster == Rnew)
                                                    onlyEmptyCluster <- 1
                                                else
                                                    onlyEmptyCluster <- 0
                                                
                                                errPseudoCentreQuant <- tmp[1,]
                                                errDaviesBouldinOrgPseudoCentres <- tmp[2,]
                                                errDaviesBouldinOrgRmsePseudoCentres <- tmp[3,]
                                                predPurityPseudoCentres <- tmp[4,]
                                                predNMIPseudoCentres <- tmp[5,]
                                                predAdjustedRandPseudoCentres <- tmp[6,]
                                                nSamplesMean <- tmp[7,]
                                                nSamplesLast <- tmp[8,]
                                                nSamplesClusterPseudoCentresMin <- tmp[9,]
                                                bestDistsMean <- tmp[10,]
                                                bestDistsMin <- tmp[11,]
                                                times <- tmp[12,]
                                                
                                                errPseudoCentreQuant <- errPseudoCentreQuant[!is.na(errPseudoCentreQuant)]
                                                errDaviesBouldinOrgPseudoCentres <- errDaviesBouldinOrgPseudoCentres[!is.na(errDaviesBouldinOrgPseudoCentres)]                                        
                                                errDaviesBouldinOrgRmsePseudoCentres <- errDaviesBouldinOrgRmsePseudoCentres[!is.na(errDaviesBouldinOrgRmsePseudoCentres)]                                        
                                                predPurityPseudoCentres <- predPurityPseudoCentres[!is.na(predPurityPseudoCentres)]
                                                predNMIPseudoCentres <- predNMIPseudoCentres[!is.na(predNMIPseudoCentres)]
                                                predAdjustedRandPseudoCentres <- predAdjustedRandPseudoCentres[!is.na(predAdjustedRandPseudoCentres)]
                                                nSamplesMean <- nSamplesMean[!is.na(nSamplesMean)]
                                                nSamplesLast <- nSamplesLast[!is.na(nSamplesLast)]
                                                nSamplesClusterPseudoCentresMin <- nSamplesClusterPseudoCentresMin[!is.na(nSamplesClusterPseudoCentresMin)]
                                                bestDistsMean <- bestDistsMean[!is.na(bestDistsMean)]
                                                bestDistsMin <- bestDistsMin[!is.na(bestDistsMin)]                                                    
                                                times <- times[!is.na(times)]
                                                
                                                if(onlyEmptyCluster)
                                                {
                                                    cat("Only have detected empty cluster for M =", M, ", K =",K,", L =",L,". Try to adjust the parameters\n")
                                                }
                                                else
                                                {
                                                    errPseudoCentreQuantMean <- mean(errPseudoCentreQuant)
                                                    errPseudoCentreQuantStd <- sd(errPseudoCentreQuant)
                                                    errDaviesBouldinOrgPseudoCentresMean <- mean(errDaviesBouldinOrgPseudoCentres)
                                                    errDaviesBouldinOrgPseudoCentresStd <- sd(errDaviesBouldinOrgPseudoCentres)                                                                                                        
                                                    errDaviesBouldinOrgRmsePseudoCentresMean <- mean(errDaviesBouldinOrgRmsePseudoCentres)
                                                    errDaviesBouldinOrgRmsePseudoCentresStd <- sd(errDaviesBouldinOrgRmsePseudoCentres)                                                                                                        
                                                    predPurityPseudoCentresMean <- mean(predPurityPseudoCentres)
                                                    predPurityPseudoCentresStd <- sd(predPurityPseudoCentres)
                                                    predNMIPseudoCentresMean <- mean(predNMIPseudoCentres)
                                                    predNMIPseudoCentresStd <- sd(predNMIPseudoCentres)
                                                    predAdjustedRandPseudoCentresMean <- mean(predAdjustedRandPseudoCentres)
                                                    predAdjustedRandPseudoCentresStd <- sd(predAdjustedRandPseudoCentres)                                                        
                                                    nSamplesMeanMean <- mean(nSamplesMean)                                                        
                                                    nSamplesLastMean <- mean(nSamplesLast)                                                    
                                                    nSamplesClusterPseudoCentresMinMean <- mean(nSamplesClusterPseudoCentresMin)
                                                    nSamplesClusterPseudoCentresMinStd <- sd(nSamplesClusterPseudoCentresMin)
                                                    nSamplesClusterPseudoCentresMinMin <- min(nSamplesClusterPseudoCentresMin)
                                                    bestDistsMeanMean <- mean(bestDistsMean)
                                                    bestDistsMinMean <- mean(bestDistsMin)
                                                    minTime <- min(times)
                                                    maxTime <- max(times)
                                                    meanTime <- mean(times)
                                                    
                                                    pkpc_results <- rbind( pkpc_results, list(kernelparamI=kernelparamI, pkpcparamI=pkpcparamI, mu=mu,M=M, 
                                                                                              labelFactor=labelFactor, labelledUnlabelledFactor=labelledUnlabelledFactor, K=K, L=L, maxDataChunks=maxDataChunks, 
                                                                                              weightFactor=weightFactor, method=method, nSuccesses=Rnew-nEmptyCluster, nSamplesAPCMean=nSamplesMeanMean, nSamplesLastMean=nSamplesLastMean, 
                                                                                              nSamplesClusterPseudoCentresMinMean=nSamplesClusterPseudoCentresMinMean, nSamplesClusterPseudoCentresMinStd=nSamplesClusterPseudoCentresMinStd, nSamplesClusterPseudoCentresMinMin=nSamplesClusterPseudoCentresMinMin, 
                                                                                              errPseudoCentreQuantMean=errPseudoCentreQuantMean, errPseudoCentreQuantStd=errPseudoCentreQuantStd,
                                                                                              errDaviesBouldinOrgPseudoCentresMean=errDaviesBouldinOrgPseudoCentresMean, errDaviesBouldinOrgPseudoCentresStd=errDaviesBouldinOrgPseudoCentresStd,
                                                                                              errDaviesBouldinOrgRmsePseudoCentresMean=errDaviesBouldinOrgRmsePseudoCentresMean, errDaviesBouldinOrgRmsePseudoCentresStd=errDaviesBouldinOrgRmsePseudoCentresStd,
                                                                                              predPurityPseudoCentresMean=predPurityPseudoCentresMean, predPurityPseudoCentresStd=predPurityPseudoCentresStd,  
                                                                                              predNMIPseudoCentresMean=predNMIPseudoCentresMean, predNMIPseudoCentresStd=predNMIPseudoCentresStd,
                                                                                              predAdjustedRandPseudoCentresMean=predAdjustedRandPseudoCentresMean, predAdjustedRandPseudoCentresStd=predAdjustedRandPseudoCentresStd,
                                                                                              bestDistsMin=bestDistsMinMean, bestDistsMean=bestDistsMeanMean, 
                                                                                              minTime=minTime, maxTime=maxTime, meanTime=meanTime) )
                                                    
                                                    write.table(pkpc_results,file=fileToSave)
                                                }
                                            } # for method
                                        } # for maxDataChunks
                                    } # for weightFactor
                                } # for L
                            } # for K
                        } # for labelledUnlabelledFactor
                    } # for labelFactor
                } # for M
            } # for mu
        } # for pkpcparam
    } # for kernelparam
    
    return( pkpc_results )
}
