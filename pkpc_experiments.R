source("pkpc.R")
source("pkpc_datasets.R")
source("pkpc_test.R")

pkpc.test.cardiotocography.start <- function(basic = TRUE, semisupervised = FALSE, evalWeight = FALSE, evalLabelledUnlabelledFactor = FALSE)
{
    # Basic parameter
    MBasicSeq <- c(530, 1060) # Basic
    K <- 2
    fuzzifier <- 1.025
    degree <- 3
    offset <- 0
    sigma <- 0
    normalizeKernel <- 1
    
    # KPC-A
    M <- 530 # KPC-A
    maxDataChunks <- 0
    a_L <- c(50, 150, 250)
    weightFactor <- 0.5
    # KPC-A Eval
    maxDataChunksEval <- 0
    maxDataChunksEvalSeq <- 0
    weightFactorEvalSeq <- c(0.3, 0.5, 0.7, 1.0)
    
    # SKC
    M_semisupervised <- 1060 # SKC
    labelledUnlabelledFactor <- 2
    labelFactorSeq <- c(0.05, 0.1, 0.3)
    # SKC Eval
    labelFactorEvalSeq <- 0.3
    labelledUnlabelledFactorEvalSeq <- c(0, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
    
    param <- list(seedSamples = 1, seedLabels = 2, seed = 3, dataset = "cardiotocography", datasetfunc = cardiotocography.preprocess, MBasicSeq = MBasicSeq, M = M, M_semisupervised = M_semisupervised, K = K, fuzzifier = fuzzifier, degree = degree, offset = offset, sigma = sigma, normalizeKernel = normalizeKernel, maxDataChunks = maxDataChunks, maxDataChunksEvalSeq = maxDataChunksEvalSeq, a_L = a_L, labelledUnlabelledFactor = labelledUnlabelledFactor, labelFactorSeq = labelFactorSeq, labelledUnlabelledFactorEvalSeq = labelledUnlabelledFactorEvalSeq, weightFactorEvalSeq = weightFactorEvalSeq, labelFactorEvalSeq = labelFactorEvalSeq, weightFactor = weightFactor, maxDataChunksEval = maxDataChunksEval)
    
    if(!semisupervised)
        pkpc.test.unsupervised.framework(param, basic = basic, evalWeight = FALSE, evalWeightTwoDataChunks = evalWeight, evalMaxDataChunks = FALSE)
    else
        pkpc.test.semisupervised.framework(param, evalLabelledUnlabelledFactor = evalLabelledUnlabelledFactor)
}

pkpc.test.gas.start <- function(basic = TRUE, semisupervised = FALSE, evalWeight = FALSE, evalWeightTwoDataChunks = FALSE, evalMaxDataChunks = FALSE, evalLabelledUnlabelledFactor = FALSE)
{
    # Basic parameter
    MBasicSeq <- c(500, 1000, 1500) # Basic
    K <- 8
    fuzzifier <- 1.025
    degree <- 3
    offset <- 0
    sigma <- 0
    normalizeKernel <- 1
    
    # KPC-A
    M <- 500 # KPC-A
    maxDataChunks <- 0
    a_L <- c(50, 150, 250)
    weightFactor <- 0.5
    # KPC-A Eval
    maxDataChunksEvalSeq <- c(5, 10, 15, 20)
    weightFactorEvalSeq <- c(0.3, 0.5, 0.7, 1.0)
    
    # SKC
    M_semisupervised <- 1000 # SKC
    labelledUnlabelledFactor <- 2
    labelFactorSeq <- c(0.05, 0.1, 0.3)
    # SKC Eval
    maxDataChunksEval <- 0
    labelFactorEvalSeq <- 0.3
    labelledUnlabelledFactorEvalSeq <- c(0, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
    
    param <- list(seedSamples = 1, seedLabels = 2, seed = 3, dataset = "gas", datasetfunc = gas.preprocess, MBasicSeq = MBasicSeq, M = M, M_semisupervised = M_semisupervised, K = K, fuzzifier = fuzzifier, degree = degree, offset = offset, sigma = sigma, normalizeKernel = normalizeKernel, maxDataChunks = maxDataChunks, maxDataChunksEvalSeq = maxDataChunksEvalSeq, a_L = a_L, labelledUnlabelledFactor = labelledUnlabelledFactor, labelFactorSeq = labelFactorSeq, labelledUnlabelledFactorEvalSeq = labelledUnlabelledFactorEvalSeq, weightFactorEvalSeq = weightFactorEvalSeq, labelFactorEvalSeq = labelFactorEvalSeq, weightFactor = weightFactor, maxDataChunksEval = maxDataChunksEval)
    
    if(!semisupervised)
        pkpc.test.unsupervised.framework(param, basic = basic, evalWeight = evalWeight, evalWeightTwoDataChunks = evalWeightTwoDataChunks, evalMaxDataChunks = evalMaxDataChunks)
    else
        pkpc.test.semisupervised.framework(param, evalLabelledUnlabelledFactor = evalLabelledUnlabelledFactor)
}

pkpc.test.pen.start <- function(basic = TRUE, semisupervised = FALSE, evalWeight = FALSE, evalWeightTwoDataChunks = FALSE, evalMaxDataChunks = FALSE, evalLabelledUnlabelledFactor = FALSE)
{
    # Basic parameter
    MBasicSeq <- c(500, 1000, 1500) # Basic
    K <- 10
    fuzzifier <- 1.1
    degree <- 0
    offset <- 0
    sigma <- 3
    normalizeKernel <- 0
    
    # KPC-A
    M <- 500 # KPC-A
    maxDataChunks <- 0
    a_L <- c(50, 150, 250)
    weightFactor <- 0.5
    # KPC-A Eval
    maxDataChunksEval <- 0
    maxDataChunksEvalSeq <- c(5, 10, 15, 20)
    weightFactorEvalSeq <- c(0.3, 0.5, 0.7, 1.0)
    
    # SKC
    M_semisupervised <- 1000 # SKC
    labelledUnlabelledFactor <- 2
    labelFactorSeq <- c(0.05, 0.1, 0.3)
    # SKC Eval
    labelFactorEvalSeq <- 0.3
    labelledUnlabelledFactorEvalSeq <- c(0, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
    
    param <- list(seedSamples = 1, seedLabels = 2, seed = 3, dataset = "pen", datasetfunc = pen.preprocess, MBasicSeq = MBasicSeq, M = M, M_semisupervised = M_semisupervised, K = K, fuzzifier = fuzzifier, degree = degree, offset = offset, sigma = sigma, normalizeKernel = normalizeKernel, maxDataChunks = maxDataChunks, maxDataChunksEvalSeq = maxDataChunksEvalSeq, a_L = a_L, labelledUnlabelledFactor = labelledUnlabelledFactor, labelFactorSeq = labelFactorSeq, labelledUnlabelledFactorEvalSeq = labelledUnlabelledFactorEvalSeq, weightFactorEvalSeq = weightFactorEvalSeq, labelFactorEvalSeq = labelFactorEvalSeq, weightFactor = weightFactor, maxDataChunksEval = maxDataChunksEval)
    
    if(!semisupervised)
        pkpc.test.unsupervised.framework(param, basic = basic, evalWeight = evalWeight, evalWeightTwoDataChunks = evalWeightTwoDataChunks, evalMaxDataChunks = evalMaxDataChunks)
    else
        pkpc.test.semisupervised.framework(param, evalLabelledUnlabelledFactor = evalLabelledUnlabelledFactor)
}

pkpc.test.miniboone.start <- function(basic = TRUE, semisupervised = FALSE, evalWeight = FALSE, evalWeightTwoDataChunks = FALSE, evalMaxDataChunks = FALSE, evalLabelledUnlabelledFactor = FALSE)
{
    # Basic parameter
    MBasicSeq <- c(500, 1000, 1500) # Basic
    K <- 2
    fuzzifier <- 1.0025
    degree <- 3
    offset <- 0
    sigma <- 0
    normalizeKernel <- 1
    
    # KPC-A
    maxDataChunks <- 50
    weightFactor <- 1
    M <- 500 # KPC-A
    a_L <- c(50, 150, 250)
    
    # KPC-A Eval
    maxDataChunksEval <- 50 # limit weightFactor-evaluation to 50 data chunks
    maxDataChunksEvalSeq <- c(10, 30, 50, 100)
    weightFactorEvalSeq <- c(0.3, 0.5, 0.7, 1.0)
    
    # SKC
    M_semisupervised <- 1000 # SKC
    labelledUnlabelledFactor <- 2
    labelFactorSeq <- c(0.05, 0.1, 0.3)
    # SKC Eval
    labelFactorEvalSeq <- 0.3
    labelledUnlabelledFactorEvalSeq <- c(0, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
    
    param <- list(seedSamples = 1, seedLabels = 2, seed = 3, dataset = "miniboone", datasetfunc = miniboone.preprocess, MBasicSeq = MBasicSeq, M = M, M_semisupervised = M_semisupervised, K = K, fuzzifier = fuzzifier, degree = degree, offset = offset, sigma = sigma, normalizeKernel = normalizeKernel, maxDataChunks = maxDataChunks, maxDataChunksEvalSeq = maxDataChunksEvalSeq, a_L = a_L, labelledUnlabelledFactor = labelledUnlabelledFactor, labelFactorSeq = labelFactorSeq, labelledUnlabelledFactorEvalSeq = labelledUnlabelledFactorEvalSeq, weightFactorEvalSeq = weightFactorEvalSeq, labelFactorEvalSeq = labelFactorEvalSeq, weightFactor = weightFactor, maxDataChunksEval = maxDataChunksEval)
    
    if(!semisupervised)
        pkpc.test.unsupervised.framework(param, basic = basic, evalWeight = evalWeight, evalWeightTwoDataChunks = evalWeightTwoDataChunks, evalMaxDataChunks = evalMaxDataChunks)
    else
        pkpc.test.semisupervised.framework(param, evalLabelledUnlabelledFactor = evalLabelledUnlabelledFactor)
}

pkpc.test.activity.start <- function(basic = TRUE, semisupervised = FALSE, evalWeight = FALSE, evalWeightTwoDataChunks = FALSE, evalMaxDataChunks = FALSE, evalLabelledUnlabelledFactor = FALSE)
{
    # Basic parameter
    MBasicSeq <- c(500, 1000, 1500) # Basic
    K <- 9
    fuzzifier <- 1.025
    degree <- 0
    offset <- 0
    sigma <- 1
    normalizeKernel <- 0
    
    # KPC-A
    M <- 500 # KPC-A
    maxDataChunks <- 50
    a_L <- c(50, 150, 250)
    weightFactor <- 0.5
    # KPC-A Eval
    maxDataChunksEval <- 50 # limit weightFactor-evaluation to 50 data chunks
    maxDataChunksEvalSeq <- c(10, 30, 50, 100)
    weightFactorEvalSeq <- c(0.3, 0.5, 0.7, 1.0)
    
    # SKC
    M_semisupervised <- 1000 # SKC
    labelledUnlabelledFactor <- 2
    labelFactorSeq <- c(0.05, 0.1, 0.3)
    # SKC Eval
    labelFactorEvalSeq <- 0.3
    labelledUnlabelledFactorEvalSeq <- c(0, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
    
    param <- list(seedSamples = 1, seedLabels = 2, seed = 3, dataset = "activity", datasetfunc = activity.preprocess, MBasicSeq = MBasicSeq, M = M, M_semisupervised = M_semisupervised, K = K, fuzzifier = fuzzifier, degree = degree, offset = offset, sigma = sigma, normalizeKernel = normalizeKernel, maxDataChunks = maxDataChunks, maxDataChunksEvalSeq = maxDataChunksEvalSeq, a_L = a_L, labelledUnlabelledFactor = labelledUnlabelledFactor, labelFactorSeq = labelFactorSeq, labelledUnlabelledFactorEvalSeq = labelledUnlabelledFactorEvalSeq, weightFactorEvalSeq = weightFactorEvalSeq, labelFactorEvalSeq = labelFactorEvalSeq, weightFactor = weightFactor, maxDataChunksEval = maxDataChunksEval)
    
    if(!semisupervised)
        pkpc.test.unsupervised.framework(param, basic = basic, evalWeight = evalWeight, evalWeightTwoDataChunks = evalWeightTwoDataChunks, evalMaxDataChunks = evalMaxDataChunks)
    else
        pkpc.test.semisupervised.framework(param, evalLabelledUnlabelledFactor = evalLabelledUnlabelledFactor)
}

pkpc.test.unsupervised <- function(datasetname, eval = 0, basic = FALSE, basicOnly = FALSE)
{
    if(datasetname == "cardiotocography")
    {
        if(basic || basicOnly)
            pkpc.test.cardiotocography.start() # standard KKM, KFCM, RNG
        
        if(basicOnly)
            return(NULL)
        
        if(eval)
            pkpc.test.cardiotocography.start(basic = FALSE, evalWeight = TRUE)
        else
            pkpc.test.cardiotocography.start(basic = FALSE)  # needs weightFactor, comes from evalWeight = TRUE
    }
    else if(datasetname == "gas")
    {
        if(basic || basicOnly)
            pkpc.test.gas.start() # standard KKM, KFCM, RNG
        
        if(basicOnly)
            return(NULL)
        
        if(eval == 1)
            pkpc.test.gas.start(basic = FALSE, evalWeight = TRUE)
        else if(eval == 2)
            pkpc.test.gas.start(basic = FALSE, evalMaxDataChunks = TRUE) # needs weightFactor, comes from evalWeight = TRUE
        else
            pkpc.test.gas.start(basic = FALSE)
    }
    else if(datasetname == "pen")
    {
        if(basic || basicOnly)
            pkpc.test.pen.start() # standard KKM, KFCM, RNG
        
        if(basicOnly)
            return(NULL)
        
        if(eval == 1)
        {
            pkpc.test.pen.start(basic = FALSE, evalWeight = TRUE)
            pkpc.test.pen.start(basic = FALSE, evalWeightTwoDataChunks = TRUE)
        }
        else if(eval == 2)
            pkpc.test.pen.start(basic = FALSE, evalMaxDataChunks = TRUE) # needs weightFactor, comes from evalWeight = TRUE
        else
            pkpc.test.pen.start(basic = FALSE)
    }
    else if(datasetname == "miniboone")
    {
        if(basic || basicOnly)
            pkpc.test.miniboone.start() # standard KKM, KFCM, RNG
        
        if(basicOnly)
            return(NULL)
        
        if(eval == 1)
        {
            pkpc.test.miniboone.start(basic = FALSE, evalWeight = TRUE)
            pkpc.test.miniboone.start(basic = FALSE, evalWeightTwoDataChunks = TRUE)
        }
        else if(eval == 2)
            pkpc.test.miniboone.start(basic = FALSE, evalMaxDataChunks = TRUE) # needs weightFactor, comes from evalWeight = TRUE
        else
            pkpc.test.miniboone.start(basic = FALSE)
    }
    else if(datasetname == "activity")
    {
        if(basic || basicOnly)
            pkpc.test.activity.start() # standard KKM, KFCM, RNG
        
        if(basicOnly)
            return(NULL)
        
        if(eval == 1)
            pkpc.test.activity.start(basic = FALSE, evalWeight = TRUE)
        else if(eval == 2)
            pkpc.test.activity.start(basic = FALSE, evalMaxDataChunks = TRUE) # needs weightFactor, comes from evalWeight = TRUE
        else
            pkpc.test.activity.start(basic = FALSE)
    }
}

pkpc.test.semisupervised <- function(datasetname, eval = FALSE)
{
    if(datasetname == "cardiotocography")
    {
        if(eval)
            pkpc.test.cardiotocography.start(semisupervised = TRUE, evalLabelledUnlabelledFactor = TRUE)
        else
            pkpc.test.cardiotocography.start(semisupervised = TRUE) # needs labelledUnlabelledFactor
    }
    else if(datasetname == "gas")
    {
        if(eval)
            pkpc.test.gas.start(semisupervised = TRUE, evalLabelledUnlabelledFactor = TRUE)
        else
            pkpc.test.gas.start(semisupervised = TRUE)
    }
    else if(datasetname == "pen")
    {
        if(eval)
            pkpc.test.pen.start(semisupervised = TRUE, evalLabelledUnlabelledFactor = TRUE)
        else
            pkpc.test.pen.start(semisupervised = TRUE)
    }
    else if(datasetname == "miniboone")
    {
        if(eval)
            pkpc.test.miniboone.start(semisupervised = TRUE, evalLabelledUnlabelledFactor = TRUE)
        else
            pkpc.test.miniboone.start(semisupervised = TRUE)
    }
    else if(datasetname == "activity")
    {
        if(eval)
            pkpc.test.activity.start(semisupervised = TRUE, evalLabelledUnlabelledFactor = TRUE)
        else
            pkpc.test.activity.start(semisupervised = TRUE)
    }
}
