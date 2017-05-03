## Performs baseline kernel clustering (Kernel K-means, Kernel Fuzzy C-means, Relational Neural Gas) for the evaluation of the clustering methods and kernel functions

source("pkpc.R")
source("pkpc_datasets.R")
source("pkpc_eval.R")

pkpc.start.eval.activity <- function()
{
    M <- 1000    
    
    a_sigma <- c(0.1, 0.5, 1, 3, 5, 7, 9)
    a_K <- getNumberClusterSeq(nClasses = 5)
    
    pkpc.start.eval(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "activity", datasetfunc = activity.preprocess, M = M, method = 5, a_K = a_K, a_sigma = a_sigma)    
}

pkpc.start.eval.activity.fuzzifier <- function()
{
    M <- c(500,1000,1500)
    a_K <- 9
    
    fuzzifierSeq <- c(1.025, 1.05, 1.1, 1.25)
    
    for(fuzzifier in fuzzifierSeq)
    {
        pkpc.eval.func(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "activity", datasetfunc = activity.preprocess, M = M, a_K = a_K, degree = 0, offset = 0, sigma = 1, normalizeKernel = 0, a_method = 4, fuzzifier = fuzzifier)    
    }    
}

pkpc.start.eval.pen <- function()
{
    M <- 1000    
    
    a_sigma <- c(0.1, 0.5, 1, 3, 5, 7, 9)
    a_K <- getNumberClusterSeq(nClasses = 10)
    
    pkpc.start.eval(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "pen", datasetfunc = pen.preprocess, M = M, method = 5, a_K = a_K, a_sigma = a_sigma)    
}

pkpc.start.eval.pen.fuzzifier <- function()
{
    M <- c(500,1000,1500)
    a_K <- 10
    
    fuzzifierSeq <- c(1.025, 1.05, 1.1, 1.25)
    
    for(fuzzifier in fuzzifierSeq)
    {
        pkpc.eval.func(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "pen", datasetfunc = pen.preprocess, M = M, a_K = a_K, degree = 0, offset = 0, sigma = 3, normalizeKernel = 0, a_method = 4, fuzzifier = fuzzifier)    
    }    
}

pkpc.start.eval.gas <- function()
{
    M <- 1000
    
    a_sigma <- c(0.1, 0.5, 1, 3, 5, 7, 9)
    a_K <- getNumberClusterSeq(nClasses = 6)
    
    pkpc.start.eval(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "gas", datasetfunc = gas.preprocess, M = M, method = 5, a_K = a_K, a_sigma = a_sigma)
}

pkpc.start.eval.gas.fuzzifier <- function()
{
    M <- c(500, 1000, 1500)
    a_K <- 8
    
    fuzzifierSeq <- c(1.025, 1.05, 1.1, 1.25)
    
    for(fuzzifier in fuzzifierSeq)
    {
        pkpc.eval.func(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "gas", datasetfunc = gas.preprocess, M = M, a_K = a_K, degree = 3, offset = 0, sigma = 0, normalizeKernel = 1, a_method = 4, fuzzifier = fuzzifier)    
    }
}

pkpc.start.eval.cardiotocography <- function()
{
    M <- 530
    
    a_sigma <- c(0.1, 0.5, 1, 3, 5, 7, 9)
    a_K <- getNumberClusterSeq(nClasses = 3)
    
    pkpc.start.eval(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "cardiotocography", datasetfunc = cardiotocography.preprocess, M = M, method = 5, a_K = a_K, a_sigma = a_sigma)    
}

pkpc.start.eval.cardiotocography.fuzzifier <- function()
{
    M <- 530
    a_K <- 2
    
    fuzzifierSeq <- c(1.025, 1.05, 1.1, 1.25)
    
    for(fuzzifier in fuzzifierSeq)
    {
        pkpc.eval.func(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "cardiotocography", datasetfunc = cardiotocography.preprocess, M = M, a_K = a_K, degree = 3, offset = 0, sigma = 0, normalizeKernel = 1, a_method = 4, fuzzifier = fuzzifier)    
    }    
}

pkpc.start.eval.miniboone <- function()
{
    M <- 1000
    
    a_sigma <- c(0.1, 0.5, 1, 3, 5, 7, 9)
    a_K <- getNumberClusterSeq(nClasses = 2)
    
    pkpc.start.eval(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "miniboone", datasetfunc = miniboone.preprocess, M = M, method = 5, a_K = a_K, a_sigma = a_sigma)    
}

pkpc.start.eval.miniboone.fuzzifier <- function()
{
    M <- c(500,1000,1500)
    a_K <- 2
    
    fuzzifierSeq <- c(1.0025, 1.01, 1.025, 1.05, 1.1, 1.25)
    
    for(fuzzifier in fuzzifierSeq)
    {
        pkpc.eval.func(seedSamples = 1, seedLabels = 2, seed = 3, datasetname = "miniboone", datasetfunc = miniboone.preprocess, M = M, a_K = a_K, degree = 3, offset = 0, sigma = 0, normalizeKernel = 1, a_method = 4, fuzzifier = fuzzifier)    
    }    
}
