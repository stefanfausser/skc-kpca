activity.preprocess <- function(removeOutliers = TRUE, removeFeatures = TRUE)
{
    x <- read.csv2("../samples/weightlifting/dataset-har-PUC-Rio-ugulino.csv")
    
    labels <- as.numeric(x[,19])
    x <- x[,-19]
    
    # remove the name of the user
    x <- x[,-1]
    
    # remove sample with NA attribute
    ind <- which(x[,17] == "-14420-11-2011 04:50:23.713")
    x <- x[-ind,]
    x[,17] <- as.numeric(as.character(x[,17]))
    
    # TODO: What to do with gender? Men = 1, Woman = 2?
    x[,1] <- as.numeric(x[,1])
    
    N <- nrow(x)
    
    return(list(N=N, x=x, labels=labels))
}

gas.preprocess <- function(removeOutliers = TRUE, removeFeatures = TRUE)
{
    # for i in $(ls *.dat); do cat "$i" | sed 's/[0-9]*://g' > "$i"-preprocessed; done
    
    printf("Gas Sensor Array Drift: Read data...")
    
    x <- NULL
    
    for(i in 1:10)
    {
        file <- paste("../samples/gas/batch", i, ".dat-preprocessed", sep="")
        
        x <- rbind(x, read.table(file))
    }
    
    printf("done\n")
    
    labels <- x[,1]
    x <- x[,-1]
    
    N <- nrow(x)
    
    return(list(N=N, x=x, labels=labels))
}

pen.preprocess <- function(removeOutliers = TRUE, removeFeatures = TRUE)
{
    x <- read.table("../samples/pen/pendigits.tra",sep=",")
    x <- rbind(x,read.table("../samples/pen/pendigits.tes",sep=","))
    
    labels <- x[,17] + 1
    x <- x[,-17]
    
    N <- nrow(x)
    
    return(list(N=N, x=x, labels=labels))
}

cardiotocography.preprocess <- function(removeOutliers = TRUE, removeFeatures = TRUE)
{
    # from CTG.xls => saved as CTG.csv
    #     x <- read.csv("CTG.csv",sep="\t",dec=",",skip=1)
    #     x <- x[1:2126,c(11:31, 44, 46)]
    
    # from cleaned-up csv:
    x <- read.csv("../samples/cardiotocography/ctg.csv",sep="\t",dec=",")
    
    # three classes (normal, suspect, pathologic), most are normal
    labels <- x[,23]
    
    x <- x[,-(22:23)]
    
    if(removeFeatures)
    {
        # almost always zero (not important)
        x <- x[,-10]
    }
    
    N <- nrow(x)
    
    return(list(N=N, x=x, labels=labels))
}

miniboone.preprocess <- function(removeOutliers = TRUE, removeFeatures = TRUE)
{
    ## MiniBooNE Dataset (UCI repos)
    
    printf("MiniBooNE: Read data...")
    x <- read.table("../samples/MiniBooNE/MiniBooNE_PID.txt",skip=1,comment.char = "")
    
    labels <- c(rep(1, 36499), rep(2, 93565))
    printf("done\n")
    
    if(removeOutliers)
    {
        # Sort out 'NA' samples (468 exactly)
        ind <- which(x[,1] == -999)
        
        # Sort out outliers (a single outlier)
        ind <- c(ind, which(x[,20] > 100000))
        
        printf("Removing following sample indices:\n")
        print(ind)
        
        x <- x[-ind,]
        labels <- labels[-ind]        
    }
    
    if(removeFeatures)
    {
        # remove the features that correlate to the labels
        
        cors <- rep(0,ncol(x))
        for(i in 1:ncol(x))
            cors[i] <- cor(x[,i],labels)    
        
        ind <- which(abs(cors) > 0.5)
        
        printf("Removing following feature indices:\n")
        print(ind)
        
        x <- x[,-ind]
    }
    
    # Choose first 128000 samples
    N <- 128000
    
    x <- x[1:N,]
    labels <- labels[1:N]
    
    return(list(N=N, x=x, labels=labels))
}

dataset.test <- function(x)
{
    for(d in 1:ncol(x))
    {
        if(is.factor(x[,d]))
            return(FALSE)
    }
    
    return(TRUE)
}
