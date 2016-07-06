## Need to locate the metadata first and look at them
## Okay that is saved out as a .csv and loads into R correctly
## Take Roshonda's output and
## Had to manually remove a " from the class level hierarchy IMPORTANT
rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging/otus")

taxaLevels <- c("otu") ##c("phylum", "class", "order", "family", "genus" )

for (taxa in taxaLevels)
    {
        inFileName <- paste(taxa, "_RawwithMetadata_otu.txt", sep = "")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t")
        numCols <- ncol(myT)
        numMetadataCols <- 6
##        myT<-as.data.frame(myT, stringsAsFactors=FALSE)
##        myColClasses <- c("character", rep("numeric", numCols - 1))
#        myT <-read.csv(inFileName, header=TRUE, colClasses=myColClasses, sep="")
#        myT[is.na(myT)]<-0

        nSamples <- dim(myT)[1] #N total number of samples
        sizeSamples <- rowSums(myT[,2:(ncol(myT) - numMetadataCols)]) #n number of sequences in a sample
        totalSize <- sum(myT[,2:(ncol(myT) - numMetadataCols)]) #x is total number of counts in the table

        ## myLogNorm <- log(myT[i,j] / sizeSamples[i] + totalSize/nSamples  + 1)

        myLogNorm <- myT

        for(i in 1:dim(myT)[1])
            {
                for(j in 2:(dim(myT)[2] - numMetadataCols))
                    {
                        myLogNorm[i,j]<-log10(myT[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
                    }
        }

        write.table(myLogNorm, file=paste(taxa,"_LogNormwithMetadata_otu.txt", sep=""), row.names = FALSE, sep="\t")
    }


