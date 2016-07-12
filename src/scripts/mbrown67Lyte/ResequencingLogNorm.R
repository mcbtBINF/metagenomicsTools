rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Resequencing/ForwardReads/")

myMetadata1 <- read.csv("../ResequencingMetadata.csv", header=TRUE)

taxaLevels <- c("phylum", "class", "order", "family", "genus" )

for (taxa in taxaLevels)
    {
        inFileName <- paste(taxa, "RawwithMetadata_R1_Pooled.txt", sep = "")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t")
        numCols <- ncol(myT)
        numMetadataCols <- 17

        nSamples <- dim(myT)[1] #N total number of samples
        sizeSamples <- rowSums(myT[,2:(ncol(myT) - numMetadataCols)]) #n number of sequences in a sample
        totalSize <- sum(myT[,2:(ncol(myT) - numMetadataCols)]) #x is total number of counts in the table


        myLogNorm <- myT
        #Exclude Metadata columns, hence the -20 (this changed)
        for(i in 1:dim(myT)[1])
            {
                for(j in 2:(dim(myT)[2] - numMetadataCols))
                    {
                        myLogNorm[i,j]<-log10(myT[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
                    }
            }

#        mergedLN<-merge(x = myLogNorm, y = myMetadata1, by = "MatchFile", all = TRUE)

        write.table(myLogNorm, file=paste(taxa,"LogNormwithMetadata_R1_Pooled.txt", sep=""), row.names = FALSE, sep="\t")
    }

