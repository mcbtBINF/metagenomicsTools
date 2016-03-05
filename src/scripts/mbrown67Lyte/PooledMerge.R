rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

myMetadata1 <- read.csv("UpdatedPooled.csv", header=TRUE)

#myMetadata2 <- read.csv("Lyte_Batch04_Run02_repeat-27418842/Lyte_Batch04_Run02_repeat_sample_sheet_WITHSEQ.csv", header=TRUE, comment.char="?")

taxaLevels <- c( "domain", "phylum", "class", "order", "family", "genus" )

for (taxa in taxaLevels)
    {
        inFileName <- paste(taxa, "File_R1.txt", sep = "")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t")
        numCols <- ncol(myT)
#        myColClasses <- c("character", rep("numeric", numCols - 1))
#        myT <-read.csv(inFileName, header=TRUE, colClasses=myColClasses, sep="")
#        myT[is.na(myT)]<-0

        nSamples <- dim(myT)[1] #N total number of samples
        sizeSamples <- rowSums(myT[,2:ncol(myT)]) #n number of sequences in a sample
        totalSize <- sum(myT[,2:ncol(myT)]) #x is total number of counts in the table
                                        # myLogNorm <- log(myT[i,j] / sizeSamples[i] + totalSize/nSamples  + 1)

#        myLogNorm <- myT

#        for(i in 1:dim(myT)[1])
#            {
#                for(j in 2:dim(myT)[2])
#                    {
#                        myLogNorm[i,j]<-log10(myT[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
#                    }
#            }

        merged<-merge(x = myT, y = myMetadata1, by = "MatchFile", all = TRUE)
 #       mergedLN<-merge(x = myLogNorm, y = myMetadata1, by = "MatchFile", all = TRUE)

        write.table(merged, file=paste(taxa,"RawwithMetadata_R1_Pooled.txt", sep=""), row.names = FALSE, sep="\t")
#        write.table(mergedLN, file=paste(taxa,"LogNormwithMetadata_R1_Pooled.txt", sep=""), row.names = FALSE, sep="\t")
    }

