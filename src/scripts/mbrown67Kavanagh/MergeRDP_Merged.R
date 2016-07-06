## Need to locate the metadata first and look at them
## Okay that is saved out as a .csv and loads into R correctly
## Take Roshonda's output and
## Had to manually remove a " from the class level hierarchy IMPORTANT
rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging")
metaColClasses<- c("character", "character", "numeric", "character", "numeric", "character")
myMetadata <- read.csv("IntAgeMetadata.csv", header=TRUE, colClasses=metaColClasses)
myMetadata <- myMetadata[-1,]
setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging/mergedReadsRDPtables")

taxaLevels <- c("phylum", "class", "order", "family", "genus" )

for (taxa in taxaLevels)
    {
        inFileName <- paste("mergedReads_RDP_",taxa, "_classified.txt", sep = "")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t")
        numCols <- ncol(myT)
        myColClasses <- c("character", rep("numeric", numCols - 1))
        myT <-read.csv(inFileName, header=TRUE, colClasses=myColClasses, sep="\t")

        myT <- t(myT)
##        myT <- as.data.frame(myT)
        colnames(myT)<-myT[1,]
        myT<-myT[-1,]
##        myT<-as.data.frame(myT, stringsAsFactors=FALSE)
##        myColClasses <- c("character", rep("numeric", numCols - 1))
#        myT <-read.csv(inFileName, header=TRUE, colClasses=myColClasses, sep="")
#        myT[is.na(myT)]<-0

##        nSamples <- dim(myT)[1] #N total number of samples
##        sizeSamples <- rowSums(myT[,2:ncol(myT)]) #n number of sequences in a sample
##        totalSize <- sum(myT[,2:ncol(myT)]) #x is total number of counts in the table

## unlist(lapply(strsplit(rownames(myT.transpose), split="_"),"[[",1))

        ## myLogNorm <- log(myT[i,j] / sizeSamples[i] + totalSize/nSamples  + 1)

#        myLogNorm <- myT

#        for(i in 1:dim(myT)[1])
#            {
#                for(j in 2:dim(myT)[2])
#                    {
#                        myLogNorm[i,j]<-log10(myT[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
#                    }
                                        #            }
        ## Break apart on Rownames and reassign.
        Sample.Name<-gsub("\\.", "-", unlist(lapply(strsplit(rownames(myT), split="_"),"[[",1)))
        ## Why does the numeric code change?
        Sample.Name<-gsub("X3754", "3612", Sample.Name)
        ## Then replace "." with "_"
        ## readDirection <- sapply(strsplit(rownames(myT), split="_"),"[[", 4)
        myT<-cbind(myT, Sample.Name)##, readDirection)
        myT<-as.data.frame(myT, stringsAsFactors=FALSE)
        ## Select only R1
        ## myT<-myT[readDirection == "R1",]
        merged<-merge(x = myT, y = myMetadata, by.y = "Sample.Name", all = TRUE)
 #       mergedLN<-merge(x = myLogNorm, y = myMetadata1, by = "MatchFile", all = TRUE)

        write.table(merged, file=paste(taxa,"_RawwithMetadata_Merged.txt", sep=""), row.names = FALSE, sep="\t")
#        write.table(mergedLN, file=paste(taxa,"LogNormwithMetadata_R1_Pooled.txt", sep=""), row.names = FALSE, sep="\t")
    }


