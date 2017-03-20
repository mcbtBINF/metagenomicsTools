rm(list=ls())

## Make sure that the Excel export didn't hurt the decimals
setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/metadata/")

myColClasses <- c(rep("character", 3), "numeric", rep("character", 4), "numeric", "character", "numeric", rep("character", 8), rep("numeric", 39))
myMetadata<-read.table("CompiledMetadataFINAL.csv", header=TRUE, comment.char="@", sep=",", colClasses = myColClasses)

colnames(myMetadata)[1] <- "MatchFile"
myMetadata$Cage <- paste(myMetadata$Sex, myMetadata$Housing)

## myMetadata$MatchFile <- as.character(myMetadata$MatchFile)


## Read in the columnclasses carefully.
setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/")

## Will later be able to go through this and work on the other datasets.
## This code can be reused to get through the qiime taxonomic levels.
myT <- read.table("LyteSharon_r01_cr.txt", sep="\t", header=TRUE, comment.char="@")

otutoTaxa <- cbind(myT$OTUID, myT$taxonomy)

myT$taxonomy <- NULL

## pivot
myT <- as.data.frame(t(myT))
colnames(myT) <- myT[1,]
myT <- myT[-1,]
## counts
myMetadata$counts <- rowSums(myT)
numMetadataCols <- dim(myMetadata)[2]
nSamples <- dim(myT)[1] ## N total number of samples
sizeSamples <- rowSums(myT)
## number of sequences in a sample
totalSize <- sum(sizeSamples) ## is total number of counts in the table

## combine
myLogNorm <- myT[,]

for(i in 1:dim(myLogNorm)[1])
{
    for(j in 1:dim(myLogNorm)[2])
    {
        myLogNorm[i,j]<-log10(myLogNorm[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
    }
}

myT <- myLogNorm
myT$MatchFile <- as.character(rownames(myT))
## myT$MatchFile <- as.factor(myT$MatchFile)
## names(myT) <- colnames(myT)
merged<-merge(x = myMetadata, y = myT, by = "MatchFile", all=TRUE)

## output to txt file
write.table(merged, file=paste("qiimeOTULogNormwithMetadata.txt", sep=""), row.names = FALSE, sep="\t")
