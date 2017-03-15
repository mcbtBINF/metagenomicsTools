rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/metadata/")

myMetadata<-read.csv("CompiledMetadataFINAL.csv", header=TRUE, comment.char="@")
colnames(myMetadata)[1] <- "MatchFile"
myMetadata$MatchFile <- as.character(myMetadata$MatchFile)
## myMetadata$Sex <- unlist(lapply(strsplit(as.character(myMetadata$Experiment..Sample.info), split=" "),"[[",2))
## myMetadata$Experiment..Sample.info <- unlist(lapply(strsplit(as.character(myMetadata$Experiment..Sample.info), split=" "),"[[",1))


setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/")
myT <- read.csv("LyteSharon_r01_cr.txt", sep="\t", header=TRUE, comment.char="@")

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
