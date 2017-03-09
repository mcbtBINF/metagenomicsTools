rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/Wastewater/")
myMetadata <- read.csv("MetaMasterAllinOne.csv", header=TRUE)

## Arbitrarily drop those with lane = 7 to get rid of duplicate tp 1 samples
myMetadata <- myMetadata[myMetadata$Lane != 7,]
rownames(myMetadata) <- myMetadata$SampleID
## May need to preprocess and read this in if problems exist.

## Makes it specific to the taxa considered
## Therefore, closed, open, denovo, metaphlan and kraken will have different merge
## scripts to do this processing.

## Change this
setwd("/Users/mbrown67/Documents/Fodor/Datasets/Wastewater/qiime/closed_reference_merged/countByTaxonomy/")

taxaLevels <- c(2, 3, 4, 5, 6, 7)
for (taxa in taxaLevels){
    ## Change this
        inFileName <- paste0("merged_otu_table_L",taxa,".txt")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t", skip=1, comment.char="@")

    ## Name the rows and fix the size
    ## This will have to be improved for the QIIME data
        rownames(myT) <- myT[,1]
        myT<- myT[,-1]
        ## tempVal <- rowSums(myT)
        ## myT$countSequences <- tempVal
        ## Rename the columns to match the metadata SampleID field
        ## colnames(myT)[1] <- "SampleID"
        colnames(myT) <- gsub("\\.", "_", colnames(myT))
        colnames(myT) <- gsub("^X", "", colnames(myT))
        ## I will need to pivot this table.
        ## numCols <- ncol(myT)
        ##        myColClasses <- c("character", rep("numeric", numCols - 1))
        ##        myT <-read.csv(inFileName, header=TRUE, colClasses=myColClasses, sep="")
        ##        myT[is.na(myT)]<-0

        ## nSamples <- dim(myT)[1] #N total number of samples
        ## sizeSamples <- rowSums(myT[,2:ncol(myT)]) #n number of sequences in a sample
        ## totalSize <- sum(myT[,2:ncol(myT)]) #x is total number of counts in the table
        myT <- t(myT)
        tempVal <- rowSums(myT)
        myT <- cbind(myT, tempVal)
        colnames(myT)[length(colnames(myT))] <- "sequenceCount"

        ## myT<- myT[-1,]
        ## myT$SampleID <- rownames(myT)

        merged<-merge(x = myT, y = myMetadata, by = "row.names", all = TRUE)
        merged <- merged[, !names(merged) %in% c("SampleID")]
        colnames(merged)[1] <- "SampleID"

        write.table(merged, file=paste("RawwithMetadata_L_",taxa,".txt", sep=""), row.names = FALSE, sep="\t")
    }


