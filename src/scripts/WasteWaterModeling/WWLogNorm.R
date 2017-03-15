rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/Wastewater/")
myMetadata <- read.csv("MetaMasterAllinOne.csv", header=TRUE)

## Arbitrarily drop those with lane = 7 to get rid of duplicate tp 1 samples
myMetadata <- myMetadata[myMetadata$Lane != 7,]
rownames(myMetadata) <- myMetadata$SampleID

numMetadataCols <- ncol(myMetadata) ## Includes the countSequences in myT
## Do a better job on the logic of the above derivation
## This will change
setwd("/Users/mbrown67/Documents/Fodor/Datasets/Wastewater/qiime/closed_reference_merged/countByTaxonomy/")

taxaLevels <- c(2, 3, 4, 5, 6, 7)

for (taxa in taxaLevels) {
    ## This will change
        inFileName <- paste0("RawwithMetadata_L_",taxa,".txt")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t")
        numCols <- ncol(myT)
    ## myT <- myT[!(c(60, 235),]
    ## myT <- myT[!60,]
    ## Drop missing rows
    ## Drop 08JUNE2016Run_Sample_63_DSA_2_3
    ## Drop Sample_MNTB_Rep1
     myT <- myT[myT$SampleID != "08JUNE2016Run_Sample_63_DSA_2_3",]
     myT <- myT[myT$SampleID != "Sample_MNTB_Rep1",]


                                        #        myColClasses <- c("character", rep("numeric", numCols - 1))
#        myT <-read.csv(inFileName, header=TRUE, colClasses=myColClasses, sep="")
#        myT[is.na(myT)]<-0

        nSamples <- dim(myT)[1] #N total number of samples
        sizeSamples <- myT$sequenceCount## rowSums(myT[,2:ncol(myT)]) #n number of sequences in a sample
        totalSize <- sum(sizeSamples) #x is total number of counts in the table
                                        # myLogNorm <- log(myT[i,j] / sizeSamples[i] + totalSize/nSamples  + 1)

#        myLogNorm <- myT

#        for(i in 1:dim(myT)[1])
#            {
#                for(j in 2:dim(myT)[2])
#                    {
#                        myLogNorm[i,j]<-log10(myT[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
#                    }
#            }
        myLogNorm <- myT
        #Exclude Metadata columns, hence the -20 (this changed)
        for(i in 1:dim(myT)[1])
            {
                for(j in 2:(dim(myT)[2] - numMetadataCols))
                    {
                        myLogNorm[i,j]<-log10(myT[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
                    }
            }

        ## mergedLN<-merge(x = myLogNorm, y = myMetadata, by = "MatchFile", all = TRUE)

        write.table(myLogNorm, file=paste("LogNormwithMetadata_L_",taxa,".txt", sep=""), row.names = FALSE, sep="\t")
    }


