rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

myMetadata1 <- read.csv("NewPooledR2.csv", header=TRUE)

#myMetadata2 <- read.csv("Lyte_Batch04_Run02_repeat-27418842/Lyte_Batch04_Run02_repeat_sample_sheet_WITHSEQ.csv", header=TRUE, comment.char="?")

taxaLevels <- c( "domain", "phylum", "class", "order", "family", "genus" )

for (taxa in taxaLevels)
    {
        inFileName <- paste(taxa, "RawwithMetadata_R2_Pooled.txt", sep = "")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t")
        numCols <- ncol(myT)
#        myColClasses <- c("character", rep("numeric", numCols - 1))
#        myT <-read.csv(inFileName, header=TRUE, colClasses=myColClasses, sep="")
#        myT[is.na(myT)]<-0

#        nSamples <- dim(myT)[1] #N total number of samples
#        sizeSamples <- rowSums(myT[,2:ncol(myT)]) #n number of sequences in a sample
#        totalSize <- sum(myT[,2:ncol(myT)]) #x is total number of counts in the table
                                        # myLogNorm <- log(myT[i,j] / sizeSamples[i] + totalSize/nSamples  + 1)

        # Have to dropped the samples that don't contribute correctly to the merged case.

        removeControls<-c( "C1", "C2", "N1", "N2", "Neg", "Pos")
        myT<-myT[!(myT$Sample_ID %in% removeControls),]
        removetrs<-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")
        myT<-myT[!(myT$Sample_ID %in% removetrs),]
        removeLow<-c("04-55_S32_L001_R2_001")
        myT<-myT[!(myT$MatchFile %in% removeLow),]

        #Have to manually drop these for some reason as they are not dropped previously
        manualDrop <- c("Neg_S40_L001_R2_001", "PCR1Neg_S65_L001_R2_001")
        myT<-myT[!(myT$MatchFile %in% manualDrop),]

        nSamples <- dim(myT)[1] #N total number of samples
        sizeSamples <- rowSums(myT[,2:(ncol(myT) - 17)]) #n number of sequences in a sample
        totalSize <- sum(myT[,2:(ncol(myT) - 17)]) #x is total number of counts in the table


        myLogNorm <- myT
        #Exclude Metadata columns, hence the -17
        for(i in 1:dim(myT)[1])
            {
                for(j in 2:(dim(myT)[2] - 17))
                    {
                        myLogNorm[i,j]<-log10(myT[i,j]/sizeSamples[i] * totalSize/nSamples  + 1)
                    }
            }

#        mergedLN<-merge(x = myLogNorm, y = myMetadata1, by = "MatchFile", all = TRUE)

        write.table(myLogNorm, file=paste(taxa,"LogNormwithMetadata_R2_Pooled_Dropped.txt", sep=""), row.names = FALSE, sep="\t")
    }

