rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

myMetadata1 <- read.csv("NewPooled.csv", header=TRUE)

#myMetadata2 <- read.csv("Lyte_Batch04_Run02_repeat-27418842/Lyte_Batch04_Run02_repeat_sample_sheet_WITHSEQ.csv", header=TRUE, comment.char="?")

taxaLevels <- c( "domain", "phylum", "class", "order", "family", "genus" )

for (taxa in taxaLevels)
    {
        inFileName <- paste(taxa, "RawwithMetadata_R1_Pooled.txt", sep = "")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t")
        numCols <- ncol(myT)
        numMetadataCols <- 20
#        myColClasses <- c("character", rep("numeric", numCols - 1))
#        myT <-read.csv(inFileName, header=TRUE, colClasses=myColClasses, sep="")
#        myT[is.na(myT)]<-0

#        nSamples <- dim(myT)[1] #N total number of samples
#        sizeSamples <- rowSums(myT[,2:ncol(myT)]) #n number of sequences in a sample
#        totalSize <- sum(myT[,2:ncol(myT)]) #x is total number of counts in the table
                                        # myLogNorm <- log(myT[i,j] / sizeSamples[i] + totalSize/nSamples  + 1)

        # Have to dropped the samples that don't contribute correctly to the merged case.
        # Dr. Fodor doesn't seem to be too worried about the way that this is currently done, so I am commenting out this region to reflect that.
        ## removeControls<-c( "C1", "C2", "N1", "N2", "Neg", "Pos")
        ## myT<-myT[!(myT$Sample_ID %in% removeControls),]
        ## removetrs<-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")
        ## myT<-myT[!(myT$Sample_ID %in% removetrs),]
        ## removeLow<-c("04-55_S32_L001_R1_001")
        ## myT<-myT[!(myT$MatchFile %in% removeLow),]

        ## #Have to manually drop these for some reason as they are not dropped previously
        ## manualDrop <- c("Neg_S40_L001_R1_001", "PCR1Neg_S65_L001_R1_001")
        ## myT<-myT[!(myT$MatchFile %in% manualDrop),]

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

