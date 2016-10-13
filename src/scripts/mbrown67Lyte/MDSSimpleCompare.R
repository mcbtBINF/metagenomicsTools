rm(list=ls())
library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" )

checkSizes <- read.csv("phylumRawwithMetadata_R1_Pooled.txt", header=TRUE, sep="", na.strings="BLAH")

checkSizesR2 <- read.csv("phylumRawwithMetadata_R2_Pooled.txt", header=TRUE, sep="", na.strings="BLAH")

simpleCounter <- read.csv("LytePooled.simple.counter", header=FALSE, sep="")
justR1 <- simpleCounter[grep("_R1_", simpleCounter[,1]),]
justR2 <- simpleCounter[grep("_R2_", simpleCounter[,1]),]
numMetadataCols <- 20

sampleSizes<-rowSums(checkSizes[,2:(ncol(checkSizes)-numMetadataCols)])
sampleSizesR2<-rowSums(checkSizesR2[,2:(ncol(checkSizesR2)-numMetadataCols)])
newCheck <- cbind(checkSizes, sampleSizes)
newCheckR2 <- cbind(checkSizesR2, sampleSizesR2)
reNewCheck<-newCheck[,c(1,dim(newCheck)[2])]
reNewCheckR2<-newCheckR2[,c(1, dim(newCheckR2)[2])]
testR1 <- merge(reNewCheck, justR1, by.x=1, by.y=1)
testR2 <- merge(reNewCheckR2, justR2, by.x=1, by.y=1)
# Try for 5000 and 10000 cutoffs


### Has 118 rows now...
for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
    SeqControl <- c("DNAC1_S69_L001_R1_001", "DNAC2_S70_L001_R1_001", "DNAN1_S71_L001_R1_001", "DNAN2_S72_L001_R1_001", "Neg_S40_L001_R1_001",
                    "PCR1Neg_S65_L001_R1_001", "PCR2Pos_S66_L001_R1_001", "Pos_S41_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% SeqControl),]

    myT <- myT[myT$MouseOrigin == "Charles River", ]

    TR <-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")

    myT <- myT[!(myT$Sample_ID %in% TR),]

    ## Dropping Duplicate sample:
    myT<- myT[!(myT$Sample_ID == "04-04_S63_L001_R1_001"),]

    myT <- myT[order(myT$Sample_ID),]

    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")
    pdf( paste(taxa, "_noControl_noTR_noHarlan_topMDS.pdf",sep=""))

    for (xrun in 1:4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }

plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=16, col=ifelse(myT$Sample_Project=="Lyte_Batch04", "red","blue"))
                    legend("topright",
c("Batch1", "Batch2", "Seq",
		    "Acute", "Chronic", "Control"),
		    pch=c(16, 17, 4, 16, 16),
		    col=c("black", "black", "black",
                    "red", "blue"))
                }
        }

    dev.off()

    write.table(myMDS$CA$u, sep="\t", file=paste("ORIGINAL_pcoa_NoSeqControl_noTR_noHarlan_", taxa, ".txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("ORIGINAL_eigenValues_NoSeqControl_noTR_noHarlan", taxa, ".txt", sep=""), sep="\t")

}

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Resequencing/ForwardReads/")

## Removes Negative and Positive Controls; Removes Harlan Labs Mice
### Here it has 110 rows...
numMetadataCols <- 19
for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT2 <-read.csv(inFileName,header=TRUE,sep="")#, na.strings="BLAH")
    numCols <- ncol(myT2)
    myT2 <- myT2[-which(is.na(myT2$MouseOrigin) == TRUE),]
    myT2 <- myT2[-which(myT2$MouseOrigin == "Harlan Labs"), ]
    myMDS <- capscale(myT2[,2:(ncol(myT2)-numMetadataCols)]~1,distance="bray")
    pdf( paste(taxa, "_noPosNeg_noHarlanMice_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
                for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
                    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""), cex=2.0,
                         pch=ifelse(myT2$Treatment=="Ctrl", 1,
                             ifelse(myT2$Condition == "Acute", 15,
                                    ifelse(myT2$Condition == "14Day", 16, ifelse(myT2$Condition == "Chronic", 17, 0)))), ##col=ifelse(myT2$Treatment=="Ctrl", "black", ifelse(myT2$Condition=="Acute", "blue", ifelse(myT2$Condition=="14Day", "orange","red")))
col=ifelse(myT2$Description == "Fodor-Lyte_Run1_2016-05-25", "red", "blue")
                         )
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun], labs=myT2$TubeNumber, #labs=myT2$Sex,
                   cex=0.7, offset=0)
                    legend("topright", inset=c(-0.25,0),
                           c("Control", "Acute", "14Day", "Chronic",
		    "Batch 1", "Batch 2"),
		    pch=c(1, 15, 16, 17, 16, 16),
		    col=c("black", "black", "black", "black",
                    "red", "blue"))
                }
        }

        dev.off()
        write.table(myMDS$CA$u, sep="\t", file=paste("Reseq_pcoa_noPosNeg_noHarlanMice_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("Reseq_eigenValues_noPosNeg_noHarlanMice_", taxa, ".txt", sep=""), sep="\t")

}

### Then I compared the sorted lists (manually) and arrived that the same samples were present in each. This does mean that the first sampling effort had everything I needed. I may want to consider dropping the low-depth read from the first sampling effort.

for(taxa in taxaLevels )
{
    setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")
    inFileName <- paste("ORIGINAL_pcoa_NoSeqControl_noTR_noHarlan_", taxa, ".txt",sep="")
    myTpcoa <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Resequencing/ForwardReads/")
    inFileName <- paste("Reseq_pcoa_noPosNeg_noHarlanMice_", taxa, ".txt",sep="")
    myT2pcoa <- read.csv(inFileName, header=TRUE, sep="")
    pdf(paste(taxa, "_CompareTopMDS.pdf", sep=""))
    for( iter in 1:10){
        plot(myTpcoa[,iter], myT2pcoa[,iter], main=paste(taxa, "_MDS_", iter, sep=""))
    }
    dev.off()
}

### The comparison of the two datasets for MDSCompare.R
## They have a unique id wrt cage number and organism number
## You are removing the pooled samples from the Harlan Labs caseThere
