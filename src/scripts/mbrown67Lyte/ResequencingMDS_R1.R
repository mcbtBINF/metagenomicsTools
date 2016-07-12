rm(list=ls())

library("vegan")
library("calibrate")

#### PURPOSE
## Check for sizes of reads
## Exclude those falling too low
## Run MDS with pch on condition (Control, 9, 14, 19) [0, 9, 14, 19]
### color on batch
### What about sex?
### What about grouping?
## 1) Full
## 2) Drop Pos and Neg Control
## 3) Drop Harlan Labs Mice
## Questions
### Can you tell if stressed? Can you tell severity of stressed?
### What about batch?
### What about sex?

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Resequencing/ForwardReads/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" )

checkSizes <- read.csv("phylumRawwithMetadata_R1_Pooled.txt", header=TRUE, sep="", na.strings="BLAH")

numMetadataCols <- 17

sampleSizes<-rowSums(checkSizes[,2:(ncol(checkSizes)-numMetadataCols)])

pdf("SampleSizesHistogram.pdf")
hist(sampleSizes)
dev.off()

pdf("SampleSizeslog10Boxplot.pdf")
## Should consider not doubly plotting the outliers for niceness's sake
boxplot(log10(sampleSizes), main="Log10 of Sequence Count", pch=NA)
stripchart(log10(sampleSizes), vertical = TRUE, method = "jitter", add = TRUE, pch=20, col = "BLUE")
dev.off()

## Full visualization; Includes Negative and Positive Controls
for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)

    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")
    pdf( paste(taxa, "topMDS.pdf",sep=""))
    for (xrun in 1:4) {
                for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
                    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""), cex=2.0,
                         pch=ifelse(myT$Treatment=="Ctrl", 1,
                             ifelse(myT$Treatment == "Exp" && myT$Condition == "Acute", 15,
                                    ifelse(myT$Treatment == "Exp" && myT$Condition == "14Day", 16, ifelse(myT$Treatment == "Exp" && myT$Condition == "Chronic", 17, 0)))), ##col=ifelse(myT$Treatment=="Ctrl", "black", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="14Day", "orange","red")))
col=ifelse(myT$Description == "Fodor-Lyte_Run1_2016-05-25", "red", "blue")
                         )
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sex,
                   cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright", inset=c(-0.25,0),
                           c("Control", "Acute", "14Day", "Chronic",
		    "Batch 1", "Batch 2"),
		    pch=c(1, 15, 16, 17, 16, 16),
		    col=c("black", "black", "black", "black",
                    "red", "blue"))
                }
        }

        dev.off()
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_", taxa, ".txt", sep=""), sep="\t")

}

## Removes Negative and Positive Controls

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="")#, na.strings="BLAH")
    numCols <- ncol(myT)
    myT <- myT[-which(is.na(myT$MouseOrigin) == TRUE),]
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")
    pdf( paste(taxa, "_noPosNeg_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
                for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
                    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""), cex=2.0,
                         pch=ifelse(myT$Treatment=="Ctrl", 1,
                             ifelse(myT$Treatment == "Exp" && myT$Condition == "Acute", 15,
                                    ifelse(myT$Treatment == "Exp" && myT$Condition == "14Day", 16, ifelse(myT$Treatment == "Exp" && myT$Condition == "Chronic", 17, 0)))), ##col=ifelse(myT$Treatment=="Ctrl", "black", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="14Day", "orange","red")))
col=ifelse(myT$Description == "Fodor-Lyte_Run1_2016-05-25", "red", "blue")
                         )
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sex,
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
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_noPosNeg_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_noPosNeg_", taxa, ".txt", sep=""), sep="\t")

}

## Removes Negative and Positive Controls; Removes Harlan Labs Mice
for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="")#, na.strings="BLAH")
    numCols <- ncol(myT)
    myT <- myT[-which(is.na(myT$MouseOrigin) == TRUE),]
    myT <- myT[-which(myT$MouseOrigin == "Harlan Labs"), ]
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")
    pdf( paste(taxa, "_noPosNeg_noHarlanMice_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
                for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
                    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""), cex=2.0,
                         pch=ifelse(myT$Treatment=="Ctrl", 1,
                             ifelse(myT$Treatment == "Exp" && myT$Condition == "Acute", 15,
                                    ifelse(myT$Treatment == "Exp" && myT$Condition == "14Day", 16, ifelse(myT$Treatment == "Exp" && myT$Condition == "Chronic", 17, 0)))), ##col=ifelse(myT$Treatment=="Ctrl", "black", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="14Day", "orange","red")))
col=ifelse(myT$Description == "Fodor-Lyte_Run1_2016-05-25", "red", "blue")
                         )
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sex,
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
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_noPosNeg_noHarlanMice_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_noPosNeg_noHarlanMice_", taxa, ".txt", sep=""), sep="\t")

}
