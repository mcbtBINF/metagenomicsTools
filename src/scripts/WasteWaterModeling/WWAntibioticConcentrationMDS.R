rm(list=ls())

library("vegan")
library("calibrate")

## Plan
## Plot for each time point
## Plot over all time points
## Get the % Variance explained as well
## Standardize MDS plots across

setwd("/Users/mbrown67/Documents/Fodor/Datasets/WasteWater/")

tp <- c(1, 2, 3)
bigT <- list()
for(tpIter in tp)
{
    inFileName <- paste( "MetaMaster_tp", tpIter, ".csv", sep ="")
    myT <-read.csv(inFileName,header=TRUE, na.strings="BLAH", stringsAsFactors=FALSE)

    ## Replace "ND" cells

    ## Replace "<LOQ" cells

    antibioticConcentration <- myT[,1:10]
    myTNoSeq <- myT[,1:31]
    myTNoRepSeqNoSeq <- myT[1:66, 1:31]
    myTAntibioticOnly <- subset(myTNoRepSeqNoSeq, !duplicated(SampleNumber))

    myT <- myTAntibioticOnly[,1:10]

    ## This replacement will have to be done more carefully later
    myT[myT=="ND"|myT=="<LOQ"] <- 0

    myT<-data.matrix(myT)

    myTAntibioticOnly[1:22, 1:10] <- myT

    myT <- myTAntibioticOnly

    numMetadataCols <- 21

    myMDS <- capscale(data.matrix(myT[,1:(ncol(myT)-numMetadataCols)])~1,distance="bray")
    pdf( paste(tpIter, "_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
        for (yrun in 2:4) {
            if(xrun == yrun){
                break
            }
            par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
            plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at tp:", tpIter, sep=""), cex=2.0,
                 pch = ifelse(grepl("Mallard", myT$Sample.Location),16,17),
                 col = ifelse(grepl("Hospital", myT$Sample.Location), "red", ifelse(grepl("Plant", myT$Sample.Location), "blue", "black"))
                 )

            textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sample.name,
                   cex=0.7, offset=0)
            ##                                 # par(xpd=TRUE)
            ## legend("topright", inset=c(-0.25,0),
            ##        c("Control", "Acute", "14Day", "Chronic",
            ##          "Batch 1", "Batch 2"),
            ##        pch=c(1, 15, 16, 17, 16, 16),
            ##        col=c("black", "black", "black", "black",
            ##            "red", "blue"))
        }
    }

    dev.off()
    write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_", tpIter, ".txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("eigenValues_", tpIter, ".txt", sep=""), sep="\t")
    bigT[[tpIter]] <- myT
}

myT <- rbind(bigT[[1]], bigT[[2]], bigT[[3]])
myMDS <- capscale(data.matrix(myT[,1:(ncol(myT)-numMetadataCols)])~1,distance="bray")
pdf( paste("ALLtimes", "_topMDS.pdf",sep=""))
for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main="PCoA across all time points", cex=2.0,
             pch = ifelse(grepl("Mallard", myT$Sample.Location),16,17),
             col = ifelse(myT$Timepoint == 1, "black", ifelse(myT$Timepoint == 2, "yellow", "green"))
             )

        textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sample.name,
               cex=0.7, offset=0)
        ##                                 # par(xpd=TRUE)
        ## legend("topright", inset=c(-0.25,0),
        ##        c("Control", "Acute", "14Day", "Chronic",
        ##          "Batch 1", "Batch 2"),
        ##        pch=c(1, 15, 16, 17, 16, 16),
        ##        col=c("black", "black", "black", "black",
        ##            "red", "blue"))
    }
}

dev.off()
write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_ALLtime", ".txt",sep=""))
write.table(myMDS$CA$eig,file=paste("eigenValues_ALLtime", ".txt", sep=""), sep="\t")

### Then I'll need to join the timepoints
