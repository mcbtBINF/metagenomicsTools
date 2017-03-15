rm(list=ls())

library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/")

inFileName <- "qiimeOTULogNormwithMetadata.txt"
myT <- read.table(inFileName,header=TRUE,sep="\t")

endMetadataIndex <- which(colnames(myT) == "counts")

## Display all 5 tissue sample types
myMDS <- capscale(myT[,(endMetadataIndex +1):ncol(myT)]~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf( paste("qiime_", "R1_topMDS.pdf",sep=""))

for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(myT$Source=="Cecal Content", "green", ifelse(myT$Source == "duo", "black", ifelse(myT$Source == "feces", "brown", ifelse(myT$Source == "ileum", "yellow", "red")))))
             ##pch=ifelse(myT[myT$Read == "R1",]$Source=="Cecal Content", "C", ifelse(myT$Source == "duo", "D", ifelse(myT[myT$Read == "R1",]$Source == "feces", "F", ifelse(myT[myT$Read == "R1",]$Source == "ileum", "I", "J")))),
             ## pch=ifelse(myT$Exp.or.Ctrl == "Exp", 17, 16)

             ##col=ifelse(myT[myT$Read == "R1",]$Exp.or.Ctrl=="Exp", "red", "blue"))
             ##col=ifelse(myT$Sex=="male", "blue", "pink"))
        ## textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],
               ##labs=ifelse(myT[myT$Read == "R1",]$Sex == "Male", "M", "F"),
               ## col=ifelse(myT$Source=="Cecal Content", "C", ifelse(myT$Source == "duo", "D", ifelse(myT$Source == "feces", "F", ifelse(myT$Source == "ileum", "I", "J")))),
               ## cex=0.7, offset=0)
        par(xpd=TRUE)

        legend("topright", inset=c(-0.35,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
    pch = 16, cex = 1.1,
                col=c("green", "blue", "brown", "yellow", "red"))
    }
}

dev.off()

tissuemyT <- rbind(myT[myT$Source == "feces",],myT[myT$Source == "Cecal Content",])
## Display all 5 tissue sample types
myMDS <- capscale(tissuemyT[,(endMetadataIndex +1):ncol(tissuemyT)]~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_cecalfecal_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_cecalfecal_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf( paste("qiime_", "cecalfecal_R1_topMDS.pdf",sep=""))

for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Source=="Cecal Content", "green", ifelse(tissuemyT$Source == "duo", "black", ifelse(tissuemyT$Source == "feces", "brown", ifelse(tissuemyT$Source == "ileum", "yellow", "red")))))
             ##pch=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "C", ifelse(tissuemyT$Source == "duo", "D", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "F", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "I", "J")))),
             ## pch=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", 17, 16),
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Exp.or.Ctrl=="Exp", "red", "blue"))
             ## col=ifelse(tissuemyT$Sex=="male", "blue", "pink"))
        ## textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],
               ##labs=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Sex == "Male", "M", "F"),
        ##       labs=ifelse(tissuemyT$Source=="Cecal Content", "C", "F"),
        ##                   cex=0.7, offset=0)
        ## par(xpd=TRUE)
                legend("topright", inset=c(-0.35,0),
               c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
    pch = 16, cex = 1.1,
                col=c("green", "black", "brown", "yellow", "red"))


    }
}

dev.off()

tissuemyT <- myT[myT$Source == "feces",]
behaviormyT <- myT[myT$Source == "feces", 20:58]
## Display all 5 tissue sample types
myMDS <- capscale(behaviormyT~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_behavior_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_behavior_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf("behavior_topMDS.pdf")

for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level", cex=2.0, ##bty="L",
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "green", ifelse(tissuemyT$Source == "duo", "blue", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "brown", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "yellow", "red")))),
             ##pch=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "C", ifelse(tissuemyT$Source == "duo", "D", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "F", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "I", "J")))),
             pch=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", 17, 16),
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Exp.or.Ctrl=="Exp", "red", "blue"))
             col=ifelse(tissuemyT$Sex=="male", "blue", "pink"))
        textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],
               ##labs=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Sex == "Male", "M", "F"),
               labs=ifelse(tissuemyT$Sex=="Cecal Content", "C", "F"),
                           cex=0.7, offset=0)
        ## par(xpd=TRUE)

        legend("topright", inset=c(-0.3,0),
               c("Experiment", "Control", "Male", "Female"),
               ## "Cecal Content (C)", "Duodenum (D)", "Feces (F)", "Ileum (I)", "Jejunum (J)"),
               pch=c(17, 16, 16, 16),
               col=c("black", "black", "blue", "pink"))
    }
}

dev.off()

epmmyT <- behaviormyT[,1:11]
## Display all 5 tissue sample types
myMDS <- capscale(epmmyT~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_epm_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_epm_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf("elevatedplusmaze_topMDS.pdf")

for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Elevated Plus Maze Test", cex=2.0, ##bty="L",
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "green", ifelse(tissuemyT$Source == "duo", "blue", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "brown", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "yellow", "red")))),
             ##pch=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "C", ifelse(tissuemyT$Source == "duo", "D", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "F", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "I", "J")))),
             pch=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", 17, 16),
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Exp.or.Ctrl=="Exp", "red", "blue"))
             col=ifelse(tissuemyT$Sex=="male", "blue", "pink"))
        textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],
               ##labs=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Sex == "Male", "M", "F"),
               labs=ifelse(tissuemyT$Source=="Cecal Content", "C", "F"),
                           cex=0.7, offset=0)
        ## par(xpd=TRUE)

        legend("topright", inset=c(-0.3,0),
               c("Experiment", "Control", "Male", "Female"),
               ## "Cecal Content (C)", "Duodenum (D)", "Feces (F)", "Ileum (I)", "Jejunum (J)"),
               pch=c(17, 16, 16, 16),
               col=c("black", "black", "blue", "pink"))
    }
}

dev.off()

ldmyT <- behaviormyT[,12:21]
## Display all 5 tissue sample types
myMDS <- capscale(ldmyT~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_lightdark_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_lightdark_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf("lightdark_topMDS.pdf")

for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Light Dark Test", cex=2.0, ##bty="L",
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "green", ifelse(tissuemyT$Source == "duo", "blue", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "brown", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "yellow", "red")))),
             ##pch=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "C", ifelse(tissuemyT$Source == "duo", "D", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "F", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "I", "J")))),
             pch=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", 17, 16),
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Exp.or.Ctrl=="Exp", "red", "blue"))
             col=ifelse(tissuemyT$Sex=="male", "blue", "pink"))
        textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],
               ##labs=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Sex == "Male", "M", "F"),
               labs=ifelse(tissuemyT$Source=="Cecal Content", "C", "F"),
                           cex=0.7, offset=0)
        ## par(xpd=TRUE)

        legend("topright", inset=c(-0.3,0),
               c("Experiment", "Control", "Male", "Female"),
               ## "Cecal Content (C)", "Duodenum (D)", "Feces (F)", "Ileum (I)", "Jejunum (J)"),
               pch=c(17, 16, 16, 16),
               col=c("black", "black", "blue", "pink"))
    }
}

dev.off()

ofmyT <- behaviormyT[,22:39]
## Display all 5 tissue sample types
myMDS <- capscale(ofmyT~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_openfield_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_openfield_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf("openfield_topMDS.pdf")

for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="Open Field Test", cex=2.0, ##bty="L",
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "green", ifelse(tissuemyT$Source == "duo", "blue", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "brown", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "yellow", "red")))),
             ##pch=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source=="Cecal Content", "C", ifelse(tissuemyT$Source == "duo", "D", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "feces", "F", ifelse(tissuemyT[tissuemyT$Read == "R1",]$Source == "ileum", "I", "J")))),
             pch=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", 17, 16),
             ##col=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Exp.or.Ctrl=="Exp", "red", "blue"))
             col=ifelse(tissuemyT$Sex=="male", "blue", "pink"))
        textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],
               ##labs=ifelse(tissuemyT[tissuemyT$Read == "R1",]$Sex == "Male", "M", "F"),
               labs=ifelse(tissuemyT$Source=="Cecal Content", "C", "F"),
                           cex=0.7, offset=0)
        ## par(xpd=TRUE)

        legend("topright", inset=c(-0.3,0),
               c("Experiment", "Control", "Male", "Female"),
               ## "Cecal Content (C)", "Duodenum (D)", "Feces (F)", "Ileum (I)", "Jejunum (J)"),
               pch=c(17, 16, 16, 16),
               col=c("black", "black", "blue", "pink"))
    }
}

dev.off()

## Compute MDS on behavioral data and output
