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
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(myT$Source=="Cecal Content", "green", ifelse(myT$Source == "duo", "black", ifelse(myT$Source == "feces", "brown", ifelse(myT$Source == "ileum", "yellow", "red")))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(myT$Sex == "male", "blue", "pink"))
                    par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Male", "Female"),
        pch = 16, cex = 1.1,
                    col=c("blue", "pink"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(myT$Exp.or.Ctrl == "Exp", "red", "black"))
                            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Experiment", "Control"),
        pch = 16, cex = 1.1,
                    col=c("red", "black"))


        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(myT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(myT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(myT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(myT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(myT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(myT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(myT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("FX1", "FC2", "FC1", "FX2", "MX1", "MC2", "MC1", "MX2"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
}

dev.off()

tissuemyT <- rbind(myT[myT$Source == "feces",],myT[myT$Source == "Cecal Content",])
## Display only 2 tissue sample types
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
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Source=="Cecal Content", "green", ifelse(tissuemyT$Source == "duo", "black", ifelse(tissuemyT$Source == "feces", "brown", ifelse(tissuemyT$Source == "ileum", "yellow", "red")))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Sex == "male", "blue", "pink"))
                    par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Male", "Female"),
        pch = 16, cex = 1.1,
                    col=c("blue", "pink"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Exp.or.Ctrl == "Exp", "red", "black"))
                            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Experiment", "Control"),
        pch = 16, cex = 1.1,
                    col=c("red", "black"))


        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(tissuemyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(tissuemyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(tissuemyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(tissuemyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(tissuemyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(tissuemyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(tissuemyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("FX1", "FC2", "FC1", "FX2", "MX1", "MC2", "MC1", "MX2"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
}

dev.off()

## tissuemyT <- myT[myT$Source == "feces",]
behaviormyT <- myT[myT$Source == "feces",]
## Display all 5 tissue sample types
myMDS <- capscale(behaviormyT[,20:58]~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_behavior_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_behavior_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf("behavior_topMDS.pdf")
for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
                    par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Male", "Female"),
        pch = 16, cex = 1.1,
                    col=c("blue", "pink"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
                            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Experiment", "Control"),
        pch = 16, cex = 1.1,
                    col=c("red", "black"))


        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("FX1", "FC2", "FC1", "FX2", "MX1", "MC2", "MC1", "MX2"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
}

dev.off()

## elevated.plus.maze
## Display all 5 tissue sample types
myMDS <- capscale(behaviormyT[,20:30]~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_epm_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_epm_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf("elevatedplusmaze_topMDS.pdf")
for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
                    par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Male", "Female"),
        pch = 16, cex = 1.1,
                    col=c("blue", "pink"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
                            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Experiment", "Control"),
        pch = 16, cex = 1.1,
                    col=c("red", "black"))


        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("FX1", "FC2", "FC1", "FX2", "MX1", "MC2", "MC1", "MX2"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
}

dev.off()

## ldmyT <- behaviormyT[,12:21]
## Display all 5 tissue sample types
myMDS <- capscale(behaviormyT[,32:40]~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_lightdark_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_lightdark_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf("lightdark_topMDS.pdf")
for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
                    par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Male", "Female"),
        pch = 16, cex = 1.1,
                    col=c("blue", "pink"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
                            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Experiment", "Control"),
        pch = 16, cex = 1.1,
                    col=c("red", "black"))


        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("FX1", "FC2", "FC1", "FX2", "MX1", "MC2", "MC1", "MX2"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
}

dev.off()

## ofmyT <- behaviormyT[,22:39]
## Display all 5 tissue sample types
myMDS <- capscale(behaviormyT[,41:58]~1,distance="bray")

write.table(myMDS$CA$u, sep="\t", file=paste("LyteSharon_r01_openfield_pcoa.txt",sep=""))
write.table(myMDS$CA$eig,file=paste("LyteSharon_r01_openfield_eigenValues.txt", sep=""), sep="\t")

percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

pdf("openfield_topMDS.pdf")
for (xrun in 1:4) {
    for (yrun in 2:4) {
        if(xrun == yrun){
            break
        }
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(2,2))
        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Sample Tissue)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Source=="Cecal Content", "green", ifelse(behaviormyT$Source == "duo", "black", ifelse(behaviormyT$Source == "feces", "brown", ifelse(behaviormyT$Source == "ileum", "yellow", "red")))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Cecal Content", "Duodenum", "Feces", "Ileum", "Jejunum"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Mouse Sex)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Sex == "male", "blue", "pink"))
                    par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Male", "Female"),
        pch = 16, cex = 1.1,
                    col=c("blue", "pink"))

        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Stress/Unstressed)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Exp.or.Ctrl == "Exp", "red", "black"))
                            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("Experiment", "Control"),
        pch = 16, cex = 1.1,
                    col=c("red", "black"))


        plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
             xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
             ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
             main="PCoA at otu level \n(Cage)", cex=2.0, ##bty="L",
             pch = 16,
             col=ifelse(behaviormyT$Cage=="female #1,2,3,4 in same cage", "green", ifelse(behaviormyT$Cage =="female #13,14,15,16 in same cage", "blue", ifelse(behaviormyT$Cage == "female #5,6,7,8 in same cage", "brown", ifelse(behaviormyT$Cage == "female #9,10,11,12 in same cage", "yellow", ifelse(behaviormyT$Cage == "male #1,2,3,4 in same cage", "red", ifelse(behaviormyT$Cage == "male #13,14,15,16 in same cage", "black", ifelse(behaviormyT$Cage == "male #5,6,7,8 in same cage", "pink", "orange"))))))))
            par(xpd=TRUE)

            legend("topright", inset=c(-0.9,0),
                   c("FX1", "FC2", "FC1", "FX2", "MX1", "MC2", "MC1", "MX2"),
        pch = 16, cex = 1.1,
                    col=c("green", "blue", "brown", "yellow", "red", "black", "pink", "orange"))
    }
}

dev.off()

## Compute MDS on behavioral data and output
