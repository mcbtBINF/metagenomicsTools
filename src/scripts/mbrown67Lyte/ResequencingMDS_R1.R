rm(list=ls())

library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Resequencing/ForwardReads/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" )

checkSizes <- read.csv("phylumRawwithMetadata_R1_Pooled.txt", header=TRUE, sep="", na.strings="BLAH")

numMetadataCols <- 18

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
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Treatment=="Ctrl", 16, ifelse(myT$Treatment =="Exp", 17, 4)), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sex,
                   cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
# 		    c("Batch01", "Batch02", "Batch03",
c("Ctrl", "Exp", "Seq",
		    "Acute", "Chronic", "Control"),
		    pch=c(16, 17, 4, 16, 16, 16),
		    col=c("black", "black", "black",
		    #"cornflowerblue", "pink"))
                    "blue", "orange", "yellow"))
                                        # Do the key
                }
        }

        dev.off()
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_", taxa, ".txt", sep=""), sep="\t")

}

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
##     SeqControl <- c("DNAC1_S69_L001_R1_001", "DNAC2_S70_L001_R1_001", "DNAN1_S71_L001_R1_001", "DNAN2_S72_L001_R1_001", "Neg_S40_L001_R1_001",
##                    "PCR1Neg_S65_L001_R1_001", "PCR2Pos_S66_L001_R1_001", "Pos_S41_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% SeqControl),]
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")
        pdf( paste(taxa, "_noSeqControl_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Treatment=="Ctrl", 16, ifelse(myT$Treatment =="Exp", 17, 4)), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sex,
                   cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
# 		    c("Batch01", "Batch02", "Batch03",
c("Ctrl", "Exp", "Seq",
		    "Acute", "Chronic", "Control"),
		    pch=c(16, 17, 4, 16, 16, 16),
		    col=c("black", "black", "black",
		    #"cornflowerblue", "pink"))
                    "blue", "orange", "yellow"))
                                        # Do the key
                }
        }

        dev.off()
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}


for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")
    pdf( paste(taxa, "_noSeqControl_checkRun_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
                                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Sample_Project=="Lyte_Batch04", 4, 16), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sex,
                   cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
# 		    c("Batch01", "Batch02", "Batch03",
c("Ctrl", "Exp", "Seq",
		    "Acute", "Chronic", "Control"),
		    pch=c(16, 17, 4, 16, 16, 16),
		    col=c("black", "black", "black",
		    #"cornflowerblue", "pink"))
                    "blue", "orange", "yellow"))
                                        # Do the key
                }
        }

        dev.off()
}

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
##    TechControl<-c("04_125", "04_125_tr", "04_101", "04_101_tr", "04_103", "04_103_tr", "04_74", "04_74_tr", "04_70", "04_70_tr", "04_40", "04_40_tr", "04_41", "04_41_tr", "04_84", "04_84_tr", "C1", "C2", "N1", "N2", "Neg", "Pos", "04_04")
    myT<-myT[myT$Sample_ID %in% TechControl,]
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")
    pdf( paste(taxa, "_noSeqControl_TechControl_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=NA, col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sample_ID,
                   cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
# 		    c("Batch01", "Batch02", "Batch03",
c("Ctrl", "Exp", "Seq",
		    "Acute", "Chronic", "Control"),
		    pch=c(16, 17, 4, 16, 16, 16),
		    col=c("black", "black", "black",
		    #"cornflowerblue", "pink"))
                    "blue", "orange", "yellow"))
                                        # Do the key
                }
        }

        dev.off()
}

