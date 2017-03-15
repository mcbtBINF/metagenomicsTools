# All
# No Control
rm(list=ls())
library("vegan")
library("calibrate")
                                        #setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Lyte_Batch04-27270381")
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

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)


#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "topMDS.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Treatment=="Ctrl", 16, ifelse(myT$Treatment =="Exp", 17, 4)), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
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
#            }
        }

        dev.off()
                                        #        }
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_", taxa, ".txt", sep=""), sep="\t")

}

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    SeqControl <- c("DNAC1_S69_L001_R1_001", "DNAC2_S70_L001_R1_001", "DNAN1_S71_L001_R1_001", "DNAN2_S72_L001_R1_001", "Neg_S40_L001_R1_001",
                    "PCR1Neg_S65_L001_R1_001", "PCR2Pos_S66_L001_R1_001", "Pos_S41_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% SeqControl),]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_noSeqControl_topMDS.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }


#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Treatment=="Ctrl", 16, ifelse(myT$Treatment =="Exp", 17, 4)), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
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
#            }
        }

        dev.off()
                                        #        }
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}


for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
#    SeqControl <- c("DNAC1_S69_L001_R1_001", "DNAC2_S70_L001_R1_001", "DNAN1_S71_L001_R1_001", "DNAN2_S72_L001_R1_001", "Neg_S40_L001_R1_001", "PCR1Neg_S65_L001_R1_001", "PCR2Pos_S66_L001_R1_001", "Pos_S41_L001_R1_001")
#    myT<-myT[!(myT$MatchFile %in% SeqControl),]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_noSeqControl_checkRun_topMDS.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }

#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Sample_Project=="Lyte_Batch04", 4, 16), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
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
#            }
        }

        dev.off()
                                        #        }
#        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
#	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    TechControl<-c("04_125", "04_125_tr", "04_101", "04_101_tr", "04_103", "04_103_tr", "04_74", "04_74_tr", "04_70", "04_70_tr", "04_40", "04_40_tr", "04_41", "04_41_tr", "04_84", "04_84_tr", "C1", "C2", "N1", "N2", "Neg", "Pos", "04_04")
    myT<-myT[myT$Sample_ID %in% TechControl,]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_noSeqControl_TechControl_topMDS.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }

#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=NA, col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
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
#            }
        }

        dev.off()
                                        #        }
#        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
#	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}
### THIS! ###
for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    SeqControl <- c("DNAC1_S69_L001_R1_001", "DNAC2_S70_L001_R1_001", "DNAN1_S71_L001_R1_001", "DNAN2_S72_L001_R1_001", "Neg_S40_L001_R1_001",
                    "PCR1Neg_S65_L001_R1_001", "PCR2Pos_S66_L001_R1_001", "Pos_S41_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% SeqControl),]

    myT <- myT[myT$MouseOrigin == "Charles River", ]

    TR <-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")

    ## Dropping Duplicate sample:
    myT<- myT[!(myT$Sample_ID == "04-04_S63_L001_R1_001"),]


##    TechControl<-c("04_125", "04_125_tr", "04_101", "04_101_tr", "04_103", "04_103_tr", "04_74", "04_74_tr", "04_70", "04_70_tr", "04_40", "04_40_tr", "04_41", "04_41_tr", "04_84", "04_84_tr", "C1", "C2", "N1", "N2", "Neg", "Pos", "04_04")
##    myT<-myT[myT$Sample_ID %in% TechControl,]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_noControl_noTR_noHarlan_topMDS.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }

#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=16, col=ifelse(myT$Sample_Project=="Lyte_Batch04", "red","blue"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
               #textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sample_ID,
               #    cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
# 		    c("Batch01", "Batch02", "Batch03",
c("Batch1", "Batch2", "Seq",
		    "Acute", "Chronic", "Control"),
		    pch=c(16, 17, 4, 16, 16),
		    col=c("black", "black", "black",
		    #"cornflowerblue", "pink"))
                    "red", "blue"))
                                        # Do the key
                }
#            }
        }

        dev.off()
                                        #        }
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl_noTR_noHarlan_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl_noTR_noHarlan", taxa, ".txt", sep=""), sep="\t")

}


#### Now for the reverse reads.

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R2_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)


#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "topMDS_R2.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Treatment=="Ctrl", 16, ifelse(myT$Treatment =="Exp", 17, 4)), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
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
#            }
        }

        dev.off()
                                        #        }
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_", taxa, ".txt", sep=""), sep="\t")

}

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R2_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    SeqControl <- c("DNAC1_S69_L001_R2_001", "DNAC2_S70_L001_R2_001", "DNAN1_S71_L001_R2_001", "DNAN2_S72_L001_R2_001", "Neg_S40_L001_R2_001",
                    "PCR1Neg_S65_L001_R2_001", "PCR2Pos_S66_L001_R2_001", "Pos_S41_L001_R2_001")
    myT<-myT[!(myT$MatchFile %in% SeqControl),]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_noSeqControl_topMDS_R2.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }


#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Treatment=="Ctrl", 16, ifelse(myT$Treatment =="Exp", 17, 4)), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
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
#            }
        }

        dev.off()
                                        #        }
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}


for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R2_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
#    SeqControl <- c("DNAC1_S69_L001_R2_001", "DNAC2_S70_L001_R2_001", "DNAN1_S71_L001_R2_001", "DNAN2_S72_L001_R2_001", "Neg_S40_L001_R2_001", "PCR1Neg_S65_L001_R2_001", "PCR2Pos_S66_L001_R2_001", "Pos_S41_L001_R2_001")
#    myT<-myT[!(myT$MatchFile %in% SeqControl),]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_noSeqControl_checkRun_topMDS_R2.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }

#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myT$Sample_Project=="Lyte_Batch04", 4, 16), col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
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
#            }
        }

        dev.off()
                                        #        }
#        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
#	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R2_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    TechControl<-c("04_125", "04_125_tr", "04_101", "04_101_tr", "04_103", "04_103_tr", "04_74", "04_74_tr", "04_70", "04_70_tr", "04_40", "04_40_tr", "04_41", "04_41_tr", "04_84", "04_84_tr", "C1", "C2", "N1", "N2", "Neg", "Pos", "04_04")
    myT<-myT[myT$Sample_ID %in% TechControl,]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_noSeqControl_TechControl_topMDS_R2.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }

#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=NA, col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow"))))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
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
#            }
        }

        dev.off()
                                        #        }
#        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
#	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R2_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    SeqControl <- c("DNAC1_S69_L001_R2_001", "DNAC2_S70_L001_R2_001", "DNAN1_S71_L001_R2_001", "DNAN2_S72_L001_R2_001", "Neg_S40_L001_R2_001",
                    "PCR1Neg_S65_L001_R2_001", "PCR2Pos_S66_L001_R2_001", "Pos_S41_L001_R2_001")
    myT<-myT[!(myT$MatchFile %in% SeqControl),]

##    TechControl<-c("04_125", "04_125_tr", "04_101", "04_101_tr", "04_103", "04_103_tr", "04_74", "04_74_tr", "04_70", "04_70_tr", "04_40", "04_40_tr", "04_41", "04_41_tr", "04_84", "04_84_tr", "C1", "C2", "N1", "N2", "Neg", "Pos", "04_04")
##    myT<-myT[myT$Sample_ID %in% TechControl,]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_checkbatch_topMDS_R2.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }

#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=16, col=ifelse(myT$Sample_Project=="Lyte_Batch04", "red","blue"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
               #textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sample_ID,
               #    cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
# 		    c("Batch01", "Batch02", "Batch03",
c("Batch1", "Batch2", "Seq",
		    "Acute", "Chronic", "Control"),
		    pch=c(16, 17, 4, 16, 16),
		    col=c("black", "black", "black",
		    #"cornflowerblue", "pink"))
                    "red", "blue"))
                                        # Do the key
                }
#            }
        }

        dev.off()
                                        #        }
#        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
#	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}

#### both check
for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    inFileNameR2 <- paste( taxa, "LogNormwithMetadata_R2_Pooled.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    myTR2 <- read.csv(inFileNameR2, header=TRUE,sep="", na.strings="BLAH")
    myT<- rbind(myT, myTR2)
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
        SeqControl <- c("DNAC1_S69_L001_R1_001", "DNAC2_S70_L001_R1_001", "DNAN1_S71_L001_R1_001", "DNAN2_S72_L001_R1_001", "Neg_S40_L001_R1_001",
                    "PCR1Neg_S65_L001_R1_001", "PCR2Pos_S66_L001_R1_001", "Pos_S41_L001_R1_001")

    SeqControlR2 <- c("DNAC1_S69_L001_R2_001", "DNAC2_S70_L001_R2_001", "DNAN1_S71_L001_R2_001", "DNAN2_S72_L001_R2_001", "Neg_S40_L001_R2_001",
                    "PCR1Neg_S65_L001_R2_001", "PCR2Pos_S66_L001_R2_001", "Pos_S41_L001_R2_001")
    myT<-myT[!(myT$MatchFile %in% SeqControl),]
    myT<-myT[!(myT$MatchFile %in% SeqControlR2),]

##    TechControl<-c("04_125", "04_125_tr", "04_101", "04_101_tr", "04_103", "04_103_tr", "04_74", "04_74_tr", "04_70", "04_70_tr", "04_40", "04_40_tr", "04_41", "04_41_tr", "04_84", "04_84_tr", "C1", "C2", "N1", "N2", "Neg", "Pos", "04_04")
##    myT<-myT[myT$Sample_ID %in% TechControl,]
                                        #    myT <- myT
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)
 #   excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91"
  #      )
# lines commented out for keeping Batch03
                                        #myT<-myT[-which(rownames(myT) %in% excludedRows),]
#    myT[is.na(myT)]<-as.character("Non")
#     print("error")
    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_checkbatch_topMDS_both.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                           for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }

#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
#                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=16, col=ifelse(myT$Sample_Project=="Lyte_Batch04", "red","blue"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
               #textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Sample_ID,
               #    cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
# 		    c("Batch01", "Batch02", "Batch03",
c("Batch1", "Batch2", "Seq",
		    "Acute", "Chronic", "Control"),
		    pch=c(16, 17, 4, 16, 16),
		    col=c("black", "black", "black",
		    #"cornflowerblue", "pink"))
                    "red", "blue"))
                                        # Do the key
                }
#            }
        }

        dev.off()
                                        #        }
#        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_NoSeqControl", taxa, ".txt",sep=""))
#	write.table(myMDS$CA$eig,file=paste("eigenValues_NoSeqControl", taxa, ".txt", sep=""), sep="\t")

}

