rm(list=ls())
library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/FodorProcessedForwardReadsNoCutoff/noLowTaxas")
#setwd("noLowTaxas/")
#setwd("C:\\MarkRatDataAug2015")
#more testing of the GitHub features

taxaLevels <- c( "p", "c", "o", "f", "g", "otu" );

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "TaxaAsColumnsLogNorm.txt.tempC", sep ="")
    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt.tempC", sep ="")
    myT <-read.table(inFileName,header=TRUE,sep="\t")
    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c("character", rep("numeric", numCols-1))
    myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,
                     colClasses=myColClasses)
    metaCols <- ncol(myMeta)
    #myMetaColClasses <- c("character", rep("numeric", metaCols-1))
    myMeta <-read.table(metadataFileName, header=TRUE, sep="\t", row.names=1)#, colClasses = myMetaColClasses)
   # rownames(myMeta)<-rownames(myT)
    excludedRows<-c("70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "50")
                 #   )
                   #,"50")
    excludedMetaRows<-c("s_70", "s_71", "s_72", "s_73", "s_74", "s_75", "s_76", "s_77", "s_78", "s_79", "s_80", "s_81", "s_82", "s_83", "s_84", "s_85", "s_86", "s_87", "s_88", "s_89", "s_90", "s_91", "s_50")

                                        # lines commented out for keeping Batch03
    myT<-myT[-which(rownames(myT) %in% excludedRows),]
    myMeta<-myMeta[-which(rownames(myMeta) %in% excludedMetaRows),]
	myMDS <- capscale(myT~1,distance="bray")


                                        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "topMDS.pdf",sep=""))
# Make a smart iterator here.
    for (xrun in 1:4) {
#            if (xrun < 4) {
                for (yrun in 1:4) {

#                    plot(myMDS$CA$u[,xrun],
#                    myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""),
#                    ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at
#                    level:", taxa,sep=""),cex=2.0,
#                    pch=ifelse(myMeta$batch == "Batch01", 16, ifelse(myMeta$batch=="Batch02",17,18)),
                                        #                    col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
                    colorPaletteFunction <- colorRampPalette(
                        colors = c("yellow", "red", "blue", "green"),
                        space = "Lab"
                        )
                    #here this should be length(unique(myMeta$cage))
                    numColors<-nlevels(myMeta$cage)
                    diamondColorColors<-colorPaletteFunction(numColors)

                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=ifelse(myMeta$acuteOrChronic=="A",16,17), col= ifelse(myMeta$treatment=="Ctrl", "red", "blue"), xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))#ifelse(myMeta$acuteOrChronic=="A","blue","orange"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myMeta$sex,
                   cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
# 		    c("Batch01", "Batch02", "Batch03",
c("Acute", "Chronic",
                                        #"Control", "Experiment"),
"Ctrl", "Exp"
                                        #  "A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"),
	#	    pch=c(16, 17, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16),
pch=c(16,17,16,16),
  col=c("black", "black",
		    #"cornflowerblue", "pink"))
      "red", "blue"
                                        # diamondColorColors[myMeta$cage]
                                        #"blue", "orange"
                          ))
                                        # Do the key
                }
#            }
        }

        dev.off()
                                        #        }
        write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("eigenValues_", taxa, ".txt", sep=""), sep="\t")

}
