rm(list=ls())
library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrolData/Carroll_Longitudinal")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" );

for(taxa in taxaLevels )
{
    inFileName <- paste( "pivoted_", taxa, "asColumnsLogNormal_R2.txt", sep ="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt", sep ="")
    myT <-read.table(inFileName,header=TRUE,sep="\t")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c("character", rep("numeric", numCols-1))
	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)

    patientB <- grep("B", rownames(myT))
    patientC <- grep("C", rownames(myT))
#    patientA <-
# lines commented out for keeping Batch03
#myT<-myT[-which(rownames(myT) %in% excludedRows),]
	myMDS <- capscale(myT~1,distance="bray")

        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "topMDS_R2.pdf",sep=""))
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
                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0,
pch=16,
                                        #pch=ifelse(myMeta$treatment=="Ctrl",16,17),
                         #Write a check for the existence of an A or B.
                         col=ifelse(rownames(myT) %in% rownames(myT)[patientB],"red",ifelse(rownames(myT) %in% rownames(myT)[patientC], "blue", "green"))
                         )
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
               textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
                   cex=0.7, offset=0)
		    # par(xpd=TRUE)
                    legend("topright",
 		    c("Patient A", "Patient B", "Patient C"),
#c("Acute", "Chronic",
#		    "Control", "Experiment"),
		    pch=c(16, 16, 16),
		    col=c("green", "red", "blue"))
                }
#            }
        }

        dev.off()
                                        #        }
        write.table(myMDS$CA$u, sep="\t", file=paste("R2_pcoa_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("R2_eigenValues_", taxa, ".txt", sep=""), sep="\t")

}
