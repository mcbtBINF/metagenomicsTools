rm(list=ls())
library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" );

countList <- list()
shannonList <- list()
index<-1

for(dropCounter in 1:25){
for(taxa in taxaLevels )
{
    rawFileName <- paste0("pivoted_", taxa, "asColumns.txt")
    myTCheck <- read.table(rawFileName, header=TRUE, sep="\t")
    myTCheck <- myTCheck[1:147,]
    rownames(myTCheck) <- myTCheck$sample
    rownames(myTCheck) <- unlist(lapply(strsplit(rownames(myTCheck),split="r1_"),"[[",2))
                                        # Remove the row labels
    myTCheck <- myTCheck[,-1]
                                        # Restricted to the forward reads now.
    myTCheck <- myTCheck[1:147,]

## Drop bad samples (will manually have to iterate this...

    if(index == 1){
        storeList<-rowSums(myTCheck)
    }
    myTCheck <- myTCheck[!(rownames(myTCheck) %in% names(sort(storeList))[1:dropCounter]),]

    mySums <- rowSums(myTCheck)
    countList[[index]]<-mySums

    Shannon <- apply(myTCheck[,1:(ncol(myTCheck))], 1, diversity)
    nzCount <- function(x) {
        return(sum(x != 0))
    }

    shannonList[[index]]<-Shannon
    index<-index+1
    nonZero <- apply(myTCheck, 1, nzCount)


    ## inFileName <- paste( "pivoted_", taxa, "asColumnsLogNormal_R2.txt", sep ="")
    inFileName <- paste(taxa, "LogNormalwithMetadataDailyR2_Edit.txt", sep="")
#    metadataFileName <- paste( taxa, "TaxaAsColumnsLogNormWithMetaData.txt", sep ="")
    myT <-read.table(inFileName,header=TRUE,sep="\t")
#    	myMeta <-read.table(metadataFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c("character", rep("numeric", numCols-1))
    myT <-read.table(inFileName,header=TRUE,sep="\t", colClasses=myColClasses)
    myT <- myT[ !is.na(myT[2]), ]
    rownames(myT)<-myT[,1]
    myT<-myT[,-1]

    myT <- myT[!(rownames(myT) %in% names(sort(storeList))[1:dropCounter]),]
    print(names(sort(storeList))[1:dropCounter])


    colors <- rep("Black", dim(myT)[1])
    colors[grep("B", rownames(myT))] <- c("Blue")
    colors[grep("C", rownames(myT))] <- c("Red")

    ## patientB <- grep("B", rownames(myT))
    ## patientC <- grep("C", rownames(myT))
#    patientA <-
# lines commented out for keeping Batch03
#myT<-myT[-which(rownames(myT) %in% excludedRows),]
	myMDS <- capscale(myT[,1:(ncol(myT)-14)]~1,distance="bray")

        # In the future add a scree plot component here.
        # Also include color options for various metadata.
        # pdf( paste(taxa, "1v2.pdf",sep=""))
                                        #        if(taxa == "p") {
        pdf( paste(taxa, "_", dropCounter, "_topMDS_R1_test.pdf",sep=""))
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
                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,"_", names(sort(storeList))[dropCounter], sep=""),cex=1.0,
pch=NA,
                                        #pch=ifelse(myMeta$treatment=="Ctrl",16,17),
                         #Write a check for the existence of an A or B.
                         col=colors
                         )
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
               text(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labels=paste0(format(Shannon, digits=2),"_",mySums), cex=0.7, offset=0, col=colors)
		    # par(xpd=TRUE)
                    legend("topright",
 		    c("Patient A", "Patient B", "Patient C"),
#c("Acute", "Chronic",
#		    "Control", "Experiment"),
		    pch=c(16, 16, 16),
		    col=c("black", "blue", "red"))
                }
#            }
        }

        dev.off()
                                        #        }
        write.table(myMDS$CA$u, sep="\t", file=paste("R1_test_pcoa_", taxa, "_", dropCounter, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("R1_test_eigenValues_", taxa, "_", dropCounter, ".txt", sep=""), sep="\t")

}
}
