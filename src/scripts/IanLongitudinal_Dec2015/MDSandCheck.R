rm(list=ls())
library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" );

countList <- list()
shannonList <- list()
index<-1
for(taxa in taxaLevels )
{
    rawFileName <- paste0("pivoted_", taxa, "asColumns.txt")
    myTCheck <- read.table(rawFileName, header=TRUE, sep="\t")

                                        # Remove the row labels

                                        # Restricted to the forward reads now.
    myTCheck <- myTCheck[1:147,]
    rownames(myTCheck) <- myTCheck$sample
    myTCheck <- myTCheck[,-1]
    rownames(myTCheck) <- unlist(lapply(strsplit(rownames(myTCheck),split="r1_"),"[[",2))


## Drop bad samples (will manually have to iterate this...
    myTCheck<-myTCheck[!(rownames(myTCheck) %in% c("37", "45", "52", "58", "7", "9"

    ## ,"B20", "C25", "C19", "C12", "B22", "46"
    ##, "B10"
 )),]

    mySums <- rowSums(myTCheck)

    Shannon <- apply(myTCheck[,1:(ncol(myTCheck))], 1, diversity)
    nzCount <- function(x) {
        return(sum(x != 0))
    }
    countList[[index]]<-mySums
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
     myT<-myT[!(rownames(myT) %in% c("37", "45", "52", "58", "7", "9"
    ## ,"B20", "C25", "C19", "C12", "B22", "46",
    ## "B10"
 )),]



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
        pdf( paste(taxa, "CURRDROP_topMDS_R1_test.pdf",sep=""))
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
                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=1.0,
pch=NA,
                                        #pch=ifelse(myMeta$treatment=="Ctrl",16,17),
                         #Write a check for the existence of an A or B.
                         col=colors
                         )
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
# plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=2.0, pch=(14 + strtoi(substr(myMeta$acuteOrChronic, 2,2))), col=ifelse(myMeta$sex=="M","cornflowerblue","pink"))
#                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=rownames(myT),
               text(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labels=paste0(format(Shannon, digits=2),"_",mySums,"_",rownames(myTCheck)), cex=0.7, offset=0, col=colors)
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
        write.table(myMDS$CA$u, sep="\t", file=paste("R1_test_pcoa_CURRDROP_", taxa, ".txt",sep=""))
	write.table(myMDS$CA$eig,file=paste("R1_test_eigenValues_CURRDROP_", taxa, ".txt", sep=""), sep="\t")

}
