rm(list=ls())
library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging/rdpClassifications")
taxaLevels <- c( "phylum", "class", "order", "family", "genus" )

checkSizes <- read.csv("phylum_RawwithMetadata_R1.txt", header=TRUE, sep="", na.strings="BLAH")

## checkSizesR2 <- read.csv("phylumRawwithMetadata_R2_Pooled.txt", header=TRUE, sep="", na.strings="BLAH")

## simpleCounter <- read.csv("LytePooled.simple.counter", header=FALSE, sep="")
## justR1 <- simpleCounter[grep("_R1_", simpleCounter[,1]),]
## justR2 <- simpleCounter[grep("_R2_", simpleCounter[,1]),]
numMetadataCols <- 6

sampleSizes<-rowSums(checkSizes[,2:(ncol(checkSizes)-numMetadataCols)])
## sampleSizesR2<-rowSums(checkSizesR2[,2:(ncol(checkSizesR2)-numMetadataCols)])
newCheck <- cbind(checkSizes, sampleSizes)
## newCheckR2 <- cbind(checkSizesR2, sampleSizesR2)
reNewCheck<-newCheck[,c(1,dim(newCheck)[2])]
## reNewCheckR2<-newCheckR2[,c(1, dim(newCheckR2)[2])]
## testR1 <- merge(reNewCheck, justR1, by.x=1, by.y=1)
## testR2 <- merge(reNewCheckR2, justR2, by.x=1, by.y=1)
# Try for 5000 and 10000 cutoffs

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa, "_LogNormwithMetadata_R1.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
#	myColClasses <- c("character", rep("numeric", numCols-1))

    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")
    pdf( paste(taxa, "topMDS.pdf",sep=""))
    for (xrun in 1:4) {
                for (yrun in 2:4) {
                    if(xrun == yrun){
                        break
                    }
                    plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),
                         ## cex=2.0,
                         cex = log(myT$Age),
                         pch=ifelse(myT$Sample.Type=="LI Lumen", 16, ifelse(myT$Sample.Type =="LI Mucosa", 17, 4)),
                         ## pch= 16 + ifelse(myT$Group == "Old", 0, 1)
                         ## col=ifelse(myT$Treatment=="Ctrl", "yellow", ifelse(myT$Condition=="Acute", "blue", ifelse(myT$Condition=="Chronic", "orange","yellow")))
##                         col = as.numeric(as.factor(myT$Animal.ID))
                         ## col = ifelse(myT$Group == "Old", "BLUE", "RED")
                         col = rainbow(19)[as.numeric(as.factor(myT$Animal.ID))]
                         )
                   textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=as.numeric(as.factor(myT$Pen.Location)),
                  cex=0.7, offset=0)
                     legend("topright",
                            c("LI Lumen", "LI Mucosa", "Feces"),
		   ##  "Acute", "Chronic", "Control"),
                            pch=c(16, 17, 4)#,
                   ##             16, 16, 16),
                   ##         col=c("black", "black", "black",
                            ##             "blue", "orange", "yellow")
                            )
                }
        }

        dev.off()

    write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_", taxa, ".txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("eigenValues_", taxa, ".txt", sep=""), sep="\t")
}
