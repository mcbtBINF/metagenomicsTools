rm(list=ls())
library("RColorBrewer")
library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging/")
KEGGLevelsLevels <- c(1, 2, 3) ## c( "phylum", "class", "order", "family", "genus" )

## checkSizes <- read.csv("phylum_RawwithMetadata_R1.txt", header=TRUE, sep="", na.strings="BLAH")

## checkSizesR2 <- read.csv("phylumRawwithMetadata_R2_Pooled.txt", header=TRUE, sep="", na.strings="BLAH")

## simpleCounter <- read.csv("LytePooled.simple.counter", header=FALSE, sep="")
## justR1 <- simpleCounter[grep("_R1_", simpleCounter[,1]),]
## justR2 <- simpleCounter[grep("_R2_", simpleCounter[,1]),]
numMetadataCols <- 6

## sampleSizes<-rowSums(checkSizes[,2:(ncol(checkSizes)-numMetadataCols)])
## sampleSizesR2<-rowSums(checkSizesR2[,2:(ncol(checkSizesR2)-numMetadataCols)])
## newCheck <- cbind(checkSizes, sampleSizes)
## newCheckR2 <- cbind(checkSizesR2, sampleSizesR2)
## reNewCheck<-newCheck[,c(1,dim(newCheck)[2])]
## reNewCheckR2<-newCheckR2[,c(1, dim(newCheckR2)[2])]
## testR1 <- merge(reNewCheck, justR1, by.x=1, by.y=1)
## testR2 <- merge(reNewCheckR2, justR2, by.x=1, by.y=1)
# Try for 5000 and 10000 cutoffs

tissueKept <- c("LI Lumen", "LI Mucosa", "Feces")
## tissueKept <- "LI Lumen"
## tissueKept <- "LI Mucosa"
## tissueKept <- "Feces"


for(KEGGLevels in KEGGLevelsLevels )
{
    inFileName <- paste("PICRUSt_", KEGGLevels, "_LogNormwithMetadata.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
    ## Restrict the tissue focus
    myT<-myT[myT$Sample.Type %in% tissueKept,]
                          colVector <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffff00", "#63dcfe", "#abd9e9", "#74add1", "#4575b4",  "#313695", "#0000ff")[as.numeric(as.factor(myT$Pen.Location))]


    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")

     pdf( paste(KEGGLevels, "_PICRUSt_", paste0(tissueKept,collapse="_"), "_topMDS.pdf",sep=""))
     for (xrun in 1:4) {
                 for (yrun in 2:4) {
                     if(xrun == yrun){
                         break
                     }
                     plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", KEGGLevels,sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                           cex=2.0,
                          ## cex = log(myT$Age),

                          pch= ifelse(myT$Group=="Young", 0, 15) + ifelse(myT$Sample.Type=="LI Lumen", 0, ifelse(myT$Sample.Type =="LI Mucosa", 1, 2)),

                          col= colVector##[as.numeric(as.factor(myT$Pen.Location))]
                          ##col = brewer.pal(12, "Set3")[as.numeric(as.factor(myT$Pen.Location))]
                          )
                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Animal.ID,
                   cex=0.7, offset=0)
                      legend("topright", c("(Young) LI Lumen", "(Young) LI Mucosa", "(Young) Feces", "(Old) LI Lumen", "(Old) LI Mucosa", "(Old) Feces"),
## 		   ##  "Acute", "Chronic", "Control"),
                             pch=c(0, 1, 2, 15, 16, 17)#,
##                    ##             16, 16, 16),
##                    ##         col=c("black", "black", "black",
##                             ##             "blue", "orange", "yellow")
                             )
                 }
         }

         dev.off()

     write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_PICRUSt_", KEGGLevels, paste0(tissueKept,collapse="_"), ".txt", sep=""))
     write.table(myMDS$CA$eig,file=paste("eigenValues_PICRUSt_", KEGGLevels, paste0(tissueKept,collapse="_"), ".txt", sep=""), sep="\t")
}


tissueKept <- "Feces"


for(KEGGLevels in KEGGLevelsLevels )
{
    inFileName <- paste("PICRUSt_", KEGGLevels, "_LogNormwithMetadata.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
    ## Restrict the tissue focus
    myT<-myT[myT$Sample.Type %in% tissueKept,]
                          colVector <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffff00", "#63dcfe", "#abd9e9", "#74add1", "#4575b4",  "#313695", "#0000ff")[as.numeric(as.factor(myT$Pen.Location))]


    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")

     pdf( paste(KEGGLevels, "_PICRUSt_", paste0(tissueKept,collapse="_"), "_feces_topMDS.pdf",sep=""))
     for (xrun in 1:4) {
                 for (yrun in 2:4) {
                     if(xrun == yrun){
                         break
                     }
                     plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", KEGGLevels,sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                           cex=2.0,
                          ## cex = log(myT$Age),

                          pch= ifelse(myT$Group=="Young", 0, 15) + ifelse(myT$Sample.Type=="LI Lumen", 0, ifelse(myT$Sample.Type =="LI Mucosa", 1, 2)),

                          col= colVector##[as.numeric(as.factor(myT$Pen.Location))]
                          ##col = brewer.pal(12, "Set3")[as.numeric(as.factor(myT$Pen.Location))]
                          )
                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Animal.ID,
                   cex=0.7, offset=0)
                      legend("topright", c("(Young) Feces", "(Old) Feces"),
## 		   ##  "Acute", "Chronic", "Control"),
                             pch=c(2, 17)#,
##                    ##             16, 16, 16),
##                    ##         col=c("black", "black", "black",
##                             ##             "blue", "orange", "yellow")
                             )
                 }
         }

         dev.off()

     write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_feces_PICRUSt_", KEGGLevels, paste0(tissueKept,collapse="_"), ".txt", sep=""))
     write.table(myMDS$CA$eig,file=paste("eigenValues_feces_PICRUSt_", KEGGLevels, paste0(tissueKept,collapse="_"), ".txt", sep=""), sep="\t")
}

tissueKept <- "LI Mucosa"


for(KEGGLevels in KEGGLevelsLevels )
{
    inFileName <- paste("PICRUSt_", KEGGLevels, "_LogNormwithMetadata.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")
    numCols <- ncol(myT)
    ## Restrict the tissue focus
    myT<-myT[myT$Sample.Type %in% tissueKept,]
                          colVector <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffff00", "#63dcfe", "#abd9e9", "#74add1", "#4575b4",  "#313695", "#0000ff")[as.numeric(as.factor(myT$Pen.Location))]


    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")

     pdf( paste(KEGGLevels, "_PICRUSt_", paste0(tissueKept,collapse="_"), "_mucosa_topMDS.pdf",sep=""))
     for (xrun in 1:4) {
                 for (yrun in 2:4) {
                     if(xrun == yrun){
                         break
                     }
                     plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", KEGGLevels,sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                           cex=2.0,
                          ## cex = log(myT$Age),

                          pch= ifelse(myT$Group=="Young", 0, 15) + ifelse(myT$Sample.Type=="LI Lumen", 0, ifelse(myT$Sample.Type =="LI Mucosa", 1, 2)),

                          col= colVector##[as.numeric(as.factor(myT$Pen.Location))]
                          ##col = brewer.pal(12, "Set3")[as.numeric(as.factor(myT$Pen.Location))]
                          )
                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Animal.ID,
                   cex=0.7, offset=0)
                      legend("topright", c("(Young) Mucosa", "(Old) Mucosa"),
## 		   ##  "Acute", "Chronic", "Control"),
                             pch=c(1, 16)#,
##                    ##             16, 16, 16),
##                    ##         col=c("black", "black", "black",
##                             ##             "blue", "orange", "yellow")
                             )
                 }
         }

         dev.off()

     write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_mucosa_PICRUSt_", KEGGLevels, paste0(tissueKept,collapse="_"), ".txt", sep=""))
     write.table(myMDS$CA$eig,file=paste("eigenValues_mucosa_PICRUSt_", KEGGLevels, paste0(tissueKept,collapse="_"), ".txt", sep=""), sep="\t")
}

tissueKept <- "LI Lumen"


for(KEGGLevels in KEGGLevelsLevels )
{
    inFileName <- paste("PICRUSt_", KEGGLevels, "_LogNormwithMetadata.txt", sep ="")
    myT <-read.csv(inFileName, header=TRUE, sep="", na.strings="BLAH")
    numCols <- ncol(myT)
    ## Restrict the tissue focus
    myT<-myT[myT$Sample.Type %in% tissueKept,]
                          colVector <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffff00", "#63dcfe", "#abd9e9", "#74add1", "#4575b4",  "#313695", "#0000ff")[as.numeric(as.factor(myT$Pen.Location))]


    myMDS <- capscale(myT[,2:(ncol(myT)-numMetadataCols)]~1,distance="bray")

     pdf( paste(KEGGLevels, "_PICRUSt_", paste0(tissueKept,collapse="_"), "_lumen_topMDS.pdf",sep=""))
     for (xrun in 1:4) {
                 for (yrun in 2:4) {
                     if(xrun == yrun){
                         break
                     }
                     plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", KEGGLevels,sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                           cex=2.0,
                          ## cex = log(myT$Age),

                          pch= ifelse(myT$Group=="Young", 0, 15) + ifelse(myT$Sample.Type=="LI Lumen", 0, ifelse(myT$Sample.Type =="LI Mucosa", 1, 2)),

                          col= colVector##[as.numeric(as.factor(myT$Pen.Location))]
                          ##col = brewer.pal(12, "Set3")[as.numeric(as.factor(myT$Pen.Location))]
                          )
                    textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Animal.ID,
                   cex=0.7, offset=0)
                      legend("topright", c("(Young) Lumen", "(Old) Lumen"),
## 		   ##  "Acute", "Chronic", "Control"),
                             pch=c(0, 15)#,
##                    ##             16, 16, 16),
##                    ##         col=c("black", "black", "black",
##                             ##             "blue", "orange", "yellow")
                             )
                 }
         }

         dev.off()

     write.table(myMDS$CA$u, sep="\t", file=paste("pcoa_lumen_PICRUSt_", KEGGLevels, paste0(tissueKept,collapse="_"), ".txt", sep=""))
     write.table(myMDS$CA$eig,file=paste("eigenValues_lumen_PICRUSt_", KEGGLevels, paste0(tissueKept,collapse="_"), ".txt", sep=""), sep="\t")
}
