rm(list=ls())

library("vegan")
library("calibrate")

baseDir <- "/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging"
setwd(baseDir)
## Assignment type
dataType <- "closedQIIMER1"
metadataDir <- paste(baseDir, "metadata", sep="/")
dataDir <- paste(baseDir, dataType, sep="/")
processedDir <- paste(dataDir, "processed", sep="/")
analysisDir <- paste(baseDir, dataType, "analysis", sep="/")
baseDataFileName <- "closed_reference_otu_table_L"

taxaLevels <- c(2:7)

for (taxa in taxaLevels){
    ## taxa <- 2
    ## Change over to data directory
    setwd(processedDir)
    inFileName <- paste(baseDataFileName, taxa, "_LogNormwithMetadata.txt", sep="")
    myT <- read.table(inFileName, header=TRUE, sep="\t")

    numCols <- ncol(myT)

    ## endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    ## The 9 could be done better...
    myColClasses <- metaColClasses<- c("character", "character", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", rep("numeric", numCols - 9))
    myT <- read.table(inFileName, header=TRUE, sep="\t", colClasses = myColClasses)

    endMetadataIndex <- which(colnames(myT) == "depthAtLevel")
    setwd(analysisDir)

    ## Get Tegan's help to fix these colors
    colVector <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4",  "#313695", "#0000ff")[as.numeric(as.factor(myT$Pen.Location))]

    ## pdf and eigenvectors and eigenvalues

    myMDS <- capscale(myT[,(endMetadataIndex +1):ncol(myT)]~1,distance="bray")
    percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

    ## Default for all tissues

    write.table(myMDS$CA$u, sep="\t", file=paste(baseDataFileName, taxa, "_pcoa.txt",sep=""), row.names = FALSE)
    write.table(myMDS$CA$eig,file=paste(baseDataFileName, taxa, "_eigenValues.txt", sep=""), sep="\t")

    pdf( paste(baseDataFileName, taxa, "_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
        for (yrun in 2:4) {
            if(xrun == yrun){
                break
            }
            par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)#, mfrow = c(1,1))

            plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
                 xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
                 ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
                 main=paste("PCoA at level:", taxa,sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                 cex=2.0,
                 pch= ifelse(myT$Group=="Young", 0, 15) + ifelse(myT$Sample.Type=="LI Lumen", 0, ifelse(myT$Sample.Type =="LI Mucosa", 1, 2)),
                 col= colVector
                 )
            textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Animal.ID,
                   cex=0.7, offset=0)
                        par(xpd=TRUE)
            legend("topright", inset=c(-0.4,0),
                   c("(Young) LI Lumen", "(Young) LI Mucosa", "(Young) Feces", "(Old) LI Lumen", "(Old) LI Mucosa", "(Old) Feces"),
                   pch=c(0, 1, 2, 15, 16, 17)#,
                   ##                    ##         col=c("black", "black", "black",
                   ##                             ##             "blue", "orange", "yellow")
                   )
        }
    }

    dev.off()



    pdf( paste(baseDataFileName, taxa, "_betterViz_topMDS.pdf",sep=""))
    for (xrun in 1:4) {
        for (yrun in 2:4) {
            if(xrun == yrun){
                break
            }
            par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)#, mfrow = c(1,1))

            plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
                 xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
                 ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
                 main=paste("PCoA at level:", taxa,sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                 cex=2.0,
                 pch= ifelse(myT$Sample.Type=="LI Lumen", 15, ifelse(myT$Sample.Type =="LI Mucosa", 16, 17)),
                 ## col= colVector
                 col=ifelse(myT$Group=="Young", "blue", "red")
                 )
            ## textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Animal.ID,
            ##        cex=0.7, offset=0)
                        par(xpd=TRUE)
            legend("topright", inset=c(-0.4,0),
                   c("Young", "Old", "LI Lumen", "LI Mucosa", "Feces"),
                   pch=c(16, 16, 15, 16, 17),
                   col=c("blue", "red", "black", "black", "black")
                   )
        }
    }

    dev.off()

pdf( paste(baseDataFileName, taxa, "_lookcageViztopMDS.pdf",sep=""))
    for (xrun in 1:4) {
        for (yrun in 2:4) {
            if(xrun == yrun){
                break
            }
            par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)#, mfrow = c(1,1))

            plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
                 xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
                 ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
                 main=paste("PCoA at level:", taxa,sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                 cex=2.0,
                 pch= ifelse(myT$Sample.Type=="LI Lumen", 15, ifelse(myT$Sample.Type =="LI Mucosa", 16, 17)),
                 ## col= colVector
                 col=ifelse(myT$Pen.Location=="315", "blue", ifelse(myT$Pen.Location =="322", "red", ifelse(myT$Pen.Location == "A10", "orange", "black")))
                 )
            ## textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myT$Animal.ID,
            ##        cex=0.7, offset=0)
                        par(xpd=TRUE)
            legend("topright", inset=c(-0.4,0),
                   c("315", "322", "A10","LI Lumen", "LI Mucosa", "Feces"),
                   pch=c(16, 16, 16, 15, 16, 17),
                   col=c("blue", "red", "orange", "black", "black", "black")
                   )
        }
    }

    dev.off()


    ## For each "Sample Type" (These are tissues)
    for (tissue in unique(myT$Sample.Type)){
        myTtissue <- myT[myT$Sample.Type == tissue,]
        myMDS <- capscale(myTtissue[,(endMetadataIndex +1):ncol(myT)]~1,distance="bray")
        percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

        pdf( paste(baseDataFileName, taxa,"_",gsub('([[:punct:]])|\\s+','_',tissue), "_topMDS.pdf",sep=""))
        write.table(myMDS$CA$u, sep="\t", file=paste(baseDataFileName, taxa, "_", gsub('([[:punct:]])|\\s+','_',tissue), "_pcoa.txt",sep=""), row.names = FALSE)
        write.table(myMDS$CA$eig,file=paste(baseDataFileName, taxa, "_", gsub('([[:punct:]])|\\s+','_',tissue), "_eigenValues.txt", sep=""), sep="\t")

        for (xrun in 1:4) {
            for (yrun in 2:4) {
                if(xrun == yrun){
                    break
                }
                plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
                     xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
                     ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
                     main=paste("PCoA at level:", taxa, " for ", tissue, sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                     cex=2.0,
                     pch= ifelse(myTtissue$Group=="Young", 0, 15) + ifelse(myTtissue$Sample.Type=="LI Lumen", 0, ifelse(myTtissue$Sample.Type =="LI Mucosa", 1, 2)),
                     col= colVector
                     )
                textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myTtissue$Animal.ID,
                       cex=0.7, offset=0)
                legend("topright", c("(Young) LI Lumen", "(Young) LI Mucosa", "(Young) Feces", "(Old) LI Lumen", "(Old) LI Mucosa", "(Old) Feces"),
                       pch=c(0, 1, 2, 15, 16, 17)#,
                       ##                    ##         col=c("black", "black", "black",
                       ##                             ##             "blue", "orange", "yellow")
                       )
            }
        }
        dev.off()

    }
    for (tissue in unique(myT$Sample.Type)){
        myTtissue <- myT[myT$Sample.Type == tissue,]
        myMDS <- capscale(myTtissue[,(endMetadataIndex +1):ncol(myT)]~1,distance="bray")
        percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

        pdf( paste(baseDataFileName, taxa,"_",gsub('([[:punct:]])|\\s+','_',tissue), "_topMDS.pdf",sep=""))
        write.table(myMDS$CA$u, sep="\t", file=paste(baseDataFileName, taxa, "_", gsub('([[:punct:]])|\\s+','_',tissue), "_pcoa.txt",sep=""), row.names = FALSE)
        write.table(myMDS$CA$eig,file=paste(baseDataFileName, taxa, "_", gsub('([[:punct:]])|\\s+','_',tissue), "_eigenValues.txt", sep=""), sep="\t")

        for (xrun in 1:4) {
            for (yrun in 2:4) {
                if(xrun == yrun){
                    break
                }
                plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
                     xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
                     ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
                     main=paste("PCoA at level:", taxa, " for ", tissue, sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                     cex=2.0,
                     pch= ifelse(myTtissue$Group=="Young", 0, 15) + ifelse(myTtissue$Sample.Type=="LI Lumen", 0, ifelse(myTtissue$Sample.Type =="LI Mucosa", 1, 2)),
                     col= colVector
                     )
                textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myTtissue$Animal.ID,
                       cex=0.7, offset=0)
                legend("topright", c("(Young) LI Lumen", "(Young) LI Mucosa", "(Young) Feces", "(Old) LI Lumen", "(Old) LI Mucosa", "(Old) Feces"),
                       pch=c(0, 1, 2, 15, 16, 17)#,
                       ##                    ##         col=c("black", "black", "black",
                       ##                             ##             "blue", "orange", "yellow")
                       )
            }
        }
        dev.off()

    }

    ## For each age grouping
    for (ageGroup in unique(myT$Group)){
        myTageGroup <- myT[myT$Group == ageGroup,]
        myMDS <- capscale(myTageGroup[,(endMetadataIndex +1):ncol(myT)]~1,distance="bray")
        percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))

        pdf( paste(baseDataFileName, taxa, "_", ageGroup, "_topMDS.pdf",sep=""))
        write.table(myMDS$CA$u, sep="\t", file=paste(baseDataFileName, taxa, "_", ageGroup, "_pcoa.txt",sep=""), row.names = FALSE)
        write.table(myMDS$CA$eig, file=paste(baseDataFileName, taxa, "_", ageGroup, "_eigenValues.txt", sep=""), sep="\t")

        for (xrun in 1:4) {
            for (yrun in 2:4) {
                if(xrun == yrun){
                    break
                }
                plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],
                     xlab=paste("MDS Axis ", xrun, ": ", format(percentVariance[xrun],digits=3), sep=""),
                     ylab=paste("MDS Axis ", yrun, ": ", format(percentVariance[yrun],digits=3), sep=""),
                     main=paste("PCoA at level:", taxa, " for ", ageGroup, sep=""), ##xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
                     cex=2.0,
                     pch= ifelse(myTageGroup$Group=="Young", 0, 15) + ifelse(myTageGroup$Sample.Type=="LI Lumen", 0, ifelse(myTageGroup$Sample.Type =="LI Mucosa", 1, 2)),
                     col= colVector
                     )
                textxy(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labs=myTageGroup$Animal.ID,
                       cex=0.7, offset=0)
                legend("topright", c("(Young) LI Lumen", "(Young) LI Mucosa", "(Young) Feces", "(Old) LI Lumen", "(Old) LI Mucosa", "(Old) Feces"),
                       pch=c(0, 1, 2, 15, 16, 17)#,
                       ##                    ##         col=c("black", "black", "black",
                       ##                             ##             "blue", "orange", "yellow")
                       )
            }
        }
        dev.off()

    }

}
