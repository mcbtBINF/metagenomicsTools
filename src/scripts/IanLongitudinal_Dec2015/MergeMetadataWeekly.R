rm(list=ls())

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrolData/Carroll_Longitudinal")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" );

for (taxa in taxaLevels)
{
    inFileName <- paste( "pivoted_", taxa, "asColumnsLogNormal_r1.txt", sep ="")
        myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c("character", rep("numeric", numCols-1))
	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)

    rownames(myT)<-sapply(strsplit(rownames(myT),split="_"),function(x) x[2])
    myT$Sample.ID<-rownames(myT)

    myMetadata <- read.table("MetadataWeeklyAll.txt", header=TRUE, sep="\t")

    merged<-merge(x = myT, y = myMetadata, by = "Sample.ID", all = TRUE)
    write.table(merged,file=paste(taxa,"LogNormalwithMetadataWeekly.txt", sep="" ),row.names=FALSE, sep="\t")
}
