rm(list=ls())
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")

library("Kendall")
library("vegan")
library("lmtest")
library("pscl")
library("nlme")
library("gtools")

taxaLevels <- c("phylum","class","order","family","genus")

index<-1
safeindex <- vector()
for(t in taxaLevels)
{

inFileName <- paste0("pivoted_",t, "asColumns.txt")
myTCheck <- read.table(inFileName, header=TRUE, sep="\t")
rownames(myTCheck) <- myTCheck$sample
# Remove the row labels
myTCheck <- myTCheck[,-1]
# Restricted to the forward reads now.
myTCheck <- myTCheck[1:147,]
mySums <- rowSums(myTCheck)
colors <- rep("Black", dim(myTCheck)[1])
## colors <- ifelse(grep("A", rownames(myTCheck)), "Black", ifelse(grep("B", rownames(myTCheck)),
colors[grep("B", rownames(myTCheck))] <- c("Blue")
colors[grep("C", rownames(myTCheck))] <- c("Red")

Shannon <- apply(myTCheck[,1:(ncol(myTCheck))], 1, diversity)

countCol <- cbind(mySums, colasFactors)
shanCol <- cbind(Shannon, colasFactors)

countCol <- as.data.frame(countCol)
shanCol <- as.data.frame(shanCol)

boxplot(mySums~colasFactors, countCol)
boxplot(Shannon~colasFactors, shanCol)

colasFactors<-as.numeric(factor(colors))



nzCount <- function(x) {
    return(sum(x != 0))
}

nonZero <- apply(myTCheck, 1, nzCount)
# Should probably drop these just to be careful...
# rownames(myTCheck[myTCheck$rowSums < 10000,])
                                        #[1] "r1_37" "r1_45" "r1_52" "r1_58" "r1_7"  "r1_9"

                                        #Get the counts:
#length(which(myTCheck[,1] > 0))
                                        #-1
                                        #-1
presence<-apply(myTCheck, 2, function (x) length(which(x > 0)))
# -1 is to remove the rowSums column.
safeindex[index]<-length(which(presence > 35)) - 1
index<-index + 1
}
