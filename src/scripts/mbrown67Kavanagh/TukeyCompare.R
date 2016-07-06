rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/KylieData/IntestinalAging/rdpClassifications")
taxaLevels <- c( "phylum", "class", "order", "family", "genus")

tissueKept <- c("LI Lumen", "LI Mucosa", "Feces")
fullSampled <- c("1077", "1079", "1087", "1092", "1127", "1150", "1164", "1223", "1266", "1316", "1342", "1365", "1406", "1474")
missingSampled <- c("1227", "1349", "1377", "1415", "1452")
LuMu <- c(fullSampled, "1349", "1415", "1452")
## tissueKept <- "LI Lumen"
## tissueKept <- "LI Mucosa"
## tissueKept <- "Feces"

##  testo<-aov(bug ~ age + sampType, data=myFrame)
## TukeyHSD(testo, which="sampType", ordered=TRUE)
## and this shows there being a significant difference between the types
## but when you restrict the analysis to only the animals with all of the samples, then you get that only Feces and Mucosa are significantly different...
## Including the other 3 samples is enough to really drive down the pvalues.


filePrefix <- paste0(c(tissueKept, "group_sampType_bycage_"), collapse="_")
mlm<- TRUE

divider <- 4

for(taxa in taxaLevels){

    inFileName <- paste( taxa, "_LogNormwithMetadata_R1.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")

    numCols <- ncol(myT)
    numMetadataCols <- 6

    names <- vector()
    pValuesSex<- vector()
    pValuesGroup<- vector()
    pValuesCage<- vector()
    iccCage <- vector()
    iccKeptCounts <- vector()
    pValuesAge<- vector()
    pValuesAnimal<- vector()
    pValuesInteraction <- vector()
    allpvals <- list()

    index <- 1

    pdf( paste(taxa, filePrefix, "boxplots.pdf", sep=""))

    myT<-myT[myT$Sample.Type %in% tissueKept,]
    myT<-myT[myT$Animal.ID %in% fullSampled,]

    for( i in 2:(ncol(myT) - numMetadataCols))
        if( sum(myT[,i] != 0 ) > nrow(myT) / divider ){
            ## Easy access names
            bug <- myT[,i]
            animal <- myT$Animal.ID
            sampType <- myT$Sample.Type
            ## sex <- myT$Sex
            cage <- myT$Pen.Location
            group <- myT$Group
            age <- myT$Age
            names[index] = names(myT)[i]
            ## Some kind of vector for the variables of interest
            myFrame <- data.frame(bug, age, group, sampType, cage, animal)
            if(mlm == TRUE){
                fullModel <- gls( bug~ group + sampType, method="REML",correlation=corCompSymm(form=~1|cage),	data = myFrame )
                reducedModel <- gls( bug~ group + sampType, method="REML", data = myFrame )

                fullModelLME <- lme(bug~ group + sampType, method="REML", random = ~1|cage, data = myFrame)

            }
            else{
                fullModelLME <- lm(bug~ group + sampType, x=TRUE)
            }
            ## Introduce the goodness of fit tests here
            ## Can streamline vectors to be automatically responsive to changes to the model
            if(mlm == TRUE){
                allNames <- rownames(anova(fullModelLME))[-1]
                allpvals[[index]] <- anova(fullModelLME)$"p-value"[-1]
            }
            else {
                allNames <- rownames(anova(fullModelLME))[-dim(anova(fullModelLME))[1]]
                allpvals[[index]]<-as.data.frame(anova(fullModelLME))[1:dim(anova(fullModelLME))[1]-1,5]
            }

            ## Random Effects and Interclass Correlation Coefficient
            if(mlm == TRUE){
                allNames<-c(allNames, "cage", "iccCage")
                allpvals[[index]]<-c(allpvals[[index]], anova(fullModelLME, reducedModel)$"p-value"[2], coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]])
                pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
                ## This is different as it is not corrected for multiple hypothesis testing.
                iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]
                ## iccCage[index] <- as.numeric(VarCorr(fullModelLME)[1,1])/ sum(as.numeric(VarCorr(fullModelLME)[,1])
            }
    ## Nice plotting
            ## Correct so that it uses corrected p-values
            if(mlm == TRUE){
                graphMain = c(names[i], paste(c(rep("p-",length(allNames)-1),""), allNames, " = ", format(allpvals[[index]], digits=3), sep=""))
            }
            else {
                graphMain = c(names[i], paste(c(rep("p-",length(allNames)),""), allNames, " = ", format(allpvals[[index]], digits=3), sep=""))

            }

            ## mlm is boolean 1/0 value. Basically, don't graph on icc
            ## need a correction for interaction terms as well
            numPlots <- length(allNames) - mlm
            par(mfrow=c(numPlots, 1),
                oma = c(1,1,0,0) + 0.1,
                mar = c(1,4,3,0) + 0.1)

            for(j in 1:numPlots){
                if(j == 1){
                    plot( bug ~ factor(get(allNames[j])), ylab = names[i], main = c(names[index],toString(graphMain[-1])) )
                    points(factor(get(allNames[j])), bug)
                }
                else{
                    testInteract <- grep(":", allNames)
                    if(length(testInteract) > 0){
                        if(j %in% testInteract){
                            intSplit <- unlist(strsplit(allNames[grep(":",allNames)[which(grep(":", allNames) == j)]],split=":"))
                            if(length(intSplit) == 2){
                                plot( bug ~ factor(c( paste( get(intSplit[1]), get(intSplit[2]),sep=""))))
                                points(factor(c( paste( get(intSplit[1]), get(intSplit[2]),sep=""))), bug)
                            }
                            else{
                                plot( bug ~ factor(c( paste( get(intSplit[1]), get(intSplit[2]), get(intSplit[3]),sep=""))))
                                points(factor(c( paste( get(intSplit[1]), get(intSplit[2]), get(intSplit[3]),sep=""))), bug, col=ifelse(group == "Old", "BLUE", "RED"))
                            }
                        }
                    }
                    else{
                        plot( bug ~ factor(get(allNames[j])), ylab = names[i])
                        points(factor(get(allNames[j])), bug, col=ifelse(group == "Old", "BLUE", "RED"))
                    }
                }
            }

            index <- index + 1
	}

    ## This should also be dynamic
    dFrame <- as.data.frame(matrix(unlist(allpvals), nrow= length(allpvals), byrow=TRUE))
    colnames(dFrame) <- unlist(allNames)
    rownames(dFrame) <- names
    orgCol <- ncol(dFrame)

    for(k in 1:(orgCol-mlm)){
        dFrame[, k + orgCol]<-p.adjust(dFrame[,k], method = "BH")
    }
    colnames(dFrame)[(orgCol + 1):length(dFrame)]<-paste0(rep("adjusted_"),colnames(dFrame)[1:(orgCol - mlm)])

    dFrame <- cbind(names, dFrame)
    ## Try and dynamically generate the name of the plot here...
    write.table(dFrame, file=paste("pValues", filePrefix, taxa, ".txt",sep=""), sep="\t",row.names=FALSE)

    ## Get the sig table working
    ## It would also be nice to have a sig-picker for the plots and a combination of all sig results into one plot.
    keepVector <- grep("adj", names(dFrame))
    sigdFrame <-dFrame[which(dFrame[,keepVector] < 0.05,arr.ind=TRUE),]
    write.table(sigdFrame, file=paste("SIGpValues", filePrefix, taxa, ".txt",sep=""), sep="\t",row.names=FALSE)

    dev.off()
}
