## Plot (standardised) residuals against fitted values to assess homogeneity.
## Make a histogram of the residuals to verify normality. You can also use a QQ-plot.
## Plot the residuals against each explanatory variable that was used in the model. If you see a pattern, you are violating the independence assumption.
## Plot the residuals against each explanatory variable not used in the model. If you see a pattern, include the omitted explanatory variable and refit the model. If the residuals patterns disappear, include the term, even if it is not significant
### IMPORTANT!
## Asses the model for influential observations. A useful tool is the Cook distance function.

## What is a sufficiently large number of replicate samples?
## abline(M5$coef[1], M5$coef[3], col="RED")
## abline(M5$coef[1] + M5$coef[2], M5$coef[3] + M5$coef[4], col="BLUE")

## plot(model)
## E <- rstandard(model)
## hist(E)
## qqnorm(E)
## plot(y = E, x=myT$Variable)
## plot(E ~ myT$Variable)
## Not sure why the next line is involved...
## abline(0,0)

## Beta <- vector()
## for( i in 2:8){
##     Mi <- summary(lm(bug ~ cage, subset = (group == i), data = myFrame))
##     Beta[i] <- Mi$coef[2,1]
## }

## Beta <- vector()
## for( i in 1:4){
##     Mi <- summary(lm(bug ~ cage, subset = (as.numeric(factor(date)) == i), data = myFrame))
##     Beta[i] <- Mi$coef[2,1]
## }


op <- par(mfrow = c(2,2), mar = c(5, 4, 1, 2))
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "Residuals", main = "")
plot(myT$var1, E, xlab = "Var1name", ylab = "Residuals")
plot(myT$va42, E, xlab = "Var2name", ylab = "Residuals")
par(op)

## Test for homogeneity
E1 <- E[myT$Var1 <= threshold]
E2 <- E[myT$Var1 > threshold]
var.test(E1, E2)

## Very sensitive to non-normality
bartlett.test(E, myT$Variable)
## Make histograms of the residuals per level
## Consider different variances per month...

## Leverage against residuals
## Measures whether any observation has extreme values of explanatory variables
## Used in the case of multiple explanatory variables
## How different it is compared to other explanatory variables.

## Cook distance is the leave-one out measure of influence
## If value is greater othatn one, kick it out.
## Influence on the estimated parameters

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/Pooled/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus")

filePrefix <- "_treatment_onebatch_onetime_"
mlm<- FALSE

dropBatch <- "Lyte_Batch04_Run01_repeat"
dropTime <- 14
dropTreatment <- "Exp"

for(taxa in taxaLevels){

    inFileName <- paste( taxa, "LogNormwithMetadata_R1_Pooled.txt", sep ="")
    myT <-read.csv(inFileName,header=TRUE,sep="", na.strings="BLAH")

    numCols <- ncol(myT)
    numMetadataCols <- 20

    ## Reprocessing for correct column data type
    myColClasses <- c("character", rep("numeric", numCols - numMetadataCols - 1), rep("character", 11), "numeric", "numeric", "numeric", "character", "character", "character", "character", "character", "character")
    myT <-read.csv(inFileName, header=TRUE, sep="", colClasses=myColClasses, na.strings="BLAH", comment.char ="@")

    ## Removing unwanted samples {controls, replicates, duplicates, poor quality samples}
    removeControls<-c( "C1", "C2", "N1", "N2", "Neg", "Pos")
    myT<-myT[!(myT$Sample_ID %in% removeControls),]

    removetrs<-c("04_125_tr", "04_101_tr", "04_103_tr", "04_74_tr", "04_70_tr", "04_40_tr", "04_41_tr", "04_84_tr")
    myT<-myT[!(myT$Sample_ID %in% removetrs),]

    removeDups<-c("04-04_S63_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% removeDups),]

    ## Controls not dropped earlier
    manualDrop <- c("Neg_S40_L001_R1_001", "PCR1Neg_S65_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% manualDrop),]

    labDrop <- c("Harlan Labs")
    myT<-myT[!(myT$MouseOrigin %in% labDrop),]

    lowSeqDrop <- c("04-55_S32_L001_R1_001")
    myT<-myT[!(myT$MatchFile %in% lowSeqDrop),]
    ## Not sure what to do about this though.
    ### myT$Condition[which(myT$Treatment == "Ctrl", arr.ind = TRUE)]<-"Control"

##    myT<-myT[!(myT$Sample_Plate == dropBatch),]
##    myT<-myT[myT$StressLength == dropTime,]

    names <- vector()
    pValuesSex<- vector()
    pValuesAcuteChronic<- vector()
    pValuesCage<- vector()
    iccCage <- vector()
    iccKeptCounts <- vector()
    pValuesBatch<- vector()
    pValuesTreatment<- vector()
    pValuesTime <- vector()
    pValuesStressLength <- vector()
    pValuesInteraction <- vector()
    allpvals <- list()

    index <- 1
    pdf( paste(taxa, filePrefix, "boxplots.pdf", sep=""))

    for( i in 2:(ncol(myT) - numMetadataCols))
        if( sum(myT[,i] != 0 ) > nrow(myT) / 4 ){
            ## Easy access names
            bug <- myT[,i]
            ac <- myT$Condition
            sex <- myT$Sex
            cage <- myT$Cage
            batch <- myT$Sample_Plate
            date <- myT$Date
            mo <- myT$MouseOrigin
            group <- myT$SpreadsheetGrouping
            treatment <- myT$Treatment
            time <- myT$StressLength
            multiWay <- paste(ac, myT$StressLength)
            multiWay[which(myT$Treatment == "Ctrl", arr.ind = TRUE)]<-"Control"
            names[index] = names(myT)[i]
            ## Some kind of vector for the variables of interest
            myFrame <- data.frame(bug, ac, time, multiWay, sex, cage, treatment, batch, date, mo, date)
            ## M1 <- lm(bug ~ sex + treatment + time + batch + group + date, data = myFrame)
            ## M2 <- lm(bug ~ sex + treatment + time + batch, data = myFrame)
            ## M3 <- lm(bug ~ sex*treatment*time*batch, data = myFrame)
            ## M4 <- lm(bug ~ treatment*batch + sex*time*batch, data = myFrame)
            ## drop1(M1, test=
            if(mlm == TRUE){
                fullModel <- gls( bug~ treatment, method="REML",correlation=corCompSymm(form=~1|factor(group)),	data = myFrame )
                reducedModel <- gls( bug~ treatment, method="REML", data = myFrame )
                fullModelLME <- lme(bug~ treatment, method="REML", random = ~1|factor(group), data = myFrame)
                ## fullModelLME <- lme(bug~ treatment + batch, method="REML", random= list(group = ~1, cage = ~1))
            }
            else{
                fullModelLME <- lm(bug~ treatment, x=TRUE)
            }
            ## Potential save time by reducing anova calls
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
                                points(factor(c( paste( get(intSplit[1]), get(intSplit[2]), get(intSplit[3]),sep=""))), bug)
                            }
                        }
                    }
                    else{
                        plot( bug ~ factor(get(allNames[j])), ylab = names[i])
                        points(factor(get(allNames[j])), bug)
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
