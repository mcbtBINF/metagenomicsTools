rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")


setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")
taxaLevels <- c( "phylum", "class", "order", "family", "genus")

filePrefix <- "_energyContent_byPatient_4_"
mlm<- FALSE

divider <- 4
## divider <- 8

numMetadataCols <- 20

for(taxa in taxaLevels){
      	inFileName <- paste(taxa,"_LogNormalwithMetadataWeekly_NearestSamplewithDay_Edit_withDailyMetadata.txt", sep="")
	myT <-read.table(inFileName,header=TRUE,sep="\t")
	numCols <- ncol(myT)
	myColClasses <- c("character", rep("numeric", numCols - 1))
	myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)
	myT <- myT[ !is.na(myT$Day), ]

	patientPValues <- vector()
	timePValues <- vector()
	interactionPValues <- vector()
	names <-vector()

        BMIPValues <- vector()
        REEPValues <- vector()
        DITPValues <- vector()
        PAEPValues <- vector()
        PAENormPValues <- vector()
        AEEPValues <- vector()
        energyPValues <- vector()
        patient <- vector()
        EnergyContentPValues <- vector()
        pValuesCage <- vector()
        iccCage <- vector()

        allpvals <- list()

        myLmpVal <- list()

	colors <- vector()
	cIndex <-1
	for ( j in 1: nrow(myT))
	{
		if( substr(myT[j,]$Sample.ID,1,1) == "B")
		{
                    colors[cIndex] <- "Blue"
                    patient[cIndex] <- "B"
		} else if (substr(myT[j,]$Sample.ID,1,1) == "C" )
		{
                    colors[cIndex] <- "Red"
                    patient[cIndex] <- "C"
		} else
		{
                    colors[cIndex] <- "Black"
                    patient[cIndex] <- "A"
		}
		cIndex = cIndex + 1

	}

	index <-1

        pdf( paste(taxa, filePrefix, "boxplots.pdf", sep=""))

##        otherTaxa <- myT[,((numCols - 2):numCols)]
##        myT <- cbind(myT[,1:(numCols - numMetadataCols)], otherTaxa, myT[,(numCols - numMetadataCols + 1):(numCols - 3)])
        ## Previous version let the Fungus and etc through.
##            for( i in 2:(ncol(myT) - (numMetadataCols - 3)))

for( i in 2:(ncol(myT) - (numMetadataCols)))
        if( sum(myT[,i] != 0 ) > nrow(myT) / divider ){
            ## Easy access names
            DIT <- myT$DIT..kcal.day.
            REE <- myT$REE..kcal.day.
            ActiveEE <- myT$Active.EE..kcal.
            DITREE <- DIT / REE

            EnergyContent <- myT$Energy.Content..cal.g
            Day <- myT$Day
            bug <- myT[,i]
            BMI <- myT$Imputed.BMI
            EnergyIntake <- myT$Energy.Intake..kcal.day.
            names[index] = names(myT)[i]
            myFrame <- data.frame(bug, patient, Day, EnergyContent, BMI, EnergyIntake)

            ## myLm <- lm( bug ~ Day*patient*EnergyContent)
            ## myLmpVal[index] <- list(summary(myLm)$coefficients[,4][-1])
            ## myAnova <- anova(myLm)
            ## names[index] = names(myT)[i]
            ## Some kind of vector for the variables of interest

            if(mlm == TRUE){
                fullModel <- gls( bug~ EnergyContent, method="REML",correlation=corCompSymm(form=~1|factor(patient)),	data = myFrame )
                reducedModel <- gls( bug~ EnergyContent, method="REML", data = myFrame )
                ## reducedModel <- lme( bug ~ treatment + batch, method="REML", data = myFrame)
                fullModelLME <- lme(bug~ EnergyContent, method="REML", random = ~1|factor(patient), data = myFrame)
                ## fullModelLME <- lme(bug~ treatment + batch, method="REML", random= list(group = ~1, cage = ~1), data = myFrame)

            }
            else{
                fullModelLME <- lm(bug~ patient*EnergyContent, x=TRUE)
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
                allNames<-c(allNames, "patient", "iccpatient")
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

    for(k in 1:(orgCol - mlm)){
        dFrame[, k + orgCol]<-p.adjust(dFrame[,k], method = "BH")
    }
    colnames(dFrame)[(orgCol + 1):length(dFrame)]<-paste0(rep("adjusted_"),colnames(dFrame)[1:(orgCol - mlm)])

    dFrame <- cbind(names, dFrame)
    ## Try and dynamically generate the name of the plot here...
    write.table(dFrame, file=paste("pValues", filePrefix, taxa, ".txt",sep=""), sep="\t",row.names=FALSE)

    ## Get the sig table working
    ## It would also be nice to have a sig-picker for the plots and a combination of all sig results into one plot.
    ## keepVector <- grep("adj", names(dFrame))
    ## sigdFrame <-dFrame[which(dFrame[,keepVector] < 0.05,arr.ind=TRUE),]
    ## write.table(sigdFrame, file=paste("SIGpValues", filePrefix, taxa, ".txt",sep=""), sep="\t",row.names=FALSE)

    dev.off()
}
