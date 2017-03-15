rm(list=ls())
library("lmtest")
library("nlme")
library("pscl")

## Addresses item 1 from the email: removing Sample 712 from the analysis

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/BombCalorimetry/")

taxaLevels <- c("phylum","class","order","family","genus")
index <- 1
ShannonP <- list()
ShannonSummary <- list()
names <- vector()
myLmpVal <- list()
pValuesPatientID <- vector()
pValuesCalorimetry <- vector()
pValuesTime <- vector()
intraclassCoefficient <- vector()

for(taxa in taxaLevels )
{
    inFileName <- paste( taxa,  "_paired_metadata.txt", sep ="")
    myT <-read.table(inFileName, header=TRUE,sep="\t")
    numCols <- ncol(myT)
    ##	myColClasses <- c(rep("character",1), rep("numeric", numCols-1))
    ##	myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)

    ##	myT <- myT[ myT$timepoint == "2" &  ! is.na(myT$calorimetryData), ]
    myT <- myT[myT$Sample != 712,]
    ##        myT <- myT[myT$Time == 1,]
    ## names <- vector()
    ## pValuesPatientID <- vector()
    ## pValuesCalorimetry <- vector()
    ## pValuesTime <- vector()
    ## meanBug <- vector()

    ## index <- 1

    pdf( paste(taxa, "_BothTimes_Absolute_Shannon.pdf", sep=""))
    myT$Shannon <- apply(myT[,3:(ncol(myT)-7)],1,diversity)
    Shannon <- myT$Shannon
    time <- factor(myT$Time)
    patientID <- myT$Sample
    calorimetry<- myT$cal.g
    myFrame <- data.frame(Shannon, time, patientID, calorimetry)
    ## Shannonlm <- lm(Shannon ~ calorimetry + time, x = TRUE)

    ## ShannonP[indexS] <- list(anova(Shannonlm)$"Pr(>F)"[1:3])
    fullModel <- gls( Shannon~  calorimetry + time,
                     method="REML",correlation=corCompSymm(form=~1|factor(patientID)),
                     data = myFrame )

    reducedModel <- gls( Shannon~  calorimetry + time, method="REML",	data = myFrame )

    fullModelLME <- lme(Shannon~  calorimetry + time, method="REML", random = ~1|factor(patientID), data = myFrame)

    pValuesCalorimetry[[index]] <- anova(fullModelLME)$"p-value"[2]
    pValuesTime[[index]] <- anova(fullModelLME)$"p-value"[3]
    pValuesPatientID[[index]] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
    intraclassCoefficient[[index]]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]
    names[[index]] <- taxa

    ## pValuesCalorimetry[index] <- anova(fullModel)$"Pr(>F)"[1]
    ## names[index] = names(myT)[i]

    ## graphMain =  paste( names(myT)[i], "pValuesTime=", format(pValuesTime[index], digits=3), "\n",
    ##     " pValuesCalorimetry= ", format(pValuesCalorimetry[index],digits=3),
    ##     " pValuesPatientID= " , format(	pValuesPatientID[index], digits=3), "\n",
    ##     " icc= " , format( intraclassCoefficient, digits=3 ), sep="")

    ## plot( bug ~ calorimetry, ylab = names[index],
    ##      main = graphMain )

    ## indexS <- indexS + 1
    index = index + 1

    ## for( i in 3:(numCols - 7))
    ##     if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
    ##         {
    ##                                     #bug <- log10( myT[,i] + 0.00001)
    ##                                     #meanBug[index] <- mean(bug)
    ##             bug <- myT[,i]
    ##             time <- factor(myT$Time)
    ##             patientID <- myT$Sample
    ##             calorimetry<- myT$cal.g

    ## myFrame <- data.frame(bug, time, patientID, calorimetry)

    ## ## fullModel <- lm( bug~  calorimetry)
    ## fullModel <- gls( bug~  calorimetry + time,
    ##                  method="REML",correlation=corCompSymm(form=~1|factor(patientID)),
    ##                  data = myFrame )

    ## reducedModel <- gls( bug~  calorimetry + time, method="REML",	data = myFrame )

    ## fullModelLME <- lme(bug~  calorimetry + time, method="REML", random = ~1|factor(patientID), data = myFrame)

    ## pValuesCalorimetry[index] <- anova(fullModelLME)$"p-value"[2]
    ## pValuesTime[index] <- anova(fullModelLME)$"p-value"[3]
    ## pValuesPatientID[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
    ## intraclassCoefficient<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]


    ## pValuesCalorimetry[index] <- anova(fullModel)$"Pr(>F)"[1]
    ## names[index] = names(myT)[i]

    graphMain =  paste( taxa)
    ## ,"_Shannon_", "pValuesTime=", format(pValuesTime[index], digits=3), "\n",
    ##     " pValuesCalorimetry= ", format(pValuesCalorimetry[index],digits=3),
    ##     " pValuesPatientID= " , format(	pValuesPatientID[index], digits=3), "\n",
    ##     " icc= " , format( intraclassCoefficient, digits=3 ), sep="")

    plot( Shannon ~ calorimetry, ylab = "Shannon",
         main = graphMain )

    ## index=index+1

##    dFrame <- data.frame( taxa, pValuesCalorimetry, pValuesTime, pValuesPatientID, intraclassCoefficient)
##    dFrame <- dFrame [order(dFrame$pValuesCalorimetry),]


##    write.table(dFrame, file=paste("pValuesFor_", taxa, "_Calorimetry_BothTimes_Patient_ABSOLUTE_Shannon.txt",sep=""), sep="\t",row.names=FALSE)
    dev.off()
}

dFrame <- data.frame(pValuesCalorimetry, pValuesTime, pValuesPatientID, intraclassCoefficient)
##    dFrame <- data.frame( myLmpVal)
    ## dFrame <- t(dFrame)
    dFramemyLm <- data.frame(names, dFrame)

    for (m in 2:(dim(dFramemyLm)[2]-1))
    {
        dFramemyLm[,dim(dFramemyLm)[2] + 1] <- p.adjust(dFramemyLm[,m], method = "BH")
        colnames(dFramemyLm)[ncol(dFramemyLm)]<-paste0("adj", colnames(dFramemyLm)[m])
    }

    write.table(dFramemyLm, file=paste("pValuesFor_", taxa, "_Calorimetry_BothTimes_Patient_ABSOLUTE_Shannon.txt",sep=""), sep="\t",row.names=FALSE)

