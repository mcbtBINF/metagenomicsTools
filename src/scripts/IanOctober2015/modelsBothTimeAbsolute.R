rm(list=ls())
library("lmtest")
library("nlme")
library("pscl")

## Addresses item 1 from the email: removing Sample 712 from the analysis

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/BombCalorimetry/")

taxaLevels <- c("phylum","class","order","family","genus")

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
    names <- vector()
    pValuesPatientID <- vector()
    pValuesCalorimetry <- vector()
    pValuesTime <- vector()
    meanBug <- vector()

    index <- 1

    pdf( paste(taxa, "_BothTimes_Absolute.pdf", sep=""))

    for( i in 3:(numCols - 7))
        if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
            {
                                        ## bug <- log10( myT[,i] + 0.00001)
                                        ## meanBug[index] <- mean(bug)
                bug <- myT[,i]
                time <- factor(myT$Time)
                patientID <- myT$Sample
                calorimetry<- myT$cal.g

                myFrame <- data.frame(bug, time, patientID, calorimetry)

                ## fullModel <- lm( bug~  calorimetry)
                fullModel <- gls( bug~  calorimetry + time,
				 method="REML",correlation=corCompSymm(form=~1|factor(patientID)),
                                 data = myFrame )

                reducedModel <- gls( bug~  calorimetry + time, method="REML",	data = myFrame )

                fullModelLME <- lme(bug~  calorimetry + time, method="REML", random = ~1|factor(patientID), data = myFrame)

                pValuesCalorimetry[index] <- anova(fullModelLME)$"p-value"[2]
                pValuesTime[index] <- anova(fullModelLME)$"p-value"[3]
                pValuesPatientID[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
                intraclassCoefficient<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]


                ## pValuesCalorimetry[index] <- anova(fullModel)$"Pr(>F)"[1]
                names[index] = names(myT)[i]

                graphMain =  paste( names(myT)[i], "pValuesTime=", format(pValuesTime[index], digits=3), "\n",
                    " pValuesCalorimetry= ", format(pValuesCalorimetry[index],digits=3),
                    " pValuesPatientID= " , format(	pValuesPatientID[index], digits=3), "\n",
                    " icc= " , format( intraclassCoefficient, digits=3 ), sep="")

                plot( bug ~ calorimetry, ylab = names[index],
                     main = graphMain )

                index=index+1

            }
    print(c(taxa, index))
    dFrame <- data.frame( names, pValuesCalorimetry, pValuesTime, pValuesPatientID, intraclassCoefficient)
    dFrame <- dFrame [order(dFrame$pValuesCalorimetry),]
    dFrame$adjustedpValuesCalorimetry <- p.adjust( dFrame$pValuesCalorimetry, method = "BH" )
    dFrame$adjustedpValuesTime <- p.adjust( dFrame$pValuesTime, method = "BH" )
    dFrame$adjustedpValuesPatientID <- p.adjust( dFrame$pValuesPatientID, method = "BH" )

    write.table(dFrame, file=paste("pValuesFor_", taxa, "_Calorimetry_BothTimes_Patient_ABSOLUTE.txt",sep=""), sep="\t",row.names=FALSE)
    dev.off()
}
