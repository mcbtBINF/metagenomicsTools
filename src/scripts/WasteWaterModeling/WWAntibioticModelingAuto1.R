rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
library("vegan")

## Testing with another comment
## Load Support Packages
## Load Data
## Visualize Residuals
## Hypothesize Models
## Will have a spatial and time component.
## conc ~ space*time
## Run Simple Models
### taxa ~ antibioticConc + (1|timePoint)

## Model 1
## SummedAntibiotic Concentration ~ Creek

## Model 2
## taxa ~ Environmental Features
>
## Model 3
## diversity ~ SummedAntibiotic Concentration

## Model 4
## MDS Axis ~ SummedAntibiotic Concentration

setwd("/Users/mbrown67/Documents/Fodor/Datasets/WasteWater/")

tp <- c(1, 2, 3)

for(tpIter in tp)
{
    inFileName <- paste( "MetaMaster_tp", tpIter, ".csv", sep ="")
    myT <-read.csv(inFileName,header=TRUE, na.strings="BLAH", stringsAsFactors=FALSE)

    ## Replace "ND" cells

    ## Replace "<LOQ" cells

    antibioticConcentration <- myT[,1:10]
    myTNoSeq <- myT[,1:31]
    myTNoRepSeqNoSeq <- myT[1:66, 1:31]
    myTAntibioticOnly <- subset(myTNoRepSeqNoSeq, !duplicated(SampleNumber))

    myT <- myTAntibioticOnly[,1:10]

    ## This replacement will have to be done more carefully later
    myT[myT=="ND"|myT=="<LOQ"] <- 0

    myT<-data.matrix(myT)

    myTAntibioticOnly[1:22, 1:10] <- myT

    myT <- myTAntibioticOnly

    numMetadataCols <- 21

    ## Derived variables
    myT$Date <- unlist(lapply(strsplit(myT$Sample.Date.Time, split = " - "), '[[', 1))
    myT$TimeofDay <- unlist(lapply(strsplit(myT$Sample.Date.Time, split = " - "), '[[', 2))
    myT$Waterway <- unlist(lapply(strsplit(myT$Sample.Location, split = " Creek"), '[[', 1))

    ## Set up vectors
    pValuesWaterway <- vector()
    pValuesDate <- vector()
    pValuesTimeofDay <- vector()
    pValuespH <- vector()
    pValuesConductivity <- vector()
    pValuesHumidity <- vector()
    pValuesAmbientTemp <- vector()
    pValuesStorageTemp <- vector()
    pValuesSampleTemp <- vector()
    pValuesConcentration <- vector()
    pValuesBatch <- vector()
    pValuesQueue <- vector()
    pValuesLatitude <- vector()
    pValuesLongitude <- vector()
    pValuesLocation <- vector()
    pValuesWaterway <- vector()
    pValuesfeIter <- vector()
    ## Build Modeling Dataframe

    index <- 1


    ## Do the modeling for each antibiotic
    for( i in 1:(ncol(myT)-numMetadataCols)){
        ## if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
        ## {
        abc <- as.numeric(myT[,i])
        batch <- myT$Batch
        queue <- myT$Queue
        waterway <- myT$Waterway
        lat <- myT$Lattitude
        long <- myT$Longitude
        ambTemp <- myT$Ambient.Temp.Celsius
        storTemp <- myT$Storage.Temp.Celsius
        sampTemp <- myT$Sample.Temp.Celsius
        conduct <- myT$Conductivity..pHmv.
        pH <- myT$pH
        humidity <- myT$Humidity
        date <- myT$Date
        timeofDay <- myT$TimeofDay

        myFrame <- data.frame(abc, batch, queue, lat, long, ambTemp, storTemp, sampTemp, conduct, pH, humidity, date, timeofDay, waterway)

        ## Single Fixed Effect Models with waterway as the random effect
        ## Do the modeling for each antibiotic for each possible single fixed effect variable
        ## the last is excluded because it is a random effect
        for(feIter in 2:(dim(myFrame)[2] - 1))
            {
                fullModel <- gls( abc~  myFrame[,feIter], method="REML",correlation=corCompSymm(form=~1|factor(waterway)),				data = myFrame )
                reducedModel <- gls( abc~  myFrame[,feIter], method="REML", data = myFrame )
                fullModelLME <- lme( abc~  myFrame[,feIter], method="REML", random = ~1|factor(waterway), data = myFrame)

                ## Model Assessment
                ## pValuesWaterway[index] <- anova(fullModelLME)$"p-value"[5]
                ## pValuesAcuteChronic[index] <- anova(fullModelLME)$"p-value"[5]
                ## pValuesSex[index] <- anova(fullModelLME)$"p-value"[3]
                ## pValuesExperiment[index] <- anova(fullModelLME)$"p-value"[5]
                ## pValuesBatch[index] <- anova(fullModelLME)$"p-value"[5]
                ## pValuesTreatment[index] <- anova(fullModelLME)$"p-value"[2]
                ## pValuesInteraction[index] <- anova(fullModelLME)$"p-value"[5]
                pValuesfeIter[index] <- anova(fullModelLME)$"p-value"[2]

                pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
                iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]

                ## Graphing
                pdf( paste(tpIter, "_",names(myT)[feIter],"_boxplots.pdf", sep=""))
                names[index] = names(myT)[i]

                graphMain =  paste( names(myT)[i], "\n",
                    names(myT)[feIter], format( pValuesfeIter[index], digits=3),
                    " pCage= " , format( pValuesCage[index], digits=3),
                    " icc= " , format( iccCage[index], digits=3 ), sep="")

                par(mfrow=c(3,1),
                    oma = c(1,1,0,0) + 0.1,
                    mar = c(1,4,2.5,0) + 0.1)

                ## plot( bug ~ factor(treatment), ylab = names[index],main = graphMain )
                ## points(factor(sex), bug)

                ## plot( bug ~ factor(c( paste( myT$Sex, myT$Treatment,sep=""))))
                ## points(factor(c( paste( myT$Sex, myT$Treatment,sep=""))), bug)

                plot( abc ~ factor(waterway), ylab=names[index])
                points(factor(waterway), abc)
                index = index + 1

                ## MHC

                dFrame <- data.frame( names, pValuesfeIter, pValuesCage, iccCage)
                dFrame$adjustedfeIter<- p.adjust( dFrame$pValuesfeIter, method = "BH" )
                dFrame$adjustedCage <- p.adjust( dFrame$pValuesCage, method = "BH" )
            }
        ## }

        ## Output p-values
        write.table(dFrame, file=paste("pValuesAntibioticConcentrationat",tpIter,"_BYcage.txt",sep=""),
                    sep="\t",row.names=FALSE)
    }
        dev.off()
}
