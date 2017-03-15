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
## taxa ~ SummedAntibiotic Concentration

## Model 2
## taxa ~ Environmental Features

## Model 3
## diversity ~ SummedAntibiotic Concentration

## Model 4
## MDS Axis ~ SummedAntibiotic Concentration

## Makes this specific to the analysis at hand
##setwd("/Users/mbrown67/Documents/Fodor/Datasets/WasteWater/closed_taxonomy_tp_1_2_3/")
setwd("/Users/mbrown67/Documents/Fodor/Datasets/Wastewater/closed4timepoints/closed_reference/convertedbyTaxonomicLevel/")

## Later will include other timepoints
tps <- c(1, 2, 3, 4)
## Give these meaningful names. Does 7 correspond to "Species"/"OTU"?
taxaLevels <- c(2, 3, 4, 5, 6, 7)

## Make this dynamic
numMetadataCols <- 40
mlm <- FALSE
## Ideally, this should pick out the variables to be analyzed
## Can parse from text entry to match up with the names in the data.frame bizness.
filePrefix <- "NamingConventionHere"

## Prescreen for largest technical replicate via the sequence count of the 2 level
## Also kick out too-small samples here...
## Fix coming in the future.

divider <- 4

for (taxa in taxaLevels){
    for (tpIter in tps){
        inFileName <- paste0("LogNormwithMetadata_R1_L_", taxa, ".txt")
        myT <- read.csv(inFileName, header=TRUE, sep = "\t")
        numCols <- ncol(myT)

        ## Restrict to timepoint tpIter
        myTtp <- myT[myT$Timepoint == tpIter,]

        factorsForThistp <- unique(factor(myTtp$Sample.Location))
        myT <- myTtp
        ## Select the largest technical replicate
        ## This can be kept as is because of the fact that sequenceCount doesn't change right now; I need to fix this in the future though
        savemyT <- myT[FALSE,]
        for(loc in factorsForThistp)){
            myTloc <- myT[myTtp$Sample.Location == loc,]
            rep<-myTloc[which(myTloc$sequenceCount == max(myTloc$sequenceCount)),]
            ## print(c(tp, loc, rep))
            ## print(rep)
            ## The order here may be backwards
            savemyT <- rbind(savemyT, rep)
        }
        myT <- savemyT


        pdf( paste(taxa, filePrefix, "boxplots.pdf", sep=""))

        for( i in 2:(ncol(myT) - numMetadataCols))
            if( sum(myT[,i] != 0 ) > nrow(myT) / divider ){
                ## Taxon abundances
                bug <- myT[,i]
                ## Diversity Metrics
                ## MDS Axes

                ## Technical
                seqCount <- myT$sequenceCount
                sampNumber <- myT$SampleNumber
                sampID <- myT$Sample.ID
                sampName <- myT$Sample.Name
                ## concBatch <- myT$Batch
                ## concQueue <- myT$Queue
                concSamp <- myT$Sample..
                anotherName <- myT$Sample
                replicateNum <- myT$Rep
                sampFCID <- myT$FCID
                sampLane <- myT$Lane
                sampRef <- myT$SampleRef
                sampIndex <- myT$Index
                sampDesc <- myT$Description
                sampControl <- myT$Control
                sampRecipe <- myT$Recipe
                sampOperator <- myT$Operator
                sampProject <- myT$SampleProject
                techInfo <- cbind(seqCount, sampNumber, sampID, sampName, concSampe, anotherName, replicateNum, sampFCID, sampLane, sampRef, sampIndex, sampDesc, sampControl, sampRecipe, sampOperator, sampProject)

                ## Abiotic
                ambTemp <- myT$Ambient.Temp.Celsius
                storeTemp <- myT$Storage.Temp.Celsius
                sampTemp <- myT$Sample.Temp.Celsius
                conduct <- myT$Conductivity..pHmv.
                pH <- myT$pH
                humidity <- myT$Humidity
                abioticInfo <- cbind(ambTemp, storeTemp, sampTemp, conduct, pH, humidity)

                ## Time-related
                dateTime <- myT$Sample.Date.Time
                compStart <- myT$Composite.Start
                compEnd <- myT$Composite.End
                sampTimepoint <- myT$Timepoint
                timeInfo <- cbind(dateTime, compStart, compEnd, sampTimepoint)

                ## Space-related
                latitude <- myT$Lattitude
                longitude <- myT$Longitude
                sampLoc <- myT$Sample.Location
                anotherLoc <- myT$Location
                majorLoc <- unlist(lapply(lapply(lapply(as.character(myT$Sample.Location), function(x) strsplit(x, split=" ")),'[[',1),'[[',1))
                spaceInfo <- cbind(latitude, longitude, sampLoc, anotherLoc, majorLoc)

                ## Antibiotic Concentrations
                ertapenem <- myT$Ertapenem..ng.L.
                amoxicillin <- myT$Amoxicillin..ng.L.
                ciprofloxacin <- myT$Ciprofloxacin..ng.L.
                doxycycline <- myT$Doxycycline..ng.L.
                azithromycin <- myT$Azithromycin..ng.L.
                clindamycin <- myT$Clindamycin..ng.L.
                sulfamethoxazole <- myT$Sulfamethoxazole..ng.L.
                cephalexin <- myT$Cephalexin..ng.L.
                trimethoprim <- myT$Trimethoprim..ng.L.
                levofloxacin <- myT$Levofloxacin..ng.L.
                abcInfo <- cbind(ertapenem, amoxicillin, ciprofloxacin, doxycycline, azithromycin, clindamycin, sulfamethoxazole, cephalexin, trimethoprim, levofloxacin)
                ## rowMeans
                ## rowSums

                ## Some kind of vector for the variables of interest
                myFrame <- data.frame(bug, techInfo, abioticInfo, timeInfo, spaceInfo, abcInfo)
                if(mlm == TRUE){
                    fullModel <- gls( bug~ abcInfo[,1], method="REML",correlation=corCompSymm(form=~1|majorLoc), data = myFrame )
                    reducedModel <- gls( bug~ abcInfo[,1], method="REML", data = myFrame )
                    ## reducedModel <- lme( bug ~ treatment + batch, method="REML", data = myFrame)
                    fullModelLME <- lme(bug~ abcInfo[,1], method="REML", random = ~1|majorLoc, data = myFrame)
                    ## fullModelLME <- lme(bug~ treatment + batch, method="REML", random= list(group = ~1, cage = ~1), data = myFrame)

                }
                else{
                    fullModelLME <- lm(bug~ abcInfo[,1], x=TRUE)
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
                    allNames<-c(allNames, "majorLocation", "iccmajorLocation")
                    allpvals[[index]]<-c(allpvals[[index]], anova(fullModelLME, reducedModel)$"p-value"[2], coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]])
                    pValuesCage[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
                    ## This is different as it is not corrected for multiple hypothesis testing.
                    iccCage[index]<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]
                    ## iccCage[index] <- as.numeric(VarCorr(fullModelLME)[1,1])/ sum(as.numeric(VarCorr(fullModelLME)[,1])
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
        keepVector <- grep("adj", names(dFrame))
        ## Try and dynamically generate the name of the plot here...
        write.table(dFrame, file=paste("pValues", filePrefix, taxa, ".txt",sep=""), sep="\t",row.names=FALSE)

        ## Get the sig table working
        ## It would also be nice to have a sig-picker for the plots and a combination of all sig results into one plot.
        index <- 1
        ## Iterate over rownames in dFrame instead
        for( i in rownames(dFrame)){
            bug <- myT[,i]
            animal <- paste0(myT$Sex,"_", myT$Animal.Nu)
            sampType <- myT$Source
            sex <- myT$Sex
            cage <- myT$Cage
            ## cage <- myT$Pen.Location
            group <- myT$Exp.or.Ctrl
            date <- myT$Experiment..Sample.info
            names[index] = names(myT)[i]
            ## Some kind of vector for the variables of interest
            myFrame <- data.frame(bug, group, sampType, sex, animal, date, cage)

            if(mlm == TRUE){
                graphMain = c(dFrame[i,1], paste(c(rep("p-adj_",length(allNames)-1),""), allNames, " = ", format(dFrame[i, keepVector], digits=3), sep=""))
            }
            else {
                graphMain = c(dFrame[i,1], paste(c(rep("p-adj_",length(allNames)),""), allNames, " = ", format(dFrame[i, keepVector], digits=3), sep=""))
            }

            ## mlm is boolean 1/0 value. Basically, don't graph on icc
            ## need a correction for interaction terms as well
            if(mlm == TRUE)
                numPlots <- length(allNames) - mlm
            else
                numPlots <- length(allNames) + 1

            par(mfrow=c(numPlots, 1),
                oma = c(1,1,0,0) + 0.1,
                mar = c(1,4,3,0) + 0.1)
            testInteract <- grep(":", allNames)
            graphMain <- graphMain[-length(graphMain)]
            for(j in 1:numPlots){
                if(j == 1){
                    plot( bug ~ factor(get(allNames[j])), ylab = "lognorm Abundance",
                         main = c(i, toString(graphMain[-1]) ))
                    points(factor(get(allNames[j])), bug)
                }
                else{
                    if(length(testInteract) > 0 && j == testInteract && j != numPlots){
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
                        if(j < numPlots){
                            plot( bug ~ factor(get(allNames[j])), ylab = names[i])
                            points(factor(get(allNames[j])), bug) ##,
                        }
                        else{
                            plot( bug ~ factor(cage), ylab = names[i], xaxt = "n")
                            points(factor(cage), bug) ##,
                            axis(1, at=1:8, labels=c("AME", "BME", "CMC", "DMC", "EFE", "FFE", "GFC", "HFC"))
                        }
                    }
                }
            }
            index <- index + 1
        }

        dev.off()
        if(dim(dFrame)[1] > maxCount){
            maxCount <- dim(dFrame)[1]
            print(maxCount)
        }
        ## This should be on unadjusted p-values
        pdf( paste(taxa, filePrefix, "p-value_histogram.pdf", sep=""))
        nonadjust <- colnames(dFrame)[-keepVector]
        nonadjust <- nonadjust[-1]
        for(phist in nonadjust){
            graphMainhist = c(taxa, filePrefix, phist)
            hist(dFrame[,phist], main=graphMainhist, xlim=c(0,1), breaks=10)
        }
        dev.off()
        ## sigdFrame <-dFrame[which(dFrame[,keepVector] < 0.05,arr.ind=TRUE),]
        ## write.table(sigdFrame, file=paste("SIGpValues", filePrefix, taxa, ".txt",sep=""), sep="\t",row.names=FALSE)

   }
}
