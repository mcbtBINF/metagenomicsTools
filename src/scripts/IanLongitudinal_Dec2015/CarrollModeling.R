rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")
# Below is what I need for the Shannon diversity.
library("vegan")
library("Kendall")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" );

for(taxa in taxaLevels )
{
        inFileName <- paste( taxa, "LogNormalwithMetadata.txt", sep ="")
        myT <-read.table(inFileName,header=TRUE,sep="\t")
# This was what was likely causing my errors yesterday
                                        #	numCols <- ncol(myT)
#	myColClasses <- c("character", rep("numeric", numCols-1))
#	myT <-read.table(inFileName,header=TRUE,sep="\t",row.names=1,colClasses=myColClasses)

    patientB <- grep("B", rownames(myT$Sample.ID))
    patientC <- grep("C", rownames(myT$Sample.ID))

                                        #Staging variables
        names <- vector()
        pValuesEnergy<- vector()
        pValuesImputedBMI<- vector()
        pValuesImputedWeight <- vector()
        pValuesDay <- vector()
        pValuesParticipant <- vector()
        iccParticipant <- vector()

        index<-1
      	pdf( paste(taxa, "boxplots.pdf", sep=""))

        #Magic number here unclear, change later
        for(i in 2:300)
        {
            if( sum(myT[,i] != 0 ) > nrow(myT) / 4 )
                {
                    bug<-myT[,i]
                    sampleID<-myT$Sample.ID
                    energy<-myT$Energy.Intake..kcal.day.
                    iBMI<-myT$Imputed.BMI
                    iWeight<-myT$Imputed.Weight..kg.
                    day<-myT$Day
                    participant<-myT$Participant

                    myFrame <- data.frame(bug, sampleID, energy, iBMI, iWeight, day, participant)

                    fullModel <- lm(bug ~ day + participant, x=TRUE)
                    fullResiduals <- sum(residuals(fullModel)^2)

                    reducedModel <- lm(bug ~ day, x=TRUE)
                    reducedResiduals <- sum(residuals(reducedModel)^2)
                    #No clue as to what should be here with respects to n
                                        #                    F = (( reducedResiduals - fullResiduals) / (3) ) / (fullResiduals / 122)
                                        #anova(fullModel)
                    names[index] = names(myT)[i]

                    graphMain =  paste( names(myT)[i], "\n",
                            " pDay=", format( pValuesDay[index], digits=3),
                            " pParticipant= " , format( pValuesParticipant[index], digits=3),
                            " icc= " , format( iccParticipant[index], digits=3 ), sep="")
                    par(mfrow=c(3,1),
                            oma = c(1,1,0,0) + 0.1,
                            mar = c(1,4,2.5,0) + 0.1)

                    plot( bug ~ factor(day), ylab = names[index],main = graphMain )
                    points(factor(day), bug)

#                        plot( bug ~ factor(c( paste( myT$Day, myT$acuteOrChronic,sep=""))))
#                        points(factor(c( paste( myT$sex, myT$acuteOrChronic,sep=""))), bug)

                    plot( bug ~ factor(participant), ylab=names[index])
                    points(factor(participant), bug)

                    index=index+1

                }

        dFrame <- data.frame( names, pValuesDay, pValuesParticipant, iccParticipant)
	dFrame <- dFrame [order(pValuesDay),]

	dFrame$adjustedDay<- p.adjust( dFrame$pValuesDay, method = "BH" )

	dFrame$adjustedParticipant <- p.adjust( dFrame$pValuesParticipant, method = "BH" )

        write.table(dFrame, file=paste("pValuesFile", taxa, ".txt",sep=""),
                    sep="\t",row.names=FALSE)

        dev.off()
}
