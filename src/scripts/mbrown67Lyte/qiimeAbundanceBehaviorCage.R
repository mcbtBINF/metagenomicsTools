rm(list=ls())

library("Kendall")
library("pscl")
library("lmtest")
library("nlme")
library("psych")
library("ggplot2")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/MarkExperiment/ArgonneSequencing/JanuaryResequencing/rawData/rg_results/")

inFileName <- "qiimeOTULogNormwithMetadata.txt"
myT <- read.table(inFileName,header=TRUE,sep="\t")

otuMapping <- read.table("LyteSharon_r01_cr_MAPPING.txt", sep="\t", header=TRUE)
otuMapping[,1] <- paste0("X",otuMapping[,1])
rownames(otuMapping)<-otuMapping[,1]

##tissues <- unique(myT$Source)
tissues <- c("feces")
for( t in tissues ) {
    subT <- myT[ myT$Source == t, ]
    ## Should ultimately be 20:58
    for ( behavior in 20:58){

        index <- 1
        names <- vector()
        pValuesGroup <- vector()
        pValuesSex <- vector()

        pValuesBehaveFromMixed <- vector()
        pValuesSexFromMixed <- vector()
        pValuesSexBehaveFromMixedInteraction <- vector()
        pValuesCage <- vector()

        pdf( paste("qiimeOTUanalysisBehavior_", t, colnames(myT)[behavior], ".pdf", sep=""))

        ## Ignores the metadata
        for ( j in c(61:3560)) {
            ## par(mfrow=c(2,2), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
            bug <- subT[,j]

            if( sum(bug != 0 ) > nrow(subT) / 4 ) {
                ## pValuesBehave[index] <- wilcox.test( subT[subT$Group=="Experimental", j],
                ##                                   subT[subT$Group=="Control", j])$p.value
                names[index] <- names(subT)[j]
                behave <- subT[,behavior]
                ## group <- factor(subT$Group)
                sex <- factor(subT$Sex)

                cage <-  factor( paste( subT$Housing, subT$Sex, sep=""), c("#5,6,7,8 in same cagefemale", "#13,14,15,16 in same cagefemale", "#1,2,3,4 in same cagefemale", "#9,10,11,12 in same cagefemale",
                                                                           "#5,6,7,8 in same cagemale", "#13,14,15,16 in same cagemale", "#1,2,3,4 in same cagemale", "#9,10,11,12 in same cagemale"))

                myFrame <- data.frame(bug, sex, behave, cage)

                ## pValuesSex[index] <- wilcox.test( subT[subT$Sex=="male", j],
                ##                                 subT[subT$Sex=="female", j])$p.value

                fullModel <-
                    gls( bug~  behave * sex, method="REML",correlation=corCompSymm(form=~1|factor(cage)),data = myFrame )

                mixedAnova <- anova(fullModel)
                pValuesBehaveFromMixed[index] <- mixedAnova$"p-value"[2]
                pValuesSexFromMixed[index] <- mixedAnova$"p-value"[3]
                pValuesSexBehaveFromMixedInteraction[index] <- mixedAnova$"p-value"[4]

                pValuesCage[index] <- anova(lm( bug ~ cage ))$"Pr(>F)"[1]

                index <- index + 1
            }
        }
        index <- 1

        namesMapping <- otuMapping[names,2]

        dFrame <- data.frame(names,namesMapping, colnames(myT)[behavior],
                             ## pValuesBehave,pValuesSex,
                             pValuesBehaveFromMixed,
                             pValuesSexFromMixed,pValuesSexBehaveFromMixedInteraction, pValuesCage)
        dFrame <- dFrame [order(dFrame$pValuesBehaveFromMixed),]
        ## dFrame$adjustedPValuesBehave <- p.adjust( dFrame$pValuesBehave, method = "BH" )
        ## dFrame$adjustedPValuesSex<- p.adjust( dFrame$pValuesSex, method = "BH" )

        dFrame$adjustedpValuesBehaveFromMixed<- p.adjust( dFrame$pValuesBehaveFromMixed, method = "BH" )
        dFrame$adjustedpValuesSexFromMixed<- p.adjust( dFrame$pValuesSexFromMixed, method = "BH" )
        dFrame$adjustedpValuesSexBehaveFromMixedInteraction<- p.adjust( dFrame$pValuesSexBehaveFromMixedInteraction, method = "BH" )
        dFrame$pValuesCageAdjust<- p.adjust( dFrame$pValuesCage, method = "BH" )
        rownames(dFrame) <- dFrame$names
        ifelse(exists("totalFrame"), totalFrame <- rbind(totalFrame, dFrame), totalFrame <- dFrame)

        write.table(dFrame, file=paste("qiimeOTUanalysisBehavior", t, colnames(myT)[behavior], ".txt",sep=""), sep="\t",row.names=FALSE)

        ## for ( j in c(61:3560)) {
        for (nameIter in rownames(dFrame)) {
            par(mfrow=c(2,2), mgp=c(3, 1.25, 0), oma = c(0, 0, 2, 0))
            bug <- subT[,nameIter]

            ## if( sum(bug != 0 ) > nrow(subT) / 4 ) {
            ## names[index] <- names(subT)[j]

            behave <- as.numeric(levels(subT[,behavior])[subT[,behavior]]
            sex <- factor(subT$Sex)

            cage <-  factor( paste( subT$Housing, subT$Sex, sep=""), c("#5,6,7,8 in same cagefemale", "#13,14,15,16 in same cagefemale", "#1,2,3,4 in same cagefemale", "#9,10,11,12 in same cagefemale",
                                                                       "#5,6,7,8 in same cagemale", "#13,14,15,16 in same cagemale", "#1,2,3,4 in same cagemale", "#9,10,11,12 in same cagemale"))

            myFrame <- data.frame(bug, sex, behave, cage)

            plot( behave, bug , ylab="Log normalized abundance",
                    main = paste("behavior p-value", format(dFrame$adjustedpValuesBehaveFromMixed[index],digits=3) ), col=ifelse(sex == "male", "blue", "pink"), pch=16)

            ##stripchart(bug ~ behave ,
            ##           data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

            boxplot( bug ~ sex,
                    main = paste("sex p-value", format(dFrame$adjustedpValuesSexFromMixed[index],digits=3) ), ylab="Log normalized abundance" )

            stripchart(bug ~ sex,
                       data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

        ##             boxplot( bug ~ behave*sex, ylab="Log normalized abundance",
        ##         main=paste("sex*behave interaction p-value", format(dFrame$adjustedpValuesSexBehaveFromMixedInteraction[index],digits=3)),
        ##         xaxt='n')
        ## axis(1, at=c(1, 2, 3, 4), labels = c("Control\nFemale", "Exp.\nFemale", "Control\nMale", "Exp.\nMale"), cex.axis=0.9)

        ## stripchart(bug ~ behave*sex,
        ##            data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index] )

            ## The below code is broken currently :I
            model1 = lm(bug ~ behave*sex, data = myFrame)
            myFrame$predicted=predict(model1)
            apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'))
            p = ggplot(myFrame, aes(x = behave, y=bug, shape=sex))+
                geom_point()+
                scale_shape_manual(values=c(1,16), name='Sex', labels=c('Female','Male'))+
geom_line(aes(x = behave, y = predicted, linetype=sex)) +
scale_linetype_discrete(name='Sex', labels=c('Female','Male'))+
labs(x = 'Behavioral Score', y = 'Log Normalized Abundance')+
apatheme
p
                    ## plot(1, type="n", axes=F, xlab="", ylab="")
            boxplot( bug ~ cage ,las=2, ylab="Log normalized abundance",
                    main=paste("cage effect p-value", format(dFrame$pValuesCageAdjust[index],digits=3)), xaxt='n')
            ##            axis(1, at=c(1, 2, 3, 4, 5, 6, 7, 8),
            ##                 labels=c("FC1", "FC2", "FX1", "FX2", "MC1", "MC2", "MX1", "MX2"),
            ##                 cex.axis=0.65)
            axis(1, at=c(1, 3, 5, 7),
                 labels=c("FC1", "FX1", "MC1", "MX1"), cex.axis=0.8)
            axis(1, at=c(2, 4, 6, 8),
                 labels=c("FC2", "FX2", "MC2", "MX2"), cex.axis=0.8)

            stripchart(bug ~ cage,
                       data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])

            mtext(paste0(otuMapping[nameIter,2],"\n",nameIter), outer=TRUE, cex = 0.7)

            index = index + 1
        }

        hist(pValuesBehaveFromMixed,breaks=20, main="")
        hist(pValuesSexFromMixed,breaks=20, main="")
        hist(pValuesSexBehaveFromMixedInteraction,breaks=20, main="")
        hist(pValuesCage,breaks=20, main="")

        dev.off()

    }
}

        write.table(totalFrame, file=paste("qiimeOTUanalysisBehaviortotalFrame", ".txt",sep=""), sep="\t",row.names=FALSE)
