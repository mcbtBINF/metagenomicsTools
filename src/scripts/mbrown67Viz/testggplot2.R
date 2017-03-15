## Purpose to explore ggplot2 in producing nice plots towards getting ready for better publication quality figures and better figures for my dissertation proposal and the ultimate dissertation

rm(list=ls())

library("ggplot2")
library("ggrepel")
## For situations where you want the graphical labels to repel each other.
## library("ggvis")
## For interactive graphics
setwd("~/Documents/Tutorials/RTutorials/ggplot2/Rgraphics/")

housing <- read.csv("dataSets/landdata-states.csv")

housing$Year <- as.numeric(substr(housing$Date, 1, 4))
housing$Qrtr <- as.numeric(substr(housing$Date, 5, 5))
housing$Date <- housing$Year + housing$Qrtr/4

## Showing basic versus ggplot2

hist(housing$Home.Value)

ggplot(housing, aes(x = Home.Value)) +
    geom_histogram()

##So it is name the data frame, choose the variable of concern, choose the layout

plot(Home.Value ~ Date,
     data=subset(housing, State == "MA"))
points(Home.Value ~ Date, col="red",
       data=subset(housing, State == "TX"))
legend(19750, 400000,
       c("MA", "TX"), title="State",
       col=c("black", "red"),
       pch=c(1, 1))
## It doesn't display the legend for some reason.


ggplot(subset(housing, State %in% c("MA", "TX")),
       aes(x=Date,
           y=Home.Value,
           color=State))+
    geom_point()
## It knows how to automatically do the legend here.
## Select the rows, select the variables of interest, decide on coloration.


hp2001Q1 <- subset(housing, Date == 2001.25)
hp2001Q1$pred.SC <- predict(lm(Structure.Cost ~ log(Land.Value), data = hp2001Q1))

p1 <- ggplot(hp2001Q1, aes(x = log(Land.Value), y = Structure.Cost))

p1 + geom_point(aes(color = Home.Value)) +
    geom_line(aes(y = pred.SC))

p1 +
  geom_point() +
  geom_text_repel(aes(label=State), size = 3)

p1 +
  geom_point(size = 2,# incorrect! 2 is not a variable
             color="red") # this is fine -- all points red

p1 +
  geom_point(aes(color=Home.Value, shape = region))

dat <- read.csv("dataSets/EconomistData.csv")
head(dat)

ggplot(dat, aes(x = CPI, y = HDI, size = HDI.Rank)) + geom_point()

p2 <- ggplot(dat, aes(x = CPI, y = HDI))

p2 +
    geom_boxplot(aes(color = Region, size = HDI.Rank))

## Some ggvis examples from a talk by Hadley at Stanford
## Almost none of the names here are currently working right
slider <- input_slider(0.1, 1, value = 0.5)
ggvis(nameofdataframe, props(~item1, ~item2)) +
    layer_line() +
    layer_smooth(span = slider, props(stroke := "red"))
