rm(list=ls())
library("Kendall")

## this function from...
## http://stackoverflow.com/questions/19096983/when-simulating-multivariate-data-for-regression-how-can-i-set-the-r-squared-e
simulate <- function(n.obs=23, beta=c(5, 3, -2), R.sq) {
    stopifnot(length(beta) == 3)
    df <- data.frame(x1=rnorm(n.obs))
    var.epsilon <- (beta[2]^2 ) * (1 - R.sq) / R.sq
    stopifnot(var.epsilon > 0)
    df$epsilon <- rnorm(n.obs, sd=sqrt(var.epsilon))
    df$y <- with(df, beta[1] + beta[2]*x1 + epsilon)
    return(df)
}
get.R.sq <- function(desired) {
    model <- lm(y ~ x1 , data=simulate(R.sq=desired))
    return(summary(model)$r.squared)
}

## confirm that we are acheiving the desired r-squareds..
df <- data.frame(desired.R.sq=seq(from=0.05, to=0.95, by=0.05))
df$actual.R.sq <- sapply(df$desired.R.sq, FUN=get.R.sq)
plot(df)
abline(a=0, b=1, col="red", lty=2)

## taxa that pass the 1/4 presence/absence cutoff (see numerical values in manuscript) multiplied by the number of psych measures (17).
## This is what will be different from the first time.
countTests <- c(6, 9, 12, 25, 50) ## For T1 and T2
## countTests <- c(7, 9, 13, 26, 51) ## For T1
## countTests <- c(7, 10, 14, 26, 48) ## For T2


## Some base level for setting R.sq
testR.sq <- 0.10
desiredPower <- 0.8

for( k in 1:1000){
    estPower <- vector()
    levelIndex <- 1
    pValuesList <- list()
    for( j in countTests){
        pValuesList[[levelIndex]] <- vector()
        for( i in 1:j)
            {
                ## This R.sq value will change 0.08, 0.2, 0.1, 0.11
                aSim <- simulate( R.sq = testR.sq)
                pValuesList[[levelIndex]][i]  = Kendall(aSim$x, aSim$y)$sl[1];
            }
        correctedPValues <- p.adjust( pValuesList[[levelIndex]] , method = "BH" )
        estPower[[levelIndex]] <- sum( correctedPValues < 0.10 )/ length( correctedPValues)
        levelIndex <- levelIndex + 1
    }
    if(FALSE %in% (desiredPower < estPower))
       {
           testR.sq <- testR.sq + 0.01
           k <- 1
       }
   }

print(testR.sq)
