rm(list=ls())

library(pwr)
powerNoMHC<-pwr.r.test(n=91, r=NULL, sig.level=0.05, power=0.80)
powerWithMHC<-pwr.r.test(n=91, r=NULL, sig.level=0.05/4522, power=0.80)
## To convert to Kendall's tau, divide the number of subjects by 1.1 according to David Clark-Carter's Quant Psych Research book.
## 91/1.1 ~ 82.7, so 82 rounded down
powerNoMHC.Kendallequiv<-pwr.r.test(n=82, r=NULL, sig.level=0.05, power=0.80)
powerWithMHC.Kendallequiv<-pwr.r.test(n=82, r=NULL, sig.level=0.05/4522, power=0.80)


powerGoal <- 0.80

## this function from...
## http://stackoverflow.com/questions/19096983/when-simulating-multivariate-data-for-regression-how-can-i-set-the-r-squared-e
## I don't know how to properly introduce variance to epsilon because I am not postulating a form of the relationship, rather I am working with a correlation
simulate <- function(n.obs = 92, R.sq) {#beta=c(5, 3, -2), R.sq) {
    ## stopifnot(length(beta) == 3)
    df <- data.frame(x1=rnorm(n.obs))

    ## I don't fully understand the purpose of this epsilon
    ## var.epsilon <- (beta[2]^2 ) * (1 - R.sq) / R.sq
    var.epsilon <- 1 / R.sq
    stopifnot(var.epsilon > 0)
    df$epsilon <- rnorm(n.obs, sd=sqrt(var.epsilon))
    df$y <- with(df, epsilon) ## beta[1] + beta[2]*x1 + epsilon)
    return(df)
}

get.R.sq <- function(desired) {
    ## model <- lm(y ~ x1 , data=simulate(R.sq=desired))
    ## return(summary(model)$r.squared)
    df <- simulate(R.sq = desired)
    y <- df$y
    x1 <- df$x1
    model <- cor(y, x1, method="kendall")
    return(model)
}

# confirm that we are acheiving the desired r-squareds..
df <- data.frame(desired.R.sq=seq(from=0.05, to=0.95, by=0.05))
df$actual.R.sq <- sapply(df$desired.R.sq, FUN=get.R.sq)
plot(df)
abline(a=0, b=1, col="red", lty=2)

Ztransform <- function(rVal, nSubjects) {
    ## from Cohen 1988
    return((arctanh(rVal) + rVal)/(2*(nSubjects - 1)))
}

## perform 1,000 simulations in order to ask
## what % of the time you would observe a significant value
## when correcting for 359 tests at a 0.05 threshold
## In our case it is 4522 tests at a 0.05 threshold
numBelow <- 0
pVals <- vector()
for( i in 1:4522)
{
	aSim <- simulate( R.sq=0.60)

	## aLm <- lm( aSim$y ~ aSim$x );

	## aPValue = anova(aLm)$"Pr(>F)"[1];
        aPValue = Kendall(aSim$y, aSim$x)$sl[1]
        pVals[i] <- aPValue
	if( aPValue < 0.05/4522)
	{
		numBelow = numBelow + 1
	}
}

# print out the power
numBelow
adjP <- p.adjust(pVals, method="BH")
