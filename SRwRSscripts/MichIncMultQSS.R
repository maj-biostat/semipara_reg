########## R script: MichIncMultQSS ##########

# For doing a multiple quantile smoothing spline
# fit to the Michigan Income Study example.

# Last changed: 18 MAY 2016

# Load required packages:

library(Ecdat) ; library(quantreg)

# Read in data:

data(Workinghours)
x <- Workinghours$age 
y <- Workinghours$income/10

# Obtain and plot multiple quantile regression fits:

ng <- 1001
xg <- seq(min(x),max(x),length=ng)
ylabString <- "other household income ('000 $US)"
mainStrings <- c("full data view","zoomed view")
lwdVal <- 2

quantVec <- c(0.01,0.05,0.25,0.5,0.75,0.95,0.99)
colVec <- c("darkmagenta","navy","blue","darkgreen",
             "gold","darkorange","red")
lambdaVal <- 3.5

fitQSS <- vector("list",length(quantVec))
fitQSSg <- vector("list",length(quantVec))
for (iq in 1:length(quantVec))
{
   fitQSS[[iq]] <- rqss(y~qss(x,lambda=lambdaVal),tau=quantVec[iq])
   xgDF <- data.frame(x=xg)
   fitQSSg[[iq]] <- predict(fitQSS[[iq]],xgDF)
}

par(mfrow=c(1,2),mai=c(0.9,0.9,0.54,0.05)) 

# First do plots with actual data:

plot(x,y,cex=0.4,bty="l",xlab="wife's age in years",
     ylab=ylabString,
     cex.lab=1.5,col="dodgerblue",main=mainStrings[1],
     cex.main=1.5,cex.axis=1.5,col.main="navy")

legend("topleft",legend=c("99% quantile","95% quantile","75% quantile",
                          "50% quantile","25% quantile","5% quantile",
                          "1% quantile"),lwd=rep(lwdVal,7),col=rev(colVec),
                          cex=1.5)

for (iq in 1:length(quantVec))
   lines(xg,fitQSSg[[iq]],col=colVec[iq],lwd=lwdVal)

plot(x,y,cex=0.4,bty="l",xlab="wife's age in years",
     ylab=ylabString,ylim=c(0,200),cex.lab=1.5,col="dodgerblue",
     main=mainStrings[2],cex.main=1.5,cex.axis=1.5,col.main="navy")

for (iq in 1:length(quantVec))
   lines(xg,fitQSSg[[iq]],col=colVec[iq],lwd=lwdVal)

############ End of MichIncMultQSS ############
