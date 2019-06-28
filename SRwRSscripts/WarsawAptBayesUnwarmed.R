########## R script: WarsawAptBayesUnwarmed ##########

# For obtaining a Bayesian penalized spline
# fit to the Warsaw apartment data using Stan
# via the package rstan. Approximate Bayesian
# inference is achieved via Markov chain
# Monte Carlo (MCMC).

# Last changed: 18 MAY 2016

# Load required packages:

library(HRW)  ;  library(rstan) 

# Load Warsaw apartment data:

data(WarsawApts)

# Standardize both predictor and response variable
# and set hyperparameter values:

xOrig <- WarsawApts$construction.date
yOrig <- WarsawApts$areaPerMzloty
mean.x <- mean(xOrig)  ; sd.x <- sd(xOrig)
mean.y <- mean(yOrig)  ; sd.y <- sd(yOrig)
x <- (xOrig - mean.x)/sd.x
y <- (yOrig - mean.y)/sd.y
sigmaBeta <- 1e5 ; Au <- 1e5 ; Aeps <- 1e5

# Obtain linear and spline basis design matrices (X and Z):

X <- cbind(1,x)
aOrig <- min(xOrig) ; bOrig <- max(xOrig)
a <- (aOrig - mean.x)/sd.x  ;  b <- (bOrig - mean.x)/sd.x
numIntKnots <- 25
intKnots <-  quantile(unique(x),seq(0,1,length=numIntKnots+2)
                      [-c(1,numIntKnots+2)])
Z <- ZOSull(x,intKnots=intKnots,range.x=c(a,b))
ncZ <- ncol(Z)

# Specify model in Stan:

npRegModel <- 
   "data
   {
      int<lower=1> n;         int<lower=1> ncZ;
      vector[n] y;            matrix[n,2] X;
      matrix[n,ncZ] Z;        real<lower=0> sigmaBeta;
      real<lower=0> Au;       real<lower=0> Aeps;
   }
   parameters 
   {
      vector[2] beta;          vector[ncZ] u;
      real<lower=0> sigmaeps;  real<lower=0> sigmau;
   }
   model 
   {
      y ~ normal(X*beta + Z*u,sigmaeps);
      u ~ normal(0,sigmau); beta ~ normal(0,sigmaBeta); 
      sigmaeps ~ cauchy(0,Aeps); sigmau ~ cauchy(0,Au);
   }"

# Store data in a list in format required by Stan:

allData <- list(n=length(x),ncZ=ncZ,y=y,X=X,Z=Z,
                sigmaBeta=sigmaBeta,Au=Au,Aeps=Aeps)

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code=npRegModel,data=allData,
                         iter=1,chains=1)

# Set MCMC sample size parameters:

nWarm <-    5       # Length of warm-up.
nKept <- 1000       # Size of the kept sample.

# Obtain MCMC samples for each parameter using Stan:

initFun <- function()
   return(list(sigmau=1,sigmaeps=0.7,beta=rep(0,2),u=rep(0,ncZ)))

stanObj <-  stan(model_code=npRegModel,data=allData,warmup=nWarm,
                 iter=(nWarm+nKept),chains=1,refresh=100,
                 fit=stanCompilObj,init=initFun,seed=3)

# Extract relevant MCMC samples:

betaMCMC <- NULL
for (j in 1:2)
{
   charVar <- paste("beta[",as.character(j),"]",sep="") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar,permuted=FALSE))
}
uMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("u[",as.character(k),"]",sep="") 
   uMCMC <- rbind(uMCMC,extract(stanObj,charVar,permuted=FALSE))
}
sigmaepsMCMC <- as.vector(extract(stanObj,"sigmaeps",permuted=FALSE))
sigmauMCMC <- as.vector(extract(stanObj,"sigmau",permuted=FALSE))

# Obtain MCMC samples of regression curves over a fine grid:

ng <- 101
xgOrig <- seq(aOrig,bOrig,length=ng)
xg <- (xgOrig - mean.x)/sd.x
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots=intKnots,range.x=c(a,b))
fhatMCMC <- Xg%*%betaMCMC + Zg%*%uMCMC

# Convert fhatMCMC matrix to original scale:

fhatMCMCOrig <- fhatMCMC*sd.y + mean.y
fhatgOrig <- apply(fhatMCMCOrig,1,mean)
credLower <- apply(fhatMCMCOrig,1,quantile,0.025)
credUpper <- apply(fhatMCMCOrig,1,quantile,0.975)

# Display the fit:

par(mai=c(1,1.1,0.1,0.1))
cex.labVal <- 2   ;   cex.axisVal <- 1.5
plot(xOrig,yOrig,type="n",xlab="construction date (year)",
     ylab="area (square meters) per million zloty",
     bty="l",xlim=range(xgOrig),ylim=range(c(credLower,credUpper,yOrig)),
     cex.lab=cex.labVal,cex.axis=cex.axisVal)
polygon(c(xgOrig,rev(xgOrig)),c(credLower,rev(credUpper)),
        col="palegreen",border=FALSE)
lines(xgOrig,fhatgOrig,col="darkgreen",lwd=2)
points(xOrig,yOrig,col="dodgerblue")
abline(v=quantile(xOrig,0.25),lty=2,col="darkorange")
abline(v=quantile(xOrig,0.50),lty=2,col="darkorange")
abline(v=quantile(xOrig,0.75),lty=2,col="darkorange")

cat("Hit Enter to continue.\n")
ans <- readline()

# Obtain samples from the posterior distribution of the
# effective degrees of freedom:

X <- cbind(rep(1,length(x)),x)
Z <- ZOSull(x,intKnots=intKnots,range.x=c(a,b))
CTC <- crossprod(cbind(X,Z)) ; Dmat <- diag(c(0,0,rep(1,ncol(Z))))
lambdaMCMC <- (sigmaepsMCMC/sigmauMCMC)^2
EDFMCMC <- rep(NA,length(lambdaMCMC))
for (i in 1:length(lambdaMCMC))
   EDFMCMC[i] <- sum(diag(solve(CTC+lambdaMCMC[i]*Dmat,CTC)))

# Convert error standard deviation MCMC sample to the orginal units:

sigmaepsOrigMCMC <- sd.y*sigmaepsMCMC

# Do some summaries and diagnostic checking of the MCMC:

indQ1 <- length(xgOrig[xgOrig<quantile(xOrig,0.25)])
indQ2 <- length(xgOrig[xgOrig<quantile(xOrig,0.50)])
indQ3 <- length(xgOrig[xgOrig<quantile(xOrig,0.75)])
fhatOrigQ1MCMC <- fhatMCMCOrig[indQ1,]
fhatOrigQ2MCMC <- fhatMCMCOrig[indQ2,]
fhatOrigQ3MCMC <- fhatMCMCOrig[indQ3,]
MCMClist <- list(cbind(EDFMCMC,sigmaepsOrigMCMC,
                 fhatOrigQ1MCMC,fhatOrigQ2MCMC,fhatOrigQ3MCMC))
parNamesVal <- list(c("effective","degrees","of freedom"),
                    c("error","standard","deviation"),
                    c("reg'n func. est.","at 1st quartile","of construc. date"),
                    c("reg'n func. est.","at 2nd quartile","of construc. date"),
                    c("reg'n func. est.","at 3rd quartile","of construc. date"))

summMCMC(MCMClist,parNames=parNamesVal)

cat("Hit Enter to continue.\n")
ans <- readline()

# Obtain chain summaries via the monitor() function:

myMCMCarray <- array(0,dim=c(length(sigmaepsMCMC),1,5))
myMCMCarray[,1,1] <- EDFMCMC
myMCMCarray[,1,2] <- sigmaepsOrigMCMC
myMCMCarray[,1,3] <- fhatOrigQ1MCMC
myMCMCarray[,1,4] <- fhatOrigQ2MCMC
myMCMCarray[,1,5] <- fhatOrigQ3MCMC
monitorAnswer <- monitor(myMCMCarray,warmup=0,print=FALSE)
dimnames(monitorAnswer)[[1]] <- c("EDF","err. st. dev.","f(Q_1)","f(Q_2)","f(Q_3)")
print(signif(monitorAnswer,4))

########## End of WarsawAptBayesUnwarmed ##########

