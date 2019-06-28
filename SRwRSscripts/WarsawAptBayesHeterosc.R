########## R script: WarsawAptBayesHeterosc ##########

# For doing hetroscedastic Bayesian penalized spline
# fitting to the Warsaw apartments running example dataset.

# Last changed: 30 MAY 2016

# Load required packages:

library(HRW) ; library(rstan)

# Obtain data:

xOrig <- WarsawApts$construction.date
yOrig <- WarsawApts$areaPerMzloty

n <- length(yOrig)

mean.x <- mean(xOrig)  ; sd.x <- sd(xOrig)
mean.y <- mean(yOrig)  ; sd.y <- sd(yOrig)

x <- (xOrig - mean.x)/sd.x
y <- (yOrig - mean.y)/sd.y

# Obtain linear and spline basis design matrices (X and Z):

X <- cbind(1,x)
aOrig <- min(xOrig) ; bOrig <- max(xOrig)
a <- (aOrig - mean.x)/sd.x  ;  b <- (bOrig - mean.x)/sd.x

numIntKnotsU <- 25
intKnotsU <-  quantile(unique(x),seq(0,1,length=numIntKnotsU+2)
                      [-c(1,numIntKnotsU+2)])
Zu <- ZOSull(x,intKnots=intKnotsU,range.x=c(a,b))
ncZu <- ncol(Zu)

numIntKnotsV <- 15
intKnotsV <-  quantile(unique(x),seq(0,1,length=numIntKnotsV+2)
                      [-c(1,numIntKnotsV+2)])
Zv <- ZOSull(x,intKnots=intKnotsV,range.x=c(a,b))
ncZv <- ncol(Zv)

# Set MCMC sample sizes:

nWarm <- 100
nKept <- 200

# Specify model in Stan:

penSplHeteroscModel <-
'data
{
   int<lower=1> n; vector[n] y;
   int<lower=1> ncZu; int<lower=1> ncZv;
   matrix[n,2] X;
   matrix[n,ncZu] Zu; matrix[n,ncZv] Zv;
   real<lower=0> sigmaBeta; real<lower=0> sigmaGamma;
   real<lower=0> Au; real<lower=0> Av;
}
parameters
{
   vector[2] beta; vector[2] gamma;
   vector[ncZu] u; vector[ncZv] v;
   real<lower=0> sigmaU; real<lower=0> sigmaV;
}
model
{
   y ~ normal(X*beta + Zu*u,exp((X*gamma + Zv*v)/2));
   u ~ normal(0,sigmaU);  v ~ normal(0,sigmaV);
   beta ~ normal(0,sigmaBeta); gamma ~ normal(0,sigmaGamma);
   sigmaU ~ cauchy(0,Au); sigmaV ~ cauchy(0,Av);
}'

# Store data in a list in format required by stan():

allData <- list(n=length(x),ncZu=ncZu,ncZv=ncZv,y=y,X=X,
                Zu=Zu,Zv=Zv,sigmaBeta=1e5,sigmaGamma=1e5,
                Au=1e5,Av=1e5)

# Compile code for model:

stanCompilObj <- stan(model_code=penSplHeteroscModel,data=allData,
                      iter=1,chains=1)

# Obtain MCMC samples for each parameter using Stan:

stanObj <- stan(model_code=penSplHeteroscModel,data=allData,
                warmup=nWarm,iter=(nWarm+nKept),
                chains=1,refresh=30,fit=stanCompilObj)

# Save Stan output:

betaMCMC <- NULL
gammaMCMC <- NULL
for (j in 1:2)
{
   charVar <- paste("beta[",as.character(j),"]",sep="") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar,permuted=FALSE))
   charVar <- paste("gamma[",as.character(j),"]",sep="") 
   gammaMCMC <- rbind(gammaMCMC,extract(stanObj,charVar,permuted=FALSE))
}   
uMCMC <- NULL
for (k in 1:ncZu)
{
   charVar <- paste("u[",as.character(k),"]",sep="") 
   uMCMC <- rbind(uMCMC,extract(stanObj,charVar,permuted=FALSE))
}
vMCMC <- NULL
for (k in 1:ncZv)
{
   charVar <- paste("v[",as.character(k),"]",sep="")
   vMCMC <- rbind(vMCMC,extract(stanObj,charVar,permuted=FALSE))
}

# Obtain MCMC samples of regression curves over a fine grid:

ng <- 101
xg <- seq(a,b,length=ng)
Xg <- cbind(rep(1,ng),xg)
Zug <- ZOSull(xg,intKnots=intKnotsU,range.x=c(a,b))
fhatMCMC <- Xg%*%betaMCMC + Zug%*%uMCMC
xgOrig <- xg*sd.x + mean.x
fhatMCMCorig <- fhatMCMC*sd.y + mean.y
fhatg <- apply(fhatMCMCorig,1,mean)
credLowerf <- apply(fhatMCMCorig,1,quantile,0.025)
credUpperf <- apply(fhatMCMCorig,1,quantile,0.975)

# Display the fitted mean and standard deviation functions:

par(mfrow=c(1,2),mai=c(1,1.1,0.8,0.2))

plot(xOrig,yOrig,type="n",xlab="construction date (year)",
     ylab="area (square meters) per million zloty",
     main="mean function estimate",col.main="navy",
     bty="l",xlim=range(xgOrig),ylim=range(c(yOrig,credLowerf,credUpperf)))
polygon(c(xgOrig,rev(xgOrig)),c(credLowerf,rev(credUpperf)),
        col="palegreen",border=FALSE)
lines(xgOrig,fhatg,col="darkgreen",lwd=2)
points(xOrig,yOrig,col="dodgerblue")

Zvg <- ZOSull(xg,intKnots=intKnotsV,range.x=c(a,b)) 
loggMCMC <- Xg%*%gammaMCMC + Zvg%*%vMCMC
sqrtgMCMC <- exp(loggMCMC/2)
sqrtgMCMCorig <- sqrtgMCMC*sd.y
sqrtghatMCMCg <- apply(sqrtgMCMCorig,1,mean)
credLowMCMCrtg <- apply(sqrtgMCMCorig,1,quantile,0.025)
credUppMCMCrtg <- apply(sqrtgMCMCorig,1,quantile,0.975)

plot(0,0,type="n",xlim=range(xgOrig),
     ylim=range(c(credLowMCMCrtg,credUppMCMCrtg)),
     bty="l",main="standard deviation function estimate",
     col.main="navy",
     xlab="construction date (year)",
     ylab="area (square meters) per million zloty")
polygon(c(xgOrig,rev(xgOrig)),c(credLowMCMCrtg,
        rev(credUppMCMCrtg)),col="palegreen",border=FALSE)
lines(xgOrig,sqrtghatMCMCg,col="darkgreen",lwd=2)

########## End of WarsawAptBayesHeterosc ##########
