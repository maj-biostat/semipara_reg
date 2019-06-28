############ R script: penSplGEV ##########

# For approximate Bayesian inference for the parameters
# in a generalized extreme value response penalized spline
# model. In this script simulated data are used.

# Last changed: 31 MAY 2019

# Load required packages:

library(HRW) ; library(rstan)

# Set sample size:

n <- 1000

# Set true function and parameter values:

fTrue <- function(x) 
   return(3*exp(-78*(x-0.38)^2)+exp(-200*(x-0.75)^2) +  2*x)

sigmaepsTrue <- 1.13 ; xiTrue <- -0.58
 
# Set hyperparemetrs:

sigmabeta <- 1e5   ;   Au <- 1e5
Aeps <- 1e5        ;   sigmaxi <- 1e5

# Generate data:

set.seed(1)
z <- -log(-log(runif(n)))

if (xiTrue!=0)
   epsilon <- sigmaepsTrue*(exp(xiTrue*z)-1)/xiTrue
if (xiTrue==0)
   epsilon <- sigmaepsTrue*z

x <- runif(n) ; y <- fTrue(x) + epsilon

# Obtain design matrices:

X <- cbind(1,x)
a <- min(x) ;   b <- max(x)
numIntKnots <- 25
intKnots <-  quantile(unique(x),seq(0,1,length = numIntKnots + 2)
                      [-c(1,numIntKnots + 2)])
Z <- ZOSull(x,intKnots = intKnots,range.x = c(a,b))
ncZ <- ncol(Z)

# Plot the data:

par(mfrow=c(1,1))
plot(x,y,col = "dodgerblue",bty="l",cex=0.5,
     main="simulated data and true regression function",
     col.main="navy")
xg <- seq(min(x),max(x),length=1001)
lines(xg,fTrue(xg),col="indianred3",lwd=2)
cat("Hit Enter to continue.\n")
ans <- readline()

# Set MCMC parameters:

nWarm <- 1000          # Length of warm-up.
nKept <- 1000          # Size of the kept sample.

# Specify model in Stan:

penSplGEVModel <- 
   "data
   {
      int<lower=1> n;            int<lower=1> ncZ;
      vector[n] y;               matrix[n,2] X;
      matrix[n,ncZ] Z;          
      real<lower=0> sigmabeta;   real<lower=0> Au;
      real<lower=0> Aeps;        real<lower=0> sigmaxi;  
   }
   parameters 
   {
      vector[2] beta;            vector[ncZ] u;
      real<lower=0> sigmaeps ;   real<lower=0> sigmau ;
      real xi ; 
   }
   transformed parameters
   {
      vector[n] arg;              vector[n] supparg;
      vector[n] logdens;
      arg = (y - X*beta - Z*u)/sigmaeps;
      supparg = 1 + xi*arg;
      for (i in 1:n)
      {
         if (xi==0)
            logdens[i] = -arg[i] - exp(-arg[i]) - log(sigmaeps);
         if (xi!=0)
         {
            if (supparg[i]>0)
            {
               logdens[i] = (-(1/xi)-1)*log(1 + xi*(arg[i]))
                           - pow((1 + xi*arg[i]),(-(1/xi)))
                           - log(sigmaeps);
            }
         }
      }
   }
   model 
   {
      for (i in 1:n) target += logdens[i];
      beta ~ normal(0,sigmabeta);   sigmaeps ~ cauchy(0,Aeps);
      sigmau ~ cauchy(0,Au);        xi ~ normal(0,sigmaxi); 
   }"

# Store data in a list in format required by Stan:

allData <- list(n=n,ncZ=ncZ,y=y,X=X,Z=Z,sigmabeta=1e5,Au=Au,Aeps=Aeps,
                sigmaxi=sigmaxi)

# Compile code for model:

stanCompilObj <- stan(model_code=penSplGEVModel,data=allData,
                      iter=1,chains=1)

# Obtain MCMC samples for each parameter using Stan:

stanObj <-  stan(model_code=penSplGEVModel,data=allData,
                 warmup=nWarm,iter=(nWarm+nKept),
                 chains=1,refresh=100,fit=stanCompilObj)

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
xiMCMC <- as.vector(extract(stanObj,"xi",permuted=FALSE))

  
# Obtain MCMC samples of regression curves over a fine grid:

ng <- 101
xg <- seq(a,b,length=ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots=intKnots,range.x=c(a,b))
fhatMCMC <- Xg%*%betaMCMC + Zg%*%uMCMC
fhatg <- apply(fhatMCMC,1,mean)
credLower <- apply(fhatMCMC,1,quantile,0.025)
credUpper <- apply(fhatMCMC,1,quantile,0.975)

# Display the fit:

par(mai=c(1,1.1,0.1,0.1))
cex.labVal <- 2   ;   cex.axisVal <- 1.5
plot(x,y,type="n",
     bty="l",xlim=range(xg),ylim=range(c(credLower,credUpper,y)),
     cex.lab=cex.labVal,cex.axis=cex.axisVal)
polygon(c(xg,rev(xg)),c(credLower,rev(credUpper)),
        col="palegreen",border=FALSE)
lines(xg,fhatg,col="darkgreen",lwd=2)
points(x,y,col="dodgerblue",cex=0.5)
lines(xg,fTrue(xg),col="indianred3")
cat("Hit Enter to continue.\n")
ans <- readline()

# Do some summaries and diagnostic checking of the MCMC:

indQ1 <- length(xg[xg<quantile(x,0.25)])
indQ2 <- length(xg[xg<quantile(x,0.50)])
indQ3 <- length(xg[xg<quantile(x,0.75)])
fhatQ1MCMC <- fhatMCMC[indQ1,]
fhatQ2MCMC <- fhatMCMC[indQ2,]
fhatQ3MCMC <- fhatMCMC[indQ3,]
MCMClist <- list(cbind(sigmaepsMCMC,xiMCMC,
                 fhatQ1MCMC,fhatQ2MCMC,fhatQ3MCMC))
parNamesVal <- list(c("error","standard","deviation"),
                    c("Gener. Extr. Val.",
                      "shape parameter",
                      expression(paste("(",xi,")"))),
                    c("reg'n func. est.","at 1st quartile","of x"),
                    c("reg'n func. est.","at 2nd quartile","of x"),
                    c("reg'n func. est.","at 3rd quartile","of x"))
summMCMC(MCMClist,parNames=parNamesVal,
        addTruthToKDE=c(sigmaepsTrue,xiTrue,fTrue(xg[indQ1]),
                        fTrue(xg[indQ2]),fTrue(xg[indQ3])))

########### End of penSplGEV ############
