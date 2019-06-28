######### R script: PIDana ##########

# For performing penalised spline-based logistic nonparametric 
# regression when the predictor is subject to missingness
# not at random for the Pima Indians diabetes data.

# Last changed: 13 AUG 2018

# Set flag for code compilation (needed if 
# running script first time in current session):

compileCode <- TRUE

# Load required packages:

library(HRW) ; library(mlbench) ; library(rstan)

# Set MCMC parameters:

nWarm <- 100        # Length of warm-up.
nKept <- 500        # Size of the kept sample.

# Obtain data:

data(PimaIndiansDiabetes)
r <- as.numeric(PimaIndiansDiabetes$mass>0)
y <- as.numeric(as.character(PimaIndiansDiabetes$diabetes)=="pos")
yxObs <- y[r==1]
yxUnobs <- y[r==0]
xObsOrig <- PimaIndiansDiabetes$mass[r==1]
mean.xObs <- mean(xObsOrig)  ; sd.xObs <- sd(xObsOrig)
xObs <- (xObsOrig - mean.xObs)/sd.xObs

# Set sample sizes:

nObs <- length(yxObs)
nUnobs <- length(yxUnobs)

# Obtain the knots based on the observed data:

ncZ <- 30
knots <- seq(min(xObs),max(xObs),length=(ncZ+2))[-c(1,ncZ+2)]
ZxObs <- outer(xObs,knots,"-")
ZxObs <- ZxObs*(ZxObs>0)

# Specify model in Stan:

logisNpRegMNARModel <- 
   'data
   {
      int<lower=1> nObs;        int<lower=1> nUnobs;
      int<lower=1> n;           int<lower=1> ncZ;
      int<lower=0,upper=1>      yxObs[nObs];
      int<lower=0,upper=1>      yxUnobs[nUnobs];
      vector[nObs] xObs;        vector[ncZ]  knots;   
      matrix[nObs,ncZ] ZxObs;    
      int<lower=0,upper=1>   r[n]; 
      real<lower=0> sigmaBeta;  real<lower=0> sigmaMu;
      real<lower=0> sigmaPhi;
      real<lower=0> Ax;         real<lower=0> Au;         
   }
   transformed data
   {
      int<lower=0,upper=1> y[n];
      for (i in 1:nObs)
         y[i] = yxObs[i]; 
      for (i in 1:nUnobs)
         y[i+nObs] = yxUnobs[i];
   }
   parameters 
   {
      vector[2] beta;            vector[ncZ] u;
      vector[2] phi;
      real muX;                  real<lower=0> sigmaX;
      real<lower=0> sigmaU;      real xUnobs[nUnobs];
   }
   transformed parameters 
   {
      matrix[n,2] X;       matrix[n,ncZ] Z;
      for (i in 1:nObs)
      {
         X[i,1] = 1   ;   X[i,2] = xObs[i] ;
         Z[i] = ZxObs[i] ;
      }
      for (i in 1:nUnobs)
      {
         X[i+nObs,1] = 1    ;   X[i+nObs,2] = xUnobs[i];
         for (k in 1:ncZ)   
            Z[i+nObs,k] = (xUnobs[i]-knots[k])*step(xUnobs[i]-knots[k]);
      }
   }
   model 
   {
      y ~ bernoulli_logit(X*beta+Z*u); 
      r ~ bernoulli_logit(X*phi);
      col(X,2) ~ normal(muX,sigmaX); 
      u ~ normal(0,sigmaU) ; beta ~ normal(0,sigmaBeta);
      muX ~ normal(0,sigmaMu); phi ~ normal(0,sigmaPhi);
      sigmaX ~ cauchy(0,Ax);   sigmaU ~ cauchy(0,Au);
   }'

# Set up input data:

allData <- list(nObs=nObs,nUnobs=nUnobs,n=(nObs+nUnobs),
                ncZ=ncZ,xObs=xObs,yxObs=yxObs,yxUnobs=yxUnobs,knots=knots,
                ZxObs=ZxObs,r=r,sigmaMu=1e5,sigmaBeta=1e5,sigmaPhi=1e5,
                sigmaX=1e5,Ax=1e5,Au=1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code=logisNpRegMNARModel,data=allData,
                         iter=1,chains=1)

# Perform MCMC:

stanObj <-  stan(model_code=logisNpRegMNARModel,data=allData,warmup=nWarm,
                 iter=(nWarm+nKept),chains=1,refresh=10,fit=stanCompilObj)

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
muXMCMC <- extract(stanObj,"muX",permuted=FALSE)
muXorigMCMC <- mean.xObs + sd.xObs*muXMCMC
sigmaXMCMC <- extract(stanObj,"sigmaX",permuted=FALSE)

sigmaXorigMCMC <- sd.xObs*sigmaXMCMC
phiMCMC <- NULL
for (j in 1:2)
{
   charVar <- paste("phi[",as.character(j),"]",sep="") 
   phiMCMC <- rbind(phiMCMC,extract(stanObj,charVar,permuted=FALSE))
}
phi1origMCMC <-  phiMCMC[2,]/sd.xObs
phi0origMCMC <-  phiMCMC[1,] - phi1origMCMC *mean.xObs
xUnobsMCMC <- NULL
for (i in 1:nUnobs)
{
   charVar <- paste("xUnobs[",as.character(i),"]",sep="") 
   xUnobsMCMC <- rbind(xUnobsMCMC,extract(stanObj,charVar,permuted=FALSE))
}

# Plot fit:

par(mai=c(1.0,0.9,0.1,0.05))
cex.labVal <- 2.1
obsCol <- "dodgerblue" ; unobsCol <- "red"
estFunCol <- "DarkGreen"; 
varBandCol <- "PaleGreen"
ng <- 101
xg <- seq(min(xObs),max(xObs),length=ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- outer(xg,knots,"-")
Zg <- Zg*(Zg>0)
etaMCMC <-   Xg%*%betaMCMC + Zg%*%uMCMC
probMCMC <- 1/(1+exp(-etaMCMC)) 
credLower <- apply(probMCMC,1,quantile,0.025)
credUpper <- apply(probMCMC,1,quantile,0.975)
probg <- apply(probMCMC,1,mean)

xgOrig <- mean.xObs + sd.xObs*xg
plot(xgOrig,probg,type="n",bty="l",xlab="body mass index",
     ylab="estimated probability of diabetes",xlim=range(xgOrig),
     ylim=c(-0.1,1.1),cex.lab=cex.labVal)
polygon(c(xgOrig,rev(xgOrig)),c(credLower,rev(credUpper)),
       col=varBandCol,border=FALSE)
lines(xgOrig,probg,lwd=2,col=estFunCol)

points(xObsOrig[yxObs==0],runif(length(xObsOrig[yxObs==0]),-0.09,-0.01),
       cex=0.5,col=obsCol)
points(xObsOrig[yxObs==1],runif(length(xObsOrig[yxObs==1]),1.01,1.09),
       cex=0.5,col=obsCol)

lines(c(min(xgOrig),max(xgOrig)),rep(0,2),lty=2,col="navy")
lines(c(min(xgOrig),max(xgOrig)),rep(1,2),,lty=2,col="navy")

# Add Bayes estimates of missing predictor data:

xUnobsOrigMCMC <- mean.xObs + sd.xObs*xUnobsMCMC
xUnobsOrigHat <- apply(xUnobsOrigMCMC,1,mean)

points(xUnobsOrigHat[yxUnobs==0],
       runif(length(xUnobsOrigHat[yxUnobs==0]),-0.09,-0.01),
       cex=0.75,col=unobsCol,lwd=2)
points(xUnobsOrigHat[yxUnobs==1],
       runif(length(xUnobsOrigHat[yxUnobs==1]),1.01,1.09),
       cex=0.75,col=unobsCol,lwd=2)

# Add legend:

legend(18,0.95,legend=c("observed data","Bayes estimates of missing data"),
       col=c(obsCol,unobsCol),pch=rep(1,2),pt.cex=0.75,pt.lwd=2,cex=1.3)

cat("Hit Enter to continue.\n")
ans <- readline()

# Display summary of MCMC samples for model parameters:

indQ1 <- length(xgOrig[xgOrig<=quantile(xObsOrig,0.25)])
probQ1MCMC <- probMCMC[indQ1,]
indQ2 <- length(xgOrig[xgOrig<=quantile(xObsOrig,0.50)])
probQ2MCMC <- probMCMC[indQ2,]
indQ3 <- length(xgOrig[xgOrig<=quantile(xObsOrig,0.75)])
probQ3MCMC <- probMCMC[indQ3,]

parms <- list(cbind(muXorigMCMC,sigmaXorigMCMC,phi0origMCMC,
                    phi1origMCMC,probQ1MCMC,probQ2MCMC,
                    probQ3MCMC,xUnobsOrigMCMC[1,]))
parNamesVal <- list(c(expression(mu[x])),c(expression(sigma[x])),
                      c(expression(phi[0])),c(expression(phi[1])),
                      c("est. probab. at","1st quant. age"),
                      c("est. probab. at","2nd quant. age"),
                      c("est. probab. at","3rd quant. age"),
                      c("1st unobserved","predictor"))

summMCMC(parms,parNames=parNamesVal,columnHeadCex=2.9,
         numerSummCex=1)

########## End of PIDana ##########
