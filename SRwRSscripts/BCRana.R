######### R script: BCRana ##########

# For performing penalised spline-based nonparametric 
# regression for classical measurement error, for
# analysis of actual data from Berry, Carroll & Ruppert 
# (Journal of the American Statistical Association, 2002).

# Last changed: 13 AUG 2018

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Load required packages:

library(HRW);  library(rstan)

# Set MCMC parameters:

nWarm <- 100            # Length of warm-up.
nKept <- 200            # Size of the kept sample.
nThin <- 1              # Thinning factor. 

# Set measurement error standard deviation:

sigmaW <- sqrt(0.35)

# Set the number of spline basis functions:

ncZ <- 30

# Load in data:

data(BCR)
w <- BCR$w   ;   y <- BCR$y
treatIndic <- as.numeric(as.character(BCR$status) == "treatment")
n <- length(y)
        
# Specify model in Stan:

BCRanaModel <- 
   'data
   {
      int<lower = 1> n;           int<lower = 1> ncZ;
      vector[n] y;              vector[n] w;        
      vector[n] treatIndic;              
      real<lower = 0> sigmaBeta;  real<lower = 0> sigmaMu;
      real<lower = 0> sigmaW;     
      real<lower = 0> Ax;         real<lower = 0> Aeps;
      real<lower = 0> Au;           
   }
   parameters 
   {
      vector[4] beta;       
      vector[ncZ] uControl;           vector[ncZ] uTreatmt;            
      vector[n] x;
      real muX;                       real<lower = 0> sigmaX;
      real<lower = 0> sigmaUcontrol;    real<lower = 0> sigmaUtreatmt;
      real<lower = 0> sigmaEps;   
   }
   transformed parameters 
   {
      matrix[n,4] X;            vector[ncZ] knots;
      matrix[n,ncZ] Z;
      for (k in 1:ncZ)
         knots[k] = ((ncZ+1-k)*min(x)+k*max(x))/(ncZ+1);
      for (i in 1:n)
      {
         X[i,1] = 1     ;   X[i,2] = treatIndic[i];
         X[i,3] = x[i]  ;   X[i,4] = treatIndic[i]*x[i];
         for (k in 1:ncZ)   
            Z[i,k] = (x[i]-knots[k])*step(x[i]-knots[k]);
      }
   }
   model 
   {
      for (i in 1:n)
         y[i] ~ normal((dot_product(beta,X[i])
                       + dot_product(uControl,((1-treatIndic[i])*Z[i]))
                       + dot_product(uTreatmt,(treatIndic[i]*Z[i]))),sigmaEps);
      x  ~ normal(muX,sigmaX); 
      w ~ normal(x,sigmaW);
      uControl ~ normal(0,sigmaUcontrol) ; uTreatmt ~ normal(0,sigmaUtreatmt); 
      beta ~ normal(0,sigmaBeta);
      muX ~ normal(0,sigmaMu); sigmaX ~ cauchy(0,Ax);       
      sigmaEps ~ cauchy(0,Aeps); 
      sigmaUcontrol ~ cauchy(0,Au); sigmaUtreatmt ~ cauchy(0,Au); 
   }'

# Set up input data:

allData <- list(n = n,ncZ = ncZ,w = w,y = y,treatIndic = treatIndic,
                sigmaW = sigmaW,sigmaMu = 1e5,
                sigmaBeta = 1e5,Ax = 1e5,Aeps = 1e5,Au = 1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = BCRanaModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = BCRanaModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 25,
                 fit = stanCompilObj)

# Extract relevant MCMC samples:

betaMCMC <- NULL
for (j in 1:4)
{
   charVar <- paste("beta[",as.character(j),"]",sep = "") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar,permuted = FALSE))
}
uControlMCMC <- NULL  ;  uTreatmtMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("uControl[",as.character(k),"]",sep = "") 
   uControlMCMC <- rbind(uControlMCMC,extract(stanObj,charVar,permuted = FALSE))
                   
   charVar <- paste("uTreatmt[",as.character(k),"]",sep = "") 
   uTreatmtMCMC <- rbind(uTreatmtMCMC,extract(stanObj,charVar,permuted = FALSE))
}
muXMCMC <- as.vector(extract(stanObj,"muX",permuted = FALSE))
sigmaXMCMC <- as.vector(extract(stanObj,"sigmaX",permuted = FALSE))
sigmaEpsMCMC <- as.vector(extract(stanObj,"sigmaEps",permuted = FALSE))
sigmaUcontrolMCMC <- as.vector(extract(stanObj,"sigmaUcontrol",permuted = FALSE))
sigmaUtreatmtMCMC <- as.vector(extract(stanObj,"sigmaUtreatmt",permuted = FALSE))
xMCMC <- NULL
for (i in 1:n)
{
   charVar <- paste("x[",as.character(i),"]",sep = "") 
   xMCMC <- rbind(xMCMC,extract(stanObj,charVar,permuted = FALSE))
}
knotsMCMC <- NULL
for (i in 1:ncZ)
{
   charVar <- paste("knots[",as.character(i),"]",sep = "") 
   knotsMCMC <- rbind(knotsMCMC,extract(stanObj,charVar,permuted = FALSE))
}

# Plot fit:

estFunCol <- "DarkGreen"; 
varBandCol <- "PaleGreen"
ng <- 201
cex.labVal <- 1.8
xLow <- -1.8 ; xUpp <- 2.5
xg <- seq(xLow,xUpp,length = ng)

nMCMC <- length(muXMCMC)
ContrastMCMC <- matrix(0,ng,nMCMC)
for (g in 1:nMCMC)
{
   ContrastMCMC[,g] <- betaMCMC[2,g]*rep(1,ng) + betaMCMC[4,g]*xg 
   for (k in 1:ncZ)
      ContrastMCMC[,g] <- (ContrastMCMC[,g] + (uTreatmtMCMC[k,g]-uControlMCMC[k,g])
                           *(xg-knotsMCMC[k,g])*(xg>knotsMCMC[k,g]))
}

credLower <- apply(ContrastMCMC,1,quantile,0.05)
credUpper <- apply(ContrastMCMC,1,quantile,0.95)
Contrastg <- apply(ContrastMCMC,1,mean)

par(mfrow = c(1,1),mai = c(1.02,0.82,0.82,0.42))
plot(0,0,type = "n",bty = "l",xlim = c(xLow,xUpp),
     ylim = c(-3,3),xlab = "true baseline score",ylab = "treatment effect",
     cex.lab = cex.labVal)
polygon(c(xg,rev(xg)),c(credLower,rev(credUpper)),col = varBandCol,border = FALSE)

lines(xg,Contrastg,lwd = 2,col = estFunCol)
abline(0,0,col = "navy")
readline("Hit Enter to continue.\n")

# Display summary of MCMC samples for model parameters:

xHat <- apply(xMCMC,1,mean)
indQ1 <- length(xg[xg <= quantile(xHat,0.25)])
ContrastQ1MCMC <- ContrastMCMC[indQ1,]
indQ2 <- length(xg[xg <= quantile(xHat,0.50)])
ContrastQ2MCMC <- ContrastMCMC[indQ2,]
indQ3 <- length(xg[xg <= quantile(xHat,0.75)])
ContrastQ3MCMC <- ContrastMCMC[indQ3,]

parms <- list(cbind(muXMCMC,sigmaXMCMC,sigmaEpsMCMC,
                    ContrastQ1MCMC,ContrastQ2MCMC,ContrastQ3MCMC))
parNamesVal <- list(c(expression(mu[x])),c(expression(sigma[x])),
                      c(expression(sigma[epsilon])),
                      c("contrast funct.","at 1st quart. of","true base. score"),
                      c("contrast funct.","at 2nd quart. of","true base. score"),
                      c("contrast funct.","at 3rd quart. of","true base. score"))

summMCMC(parms,parNames = parNamesVal,numerSummCex = 1,columnHeadCex = 2.9,
         KDEvertLine = FALSE)

readline("Hit Enter to continue.\n")

# Display summary for some of the unobserved x values:

sampInds <-  c(10,18,27,44,59)
parms <- list(t(xMCMC[sampInds,]))
parNamesVal <- list(c(expression(x[10])),
                      c(expression(x[18])),
                      c(expression(x[27])),
                      c(expression(x[44])),
                      c(expression(x[59])))

summMCMC(parms,parNames = parNamesVal,columnHeadCex = 2.9,KDEvertLine = FALSE)

############ End of BCRana ############

