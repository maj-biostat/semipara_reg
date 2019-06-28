########## R script: penSplinesDetails ##########

# For providing details on penalized splines.

# Last changed: 16 MAY 2016 by M.P. Wand.

# Load required packages:

library(HRW) ; library(nlme) ; library(splines)

# Define wait() function:

wait <- function()
{
   cat("Hit Enter to continue\n")
   ans <- readline()
   invisible()
}

# Load illustrative dataset:

data(WarsawApts)

x <- WarsawApts$construction.date
y <- WarsawApts$areaPerMzloty

ng <- 1001
a <- min(x) ; b <- max(x)
range.x <- c(a,b)
xg <- seq(a,b,length=ng)

X <- cbind(rep(1,length(x)),x)
numIntKnots <- 20
intKnots <- quantile(unique(x),seq(0,1,length=
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])

# Define empty plot function:

empty.plot <- function(bty="n") 
{
   par(mfrow=c(1,1))
   plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n", 
        yaxt="n",xlab="",ylab="",bty="o")
   invisible()
}

# Explain point of this script:

empty.plot()
text(0.1,0.9,adj=0,"This R script explains how penalized splines",col="black",cex=1.2)
text(0.1,0.7,adj=0,"based on truncated lines can be improved",col="black",cex=1.2)
text(0.1,0.5,adj=0,"via the use of better basis functions.",col="black",cex=1.2)
text(0.1,0.3,adj=0,"It works best with this graphics window made large.",col="black",cex=1.2)
text(0.1,0.1,adj=0,"and the R session window made tall and thin.",col="black",cex=1.2)

wait()

# Introduce data:

empty.plot()
text(0.1,0.6,adj=0,"The next screen shows a scatterplot",col="black",cex=1.2)
text(0.1,0.4,adj=0,"of Warsaw apartments data that we wish to smooth.",col="black",cex=1.2)
wait()

par(mfrow=c(1,1))
plot(x,y,bty="l",type="n",xlab="construction date (years)",
      ylab="area (square metres) per million zloty",xlim=range.x)
points(x,y,col="dodgerblue",lwd=2)
wait()

# Explain truncated line-based smoother:

empty.plot()
text(0.1,0.95,adj=0,"The truncated line penalized spline takes the form",col="black",cex=1.2)
text(0.1,0.80,adj=0,expression(paste(hat(y),"=C(",C^{T},"C+",lambda,D,")",phantom()^{-1},C^{T},"y")),
    col="black",cex=1.2)
text(0.1,0.65,adj=0,"where C contains truncated line basis functions",cex=1.2) 
text(0.1,0.50,adj=0,"(plus the straight line basis functions: 1 and x)",col="black",cex=1.2)
text(0.1,0.35,adj=0,"and D=diag(0,0,1,...,1). This diagonal form of the penalty",
     col="black",cex=1.2)
text(0.1,0.2,adj=0," allows easy implementation in mixed model software (Sessions 2 and 3).",col="black",cex=1.2)
text(0.1,0.05,adj=0,"The next screen shows the fit and basis functions.",col="black",cex=1.2)
    
wait()

# Do truncated line fit:

Z <- outer(x,intKnots,"-")
Z <- Z*(Z>0)
group <- rep(1,length(x))
data.fr <- groupedData(y~x|group,data=data.frame(x,y))
fit <- lme(y~-1+X,random=pdIdent(~-1+Z),data=data.fr)
betaHat <- fit$coef$fixed
uHat <- unlist(fit$coef$random)

Xg <- cbind(rep(1,ng),xg)
Zg <- outer(xg,intKnots,"-")
Zg <- Zg*(Zg>0)
fHatg <- Xg%*%betaHat + Zg%*%uHat

par(mfrow=c(2,1))
plot(x,y,bty="l",xlab="construction date (years)",
     ylab="area (square metres) per million zloty",xlim=range.x)
points(x,y,col="dodgerblue",lwd=2)
lines(xg,fHatg,col="darkgreen",lwd=2)
lines(c(min(xg),max(xg)),rep(min(y),2),col="slateblue")
for (k in 1:numIntKnots)
   points(intKnots[k],min(y),pch=18,cex=2,col="darkmagenta")

plot(0,0,type="l",xlim=range(xg),ylim=c(-0.1,1.2),bty="l",
     xlab="construction date (years)",ylab="basis function",
     main="spline basis functions")

for (k in 1:numIntKnots)
{
   Zkg <- outer(xg,intKnots[k],"-")
   Zkg <- 0.2*Zkg*(Zkg>0)
   lines(xg,Zkg,col=k,lwd=2)
}
lines(c(min(xg),max(xg)),rep(0,2),col="slateblue")
for (k in 1:numIntKnots)
   points(intKnots[k],0,pch=18,cex=2,col="darkmagenta")
wait()

# Comment on the short-comings of truncated line-based smoothing:

empty.plot()
text(0.1,0.9,adj=0,"Although they are nice and simple, truncated line basis",col="black",cex=1.2)
text(0.1,0.7,adj=0,"functions have room for improvement. They are not very smooth,",
     col="black",cex=1.2)
text(0.1,0.5,adj=0,"leading to a kinkiness in the penalized spline (unless many more knots",
    col="black",cex=1.2)
text(0.1,0.3,adj=0,"are added). In addition, truncated lines are unbounded and far from",
    col="black",cex=1.2)
text(0.1,0.1,adj=0,"orthogonal. This can lead to computational problems.",
    col="black",cex=1.2)
wait()

# Introduce O'Sullivan penalized splines:

empty.plot()
text(0.1,0.9,adj=0,"O'Sullivan penalized splines are of the form",col="black",cex=1.2)
text(0.1,0.7,adj=0,expression(paste(hat(y),"=B(",B^{T},"B+",lambda,Omega,")",phantom()^{-1},B^{T},"y")),
    col="black",cex=1.2)
text(0.1,0.5,adj=0,"where the B matrix contains cubic B-splines and the (i,j)",col="black",cex=1.2)
text(0.1,0.4,adj=0,"entry of the penalty matrix is",col="black",cex=1.2)
text(0.1,0.3,adj=0,expression(paste(Omega[ij],"=",integral(g(x)*dx,a,b))),
     col="black",cex=1.2)
text(0.1,0.1,adj=0,expression(paste("and g(x)=",B[i]^{(2)},"(x)",B[j]^{(2)},"(x).")),
     col="black",cex=1.2)
wait()

empty.plot()
text(0.1,0.9,adj=0,"The next screen shows the O'Sullivan penalized spline fit",
     col="black",cex=1.2)
text(0.1,0.7,adj=0,"and the cubic B-spline basis functions.",col="black",cex=1.2)
text(0.1,0.5,adj=0,"They are smoother, bounded and closer to being orthogonal",col="black",cex=1.2)
text(0.1,0.3,adj=0,"which makes them more preferrable in practice.",col="black",cex=1.2)
text(0.1,0.1,adj=0,"They also have very good boundary and extrapolation properties.",
     col="black",cex=1.2)

wait()

Z <- ZOSull(x,range.x,intKnots)

fit <- lme(y~-1+X,random=pdIdent(~-1+Z),data=data.fr)
betaHat <- fit$coef$fixed
uHat <- unlist(fit$coef$random)
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,range.x,intKnots)
fHatg <- Xg%*%betaHat + Zg%*%uHat

par(mfrow=c(2,1))
plot(x,y,bty="l",xlab="construction date (years)",
     ylab="area (square meters) per million zloty",xlim=range.x)
points(x,y,col="dodgerblue",lwd=2)
lines(xg,fHatg,col="darkgreen",lwd=2)
lines(c(min(xg),max(xg)),rep(min(y),2),col="slateblue")
for (k in 1:numIntKnots)
   points(intKnots[k],min(y),pch=18,cex=2,col="darkmagenta")

plot(0,0,type="l",xlim=range(xg),ylim=c(-0.1,1.2),bty="l",
     xlab="construction date (years)",
     ylab="area (square meters) per million zloty",
     main="spline basis functions")

Bgmat <- bs(xg,knots=intKnots,degree=3,
            Boundary.knots=c(a,b),intercept=TRUE)
for (k in 1:ncol(Bgmat))
   lines(xg,Bgmat[,k],col=k,lwd=2)
lines(c(min(xg),max(xg)),rep(0,2),col="slateblue")
for (k in 1:numIntKnots)
   points(intKnots[k],0,pch=18,cex=2,col="darkmagenta")
wait()

# Explain the need to transform the basis for mixed model representation:

empty.plot()
text(0.1,0.95,adj=0,"Since the penalty of O'Sullivan splines is non-diagonal",
     col="black",cex=1.2)
text(0.1,0.80,adj=0,"implementation in mixed model software (Sessions 2 and 3)",
     col="black",cex=1.2)
text(0.1,0.65,adj=0,"is not straightforward. To overcome this, we transform",
     col="black",cex=1.2)
text(0.1,0.50,adj=0,"the cubic B-spline basis functions so that the penalty",
     col="black",cex=1.2)
text(0.1,0.35,adj=0,"is of the form D=diag(0,0,1,...,1). The next screen shows",
     col="black",cex=1.2)
text(0.1,0.20,adj=0,"the spline basis functions after undergoing this transformation.",
     col="black",cex=1.2)
    
wait()

# Show the transformed basis functions:

par(mfrow=c(2,1))
plot(x,y,bty="l",xlab="construction date (years)",
     ylab="area (square meters) per million zloty",xlim=range.x)
points(x,y,col="dodgerblue",lwd=2)
lines(xg,fHatg,col="darkgreen",lwd=2)
lines(c(min(xg),max(xg)),rep(min(y),2),col="slateblue")
for (k in 1:numIntKnots)
   points(intKnots[k],min(y),pch=18,cex=2,col="darkmagenta")

plot(0,0,type="l",xlim=range(xg),ylim=range(Zg),bty="l",
     xlab="construction date (years)",ylab="basis function",
     main="spline basis functions")

for (k in 1:ncol(Zg))
   lines(xg,Zg[,k],col=k,lwd=2)
lines(c(min(xg),max(xg)),rep(0,2),col="slateblue")
for (k in 1:numIntKnots)
   points(intKnots[k],0,pch=18,cex=2,col="darkmagenta")
wait()

# Describe nature of new basis functions:

empty.plot()
text(0.1,0.95,adj=0,"Note that the transformed functions have an interesting",
     col="black",cex=1.2)
text(0.1,0.83,adj=0,"multi-resolution nature. This is easier to see in the",
     col="black",cex=1.2)
text(0.1,0.71,adj=0,"next screen, which blows up the higher frequency basis functions.",col="black",cex=1.2)
text(0.1,0.59,adj=0,"This functions have an eigenvector interpretation based",
     col="black",cex=1.2)
text(0.1,0.47,adj=0,expression(paste(" on the spectral decomposition of the ",Omega," matrix.")),
     col="black",cex=1.2)
text(0.1,0.35,adj=0,"Most importantly, the penalized spline is of the form",col="black",cex=1.2)
text(0.1,0.23,adj=0,expression(paste(hat(y),"=C(",C^{T},"C+",lambda,D,")",phantom()^{-1},C^{T},"y,    D=diag(0,0,1,...,1)")),
    col="black",cex=1.2)
text(0.1,0.11,adj=0,"which facilitates direct implementation into mixed model and Bayesian models.",
    col="black",cex=1.2)

wait()

# Show the blown up version of the basis functions:

par(mfrow=c(2,1))
plot(x,y,bty="l",xlab="construction date (years)",
     ylab="area (square meters) per million zloty",xlim=range.x)
points(x,y,col="dodgerblue",lwd=2)
lines(xg,fHatg,col="darkgreen",lwd=2)
lines(c(min(xg),max(xg)),rep(min(y),2),col="slateblue")
for (k in 1:numIntKnots)
   points(intKnots[k],min(y),pch=18,cex=2,col="darkmagenta")

plot(0,0,type="l",xlim=range(xg),ylim=c(-2,2),bty="l",
     xlab="age",ylab="basis function",main="spline basis functions (vertically magnified)")

for (k in 1:ncol(Zg))
   lines(xg,Zg[,k],col=k,lwd=2)
lines(c(min(xg),max(xg)),rep(0,2),col="slateblue")
for (k in 1:numIntKnots)
   points(intKnots[k],0,pch=18,cex=2,col="darkmagenta")

########## End of penSplinesDetails ##########
