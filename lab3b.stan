data
{
  int<lower=1> n;         
  int<lower=1> ncZ;
  vector[n] y;            
  matrix[n,2] X;
  matrix[n,ncZ] Z;        
  real<lower=0> sigmaBeta;
  real<lower=0> Au;       
  real<lower=0> Aeps;
}
parameters 
{
  vector[2] beta;          
  vector[ncZ] u;
  real<lower=0> sigmaeps;  
  real<lower=0> sigmau;
}
model 
{
  y ~ normal(X*beta + Z*u,sigmaeps);
  u ~ normal(0,sigmau); 
  beta ~ normal(0,sigmaBeta); 
  sigmaeps ~ cauchy(0,Aeps); 
  sigmau ~ cauchy(0,Au);
}

