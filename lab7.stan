data
{
  int<lower=1> n;            
  int<lower=1> ncZ;
  vector[n] y;               
  matrix[n,2] X;
  matrix[n,ncZ] Z;          
  real<lower=0> sigmabeta;   
  real<lower=0> Au;
  real<lower=0> Aeps;        
  real<lower=0> sigmaxi;  
}
parameters 
{
  vector[2] beta;            
  vector[ncZ] u;
  real<lower=0> sigmaeps ;   
  real<lower=0> sigmau ;
  real xi ; 
}
transformed parameters
{
  vector[n] arg;              
  vector[n] supparg;
  vector[n] logdens;
  arg = (y - X*beta - Z*u)/sigmaeps;
  supparg = 1 + xi*arg;
  for (i in 1:n)
  {
     if (xi==0){
       logdens[i] = -arg[i] - exp(-arg[i]) - log(sigmaeps);
     }
        
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
  for (i in 1:n) {
    target += logdens[i];
  }
  
  beta ~ normal(0,sigmabeta);   
  sigmaeps ~ cauchy(0,Aeps);
  sigmau ~ cauchy(0,Au);        
  xi ~ normal(0,sigmaxi); 
}
