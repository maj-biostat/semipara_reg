data {                          
  int<lower=0> x[100]; 
}
parameters {
  real<lower=0, upper=1> p;
}
// transformed parameters {
//   vector[N] eta;
//   vector[N_pred] eta_new;
//   eta = alpha +  X * b;
//   eta_new = alpha +  X_new * b;
// 
// }
model {

  // target += uniform_lpdf(p | 0, 1);
  // 
  // // likelihood 
  // target += bernoulli_logit_lpmf(y | eta);
  
  target += uniform_lpdf(p | 0, 1);
  target += binomial_lpmf(x | 1, p);
  


}
// generated quantities { 
//   
//   /*
//   predict the values for those unenrolled
//   */
//   vector[N_pred] y_rep;
// 
//   for (i in 1:N_pred) {
//     y_rep[i] = bernoulli_rng(inv_logit(eta_new[i]));
//   }
// 
// }
