data{
  int<lower=1> N;                      // Number of observations
  int<lower=0> N_covariates;           // Number of covariates
  vector[N] y;                         // Observations
  matrix[N,N] W;                       // Adjacency matrix
  matrix[N,N_covariates] B;            // Design matrix
}

parameters{
  real<lower=0> sigma;           // Standard deviation of NIG noise
  real<lower=0, upper=1> rho;    // Spatial range parameter
  vector[N_covariates] beta;     // Regression coefficients  
}

transformed parameters{
  vector[N] x = (y - B*beta)/sigma;                   // Spatial effects
}

model{
  matrix[N,N] D = add_diag(-rho*W, 1);                             // D = I - rho W;
  x ~ multi_normal_prec(rep_vector(0,N), D'*D);

  //prior layer---------------------------------
  //prior for intercept
  beta[1]     ~ normal(0.0, 1000);
  //prior for regression coefficients
  beta[2]     ~ normal(0.0, 10);
  beta[3]     ~ normal(0.0, 10);
  //prior for kappa
  rho     ~ uniform(0,1);
  //prior for sigma
  sigma   ~ normal(0,5);

}

generated quantities { 
  vector[N] y_fit = B*beta + sigma*x;
}


