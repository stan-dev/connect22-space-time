functions {

  /**
   * Returns the log-likelihood of the NIG distribution with variance 
   * correction at x (Var[x] = sigma^2*h)
   *
   * @param x 
   * @param sigma Standard deviation 
   * @param eta   First flexibility parameter (controls heaviness of tails)
   * @param zeta  Second flexibility parameter (controls asymmetry of tails)
   * @param h     Distance between locations or area of the basis function
   * @return      Log-likelihood of the NIG distribution with variance correction
   */
  real NIG_var_correction_lpdf(real x, real eta, real zeta, real h){
    real sigmas     = 1/sqrt(1+zeta^2*eta); //variance correction
    real hyp_alpha  = sqrt(1/eta+zeta^2)/sigmas;
    real hyp_beta   = zeta/sigmas;
    real hyp_delta  = sigmas*sqrt(1/eta)*h;
    real hyp_mu     = -sigmas*zeta*h;
    return (sqrt(hyp_alpha^2 - hyp_beta^2)*hyp_delta + hyp_beta*(x - hyp_mu) + log(hyp_alpha) + 
            log(hyp_delta) - log(pi()) - 0.5*log(hyp_delta^2 + (x - hyp_mu)^2) + 
            log(modified_bessel_second_kind(1, hyp_alpha*sqrt(hyp_delta^2 + (x - hyp_mu)^2)))); 
  }


  /**
   * Returns the partial sums (the log-likelihood terms from 'start' to 'end' (inclusive) )
   * of the  log-likelihood of the vector X defined by DX=Lambda, where D is a matrix and Lambda
   * is independent NIG noise: log dens_X(x) = log determinant(D) + sum_{i=1}^N log dens_NIG([Dx]_i).
   * This partial sum function is the first argument of the reduce_sum() function which
   * automatically chooses partial sums partitioning based on a dynamic scheduling algorithm,
   * allowing for within-chain parallelization
   *
   * @param slice_Dx The subset of DX for which this partial sum is responsible (slice_DX = DX[start:end])
   * @param start    First term in the partial sum
   * @param end      Last term in the partial sum (inclusive)
   * @param sigma    Standard deviation 
   * @param eta      First flexibility parameter (controls heaviness of tails)
   * @param zeta     Second flexibility parameter (controls asymmetry of tails)
   * @param h        Vector containing the distance between locations or area of the basis functions
   * @return         Returns the partial sum of the log likelihood of X from start to end
   */
  real partial_sum_lpdf(real[] slice_DX, int start, int end, real eta, real zeta, vector h){
    int N = end - start + 1; 
    real llik_partial;
    vector[N] ll;
    for(i in 1:N){
      ll[i] = NIG_var_correction_lpdf( slice_DX[i] | eta, zeta, h[start + i - 1]);}
    llik_partial = sum(ll);
    return llik_partial;
  }

  /**
   * Returns the log-likelihood of X where DX=Lambda, D is a matrix and 
   * Lambda is a vector of independent NIG noise.
   *
   * @param X
   * @param D               Matrix D
   * @param Dv              Array of integer column indices for the non-zero values in D
   * @param Du              Array of integer indices indicating where a given row’s values start
   * @param sizes           Array of integers containing the number of rows, collums and number non-zero elements of D
   * @param sigma           Standard deviation 
   * @param etas            First  flexibility parameter (controls heaviness of tails) 
   * @param mus             Second flexibility parameter (controls asymmetry of tails) 
   * @param h               Vector containing the distance between locations or area of the basis functions
   * @param compute_det     Compute the log determinant of D if D depends on parameter (compute_det=1). Otherwise set compute_det=0.
   * @return                Log-likelihood of X
   */
  real nig_model_lpdf(vector X, matrix D, real etas, real zetas, vector h, int compute_det){

    //convert tail correction parameters (etas, zetas) to original parameterization (eta, zeta) 
    real eta = etas*(1+zetas^2-fabs(zetas)*sqrt(1+zetas^2))^2;
    real zeta  = zetas/sqrt(eta);

    //log-likelihood of the random vector X
    real llik_total = reduce_sum(partial_sum_lpdf, to_array_1d(D*X), 1, eta, zeta, h);

    //If D depends on parameter add log determinant
    if(compute_det==1){
      llik_total  += sum(log(diagonal(cholesky_decompose(D'*D))));    //Only works for spd precision matrices!!! 
    }
     
    return(llik_total);
   }

}

data{
  int<lower=1> N;                      // Number of observations
  int<lower=0> N_covariates;           // Number of covariates
  vector[N] y;                         // Observations
  matrix[N,N] W;                       // Adjacency matrix
  matrix[N,N_covariates] B;            // Design matrix

  real<lower=0> thetaetas;             // Rate constant for PC prior of nus
  real<lower=0> thetazetas;            // Rate constant for PC prior of zetas
}

transformed data{
  vector[N] h = rep_vector(1,N);          // In this problem h is a vector of ones
}

parameters{
  real<lower=0> sigma;           // Standard deviation of NIG noise
  real<lower=0, upper=1> rho;    // Spatial range parameter
  real<lower=0> etas;            // First flexibility parameter (controls heaviness of tails) 
  real zetas;                    // Second flexibility parameter (controls asymmetry of tails) 
  vector[N_covariates] beta;     // Regression coefficients
}

transformed parameters{
  vector[N] x = (y - B*beta)/sigma;                       // Spatial effects
}

model{

  matrix[N,N] D = add_diag(-rho*W, 1);
  x ~ nig_model(D, etas, zetas, h, 1);

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
  //prior for etas
  target += -thetaetas*etas;
  //prior for zetas
  target += -thetazetas*fabs(zetas);

}
