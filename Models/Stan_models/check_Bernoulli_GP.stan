functions{
  
  matrix covariance_GP2(matrix  x, real zeta, real rho, real delta){
    int N = dims(x)[1];
    matrix[N,N] K;
    for (i in 1:(N-1)){
      K[i,i] = zeta + delta;
      for (j in (i+1):N){
        K[i,j] = zeta * exp(-rho * x[i,j]^2);
        K[j,i] = K[i,j];
      }
    }
    K[N,N] = zeta + delta;
    return K;
  }
}

////////// Declare the data
data {
///// Data dimension
  // Model
  int<lower=0> N_studies;
  int<lower=0> N_sites;                          // number of sites
  int<lower=0> N_households;                     // number of households (hh)
  int<lower=0> N_recalls;                        // number of recalls (can be more than 1/hh)
  int<lower=0> N_education_levels;                // household data - number of education levels - random factor
  // Imputation
  int<lower=0> N_missing_AME;
  
///// Data
  int<lower=0,upper=1> consumption[N_recalls];                    // household data - consumption yes/no
  
///// Variables
  int<lower=1,upper=N_studies> study[N_recalls];        // covariates - site/site identifier
  int<lower=1,upper=N_sites> site[N_recalls];        // covariates - site/site identifier
  int<lower=1,upper=N_households> hh[N_recalls];        // covariates - site/site identifier
  vector[N_recalls] HPD;                      // covariates - forest cover
  vector[N_recalls] HDI;
  vector[N_recalls] accessibility;
  vector[N_recalls] FCI;
  int<lower=0,upper=N_education_levels> education[N_recalls];                    // household data - consumption yes/no
  vector[N_recalls] AME;
  int AME_missidx[N_missing_AME];
  int<lower=-1> days[N_recalls] ;              // sampling effort in days
  matrix[N_sites,N_sites] distance_matrix;
}

transformed data {
  real delta = 1e-9;
}

////////// Declare the parameters to be estiamted
parameters {
  //linear model of probability of consumption , pi
  real a0[N_studies];                          // intercept
  real a1;                              // slope - HPD
  real a2;                              // slope - HDI
  real a3;                              // slope - accessibility
  real a4;                              // slope - distance to intact forest
  real a5[N_education_levels];                   // random factor, by education level
  real<lower=0> a7;
  
  real<lower=0> zeta;                            // gaussian process
  real<lower=0> rho;                             // gaussian process
  vector [N_sites] eta;                          // gaussian process
}

////////// Declare the linear models
transformed parameters{
  vector<lower=0,upper=1>[N_recalls] p;                           // mean propability of consumtpion
  matrix[N_sites, N_sites] K;
  matrix[N_sites, N_sites] L_K;
  vector [N_sites] f;
  
  // Gaussian -process - spatial autocorrelation
  K = covariance_GP2(distance_matrix, zeta, rho, delta);
  L_K = cholesky_decompose(K);
  f = L_K * eta;

  for (i in 1: N_recalls){                       // parameters are modelled with covariates
    
    p[i] = inv_logit(a0[study[i]] + a1 * HPD[i] +  a2 * HDI[i] + a3 * accessibility[i] + a4 * FCI[i]  + a5[education[i]] + a7 * days[i] + f[site[i]]);
    
  }
   
}
////////// Declare the model
model {
 
 // Priors
 a0 ~ normal(0,1.4);
 a1 ~ normal(0,0.5);
 a2 ~ normal(0,0.5);
 a3 ~ normal(0,0.5);
 a4 ~ normal(0,0.5);
 a5 ~ normal(0,1.4);
 a7 ~ exponential(1);
 
 
 zeta ~ normal(0,1);
 rho ~ inv_gamma(5, 5);
 eta ~ normal(0,1);
 
 for (i in 1: N_recalls){                                            // yes/no data are modelled with a Bernoulli distribution
   
  consumption[i] ~ bernoulli(p[i]);
  
 }
  
}
////////// Make predictions and calculations
generated quantities{
  
  /////Define genrated quantities
  real log_lik [N_recalls];//predicted mean density
  
/////Compute log-likelihood
   for (i in 1:N_recalls){
     
     log_lik[i] =  bernoulli_lpmf(consumption[i] | p[i]); //occupancy rate is only conditional on the probability of detection p
      
    }

   
}

