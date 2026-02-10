////////// Declare functions to be used

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

data {
  ///// Data dimension
  // Model
  int<lower=0> N_studies;
  int<lower=0> N_sites;                          // number of sites
  int<lower=0> N_households;                     // number of households (hh)
  int<lower=0> N_recalls;                        // number of recalls (can be more than 1/hh)
  int<lower=0> N_education_levels;                // household data - number of education levels - random factor
  int<lower=0> N_ethnicity;                      // household data - number of ethnicity types - random factor
  int<lower=0> N_years;
 
///// Data
  vector[N_households] frequency;                   // household data - frequency of consumption
  
///// Variables
  int<lower=1,upper=N_studies> study[N_households];        // covariates - site/site identifier
  int<lower=1,upper=N_sites> site[N_households];        // covariates - site/site identifier
  int<lower=1,upper=N_households> hh[N_households];        // covariates - site/site identifier
  vector[N_households] HPD;                      // covariates - forest cover
  vector[N_households] HDI;
  vector[N_households] accessibility;
  vector[N_households] FCI;
  int education[N_households];                    // household data - consumption yes/no
  int<lower=-1> mdays[N_households] ;
  matrix[N_sites,N_sites] distance_matrix;

}

transformed data{
  
  real delta = 1e-9;
 }

parameters {
  // other parameters
  real b0[N_studies];                          // intercept
  real b1;                              // slope - HPD
  real b2;                              // slope - HDI
  real b3;                              // slope - accessibility
  real b4;                              // slope - distance to intact forest
  real b5[N_education_levels];                   // random factor, by education level
  
  //vector[N_households] phi;
  real<lower=0> kappa;                           // overdispersion parameter, beta distribution
  real<lower=0> sigma;                                    // uncertainty if a whole year was recalls
  real<lower=0> xi;                            // gaussian process
  real<lower=0> omicron;                             // gaussian process
  vector [N_sites] omega;                          // gaussian process
  vector[N_households] tau;
}

transformed parameters {
  vector<lower=0>[N_households] sd_frequency;
  vector[N_households] mu_phi;
  vector<lower=0,upper=1>[N_households] phi;
  
  matrix[N_sites, N_sites] Y;
  matrix[N_sites, N_sites] L_Y;
  vector [N_sites] lambda;
  // Gaussian -process - spatial autocorrelation - frequency of consumption
  Y = covariance_GP2(distance_matrix, xi, omicron, delta);
  L_Y = cholesky_decompose(Y);
  lambda = L_Y * omega;
  
  for (i in 1:N_households){
   sd_frequency[i] = sigma * (365-mdays[i]);
  }
  
  for (i in 1:N_households){
    mu_phi[i] = b0[study[i]] + b1 * HPD[i] + b2 * HDI[i] + b3 * accessibility[i] + b4 * FCI[i] + b5[education[i]] + lambda[site[i]];
  }
  
   for (i in 1:N_households){
    phi[i]= inv_logit(mu_phi[i]+sd_frequency[i]*tau[i]);
  }
  
}
  
model{
   
   kappa ~ exponential(1);
   sigma ~ exponential(5);
   b0 ~ normal(0,1.4);
   b1 ~ normal(0,0.5);
   b2 ~ normal(0,0.5);
   b3 ~ normal(0,0.5);
   b4 ~ normal(0,0.5);
   b5 ~ normal(0,1.4);

   xi ~ normal(0,1);
   omicron ~ inv_gamma(5, 5);
   omega ~ normal(0,1);
   
   tau~ normal(0,1);

   for (i in 1:N_households){
     if(frequency[i]>0){
       
       frequency[i] ~ beta(phi[i] * kappa, (1.0 - phi[i]) * kappa);     //observed frequencies on the natural scale are modelled with a Beta distribution- the uncertainty estiamted above is thus propagated in the estimated parameters
       
     }
   }     

}

generated quantities{
  
   /////Define genrated quantities
  real log_lik [N_households];//predicted mean density
  
/////Compute log-likelihood
   for (i in 1:N_households){
        log_lik[i] =  beta_lpdf(frequency[i] | phi[i] * kappa, (1.0 - phi[i]) * kappa);
     }

   
}


