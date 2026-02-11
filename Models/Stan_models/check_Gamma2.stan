///// WILDMEAT CONSUMPTION INTEGRATED ANALYSIS /////
// Mattia Bessone //

////////// Declare functions to be used

functions{
  
  vector merge_missing(int[] miss_indexes, vector x_obs, vector x_miss){
    int N = dims(x_obs)[1];
    int N_miss = dims(x_miss)[1];
    vector[N] merged;
    merged = x_obs;
    for (i in 1:N_miss)
      merged[miss_indexes[i]] = x_miss[i];
    return merged;
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
  vector[N_recalls] quantity;                    // household data - qunatity consumed
  
///// Variables
  int<lower=1,upper=N_studies> study[N_recalls];        // covariates - site/site identifier
  int<lower=1,upper=N_sites> site[N_recalls];        // covariates - site/site identifier
  int<lower=1,upper=N_households> hh[N_recalls];        // covariates - site/site identifier
  vector[N_recalls] HPD;                      // covariates - forest cover
  vector[N_recalls] HDI;
  vector[N_recalls] accessibility;
  vector[N_recalls] FCI;
  int education[N_recalls];                    // household data - consumption yes/no
  vector[N_recalls] AME;
  int AME_missidx[N_missing_AME];
  
}

////////// Declare the parameters to be estiamted
parameters {

  //linear model of consumption quantities, mu
  real c0[N_studies];                          // intercept
  real c1;                                        // spatial autocorrelation
  real c2[N_education_levels];
  real c3[N_households];
  
  real<lower=0> theta;                           // overdispersion parameter, gamma distribution
  real nu;                                       // AME imputation
  vector<lower=0> [N_missing_AME] AME_imputed;            // AME imputation
  real<lower=0> sigma_AME;                       // AME imputation
  
}

////////// Declare the linear models
transformed parameters{
  vector<lower=0>[N_recalls] mu;                          // mean quantity consumed
  vector<lower=0>[N_recalls] AME_merged;
  
  AME_merged = merge_missing(AME_missidx,to_vector(AME),AME_imputed);
  
  for (i in 1: N_recalls){                       // parameters are modelled with covariates
    
    mu[i] = exp(c0[study[i]] + c1 * AME_merged[i] + c2[education[i]] + c3[hh[i]]);
  
  }
   
}
////////// Declare the model
model {
 
 // Priors
 c0 ~ normal(0,5);
 c1 ~ normal(0,0.5);
 c2 ~ normal(0,5);
 c3 ~ normal(0,5);
 
 theta ~ exponential(1);
 nu ~ normal(6,1.5);
 sigma_AME ~ exponential(1);
 
 AME_merged ~ normal(nu, sigma_AME);
 
  
 for (i in 1: N_recalls){                                            // yes/no data are modelled with a Bernoulli distribution
 
  if (quantity[i] > 0){ // if quantity data are available
      
      quantity[i] / (AME_merged[i]) ~ gamma(mu[i] * theta, theta);   //observed quantities consumed are modelled with a Gamma distribution
      
    }
 }
  
}
////////// Make predictions and calculations
generated quantities{
  
  /////Define genrated quantities
 // real grams_consumed_per_AME [N_recalls];
  real log_lik [N_recalls];//predicted mean density


/////Compute log-likelihood
   for (i in 1:N_recalls){
     log_lik[i] =  gamma_lpdf(quantity[i] / (AME_merged[i]) | mu[i]*theta, theta); //occupancy rate is only conditional on the probability of detection p
    }
   
}


