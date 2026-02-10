///// WILDMEAT CONSUMPTION INTEGRATED ANALYSIS

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

// Declare the data

data {
/////
  int<lower=0> N_recalls;                        // number of households (hh) sampled: at least for consumption = yes/no
  int<lower=0> N_region;                         // number of regions/sites
  int<lower=0> N_sites;
  int<lower=0> N_ethnicity;                      // household data - number of ethnicity types - random factor
  int<lower=0> N_missing_AME;
///// Data
  int consumption[N_recalls];                    // household data - consumption yes/no
  vector[N_recalls] frequency;                   // household data - frequency of consumption
  vector[N_recalls] quantity;                    // household data - qunatity consumed
  
///// Variables
  int<lower=1,upper=3> region[N_recalls];        // covariates - region/site identifier
  vector[N_recalls] forest;                      // covariates - forest cover
  int<lower=1,upper=2> ethnicity[N_recalls];     // covariates - ethnicity
  vector[N_recalls] accessibility;                        // covariates - mean body mass index of consumed species
  vector<lower=-1>[N_recalls] days;              // sampling effort in days
  int AME[N_recalls];                // number of people in each sampled household - e.g. adult male equivalent
  int AME_missidx[N_missing_AME];
  matrix[N_sites,N_sites] distance_matrix;
  int<lower=1> site [N_recalls];
///// Prediction
  int sampled_hh [N_region];                     // number of sampled households
  int pred_region[N_sites];                      // region identifier of sites/cells for prediction
  int<lower=0> population[N_sites];             // population in each region
  real pred_forest [N_sites];                    // forest cover in each site/cell
  real pred_accessibility[N_sites];

}

transformed data{
  real delta = 1e-9;
}


parameters {
  //linear model of probability of consumption , pi
  real a1[N_region];                             //intercept
  real a2[N_region];                             //slope of forest cover, by region
  real a3[N_ethnicity];                          //random factor, by ethnicity
  
  //linear model of consumption frequency, phi
  real b1[N_region];                             //intercept
  real b2[N_region];                             //slope of accessibility, by region
  real b3[N_ethnicity];                          //random factor, by ethnicity
 
  //linear model of consumption quantities, mu
  real c1[N_region];                             //intercept
  real c2[N_region];                             //slope of accessibility, by region
  real c3;
  
  // other parameters
  real<lower=0> kappa;                           // overdispersion parameter, beta distribution
  real<lower=0> sigma;                                    // uncertainty if a whole year was sampled
  real<lower=0> theta;                           // overdispersion parameter, gamma distribution
  vector[N_recalls] tau;                      
  
  // imputation
  real nu;                                       // AME imputation
  vector<lower=0> [N_missing_AME] AME_imputed;            // AME imputation
  real<lower=0> sigma_AME;                       // AME imputation
  
  //// Gaussian process
  // Consumption probability
  real<lower=0> zeta;                            
  real<lower=0> rho;
  vector [N_sites] eta;
  // Frequency of consumtpion
  real<lower=0> xi;
  real<lower=0> omicron;
  vector [N_sites] omega;
  
}


transformed parameters{
  vector<lower=0,upper=1>[N_recalls] p;                           // mean propability of consumtpion
  vector<lower=0>[N_recalls] sd_frequency;       // uncertainty of observed frequency of consumption
  vector[N_recalls] mu_phi;                      // mean observed frequency of consumption on the logit scale
  vector<lower=0,upper=1>[N_recalls] phi;                      // mean observed frequency of consumption on the logit scale
  vector<lower=0>[N_recalls] mu;                          // mean quantity consumed
  vector<lower=0>[N_recalls] AME_merged;
  matrix[N_sites, N_sites] X;
  matrix[N_sites, N_sites] L_X;
  vector [N_sites] epsilon;
  matrix[N_sites, N_sites] Y;
  matrix[N_sites, N_sites] L_Y;
  vector [N_sites] lambda;
  
  // Gaussian -process - spatial autocorrelation - consumption probability
  X = covariance_GP2(distance_matrix, zeta, rho, delta);
  L_X = cholesky_decompose(X);
  epsilon = L_X * eta;
  
  // Gaussian -process - spatial autocorrelation - frequency of consumption
  Y = covariance_GP2(distance_matrix, xi, omicron, delta);
  L_Y = cholesky_decompose(Y);
  lambda = L_Y * omega;
  
  // AME Imputation
  AME_merged = merge_missing(AME_missidx,to_vector(AME),AME_imputed);
 
  
    
  for (i in 1: N_recalls){                       // parameters are modelled with covariates
    
    p[i] = inv_logit(a1[region[i]] + a2[region[i]] * forest[i] + a3[ethnicity[i]] + epsilon[site[i]]);
    
    sd_frequency[i] = sigma * (365 - days[i]); // uncertainty of observed frequency of consumption depends on sampling effort - the longer the less uncertain
    mu_phi[i] = b1[region[i]] + b2[region[i]] * accessibility[i] + b3[ethnicity[i]] + lambda[site[i]];
    
    mu[i] = exp(c1[region[i]] + c2[region[i]] * accessibility[i] + c3 * AME_merged[i]);
  
  }
  
  for (i in 1: N_recalls){ 
  phi[i] = inv_logit(mu_phi[i]+sd_frequency[i]*tau[i]);
  }

}

model {
 // Priors
 a1 ~ normal(0,1.4);
 a2 ~ normal(0,0.5);
 a3 ~ normal(0,1.4);
 b1 ~ normal(0,1.4);
 b2 ~ normal(0,0.5);
 b3 ~ normal(0,1.4);
 c1 ~ normal(0,10);
 c2 ~ normal(0,0.5);
 c3 ~ normal(0,0.5);
 
 tau ~ std_normal();
 kappa ~ exponential(1);
 sigma ~ exponential(10);
 theta ~ exponential(1);
 
 nu ~ normal(0,2);
 sigma_AME ~ exponential(1);
 
 zeta ~ normal(0,1);
 rho ~ inv_gamma(5, 5);
 eta ~ normal(0,1);
 xi ~ normal(0,1);
 omicron ~ inv_gamma(5, 5);
 omega ~ normal(0,1);
 
 AME_merged ~ normal(nu, sigma_AME);
 
 //phi ~ normal(mu_phi,sd_frequency);
 
 for (i in 1: N_recalls){                                            // yes/no data are modelled with a Bernoulli distribution
   
  consumption[i] ~ bernoulli(p[i]);
 }
 
 for (i in 1:N_recalls){
    if (frequency[i] > 0 ){                                            // if frequency data are available
     frequency[i] ~ beta(phi[i] * kappa, (1.0 - phi[i]) * kappa);     //observed frequencies on the natural scale are modelled with a Beta distribution- the uncertainty estiamted above is thus propagated in the estimated parameters
    }
 }
  
 for (i in 1:N_recalls){
    
  if (quantity[i] > 0){                                            // if quantity data are available
    
    quantity[i] / AME_merged[i] ~ gamma(mu[i] * theta, theta);     //observed quantities consumed are modelled with a Gamma distribution
    
  } 
 }
  
}

generated quantities{
  
  real pred_p [N_sites];                       // compute mean probability of consumption
  int  pred_consume_occurrence [N_sites];        // predicted probability consumption in each site
  real pred_phi [N_sites];
  real pred_frequency [N_sites];
  real mean_quantity;                 // compute mean quantity consumed per day for each region
  
  real grams_consumed[N_sites];           // predicted probability consumption in each site
  real proportion_sites_consuming;                    // predicted bushmeat consumed in each site
  real mean_frequency_consumption;
  real total_grams_consumed;
 
  // Estiamate p, phi and mu for each prediction grid cell
  for (i in 1:N_sites){
    
    pred_p[i] = inv_logit(a1[pred_region[i]] + a2[pred_region[i]] * pred_forest[i] + mean(a3));
    pred_phi[i] = inv_logit(b1[pred_region[i]] + b2[pred_region[i]] * pred_accessibility[i] + mean(b3));

  }
  
  for (i in 1:N_sites){
      
     pred_consume_occurrence[i] = bernoulli_rng(pred_p[i]);
     pred_frequency[i] = beta_rng(pred_phi[i] * kappa, (1.0 - pred_phi[i]) * kappa);
  }
  
  mean_quantity = mean(mu);
   
  // Calculate amount consumed in each site, for each year
  for (i in 1:N_sites){
    
    grams_consumed[i] =  pred_consume_occurrence[i] * pred_frequency[i] * mean_quantity * population[i];
}
  
   // Calculate number of sites where wildmeat is consumed & total amount consumed in the study area
    proportion_sites_consuming = sum(pred_consume_occurrence) * 1.0 / N_sites;
    mean_frequency_consumption = mean(pred_frequency) * mean(pred_p);
    total_grams_consumed = sum(grams_consumed);

}

