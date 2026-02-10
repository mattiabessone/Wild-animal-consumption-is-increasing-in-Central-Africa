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
  int<lower=0> N_studies;                                            // number of studies
  int<lower=0> N_sites;                                              // number of sites
  int<lower=0> N_households;                                         // number of households (hh)
  int<lower=0> N_recalls;                                            // number of recalls (can be more than 1/hh)
  int<lower=0> N_education_levels;                                   // number of education levels - random factor
  int<lower=0> N_survey_types;                                       // number of survey types - random factor
  int<lower=0> N_location_types;                                     // number of location types - random factor
  int<lower=0> N_periods;                                            // number of periods
   // Imputation
  int<lower=0> N_missing_AME;                                        // number of data ppoints with missing AME information
  // Prediction
  int<lower=0> N_scenarios;                                          // number of scenarios used for prediction
  int<lower=0> N_sites_pred;                                         // number of cells of prediction grid
  
///// Data
  int consumption[N_recalls];                                        // household data - consumption yes/no
  vector[N_households] frequency;                                    // household data - frequency of consumption
  vector[N_recalls] quantity;                                        // household data - qunatity consumed
  
///// Variables
  int<lower=1,upper=N_studies> study[N_recalls];                     // random factor - study identifier
  int<lower=1,upper=N_sites> site[N_recalls];                        // random factor - site identifier
  int<lower=1,upper=N_households> household[N_recalls];              // random factor - household identifier
  vector[N_recalls] HPD;                                             // covariates - human population density
  vector[N_recalls] HDI;                                             // covariates - human development index
  vector[N_recalls] REM;                                             // covariates - remoteness
  vector[N_recalls] FCI;                                             // covariates - forest condition index
  int education[N_recalls];                                          // factor - education level
  int location_type[N_recalls];                                      // factor - location type
  int survey_type[N_recalls];                                        // factor - survey type
  int period[N_recalls];                                             // factor - period
  int<lower=-1> days[N_recalls] ;                                    // recall length in days
  //// Variables aggregated by household for frequency model  
  int<lower=1,upper=N_studies> study_fr[N_households];               // random factor - study identifier
  int<lower=1,upper=N_sites> site_fr[N_households];                  // random factor - site identifier
  vector[N_households] HPD_fr;                                       // covariates - human population density
  vector[N_households] HDI_fr;                                       // covariates - human developmnet index
  vector[N_households] REM_fr;                                       // covariates - remoteness
  vector[N_households] FCI_fr;                                       // covariates - forest condition index
  vector[N_recalls] AME;                                             // covarites - AME
  int education_fr[N_households];                                    // factor - education level
  int location_type_fr[N_households];                                // factor - location type
  int period_fr[N_households];                                       // factor - period
  int<lower=-1> mdays[N_households] ;                                // recall length in days
  //// Spatial autocorrelation & AME imputation  
  int AME_missidx[N_missing_AME];                                    // position of missing AME value in the AME vector
  real mean_AME;                                                     // mean of available AME value - used to centre prior
  matrix[N_sites,N_sites] distance_matrix;                           // distance matrix between sites

///// Prediction
  matrix[N_sites_pred,N_scenarios] AME_pred;                        // AME value of each cell, based on population in prediction cell each cell
  matrix[N_sites_pred,N_scenarios] HPD_pred;                        // covariates - for each scenario, human population density in each cell
  matrix[N_sites_pred,N_scenarios] HDI_pred;                        // covariates - for each scenario, human development index in each cell
  vector[N_sites_pred] REM_pred;                                    // covariates - remoteness in each cell (fixed)
  vector[N_sites_pred] FCI_pred;                                    // covariates - forest condition index in each cell (fixed) 
  matrix[N_education_levels-1,N_scenarios] ED_pred [N_sites_pred];  // factor - for each scenario, proportion of people with education > primary AND  proportion of people with education < secondary, in each cell
  matrix[N_location_types,N_scenarios] LT_pred [N_sites_pred];      // factor - for each scenario, proportion of people lving in villages, towns and cities in each cell
  
}

transformed data {
  real delta = 1e-9;                                                 // gaussian process
}

////////// Declare the parameters to be estiamted
parameters {
  //linear model of probability of consumption , pi
  real a0[N_studies];                   // intercept, varying by study
  real a1[N_location_types];            // slope - HPD, interaction with settlement type
  real a2;                              // slope - HDI
  real a3;                              // slope - REM
  real a4;                              // slope - FCI
  real a5[N_education_levels];          // factor - education
  real a6[N_households];                // random factor - household
  real<lower=0> a7;                     // parameter defining the increment in π as a function of recall duration 
  
  //linear model of consumption frequency, phi
  real b0[N_studies];                   // intercept, varying by study
  real b1[N_location_types];            // slope - HPD, interaction with settlement type
  real b2;                              // slope - HDI
  real b3;                              // slope - REM
  real b4;                              // slope - FCI
  real b5[N_education_levels];          // factor - education
  
  //linear model of consumption quantities, mu
  real c0[N_studies];                          // intercept, varying by study
  real c1;                                     // slope - AME
  real c2[N_education_levels];                 // factor - education
  real c3[N_location_types];                   // factor - location type
  real c4[N_survey_types];                     // factor - survey type
  real c5[N_households];                       // random factor - household
  
  // other parameters
  real<lower=0> kappa;                           // overdispersion parameter, beta distribution
  real<lower=0> sigma;                           // uncertainty if a whole period was recalls
  vector [N_households] tau;                     // 
  
  real<lower=0> theta;                           // overdispersion parameter, gamma distribution
  
  real<lower=0> zeta;                            // gaussian process
  real<lower=0> rho;                             // gaussian process
  vector [N_sites] eta;                          // gaussian process
   
  real<lower=0> xi;                              // gaussian process
  real<lower=0> omicron;                         // gaussian process
  vector [N_sites] omega;                        // gaussian process
  
  real nu;                                       // AME imputation
  vector<lower=0> [N_missing_AME] AME_imputed;   // AME imputation
  real<lower=0> psi;                             // AME imputation

}

////////// Declare the model
model {
 vector [N_sites] epsilon;
 real p;
 vector [N_sites] lambda;
 real mu_phi;
 real sd_frequency;
 vector [N_recalls] AME_merged;
 real mu;
 
   {
    matrix[N_sites, N_sites] L_X;
    matrix[N_sites, N_sites] X = covariance_GP2(distance_matrix, zeta, rho, delta);
    
    matrix[N_sites, N_sites] L_Y;
    matrix[N_sites, N_sites] Y = covariance_GP2(distance_matrix, xi, omicron, delta);

    L_X = cholesky_decompose(X);
    epsilon = L_X * eta;
    
    L_Y = cholesky_decompose(Y);
    lambda = L_Y * omega;
    
    AME_merged = merge_missing(AME_missidx,to_vector(AME),AME_imputed);

  }
 // Priors
 a0 ~ normal(0,1.4);
 a1 ~ normal(0,0.5);
 a2 ~ normal(0,0.5);
 a3 ~ normal(0,0.5);
 a4 ~ normal(0,0.5);
 a5 ~ normal(0,1.4);
 a6 ~ normal(0,1.4);
 a7 ~ exponential(1);
 
 b0 ~ normal(0,1.4);
 b1 ~ normal(0,0.5);
 b2 ~ normal(0,0.5);
 b3 ~ normal(0,0.5);
 b4 ~ normal(0,0.5);
 b5 ~ normal(0,1.4);
 
 c0 ~ normal(0,1);
 c1 ~ normal(0,0.5);
 c2 ~ normal(0,1);
 c3 ~ normal(0,1);
 c4 ~ normal(0,1);
 c5 ~ normal(0,1);

 kappa ~ exponential(0.1);
 sigma ~ exponential(10);
 tau ~ normal(0,0.5);
 theta ~ exponential(1);
 
 zeta ~ normal(0,1);
 rho ~ inv_gamma(5, 5);
 eta ~ normal(0,1);
 xi ~ normal(0,1);
 omicron ~ inv_gamma(5, 5);
 omega ~ normal(0,1);
 
 nu ~ normal(mean_AME,1);
 psi ~ exponential(1);
 AME_merged ~ normal(nu, psi);
 
 for (r in 1:N_recalls){                                             // yes/no data are modelled with a Bernoulli distribution
  {
    p = inv_logit(a0[study[r]] + a1[location_type[r]] * HPD[r] + a2 * REM[r] + a3 * HDI[r] + a4 * FCI[r] + a5[education[r]] + a6[household[r]] + a7 * days[r] + epsilon[site[r]]);
    consumption[r] ~ bernoulli(p);
  }
 }
 
 for (h in 1:N_households){
   {
     sd_frequency = sigma * (365 - mdays[h]);
     mu_phi = b0[study_fr[h]] + b1[location_type_fr[h]] * HPD_fr[h] + b2 * REM_fr[h] + b3 * HDI_fr[h] + b4 * FCI_fr[h] + b5[education_fr[h]] + lambda[site_fr[h]];
     
   if (frequency[h] > 0 ){                                           // if frequency data are available
    frequency[h] ~ beta(inv_logit(mu_phi + sd_frequency * tau[h]) * kappa,
    (1.0 - inv_logit(mu_phi + sd_frequency * tau[h])) * kappa);      //observed frequencies on the natural scale are modelled with a Beta distribution- the uncertainty estimated above is thus propagated in the estimated parameters
   }
   }
 }
 
 for (r in 1:N_recalls){                                             // yes/no data are modelled with a Bernoulli distribution
 
   {
     mu = exp(c0[study[r]] + c1 * AME_merged[r] + c2[education[r]] + c3[location_type[r]] + c4[survey_type[r]] + c5[household[r]]);
  
   if (quantity[r] > 0 ){                                            // if frequency data are available
      (quantity[r] / AME_merged[r]) ~ gamma(mu * theta, theta);      //observed quantities consumed are modelled with a Gamma distribution
   }
   }
 }
  
}
////////// Make predictions and calculations
generated quantities{
  
  // predicted factorial coefficients (education level and human population density)
  
  real pred_a1 [N_sites_pred,N_scenarios];
  real pred_b1 [N_sites_pred,N_scenarios];
  real pred_a5 [N_sites_pred,N_scenarios];
  real pred_b5 [N_sites_pred,N_scenarios];
  real pred_c2 [N_sites_pred,N_scenarios];
  real pred_c3 [N_sites_pred,N_scenarios];
 
  real pred_consumption_probability [N_sites_pred,N_scenarios];      // predicted probability consumption in each cell, for each scenario
  real pred_frequency [N_sites_pred,N_scenarios];                    // predicted frequency of consumption in each cell, for each scenario
  real pred_quantity [N_sites_pred,N_scenarios];                     // predicted qunatity consumed in each cell, for each scenario
  
  real consumption_rates [N_sites_pred,N_scenarios];                 // predicted cosumption rates in each cell, for each scenario
  real tonnes_consumed[N_sites_pred,N_scenarios];                    // predicted tonnes consumed in each cell, for each scenario
  real total_tonnes_consumed[N_scenarios];                           // total predicted bushmeat consumed in tonnes
  
  
  for (j in 1:N_sites_pred){
    
    pred_a1[j,1] = a1[1] * LT_pred[j,1,1] + a1[2] * LT_pred[j,2,1] + a1[3] * LT_pred[j,3,1];
    pred_a1[j,2] = a1[1] * LT_pred[j,1,2] + a1[2] * LT_pred[j,2,2] + a1[3] * LT_pred[j,3,2];
    pred_a1[j,3] = a1[1] * LT_pred[j,1,3] + a1[2] * LT_pred[j,2,3] + a1[3] * LT_pred[j,3,3];
    
    pred_b1[j,1] = b1[1] * LT_pred[j,1,1] + b1[2] * LT_pred[j,2,1] + b1[3] * LT_pred[j,3,1];
    pred_b1[j,2] = b1[1] * LT_pred[j,1,2] + b1[2] * LT_pred[j,2,2] + b1[3] * LT_pred[j,3,2];
    pred_b1[j,3] = b1[1] * LT_pred[j,1,3] + b1[2] * LT_pred[j,2,3] + b1[3] * LT_pred[j,3,3];
    
    pred_c3[j,1] = c3[1] * LT_pred[j,1,1] + c3[2] * LT_pred[j,2,1] + c3[3] * LT_pred[j,3,1];
    pred_c3[j,2] = c3[1] * LT_pred[j,1,2] + c3[2] * LT_pred[j,2,2] + c3[3] * LT_pred[j,3,2];
    pred_c3[j,3] = c3[1] * LT_pred[j,1,3] + c3[2] * LT_pred[j,2,3] + c3[3] * LT_pred[j,3,3];
    
    pred_a5[j,1] = a5[1] * (1 - ED_pred[j,1,1]) + a5[2] * (1 - ED_pred[j,2,1]) + a5[3] * (1 - ED_pred[j,3,1]) + a5[4] * ED_pred[j,4,1] + a5[5] * ED_pred[j,5,1] + a5[6] * ED_pred[j,6,1];
    pred_a5[j,2] = a5[1] * (1 - ED_pred[j,1,2]) + a5[2] * (1 - ED_pred[j,2,2]) + a5[3] * (1 - ED_pred[j,3,2]) + a5[4] * ED_pred[j,4,2] + a5[5] * ED_pred[j,5,2] + a5[6] * ED_pred[j,6,2];
    pred_a5[j,3] = a5[1] * (1 - ED_pred[j,1,3]) + a5[2] * (1 - ED_pred[j,2,3]) + a5[3] * (1 - ED_pred[j,3,3]) + a5[4] * ED_pred[j,4,3] + a5[5] * ED_pred[j,5,3] + a5[6] * ED_pred[j,6,3];
    
    pred_b5[j,1] = b5[1] * (1 - ED_pred[j,1,1]) + b5[2] * (1 - ED_pred[j,2,1]) + b5[3] * (1 - ED_pred[j,3,1]) + b5[4] * ED_pred[j,4,1] + b5[5] * ED_pred[j,5,1] + b5[6] * ED_pred[j,6,1];
    pred_b5[j,2] = b5[1] * (1 - ED_pred[j,1,2]) + b5[2] * (1 - ED_pred[j,2,2]) + b5[3] * (1 - ED_pred[j,3,2]) + b5[4] * ED_pred[j,4,2] + b5[5] * ED_pred[j,5,2] + b5[6] * ED_pred[j,6,2];
    pred_b5[j,3] = b5[1] * (1 - ED_pred[j,1,3]) + b5[2] * (1 - ED_pred[j,2,3]) + b5[3] * (1 - ED_pred[j,3,3]) + b5[4] * ED_pred[j,4,3] + b5[5] * ED_pred[j,5,3] + b5[6] * ED_pred[j,6,3];
    
    pred_c2[j,1] = c2[1] * (1 - ED_pred[j,1,1]) + c2[2] * (1 - ED_pred[j,2,1]) + c2[3] * (1 - ED_pred[j,3,1]) + c2[4] * ED_pred[j,4,1] + c2[5] * ED_pred[j,5,1] + c2[6] * ED_pred[j,6,1];
    pred_c2[j,2] = c2[1] * (1 - ED_pred[j,1,2]) + c2[2] * (1 - ED_pred[j,2,2]) + c2[3] * (1 - ED_pred[j,3,2]) + c2[4] * ED_pred[j,4,2] + c2[5] * ED_pred[j,5,2] + c2[6] * ED_pred[j,6,2];
    pred_c2[j,3] = c2[1] * (1 - ED_pred[j,1,3]) + c2[2] * (1 - ED_pred[j,2,3]) + c2[3] * (1 - ED_pred[j,3,3]) + c2[4] * ED_pred[j,4,3] + c2[5] * ED_pred[j,5,3] + c2[6] * ED_pred[j,6,3];
  
       
  }
  
  // Estimate p, phi and mu for each prediction grid cell
  for (j in 1:N_sites_pred){
    
    pred_consumption_probability[j,1] = inv_logit(mean(a0) + pred_a1[j,1] * HPD_pred[j,1] + a2 * REM_pred[j] + a3 * HDI_pred[j,1] + a4 * FCI_pred[j] + pred_a5[j,1] + a7);
    pred_consumption_probability[j,2] = inv_logit(mean(a0) + pred_a1[j,2] * HPD_pred[j,2] + a2 * REM_pred[j] + a3 * HDI_pred[j,2] + a4 * FCI_pred[j] + pred_a5[j,2] + a7);
    pred_consumption_probability[j,3] = inv_logit(mean(a0) + pred_a1[j,3] * HPD_pred[j,3] + a2 * REM_pred[j] + a3 * HDI_pred[j,3] + a4 * FCI_pred[j] + pred_a5[j,3] + a7);

  
    pred_frequency[j,1] = beta_rng(inv_logit(mean(b0) + pred_b1[j,1] * HPD_pred[j,1] + b2 * REM_pred[j] + b3 * HDI_pred[j,1] + b4 * FCI_pred[j] + pred_b5[j,1]) * kappa, 
    (1.0 - inv_logit(mean(b0) + pred_b1[j,1] * HPD_pred[j,1] + b2 * REM_pred[j] + b3 * HDI_pred[j,1] + b4 * FCI_pred[j] + pred_b5[j,1])) * kappa);
    
    pred_frequency[j,2] = beta_rng(inv_logit(mean(b0) + pred_b1[j,2] * HPD_pred[j,2] + b2 * REM_pred[j] + b3 * HDI_pred[j,2] + b4 * FCI_pred[j] + pred_b5[j,2]) * kappa, 
    (1.0 - inv_logit(mean(b0) + pred_b1[j,2] * HPD_pred[j,2] + b2 * REM_pred[j] + b3 * HDI_pred[j,2] + b4 * FCI_pred[j] + pred_b5[j,2])) * kappa);
    
    pred_frequency[j,3] = beta_rng(inv_logit(mean(b0) + pred_b1[j,3] * HPD_pred[j,3] + b2 * REM_pred[j] + b3 * HDI_pred[j,3] + b4 * FCI_pred[j] + pred_b5[j,3]) * kappa, 
    (1.0 - inv_logit(mean(b0) + pred_b1[j,3] * HPD_pred[j,3] + b2 * REM_pred[j] + b3 * HDI_pred[j,3] + b4 * FCI_pred[j] + pred_b5[j,3])) * kappa);
    
    pred_quantity[j,1] = gamma_rng(exp(mean(c0) + c1 * mean_AME + pred_c2[j,1] + pred_c3[j,1] + mean(c4)) * theta, theta);
    pred_quantity[j,2] = gamma_rng(exp(mean(c0) + c1 * mean_AME + pred_c2[j,2] + pred_c3[j,2] + mean(c4)) * theta, theta);
    pred_quantity[j,3] = gamma_rng(exp(mean(c0) + c1 * mean_AME + pred_c2[j,3] + pred_c3[j,3] + mean(c4)) * theta, theta);
  
  }
  
  // Calculate amount consumed in each period
  for (j in 1:N_sites_pred){
    for (z in 1:N_scenarios){
      
      consumption_rates[j,z] =  pred_consumption_probability[j,z] * pred_frequency[j,z] * pred_quantity[j,z];
      tonnes_consumed[j,z] =  (consumption_rates[j,z] * AME_pred [j,z] * 365) / 1000;
    
    }
  }
  
   // Calculate number of sites where wildmeat is consumed & total amount consumed in the study area
    for (z in 1:N_scenarios){
      
      total_tonnes_consumed[z] = sum(tonnes_consumed[1:N_sites_pred,z])/1000000; // expressed in million tonnes
    
    }

}

