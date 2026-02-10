### WILDMEAT CONSUMPTION SIMULATION STUDY ###
### Mattia Bessone, PhD ###
#==========================================================================================================#
# We evaluate the model's performances retrieving the parameters defined in
# the simulation under 3 scenarios of sampling coverage ##

# Load packages
library (rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

n_rep<-100
n_iter<-1000
n_draws<-n_iter/2
rows <- 30                   # Set side of square region, number of cells = 900
                             # our study area = 1121
# State parameters

a1 <- c(0.8,-1,-3.5)         # Probability of consumption - intercept, 3 levels
a2 <- c(0.4,0.1,0.3)         # Probability of consumption - slope, 3 levels
a3 <- c(2,0.3)               # Probability of consumption - random factor, 2 levels
zeta <- 0.8                  # Gaussian process, i.e. spatial autocorrelation
rho<-2                       # Gaussian process, i.e. spatial autocorrelation               

b1<-c(0.9,0.5,-0.5)          # Frequency of consumption - intercept, 3 levels
b2<-c(0.15,0.08,0.2)         # Frequency of consumption - slope, 3 levels  
b3<-c(1,-0.9)                # Frequency of consumption - random factor, 2 levels
kappa <- 15                  # Frequency of consumption - overdispersion parameter
sigma <- 0.03                # Standard deviation of observed frequency of consumption if a full year is monitored
xi<-0.3                      # Gaussian process, i.e. spatial autocorrelation
omicron<-1.2                 # Gaussian process, i.e. spatial autocorrelation

c1<-c(5.8,5.5,5)             # Quantity consumed - intercept, 3 levels
c2<-c(0.01,-0.02,0.06)       # Quantity consumed - slope, 3 levels
c3<- -0.07                   # Quantity consumed - slope AME - 1 levels
theta<-10                    # Quantity consumed - over-dispersion parameter
prob_missing_AME<-0.2        # Probability that the number of AME was not recorded

##### SCENARIO 1 #####
# Define coverage first scenario
coverage <- 0.05  # Set percentage of cells surveyed - this will vary
# Simulate datasets and run models
source("Code/run_simulation_coverage_1.R")
# Store results
source("Code/save_objects_coverage_1.R")
# Remove objects to save memory
rm(est_tot_consumption_1,est_mean_frequency_1,est_mean_f1_1,est_mean_f2_1,
   est_mean_f3_1,est_mean_consprob_1,est_a2_1_1,est_a2_2_1,est_a2_3_1,est_a3_1_1,
   est_a3_2_1,est_b2_1_1,est_b2_2_1,est_b2_3_1,est_b3_1_1,est_b3_2_1,est_c2_1_1,
   est_c2_2_1,est_c2_3_1,sim_tot_consumption_1,sim_mean_frequency_1,sim_mean_f1_1,
   sim_mean_f2_1,sim_mean_f3_1,sim_mean_consprob_1)

##### SCENARIO 2 #####
# Define coverage second scenario
coverage <- 0.10  # Set percentage of cells surveyed - this will vary
# Simulate datasets and run models
source("Code/run_simulation_coverage_2.R")
# Store results
source("Code/save_objects_coverage_2.R")
# Remove objects to save memory
rm(est_tot_consumption_2,est_mean_frequency_2,est_mean_f1_2,est_mean_f2_2,
   est_mean_f3_2,est_mean_consprob_2,est_a2_1_2,est_a2_2_2,est_a2_3_2,est_a3_1_2,
   est_a3_2_2,est_b2_1_2,est_b2_2_2,est_b2_3_2,est_b3_1_2,est_b3_2_2,est_c2_1_2,
   est_c2_2_2,est_c2_3_2,sim_tot_consumption_2,sim_mean_frequency_2,sim_mean_f1_2,
   sim_mean_f2_2,sim_mean_f3_2,sim_mean_consprob_2)

##### SCENARIO 3 #####
# Define coverage third scenario
coverage <- 0.15  # Set percentage of cells surveyed - this will vary
# Simulate datasets and run models
source("Code/run_simulation_coverage_3.R")
# Store results
source("Code/save_objects_coverage_3.R")
# Remove objects to save memory
rm(est_tot_consumption_3,est_mean_frequency_3,est_mean_f1_3,est_mean_f2_3,
   est_mean_f3_3,est_mean_consprob_3,est_a2_1_3,est_a2_2_3,est_a2_3_3,est_a3_1_3,
   est_a3_2_3,est_b2_1_3,est_b2_2_3,est_b2_3_3,est_b3_1_3,est_b3_2_3,est_c2_1_3,
   est_c2_2_3,est_c2_3_3,sim_tot_consumption_3,sim_mean_frequency_3,sim_mean_f1_3,
   sim_mean_f2_3,sim_mean_f3_3,sim_mean_consprob_3)