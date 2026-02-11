# Bessone et al. (2025) - Wild animal consumption is increasing across Central Africa #

#### 1) Run full models described in the manuscript ####
# Load packages
library(datawizard)
library(rstan)
options(mc.cores = parallel::detectCores())

##### 1.1) Model used for final predictions (3 location types) ####
# Organise data list
source("R_code/list_data_m1.R")
# Run model
m1<-stan(data=data,"Stan_models/consumption_model.stan",chains=4,iter=2000,warmup=1000,init="0",save_warmup=FALSE,control=list(max_treedepth=12))

##### 1.2) Model investigating the interaction between location type and education (7 levels) ####
# Organise data list
source("R_code/list_data_m2.R")
# Run model
m2<-stan(data=data,"Stan_models/consumption_model_EDint.stan",chains=4,iter=2000,warmup=1000,init="0",save_warmup=FALSE,control=list(max_treedepth=12))

##### 1.3) Model investigating only 2 location types (rural vs. urban) ####
# Organise data list
source("R_code/list_data_m3.R")
# Run model
m3<-stan(data=data,"Stan_models/consumption_model_2LT.stan",chains=4,iter=2000,warmup=1000,init="0",save_warmup=FALSE,control=list(max_treedepth=12))

#### 2) Run models on data subset (4% of the full data set) ####
##### 2.1) Model used for final predictions (3 location types) ####
# Organise data list
source("R_code/list_data_m1s.R")
# Run model
m1s<-stan(data=data,"Stan_models/consumption_model.stan",chains=4,iter=2000,warmup=1000,init="0",save_warmup=FALSE,control=list(max_treedepth=12))

##### 2.2) Model investigating the interaction between location type and education (7 levels) ####
# Organise data list
source("R_code/list_data_m2s.R")
# Run model
m2s<-stan(data=data,"Stan_models/consumption_model_EDint.stan",chains=4,iter=2000,warmup=1000,init="0",save_warmup=FALSE,control=list(max_treedepth=12))

##### 2.3) Model investigating only 2 location types (rural vs. urban) ####
# Organise data list
source("R_code/list_data_m3s.R")
# Run model
m3s<-stan(data=data,"Stan_models/consumption_model_2LT.stan",chains=4,iter=2000,warmup=1000,init="0",save_warmup=FALSE,control=list(max_treedepth=12))
