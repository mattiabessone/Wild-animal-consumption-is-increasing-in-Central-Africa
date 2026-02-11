# Bessone et al. (2025) - Wild animal consumption is increasing across Central Africa #

#### Run model selection precedure described in the manuscript ####

# Load packages and data
library(rstan)
library(loo)
# Import demo dataset
d<-read.csv("Data/demo_dataset.csv") # demo consumption data
d_mat<-read.csv("Data/demo_distance_matrix.csv") # distance matrix used to assess the models
# To request access to full dataset see "Data/Data_availbility.md" or Data availbility statement on the paper

#### 1) Sub-model consumption probability ####
# Prepare data
d_mat[,1]<-NULL
colnames(d_mat)<-NULL
d$study_char_id<-match(d$study_id,unique(d$study_id))
d$site_char_id<-match(d$location_id,unique(d$location_id))
d$household_char_id<-match(d$household_id,unique(d$household_id))
d$education_char_id<-match(d$education,unique(d$education))
AME_missidx<-which(is.na(d$num_ame))
N_missing_AME<-nlevels(as.factor(AME_missidx))
d$AME<-d$num_ame
d$AME[is.na(d$num_ame)] <- -10
N_studies<-nlevels(as.factor(d$study_id))
N_sites<-nlevels(as.factor(d$location_id))
N_households<-nlevels(as.factor(d$household_id))
N_recalls<-length(d$recall_id)
N_education_levels<-nlevels(as.factor(d$education))

# Wrap data as list
data<-list(N_studies=N_studies,N_sites=N_sites,N_households=N_households,N_recalls=N_recalls
           ,N_education_levels=N_education_levels,N_missing_AME=N_missing_AME
           ,consumption=d$consumption,frequency=d$frequency,quantity=d$quantity
           ,study=d$study_char_id,site=d$site_char_id,hh=d$household_char_id
           ,HPD=scale(d$HPD)[,1],HDI=scale(d$HDI)[,1],accessibility=scale(d$REM)[,1],FCI=scale(d$FCI)[,1]
           ,education=d$education_char_id,days=d$days
           ,AME=d$AME,AME_missidx=AME_missidx,mdays=d$mdays,distance_matrix=d_mat)

##### 1.1) Test full model without spatial autocorrelation ####
m1_full<-stan(data=data,"Stan_models/check_Bernoulli.stan",chains=2,cores=2,iter=2000)
pairs(m1_full,pars=c("a1","a2","a3","a4"))
loo_m1_full<-loo(m1_full)
##### 1.2) Test full model with spatial autocorrelation ####
m1_full_GP<-stan(data=data,"Stan_models/check_Bernoulli_GP.stan",chains=2,cores=2,iter=2000)
pairs(m1_full_GP,pars=c("a1","a2","a3","a4"))
loo_m1_full_GP<-loo(m1_full_GP)
##### 1.3) Compare ELPD ####
loo_compare(loo_m1_full,loo_m1_full_GP)

#### 2) Sub-model frequency of consumption ####
# We only make use of data-points with frequency > 0
d_fr<-d[d$frequency>0,]
# Then we aggregate the data by HH
d_fr<-aggregate(d_fr,list(d_fr$household_id),mean)
# We remove sites with no entries with frequency > 0
d_mat_fr<-d_mat[-c(101:127),-c(101:127)] # ATTENTION: ignore this line when runnong on full dataset
# Prepare data
d_fr$study_char_id<-match(d_fr$study_id,unique(d_fr$study_id))
d_fr$site_char_id<-match(d_fr$location_id,unique(d_fr$location_id))
d_fr$household_char_id<-match(d_fr$household_id,unique(d_fr$household_id))
d_fr$education_char_id<-match(d_fr$education,unique(d_fr$education))
AME_missidx<-which(is.na(d_fr$num_ame))
N_missing_AME<-nlevels(as.factor(AME_missidx))
d_fr$AME<-d_fr$num_ame
d_fr$AME[is.na(d_fr$num_ame)] <- -10
N_studies<-nlevels(as.factor(d_fr$study_id))
N_sites<-nlevels(as.factor(d_fr$location_id))
N_households<-nlevels(as.factor(d_fr$household_id))
N_recalls<-length(d_fr$recall_id)
N_education_levels<-nlevels(as.factor(d_fr$education))

# Wrap data as list
data_fr<-list(N_studies=N_studies,N_sites=N_sites,N_households=N_households,N_recalls=N_recalls
              ,N_education_levels=N_education_levels,N_missing_AME=N_missing_AME
              ,consumption=d_fr$consumption,frequency=d_fr$frequency,quantity=d_fr$quantity
              ,study=d_fr$study_char_id,site=d_fr$site_char_id,hh=d_fr$household_char_id
              ,HPD=scale(d_fr$HPD)[,1],HDI=scale(d_fr$HDI)[,1],accessibility=scale(d_fr$REM)[,1],FCI=scale(d_fr$FCI)[,1]
              ,education=d_fr$education_char_id,days=d_fr$days
              ,AME=d_fr$AME,AME_missidx=AME_missidx,mdays=d_fr$mdays,distance_matrix=d_mat_fr)
##### 2.1) Test full model without spatial autocorrelation ####
m2_full<-stan(data=data_fr,"Stan_models/check_Beta_hh.stan",chains=2,cores=2,iter=2000,init="0")
pairs(m2_full,pars=c("b1","b2","b3","b4"))
loo_m2_full<-loo(m2_full)
##### 2.2) Test full model with spatial autocorrelation ####
m2_full_GP<-stan(data=data_fr,"Stan_models/check_Beta_hh_GP.stan",chains=2,cores=2,iter=2000,init="0")
pairs(m2_full_GP,pars=c("b1","b2","b3","b4"))
loo_m2_full_GP<-loo(m2_full_GP)
##### 2.3) Compare ELPD ####
loo_compare(loo_m2_full,loo_m2_full_GP)

#### 3) Sub-model quantity consumed ####
# We only make use of data-points with frequency > 0
d_q<-d[d$quantity>0,]
# Prepare data
d_q$study_char_id<-match(d_q$study_id,unique(d_q$study_id))
d_q$site_char_id<-match(d_q$location_id,unique(d_q$location_id))
d_q$household_char_id<-match(d_q$household_id,unique(d_q$household_id))
d_q$education_char_id<-match(d_q$education,unique(d_q$education))
AME_missidx<-which(is.na(d_q$num_ame))
N_missing_AME<-nlevels(as.factor(AME_missidx))
d_q$AME<-d_q$num_ame
d_q$AME[is.na(d_q$num_ame)] <- -10
N_studies<-nlevels(as.factor(d_q$study_id))
N_sites<-nlevels(as.factor(d_q$location_id))
N_households<-nlevels(as.factor(d_q$household_id))
N_recalls<-length(d_q$recall_id)
N_education_levels<-nlevels(as.factor(d_q$education))

# Wrap data as list
data_q<-list(N_studies=N_studies,N_sites=N_sites,N_households=N_households,N_recalls=N_recalls
              ,N_education_levels=N_education_levels,N_missing_AME=N_missing_AME
              ,consumption=d_q$consumption,frequency=d_q$frequency,quantity=d_q$quantity
              ,study=d_q$study_char_id,site=d_q$site_char_id,hh=d_q$household_char_id
              ,HPD=scale(d_q$HPD)[,1],HDI=scale(d_q$HDI)[,1],accessibility=scale(d_q$REM)[,1],FCI=scale(d_q$FCI)[,1]
              ,education=d_q$education_char_id,days=d_q$days
              ,AME=d_q$AME,AME_missidx=AME_missidx,mdays=d_q$mdays,distance_matrix=d_mat)
##### 3.1) Test full model without spatial autocorrelation ####
m3_full<-stan(data=data_q,"Stan_models/check_Gamma.stan",chains=2,cores=2,iter=2000)
loo_m3_full<-loo(m3_full)
pairs(m3_full,pars=c("c1","c2","c3","c4","c5"))
##### 3.2) Test full model without continuous covariates ####
m3_null<-stan(data=data_q,"Stan_models/check_Gamma2.stan",chains=2,cores=2,iter=2000)
loo_m3_null<-loo(m3_null)
##### 2.3) Compare ELPD ####
loo_compare(loo_m3_full,loo_m3_null)
