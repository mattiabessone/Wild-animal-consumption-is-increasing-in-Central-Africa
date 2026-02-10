# Load pacakges and data
library(rstan)
library(loo)
d<-read.csv("~/Data/consumption_data.csv") # consumption data
d_mat<-read.csv("~/Data/distance_matirx_subset.csv") #distance matrix used to assess the models
# To request access see "Data/Data_availbility.md" or Data availbility statemnet on the paper


##Prepare data####
d<-subset(d,d$study_id>9 & d$study_id<15) # subset of data used to assess models
d_mat[,1]<-NULL
colnames(d_mat)<-NULL
d$study_char_id<-match(d$study_id,unique(d$study_id))
d$site_char_id<-match(d$site_id,unique(d$site_id))
d$household_char_id<-match(d$household_id,unique(d$household_id))
d$education_char_id<-match(d$education,unique(d$education))
AME_missidx<-which(is.na(d$num_ame))
N_missing_AME<-nlevels(as.factor(AME_missidx))
d$AME<-d$num_ame
d$AME[is.na(d$num_ame)] <- -10
N_studies<-nlevels(as.factor(d$study_id))
N_sites<-nlevels(as.factor(d$site_id))
N_households<-nlevels(as.factor(d$household_id))
N_recalls<-length(d$recall_id)
N_education_levels<-nlevels(as.factor(d$education))

# Wrap data as list
data<-list(N_studies=N_studies,N_sites=N_sites,N_households=N_households,N_recalls=N_recalls
           ,N_education_levels=N_education_levels,N_missing_AME=N_missing_AME
           ,consumption=d$consumption,frequency=d$frequency,quantity=d$quantity+0.000001
           ,study=d$study_char_id,site=d$site_char_id,hh=d$household_char_id
           ,HPD=scale(d$HPD)[,1],HDI=scale(d$HDI)[,1],accessibility=scale(d$REM)[,1],FCI=scale(d$FCI)[,1]
           ,education=d$education_char_id,days=d$days
           ,AME=d$AME,AME_missidx=AME_missidx,mdays=d$mdays,distance_matrix=d_mat)

# Check consumption probability
m1_full<-stan(data=data,"check_Bernoulli.stan",chains=2,cores=2,iter=2000) #Test full model
pairs(m1_full,pars=c("a1","a2","a3","a4"))
loo_m1_full<-loo(m1_full)

m1_full_GP<-stan(data=data,"check_Bernoulli_GP.stan",chains=2,cores=2,iter=2000) #Test full model
pairs(m1_full_GP,pars=c("a1","a2","a3","a4"))
loo_m1_full_GP<-loo(m1_full_GP)

loo_compare(loo_m1_full,loo_m1_full_GP)

# Check quantity
m3_full<-stan(data=data,"check_Gamma2.stan",chains=2,cores=2,iter=2000) #Test full model
loo_m3_full<-loo(m3_full)
pairs(m3_full,pars=c("c1","c2","c3","c4","c5","c6"))


# Check frequency
# First we must aggregate the data by HH
d_fr<-aggregate(d,list(d$household_id),mean)
# Wrap data as list
data_fr<-list(N_studies=N_studies,N_sites=N_sites,N_households=N_households,N_recalls=N_recalls
           ,N_education_levels=N_education_levels,N_missing_AME=N_missing_AME
           ,consumption=d_fr$consumption,frequency=d_fr$frequency+0.0000001,quantity=d_fr$quantity
           ,study=d_fr$study_char_id,site=d_fr$site_char_id,hh=d_fr$household_char_id
           ,HPD=scale(d_fr$HPD)[,1],HDI=scale(d_fr$HDI)[,1],accessibility=scale(d_fr$REM)[,1],FCI=scale(d_fr$FCI)[,1]
           ,education=d_fr$education_char_id,days=d_fr$days
           ,AME=d_fr$AME,AME_missidx=AME_missidx,mdays=d_fr$mdays,distance_matrix=d_mat)

m2_full<-stan(data=data_fr,"check_Beta_hh.stan",chains=2,cores=2,iter=2000,init="0") #Test full model
pairs(m2_full,pars=c("b1","b2","b3","b4"))
loo_m2_full<-loo(m2_full)

m2_full_GP<-stan(data=data_fr,"check_Beta_hh_GP.stan",chains=2,cores=2,iter=2000,init="0") #Test full model
pairs(m2_full_GP,pars=c("b1","b2","b3","b4"))
loo_m2_full_GP<-loo(m2_full_GP)

loo_compare(loo_m2_full,loo_m2_full_GP)
