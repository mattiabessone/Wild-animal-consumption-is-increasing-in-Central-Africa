# Define vectors to store draws
est_tot_consumption_3<-vector() # Estimated consumed wildmeat (n_draws * n_rep = 1000 draws for each dataset)
est_mean_frequency_3<-vector()
est_mean_consprob_3<-vector()
est_a2_1_3<-vector()
est_a2_2_3<-vector()
est_a2_3_3<-vector()
est_a3_1_3<-vector()
est_a3_2_3<-vector()
est_b2_1_3<-vector()
est_b2_2_3<-vector()
est_b2_3_3<-vector()
est_b3_1_3<-vector()
est_b3_2_3<-vector()
est_c2_1_3<-vector()
est_c2_2_3<-vector()
est_c2_3_3<-vector()
est_c3_3<-vector()
est_sigma_3<-vector()
LMs_uncertainty_3<-matrix(NA,n_rep,2)
sd_phi<-vector()

sim_tot_consumption_3<-vector() 
sim_mean_frequency_3<-vector()
sim_mean_consprob_3<-vector()


source("Code/sim_consumption.R")
dat<-list(N_recalls=sum(surveyed),N_region=n_regions,N_ethnicity=2L,N_missing_AME=N_missing_AME,N_sites=n_sites,
          sampled_region=c(length(region[region==1]),length(region[region==2]),length(region[region==3])),
          consumption=df$consumption,frequency=df$observed_frequency,quantity=df$quantity,
          region=df$region,forest=df$forest,AME=df$AME,AME_missidx=AME_missing,
          ethnicity=df$ethnicity,accessibility=df$access,
          N_sites=n_sites,pred_region=as.vector(study_area),pred_forest=as.vector(standardize(forest_cov)),
          pred_accessibility=as.vector(standardize(accessibility)),
          sampled_hh=c(length(df$region[df$region==1]),length(df$region[df$region==2]),length(df$region[df$region==3])),
          population=pop_by_site,days=df$days,distance_matrix=D,site=df$site)
sim_tot_consumption_3<-tot_grams_consumed
sim_mean_frequency_3<-average_frequency
sim_mean_consprob_3<-probability_consumption

detach(package:rethinking,unload=TRUE)
detach(package:cmdstanr,unload=TRUE)
fit<-stan(data=dat,"Model/simulation_model.stan",chains=1,iter=n_iter,init="0")
est_tot_consumption_3<-as.vector(extract(fit,pars="total_grams_consumed")$total_grams_consumed)
est_mean_frequency_3<-as.vector(extract(fit,pars="mean_frequency_consumption")$mean_frequency_consumption)
est_mean_consprob_3<-as.vector(extract(fit,pars="proportion_sites_consuming")$proportion_sites_consuming)
est_a2_1_3<-as.vector(extract(fit,pars="a2")$a2[,1])
est_a2_2_3<-as.vector(extract(fit,pars="a2")$a2[,2])
est_a2_3_3<-as.vector(extract(fit,pars="a2")$a2[,3])
est_a3_1_3<-as.vector(extract(fit,pars="a3")$a3[,1])
est_a3_2_3<-as.vector(extract(fit,pars="a3")$a3[,2])
est_b2_1_3<-as.vector(extract(fit,pars="b2")$b2[,1])
est_b2_2_3<-as.vector(extract(fit,pars="b2")$b2[,2])
est_b2_3_3<-as.vector(extract(fit,pars="b2")$b2[,3])
est_b3_1_3<-as.vector(extract(fit,pars="b3")$b3[,1])
est_b3_2_3<-as.vector(extract(fit,pars="b3")$b3[,2])
est_c2_1_3<-as.vector(extract(fit,pars="c2")$c2[,1])
est_c2_2_3<-as.vector(extract(fit,pars="c2")$c2[,2])
est_c2_3_3<-as.vector(extract(fit,pars="c2")$c2[,3])
est_c3_3<-as.vector(extract(fit,pars="c3")$c3)
est_sigma_3<-as.vector(extract(fit,pars="sigma")$sigma)

for (i in 1:dat$N_recalls){
  sd_phi[i]<-sd(extract(fit,pars="phi")$phi[,i])
}
sd_check<-as.data.frame(cbind(days=dat$days,sd_phi=sd_phi))
sd_check<-subset(sd_check,sd_check$days>0)
lm<-lm(sd_check$sd_phi~sd_check$days)
LMs_uncertainty3[1,1]<-lm$coefficients[1]
LMs_uncertainty_3[1,2]<-lm$coefficients[2]

for(i in 2:n_rep){
  source("Code/sim_consumption.R")
  dat<-list(N_recalls=sum(surveyed),N_region=n_regions,N_ethnicity=2L,N_missing_AME=N_missing_AME,N_sites=n_sites,
            sampled_region=c(length(region[region==1]),length(region[region==2]),length(region[region==3])),
            consumption=df$consumption,frequency=df$observed_frequency,quantity=df$quantity,
            region=df$region,forest=df$forest,AME=df$AME,AME_missidx=AME_missing,
            ethnicity=df$ethnicity,accessibility=df$access,
            N_sites=n_sites,pred_region=as.vector(study_area),pred_forest=as.vector(standardize(forest_cov)),
            pred_accessibility=as.vector(standardize(accessibility)),
            sampled_hh=c(length(df$region[df$region==1]),length(df$region[df$region==2]),length(df$region[df$region==3])),
            population=pop_by_site,days=df$days,distance_matrix=D,site=df$site)
  sim_tot_consumption_3<-append(sim_tot_consumption_3,tot_grams_consumed)
  sim_mean_frequency_3<-append( sim_mean_frequency_3,average_frequency)
  sim_mean_consprob_3<-append(sim_mean_consprob_3,probability_consumption)
  
  detach(package:rethinking,unload=TRUE)
  detach(package:cmdstanr,unload=TRUE)
  fit<-stan(data=dat,"Model/simulation_model.stan",chains=1,iter=n_iter,init="0")
  est_tot_consumption_3<-append(est_tot_consumption_3,as.vector(extract(fit,pars="total_grams_consumed")$total_grams_consumed))
  est_mean_frequency_3<-append(est_mean_frequency_3,as.vector(extract(fit,pars="mean_frequency_consumption")$mean_frequency_consumption))
  est_mean_consprob_3<-append(est_mean_consprob_3,as.vector(extract(fit,pars="proportion_sites_consuming")$proportion_sites_consuming))
  est_a2_1_3<-append(est_a2_1_3,as.vector(extract(fit,pars="a2")$a2[,1]))
  est_a2_2_3<-append(est_a2_2_3,as.vector(extract(fit,pars="a2")$a2[,2]))
  est_a2_3_3<-append(est_a2_3_3,as.vector(extract(fit,pars="a2")$a2[,3]))
  est_a3_1_3<-append(est_a3_1_3,as.vector(extract(fit,pars="a3")$a3[,1]))
  est_a3_2_3<-append(est_a3_2_3,as.vector(extract(fit,pars="a3")$a3[,2]))
  est_b2_1_3<-append(est_b2_1_3,as.vector(extract(fit,pars="b2")$b2[,1]))
  est_b2_2_3<-append(est_b2_2_3,as.vector(extract(fit,pars="b2")$b2[,2]))
  est_b2_3_3<-append(est_b2_3_3,as.vector(extract(fit,pars="b2")$b2[,3]))
  est_b3_1_3<-append(est_b3_1_3,as.vector(extract(fit,pars="b3")$b3[,1]))
  est_b3_2_3<-append(est_b3_2_3,as.vector(extract(fit,pars="b3")$b3[,2]))
  est_c2_1_3<-append(est_c2_1_3,as.vector(extract(fit,pars="c2")$c2[,1]))
  est_c2_2_3<-append(est_c2_2_3,as.vector(extract(fit,pars="c2")$c2[,2]))
  est_c2_3_3<-append(est_c2_3_3,as.vector(extract(fit,pars="c2")$c2[,3]))
  est_c3_3<-append(est_c3_3,as.vector(extract(fit,pars="c3")$c3))
  est_sigma_3<-append(est_sigma_1,as.vector(extract(fit,pars="sigma")$sigma))
  for (j in 1:dat$N_recalls){
    sd_phi[j]<-sd(extract(fit,pars="phi")$phi[,j])
  }
  sd_phi<-subset(sd_phi,dat$days>0)
  sd_check<-as.data.frame(cbind(days=dat$days[dat$days>0],sd_phi=sd_phi))
  lm<-lm(sd_check$sd_phi~sd_check$days)
  for(i in 2:n_rep){
    LMs_uncertainty_1[i,1]<-as.numeric(lm$coefficients[1])
    LMs_uncertainty_1[i,2]<-as.numeric(lm$coefficients[2])
  }
}
