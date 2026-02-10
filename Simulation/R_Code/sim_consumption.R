#### SIMULATION STUDY - Wild meat consumption data generation - Mattia Bessone ####
library (rethinking)
library(LaplacesDemon)
###### Define study area
study_area <- matrix (0,rows,rows)                                              # matrix, representing study area (each site = 25 km2)
n_sites <- rows^2                                                               # number of sites
site_id<-as.vector(study_area)
site_id<-seq(1,n_sites,by=1) 
n_regions <- 3                                                                  # define 3 regions
regions <- c(1,2,3)
study_area [1:rows,1:(rows*0.4)] <- 1                                           # define region 1
study_area [1:rows,((rows*0.4)+1):(rows*0.8)] <- 2                              # define region 2
study_area<-ifelse(study_area>0,study_area,3)                                   # define region 3

###### Create distance matrix
# First we define coordinates of the center of each cell
x<-matrix(NA,rows,rows)
x[,1]<-rep(0.5,rows)
for(j in 2:rows){
  x[1:rows,j]<- rep(0.5+(j-1),rows)
}
y<-matrix(NA,rows,rows)
y[1,]<-rep(rows-0.5,rows)

for(i in 2:rows){
  y[i,1:rows]<- rep(rows-(-0.5+i),rows)
}

lat<-as.vector(x)
long<-as.vector(y)

# Then we generate a random location around that point within the cell
for(i in 1:rows^2){
  lat[i]<-runif(1,lat[i]-0.5,lat[i]+0.5)
  long[i]<-runif(1,long[i]-0.5,long[i]+0.5)
}
coords<-as.data.frame(cbind(long,lat))
#And we calculate distances between sites
D<-as.matrix(dist(coords))

###### Simulate 2 continuous variables, e.g. forest cover and accessibility
forest_cov <- matrix(0,rows,rows)
alpha1 <- c(0.75,0.40,0.20)                                                      # define a for each region
beta1 <- c(0.17,0.31,0.61)                                                       # define b for each region

accessibility <- matrix(0,rows,rows)
alpha2 <- c(25,60,120)                                                      # define a for each region
beta2 <- 20                                                       # define b for each region

for (i in 1:n_sites){                                                           # then we populate the matrix
  if(study_area[i]==1) forest_cov[i] <- rbeta(1,alpha1[1],beta1[1])               # populate region 1
  if(study_area[i]==2) forest_cov[i] <- rbeta(1,alpha1[2],beta1[2])               # populate region 2
  if(study_area[i]==3) forest_cov[i] <- rbeta(1,alpha1[3],beta1[3])               # populate region 3
}

for (i in 1:n_sites){                                                           # then we populate the matrix
  if(study_area[i]==1) accessibility[i] <- rgamma2(1,alpha2[1],beta2)               # populate region 1
  if(study_area[i]==2) accessibility[i] <- rgamma2(1,alpha2[2],beta2)               # populate region 2
  if(study_area[i]==3) accessibility[i] <- rgamma2(1,alpha2[3],beta2)               # populate region 3
}

accessibility <- standardize(accessibility)

###### Simulate number of household per site
hh <- matrix(0,rows,rows)                                                       # create a matrix specifying number of households in each site
hh_dens<- c(20,30,50)                                                           # set mean number of HH per region

for (i in 1:n_sites){                                                           # populate matrix
  if(study_area[i] == 1) hh[i] <- ceiling(rgamma2(1,hh_dens[1],1))
  if(study_area[i] == 2) hh[i] <- ceiling(rgamma2(1,hh_dens[2],1))
  if(study_area[i] == 3) hh[i] <- ceiling(rgamma2(1,hh_dens[3],1))
}

N_household<-sum(hh)                                                  # Specify number of HHs
hh_id<-seq(1,N_household,by=1)                                                  # ID the households

##### Simulate number of people per household (Adult Male Equivalent - AME)
mean_n_hh <- 5                                                                  # set average number AME per household

# Simulate components of each households
hh_people <- c()                                                                # create empty vector

for (i in 1:N_household){                                                       # populate vector
  hh_people[i] <- ceiling(rgamma2(1,mean_n_hh,2))
}

#Check total population (AME) of study area
tot_pop<-sum(hh_people)
tot_pop

##### Simulate sampled sites (yes/no)
sampled <- matrix(NA,rows,rows)                                                 # create matrix specifying if a site was surveyed or not
frequency<-matrix(0,rows,rows)
quantity<-matrix(0,rows,rows)
for(i in 1:n_sites) sampled[i] <- rbern(1,coverage)                             # populate matrix according to coverage

# Simulate surveys recording frequency
prob_frequency <- 0.8
for(i in 1:n_sites) if(sampled[i] == 1) frequency[i] <- rbern(1,prob_frequency)

# Simulate surveys recording quantities
prob_quantity <- 0.5
for(i in 1:n_sites) if(frequency[i] == 1) quantity[i] <- rbern(1,prob_quantity)

# Create dataframe of study area
df<-as.data.frame(cbind(site_id,as.vector(study_area),as.vector(forest_cov),as.vector(accessibility),as.vector(hh),as.vector(sampled),as.vector(frequency),as.vector(quantity)))
names(df)<-c("site_id","region","forest_cov","accessibility","n_hh","sampled","frequency","quantity")


# Create vectors replicating the info contained in the first site per the number of HH in the first site
source("make_vector_fun.R")

hh_site<-make_vector(df$site_id[1],df$n_hh[1],df$site_id,df$n_hh,n_sites)
hh_region<-make_vector(df$region[1],df$n_hh[1],df$region,df$n_hh,n_sites)
hh_forest<-make_vector(df$forest_cov[1],df$n_hh[1],df$forest_cov,df$n_hh,n_sites)
hh_accessibility<-make_vector(df$accessibility[1],df$n_hh[1],df$accessibility,df$n_hh,n_sites)
hh_sampled<-make_vector(df$sampled[1],df$n_hh[1],df$sampled,df$n_hh,n_sites)
hh_frequency<-make_vector(df$frequency[1],df$n_hh[1],df$frequency,df$n_hh,n_sites)
hh_quantity<-make_vector(df$quantity[1],df$n_hh[1],df$quantity,df$n_hh,n_sites)


#Aggregate population by site
pop_by_site<-aggregate(hh_people,by=list(hh_site),FUN=sum)
pop_by_site<-pop_by_site[,2]

###### Simulate 1 factorial variables at the HH level, e.g., ethnicity (2 levels)
ethnicity_prob <- c(0.3,0.7)                                                    # fix for all regions

hh_ethnicity <- vector(length=N_household)                                      # create empty vector

for (i in 1:N_household){                                                       # populate vector
  hh_ethnicity[i] <- rcategorical(1,ethnicity_prob)
}

###### Simulate survey of HH in sampled sites
sampling_effort <- c(0.8,0.7,0.5)
surveyed <- vector(length=N_household)

for (i in 1:N_household){
  if (hh_sampled[i]==1) surveyed[i]<-rbern(1,sampling_effort[1])
}


##### LEVEL 1 - simulate if a hh consumes bushmeat or not, binary, conditional on region and forest coverage

pi <- vector(length=N_household)                                                # we create a vector of consumption (yes/no)  for each household, based on mean_pi
mean_pi<-vector(length=N_household)

# And simulate spatial autocorrelation between sites
delta<-1e-9                  # Fix delta, to avoid 0s and improve sampling

X<-matrix(NA,n_sites,n_sites)
for (i in 1:(n_sites-1)){
  X[i,i] <- zeta + delta;
  for (j in (i+1):n_sites){
    X[i,j] <- zeta * exp(-rho * D[i,j]^2);
    X[j,i] <- X[i,j];
  }
}
X[n_sites,n_sites]<-zeta+delta

eta<-rnorm(n_sites)
epsilon<- c(t(chol(X))%*%eta)

# Then we model the mean consumption probability conditional on covariates and spatial autocorrelation

for (i in 1:N_household){
  mean_pi[i] <- inv_logit(a1[hh_region[i]] + a2[hh_region[i]] * hh_forest[i] + a3[hh_ethnicity[i]] + epsilon[hh_site[i]])
}

for (i in 1:N_household){
  pi[i] <- rbern(1,mean_pi[i])
}

# LEVEL 2 - simulate frequency of consumption, days / week, in sites that consume bushmeat 
# We model frequency as Beta distributed (pars = a & b), with a == mean * kappa and b = (1 - mean) * kappa
mean_phi<-rep(-9,N_household)
mu_phi<-rep(-9,N_household)

# And simulate spatial autocorrelation between sites
Y<-matrix(NA,n_sites,n_sites)
for (i in 1:(n_sites-1)){
  Y[i,i] <- xi + delta;
  for (j in (i+1):n_sites){
    Y[i,j] <- xi * exp(-omicron * D[i,j]^2);
    Y[j,i] <- Y[i,j];
  }
}
Y[n_sites,n_sites]<-xi+delta

omega<-rnorm(n_sites)
lambda<- c(t(chol(Y))%*%omega)

# Then We model mean frequency phi, conditional on ethnicity, education and accessibility and spatial autocorrelation


for (i in 1:N_household){
  if(pi[i]==1){
    mu_phi[i] <- b1[hh_region[i]] + b2[hh_region[i]] * hh_accessibility[i] + b3[hh_ethnicity[i]] + lambda[hh_site[i]] 
  }
}

phi <- rep(0,N_household)

for (i in 1:N_household){
  if(pi[i]==1){
    phi[i] <- rbeta(1,inv_logit(mu_phi[i])*kappa,(1-inv_logit(mu_phi[i]))*kappa)
  }
}

# LEVEL 3 - simulate amounts consumed per day in households (not the individual) that consume bushmeat at least once per week, i.e. frequency  > 0
# We model amounts consumed as Gamma distributed, pars = mu (mean) and theta, conditional on region, education, educated people consume less frequently but higher, we also include MBMI, and define it n.s.
mean_mu<-rep(0,N_household)
mu <- rep(0,N_household)

for (i in 1:N_household){
  if(phi[i]>0){
    mean_mu[i] <- exp(c1[hh_region[i]] + c2[hh_region[i]] * hh_accessibility[i] + c3 * hh_people[i])
  }
}

# Then, when we generate the data, we divide the daily consumption in the household by the numner of inhabitant

for (i in 1:N_household){
  if(mean_mu[i]>0){
    mu[i] <- rgamma2(1,mean_mu[i],theta)
  }
}

##### Calculate parameters to be retrieved by the model
consumed<-vector(length=N_household)
for (i in 1:N_household) consumed[i] <- pi[i] * phi[i] * mu[i] * hh_people[i]
tot_grams_consumed <- sum(consumed)
probability_consumption <- mean(pi)
average_frequency <- mean(phi)


##### Create data.frame containing data at the household level
#extract only surveyed hh within surveyed sites
hh_consum<- ifelse(surveyed==1,pi,-1)                                           # this vector contains sampled hh, with observed consumption (yes/no). Not sampled hh are coded as -1
temp<-c()
for (i in 1:N_household) temp[i]<-ifelse(hh_consum[i]>-1, rbern (1,sampling_effort[2]), -1)
hh_freq<-ifelse(temp>0, phi, -1)

temp<-c()
for (i in 1:N_household) temp[i]<-ifelse(hh_freq[i]>-1, rbern (1,sampling_effort[3]), -1)
hh_quant<-ifelse(temp>0, mu, -1)

#### Prepare data for import in Stan
site<-subset(hh_site,surveyed==1)
consumption<-subset(hh_consum,surveyed==1)
people<-subset(hh_people,surveyed==1)
frequency<-subset(hh_freq,surveyed==1) 
frequency<-ifelse(frequency>0.999,0.99,frequency)                               # Stan would not sample from a Beta distribution if input value = 1 or = 0
days<-rep(-1,length(frequency))
sd_freq<-rep(-1,N_household)

observed_frequency<-rep(-1,length(frequency))

for (i in 1:length(frequency)){
  if(frequency[i]>0){
    days[i]<-round(runif(1,1,360),0)
    sd_freq[i]<-sigma*(365-days[i])
    observed_frequency[i]=inv_logit(rnorm(1,logit(frequency[i]),sd_freq[i]))
  }
}

quantity<-subset(hh_quant,surveyed==1) * people # CAREFUL HERE!
region<-subset(hh_region,surveyed==1)
forest<-standardize(subset(hh_forest,surveyed==1))

mean_AME<-mean(people)
impute<-vector()
AME<-vector()
for (i in 1:length(people)){                                                       # populate vector
  impute[i]<- rbern(1,prob_missing_AME)
}
AME<-ifelse(impute==1,-10,people)
AME_missing<-which(AME<0)
N_missing_AME<-length(AME_missing)

ethnicity<-subset(hh_ethnicity, surveyed==1)
filter<-subset(hh_freq,hh_freq>-1 | surveyed==1)
access<-subset(hh_accessibility,hh_freq>-1 | surveyed==1)
access<-ifelse(filter>-1,hh_accessibility,0)

df<-as.data.frame(cbind(consumption,observed_frequency,quantity,region,forest,AME,ethnicity,access,days,site))
df<-df[order(df$region,decreasing = FALSE), ] 

