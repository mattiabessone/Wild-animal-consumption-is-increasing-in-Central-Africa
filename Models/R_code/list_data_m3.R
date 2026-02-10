# Import again consumption data
d <- read.csv("Data/consumption_data.csv") 
d_mat<-read.csv("Data/distance_matrix.csv")/1000 # Distance matrix for autocorrelation in km
d_pred<-read.csv("Data/predictions_2LT.csv") # Prediction data (874 cells)

# Define 2 location type levels: 1="rural", 2="urban"
d$location_type<-ifelse(d$location_type>1,2,1)
# Define 3 education levels: 1="<secondary",2=">primary",3="unknown"
d$education<-ifelse(d$education<3,1,ifelse(d$education>6,3,2))
# Define 3 study types: 1="short-recall",2="long-recall",3="cooking-pot"
d$survey_type<-ifelse(d$survey_type<4,1,ifelse(d$survey_type==7,3,2))
# Define number of levels for all categorical variables
N_studies<-nlevels(as.factor(d$study_id))
N_sites<-nlevels(as.factor(d$location_id))
N_households<-nlevels(as.factor(d$household_id))
N_recalls<-length(d$recall_id)
N_education_levels<-nlevels(as.factor(d$education))
N_location_types<-nlevels(as.factor(d$location_type))
N_survey_types<-3L
N_periods<-nlevels(as.factor(d$period))

# Define location of Adult Male Equivalent (AME) missing data points
AME_missidx<-which(is.na(d$num_ame))
# And number of AME missing data points
N_missing_AME<-nlevels(as.factor(AME_missidx))
# Assign a dummy (negative) value to missing AME
d$num_ame[is.na(d$num_ame)] <- -10
# And calculate mean AME to center relevant prior
mean_AME<-mean(subset(d$num_ame,d$num_ame>0))

# Ensure there's no mistakes and order study, location and household ID. 
d$study_id <-match(d$study_id,unique(d$study_id))
d$location_id <-match(d$location_id,unique(d$location_id))
d$household_id <-match(d$household_id,unique(d$household_id))

#### Create frequency data by aggregating recalls by households
d_f<-aggregate(d,list(d$household_id),mean)

# Remove first column and column names from distance matrix
d_mat[,1]<-NULL
colnames(d_mat)<-NULL

# Organise prediction data
# Define number of scenarios and number of cells
N_scenarios<-3L
N_sites_pred<-length(d_pred$cell_ID)
# Wrap dynamic layers into matrices and fixed layer into vectors
HPD_pred<-as.matrix(cbind(standardize(d_pred$HPD_past),standardize(d_pred$HPD_recent),standardize(d_pred$HPD_present)))
REM_pred<-as.vector(standardize(d_pred$REM))
HDI_pred<-as.matrix(cbind(standardize(d_pred$HDI_past),standardize(d_pred$HDI_recent),standardize(d_pred$HDI_present)))
FCI_pred<-as.vector(standardize(d_pred$FCI))
ED_pred<-as.matrix(cbind(d_pred$ED_past,d_pred$ED_recent,d_pred$ED_present))
LT_pred<-array(NA,dim=c(N_sites_pred,N_location_types,N_scenarios))
LT_pred[,,1]<-as.matrix(cbind(d_pred$LT1_past,d_pred$LT2_past,d_pred$LT3_past))
LT_pred[,,2]<-as.matrix(cbind(d_pred$LT1_recent,d_pred$LT2_recent,d_pred$LT3_recent))
LT_pred[,,3]<-as.matrix(cbind(d_pred$LT1_present,d_pred$LT2_present,d_pred$LT3_present))
AME_pred<-as.matrix(cbind(d_pred$AME_past,d_pred$AME_recent,d_pred$AME_present))

# List data
data<-list(N_studies=N_studies,N_sites=N_sites,N_households=N_households,N_recalls=N_recalls
           ,N_education_levels=N_education_levels,N_periods=N_periods
           ,N_location_types=N_location_types,N_survey_types=N_survey_types,N_missing_AME=N_missing_AME
           ,consumption=d$consumption,frequency=d_f$frequency,quantity=d$quantity
           ,study=d$study_id,site=d$location_id,household=d$household_id
           ,HPD=standardize(d$HPD),HDI=standardize(d$HDI),REM=standardize(d$REM),FCI=standardize(d$FCI)
           ,education=d$education,period=d$period,location_type=d$location_type,survey_type=d$survey_type
           ,study_fr=d_f$study_id,site_fr=d_f$location_id,survey_type_fr=d_f$survey_type
           ,HPD_fr=standardize(d_f$HPD),HDI_fr=standardize(d_f$HDI),REM_fr=standardize(d_f$REM),FCI_fr=standardize(d_f$FCI)
           ,education_fr=d_f$education,period_fr=d_f$period,location_type_fr=d_f$location_type
           ,AME=d$num_ame,AME_missidx=AME_missidx,mean_AME=mean_AME,mdays=d_f$mdays,days=d$days,distance_matrix=d_mat
           ,N_scenarios=N_scenarios,N_sites_pred=N_sites_pred,HPD_pred=HPD_pred,REM_pred=REM_pred,HDI_pred=HDI_pred,FCI_pred=FCI_pred
           ,ED_pred=ED_pred,LT_pred=LT_pred,AME_pred=AME_pred)
