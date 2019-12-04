# data processing
##################### 
# DATA PREP
library(haven) # for reading .sas data files
library(tidyverse) #for data manipulation
select <- dplyr::select
d <- read_sas('Data/Triple RCT COPD/exacevents.sas7bdat')

copd <- d %>% 
  select(-matches(paste0(c('LAMA','LABA','ICS','statin','azithromycin'),collapse='|')), # covariates used for randomization, so would probably have no predictive power
         -antib1yr, # remove rescaled variables
         -bmi10,
         -fev1pp100,
         -sgrq10,
         -age10,
         -tte0) %>% # time to event; we need gaptime
  select(-hosp1yr) %>%  # non-avaialble for the optimal cohort
  mutate(fev1fvcratio = fev1/fvc) %>% 
  filter(!(STATCOPE==1)) %>% # dont use STATCOPE; cannot find deaths
  select(-STATCOPE,
         -MACRO, # use the optimal indicator instead "OPTIMAL"
         -trial) %>%  # trial is equivalent to OPTIMAL 
  arrange(ID,event_date)

length(unique(copd$ID)) #1566
nrow(copd) #3822

# find patient with any missing data
na.IDs <- apply(copd,1,anyNA)
na.IDs <- copd$ID[which(na.IDs==T)] %>% unique()

copd <- copd %>% 
  filter(!(ID %in% na.IDs))

length(na.IDs) # 43 patients removed
nrow(copd)

# remove duplicates
copd <- copd %>% distinct()
nrow(copd) #no duplicates

# recurrent event indicator
copd <- copd %>% 
  mutate(recurrent_event = (event==2))

# max.period
max.period <- copd %>% 
  group_by(ID) %>% 
  summarise(max_period = max(period)) %>% 
  select(ID,max_period)

copd <- left_join(copd,max.period,by='ID')

# create gap time response variables
copd <- copd %>% 
  filter(!(period < max_period & recurrent_event==0 )) %>% 
  mutate(gaptime = ifelse(ID != lag(ID),event_date-Rand_Date,
                          event_date - lag(event_date)))

copd$gaptime[1] <- copd$event_date[1]-copd$Rand_Date[1]

# max.period
max.period <- copd %>% 
  group_by(ID) %>% 
  summarise(max_period = n()) %>% 
  select(ID,max_period)

copd <- left_join(copd %>% select(-max_period),max.period,by='ID')

# re-define period
tmp_period <- rep(0,nrow(copd))
tmp_period[1] <- 1
copd.id <- copd$ID
for (i in 2:nrow(copd)){
 if(copd.id[i] != copd.id[i-1]){
   tmp_period[i] <- 1
 } else{
   tmp_period[i] <- tmp_period[i-1]+1
 }
}
copd$period <- tmp_period
View(copd)

# double check
copd.re.sum <- copd %>% 
  group_by(ID) %>% 
  summarize(n=n(),
            sum_re = sum(recurrent_event),
            max_period=max(max_period)) %>% 
  mutate(sanity_check= (n-sum_re)==1 )
sum(copd.re.sum$sanity_check==F)

# check which variables are time-dependent
copd %>%
  group_by(ID) %>%
  summarise_at(vars(trtgroup:fev1fvcratio),list(sd=sd)) -> tmp

time_independent <- apply(tmp[,-1],2, function(x){ all(x==0 || is.na(x))})
time_independent[time_independent==T]    
time_independent[time_independent==F]   

# baseline covariates
baseline.covariates <- c('age','gender','packyears','bmi',
                         'fev1','fev1pp', 'fev1fvcratio','fvc', 'stgtotalc00',
                         'nowsmk','oxygen','sgrq',
                         'OPTIMAL',
                         'rand_month','rand_season','rand_seasonality')

# factors
factor.covariates <- c('ID','trtgroup',
                       'oxygen','nowsmk',
                       'gender','OPTIMAL',
                       'rand_month','rand_season','rand_seasonality',
                       'event_month','event_season','event_seasonality')

# change them to factor variables
copd[,factor.covariates] <- lapply(copd[,factor.covariates],as.factor)

# compute the death indicator
optimal <- read_csv("Data/OPTIMAL/handydata.csv")
optimal <- optimal %>%
  filter(died==1) %>% 
  select(subno,Followuptime,deathday,died) %>% 
  mutate(sanity_check = Followuptime==deathday) %>% 
  mutate(ID = paste0(subno,'-','OPTIM')) %>% 
  select(-subno)
optimal.death.id <- optimal$ID 
# there are two the last follow up times do not match the death time;
# died 12 and 6 days after the last follow up times
optimal.problem <- optimal %>%
  filter(sanity_check==FALSE)

# MACRO
macro <- read_sas("Data/MACRO/Macro.sas7bdat")
# look <- macro %>% 
#   filter(ID %in% 'A100008')
macro <- macro %>% 
  select(ID,dayofdeath) %>% 
  na.omit() 
macro <- macro %>% 
  mutate(ID = paste0(ID,'-','MACRO'))
macro.death.id <- macro$ID

copd <- copd %>% 
  mutate(death = ID %in% c(optimal.death.id,macro.death.id)) %>% 
  arrange(ID,event_date) %>% 
  mutate(event_time = event_date-Rand_Date)

# ASSUMPTION: cannot collect more data after a patient dies
# Hence, if a patient died, the max. period is the period in which the patient died

copd <- copd %>% 
  mutate(death=ifelse(period != max_period,FALSE,death)) 

# %>% 
#   mutate(recurrent_event = ifelse(period!=max_period,TRUE,FALSE))

copd <- copd %>% 
  arrange(ID,period) %>% 
  mutate(gaptime = ifelse(gaptime==0,event_time,gaptime)) %>% # sanity check
  mutate(death = as.numeric(death),
         recurrent_event=as.numeric(recurrent_event))

copd %>% 
  filter(death==1) %>% View()
# create the sum of the previous recurrent events
re.event.sum <- rep(NA,nrow(copd))
IDs <- copd$ID
IDs.lag <- lag(IDs)
ID.switch <- (IDs != IDs.lag)
ID.switch[1] <- TRUE
copd.re <- copd$recurrent_event
re.event.sum[1] <- 1
for(i in 2:nrow(copd)){
  if(ID.switch[i]){
    re.event.sum[i] <- ifelse(ID.switch[i],0,re.event.sum) # reset the counter for a new ID
    re.event.sum[i] <- ifelse(copd.re[i]==1,1,0)
  } else{
    re.event.sum[i] <- ifelse(copd.re[i]==1,re.event.sum[i-1]+1,re.event.sum[i-1])
  }
}
tmp <- copd %>% select(ID,recurrent_event)
tmp$num_re <- re.event.sum
tmp <- tmp %>% 
  mutate(num_re = ifelse(recurrent_event,num_re-1,num_re))
copd$num_recurrent_events <- tmp$num_re
copd <- copd %>% 
  select(ID,Days_In_Study,OPTIMAL,trtgroup,gender,age,BMI,
         nowsmk:fev1,fvc,fev1fvcratio,sgrq,ster1yr,packyears,
         Rand_Date,rand_month:rand_seasonality,
         period,event,event_date:event_seasonality,
         event_time,gaptime,recurrent_event,num_recurrent_events,death,max_period) 
# %>% 
  # mutate(death = ifelse(death==1 & event_time==365, 0, death))

# save the dataset
write_rds(copd,'SPPH501_copd_corrected.rds')


# visualization
####################################################
library(tidyverse)
copd_visual <- read_rds('SPPH501_copd_corrected.rds')
copd_summary_stats <- copd_visual
copd_summary_stats[,c(22,23,30:32)]  <- lapply(copd_summary_stats[,c(22,23,30:32)],as.factor)
summary(copd_summary_stats[,1:7] %>% distinct())[,-1] # basic
summary(copd_summary_stats[,c(1,8:17)] %>% distinct())[,-1] # lung
summary(copd_summary_stats[,c(1,18:21)] %>% distinct())[,-1] # randomization date
summary(copd_summary_stats[,22:27]) # event
summary(copd_summary_stats[,28:32])
table(copd_summary_stats$max_period)

library(reReg)
reObj <- with(copd_visual, Recur(time=event_time, id=ID, event=recurrent_event,death))
plot(reObj,cex=3,xlab= 'Time in days', ylab='Patients',
     main= '',
     terminal.name='Death',
     recurrent.name = 'Severe COPD Exacerbation')

plotEvents(reObj ~ trtgroup, data = copd_visual,
           mian='',
           cex=3,xlab= 'Time in days', ylab='Patients',
           terminal.name='Death',
           recurrent.name = 'Severe COPD Exacerbation')

plotEvents(reObj ~ gender, data = copd_visual,
           mian='',
           cex=3,xlab= 'Time in days', ylab='Patients',
           terminal.name='Death',
           recurrent.name = 'Severe COPD Exacerbation')

plotCSM(reObj~ gender +trtgroup, data = copd_visual, main = "",
        xlab='Time in days',
        ylab='Cumulative sample mean')

plotCSM(reObj~ age, data = copd_visual, main = "",
        xlab='Time in days',
        ylab='Cumulative sample mean')


library(ggplot2)
ggplot(data=copd_visual %>% filter(recurrent_event==1) %>%
         mutate(num_recurrent_events=ifelse(num_recurrent_events>3,3,num_recurrent_events)) %>% 
         mutate(num_recurrent_events=as.factor(num_recurrent_events)),aes(x=gaptime,colour=num_recurrent_events)) +
  geom_density() +
  xlab('time to get a subsequent severe COPD exacerbation (days)')

ggplot(data=copd_visual %>% filter(recurrent_event==1) %>%
         mutate(num_recurrent_events=ifelse(num_recurrent_events>3,3,num_recurrent_events)) %>% 
         mutate(num_recurrent_events=as.factor(num_recurrent_events)),aes(y=gaptime,x=fev1pp,colour=num_recurrent_events)) +
  geom_line() +
  xlab('time to get a subsequent severe COPD exacerbation (days)')

# second data processing for modeling
############################
# modeling
library(tidyverse)
library(frailtypack)
select <- dplyr::select
copd <- read_rds('SPPH501_copd_corrected.rds')
baseline.covariates <- c('trtgroup','age','gender','packyears','BMI',
                         'fev1','fev1pp', 'fev1fvcratio','fvc', 'stgtotalc00',
                         'nowsmk','oxygen','sgrq',
                         'OPTIMAL',
                         'rand_month','rand_season','rand_seasonality')

# perform some preprocessing
copd <- copd %>% 
  mutate(stgtotalc00 = stgtotalc00/100,
         fev1pp = fev1pp/100,
         sgrq = sgrq / 100) %>% 
  mutate(packyears = ifelse(packyears > age, packyears/12,packyears)) %>% 
  mutate(packyears= ifelse(age-packyears < 18, packyears/12,packyears)) %>% 
  mutate(packyears = packyears/100) %>% 
  mutate(age = age/100,
         gaptime = gaptime/365,
         BMI= BMI/100)

# time-varying covariates
timevarying.covariates <- c('num_recurrent_events','event_season','event_seasonality')

paste(c(baseline.covariates,timevarying.covariates),collapse='+')

write_rds(copd,'copd_modeling_corrected.rds')

# third data processing
######################
library(tidyverse)
library(frailtypack)
select <-dplyr::select

copd <- read_rds('copd_modeling_corrected.rds') %>% 
  filter(!(ID %in%  c('C100420-MACRO',' 2005-OPTIM'))) %>% 
  # mutate(death = ifelse(death==1 & event_time==365, 0, death)) %>%
  # mutate(recurrent_event = ifelse(event==0,0,recurrent_event)) %>% 
  mutate(num_recurrent_events = ifelse(num_recurrent_events>3,3,num_recurrent_events)) %>% 
  mutate(num_recurrent_events = as.factor(num_recurrent_events),
         num_re =num_recurrent_events) %>% 
  mutate(male = gender) %>% 
  select(-gender) %>% 
  mutate(trt = trtgroup)
# Salmeterol + Fluticasone
copd$trtgroup <- plyr::revalue(copd$trtgroup,c("1"="placebo","2"="Sal","3"="Sal_Flu"))
copd$trt <- copd$trtgroup

write_rds(copd,'copd_modeling_final_corrected.rds')
# garbage
######################
# copd %>%
#   filter(event_time==365 & death == 1) %>% View()
# 
# copd %>% 
#   filter(event!=0 & death == 1) 

# copd %>% dim()
# copd %>% filter(OPTIMAL==1) %>% dim()
# 
# rand.ID <- sample(unique(copd$ID),800)
# copd %>% 
#   filter(ID %in% rand.ID & death==1)
# copd %>% 
#   filter(ID %in% rand.ID) %>% 
#   dim()
# 

# copd.macro <- copd %>% 
#   filter(OPTIMAL != 1)
# 
# copd.optimal <- copd %>%
#   filter(OPTIMAL==1) %>%
#   dplyr::select(-OPTIMAL)


# use a shared frailty model on each of the processes to obtain a reasonable estimate
# of the smoothing parameters
# macro.ids <- unique(copd.macro$ID)
# max.period
# set.seed(2)
# 
# tmp <- copd.macro %>%
#   filter(ID %in% sample(macro.ids,400)) 
# summary(tmp)
# table(tmp$period)
# 
# tmp.ids <- unique(tmp$ID)
# 
# save.image('tmp.RData')
# load('tmp.RData')

# plotting 

# fitting frailty models
# name generation for SAS...
#########################
library(tidyverse)
library(frailtypack)
select <-dplyr::select
copd <- read_rds('copd_modeling_final_corrected.rds')
sha.gap <- frailtyPenal(Surv(gaptime,recurrent_event) ~ cluster(ID)+ 
                          trt +
                          age +
                          male +
                          # packyears +
                          BMI +
                          # fev1+
                          fev1pp+
                          fev1fvcratio+
                          # fvc+
                          # stgtotalc00+
                          nowsmk+
                          oxygen+
                          # sgrq+
                          OPTIMAL+
                          # rand_month+
                          rand_season+
                          # rand_seasonality +
                          num_re,
                        # event_seasonality,
                        data=copd,
                        n.knots = 10,
                        kappa=  1,
                        cross.validation = T,
                        maxit = 1000)
print(sha.gap)
plot(sha.gap)

sha.death <- frailtyPenal(Surv(time = gaptime,event=death) ~ cluster(ID)+ 
                            trt +
                            age +
                            male +
                            packyears +
                            BMI +
                            # fev1+
                            fev1pp+
                            fev1fvcratio+
                            # fvc+
                            nowsmk+
                            oxygen+
                            OPTIMAL+
                            num_re,
                          data=copd,
                          n.knots = 10,
                          kappa= 1,
                          cross.validation = T,
                          maxit = 1000)

print(sha.death)
plot(sha.death)

# prepping data for SAS
#########################
original.ID.values <- copd$ID %>% levels()
copd$ID <- plyr::mapvalues(copd$ID,from=original.ID.values,to=c(1:length(original.ID.values)))

sha.gap.model.matrix <- model.matrix(gaptime~recurrent_event+ 
                                       trt +
                                       age +
                                       male +
                                       packyears +
                                       BMI +
                                       # fev1+
                                       fev1pp+
                                       fev1fvcratio+
                                       # fvc+
                                       # stgtotalc00+
                                       nowsmk+
                                       oxygen+
                                       # sgrq+
                                       OPTIMAL+
                                       # rand_month+
                                       rand_season+
                                       # rand_seasonality +
                                       # tte0+
                                       num_re+
                                       event_season +
                                       ID,data=copd %>% 
                                       mutate(ID = as.numeric(ID)))


sha.death.model.matrix <- model.matrix(gaptime~death+
                                         # trt +
                                         age +
                                         male +
                                         packyears +
                                         BMI +
                                         # fev1+
                                         # fev1pp+
                                         # fev1fvcratio+
                                         # fvc+
                                         # stgtotalc00+
                                         nowsmk+
                                         oxygen+
                                         # sgrq+
                                         OPTIMAL+
                                         # rand_month+
                                         rand_season+
                                         # rand_seasonality +
                                         # tte0+
                                         num_re+
                                         # event_season +
                                         ID, data=copd %>% 
                                         mutate(ID = as.numeric(ID)))

match(colnames(sha.gap.model.matrix),colnames(sha.death.model.matrix))

match(colnames(sha.death.model.matrix),colnames(sha.gap.model.matrix)) #only death missing

sha.joint.model.matrix <- model.matrix(fev1~gaptime+recurrent_event+death+ 
                                         trt +
                                         age +
                                         male +
                                         packyears +
                                         BMI +
                                         fev1pp+
                                         fev1fvcratio+
                                         nowsmk+
                                         oxygen+
                                         OPTIMAL+
                                         rand_season+
                                         event_season+
                                         num_re+
                                         ID,data=copd %>% 
                                         mutate(ID = as.numeric(ID)))


sha.joint.model.matrix <- sha.joint.model.matrix %>% 
  as.data.frame() %>% 
  mutate(male = male1,
         nowsmk = nowsmk1,
         oxygen =oxygen1,
         OPTIMAL = OPTIMAL1) %>% 
  select(-c(male1,nowsmk1,oxygen1,OPTIMAL1)) %>% 
  select(-1)

write_csv(sha.joint.model.matrix,'copd_sas_corrected.csv')

d <- read_csv('copd_sas.csv')
cn <- colnames(d)
cn <- cn[-c(1:3,23)]
cn_d <- cn[-c(17:19)]

sha.gap.coef <- sha.gap$coef
names(sha.gap.coef)[c(4,9,10,11)] <- gsub('.{1}$', '', names(sha.gap.coef)[c(4,9,10,11)])
match(names(sha.gap.coef),cn)

sha.death.coef <- sha.death$coef
names(sha.death.coef)[c(4,9,10,11)] <- gsub('.{1}$', '', names(sha.death.coef)[c(4,9,10,11)])
match(names(sha.death.coef),cn_d)

cat(paste(paste0('b_',cn,'=',round(sha.gap.coef[match(cn,names(sha.gap.coef))],4)),collapse='\n'))

cat(paste(paste0('b_d_',cn,'=',round(sha.death.coef[match(cn_d,names(sha.death.coef))],4)),collapse='\n'))

paste(paste0('b_d_',cn),collapse='= ')

cat(paste(paste0(paste0('b_',cn),'*',cn),collapse='+\n'))

cat(paste(paste0(paste0('b_d_',cn),'*',cn),collapse='+\n'))

cat(paste(paste0('ESTIMATE ',"'re_",cn,"' ","b_",cn,'*','EXP(ln_gamma);'),collapse='\n'))
cat(paste(paste0('ESTIMATE ',"'d_",cn_d,"' ","b_d_",cn_d,'*','EXP(ln_gamma_d);'),collapse='\n'))


## parametric frailty model
library(parfm)
par.sha.gap <- parfm(Surv(gaptime,recurrent_event)~ 
                       trt +
                       age +
                       male +
                       packyears +
                       BMI +
                       # fev1+
                       fev1pp+
                       fev1fvcratio+
                       # fvc+
                       # stgtotalc00+
                       nowsmk+
                       oxygen+
                       # sgrq+
                       OPTIMAL+
                       # rand_month+
                       rand_season+
                       # rand_seasonality +
                       tte0+
                       num_re+
                       event_season,
                     # event_seasonality,
                     cluster="ID",frailty='gamma',dist='weibull',data=copd)


par.sha.death <- parfm(Surv(gaptime,death)~   trt +
                         age + 
                         male +
                         packyears +
                         BMI +
                         # fev1+
                         fev1pp+
                         fev1fvcratio+
                         # fvc+
                         # stgtotalc00+
                         nowsmk+
                         oxygen+
                         # sgrq+
                         OPTIMAL+
                         # rand_month+
                         rand_season+
                         # rand_seasonality +
                         tte0+
                         num_re,
                       cluster="ID",frailty='gamma',dist='weibull',data=copd)
##############################

# joint modeling in R


# garbage
#######################
# copd.ids <- unique(copd$ID)
# death.ids <- copd %>% filter(death==1) %>% dplyr::select(ID) %>% unlist() %>% as.vector()
# set.seed(2)
# tmp <- copd %>%
#   filter(max_period !=1) %>% 
#   filter(ID %in% c(sample(copd.ids,500),death.ids)) %>% 
#   mutate(id=ID)

# data("readmission")
# glimpse(readmission) %>% View()
# readmission %>% 
#   group_by(id) %>% 
#   tally() -> look
# table(look$n)
# 
# tmp %>% dplyr::select(ID,period,event_time,gaptime,death,recurrent_event) %>% View()
# library(haven)
# macro <- read_sas("Data/MACRO/Macro.sas7bdat")
# macro %>% 
#   filter(ID %in% 'I101279')


###################
library(tidyverse)
library(frailtypack)
select <- dplyr::select
copd <- read_rds('copd_modeling_final_corrected.rds')

sha.gap <- frailtyPenal(Surv(gaptime,recurrent_event) ~ cluster(ID)+ 
                          trt +
                          # age +
                          male +
                          # packyears +
                          # BMI +
                          fev1pp+
                          # fev1fvcratio+
                          # nowsmk+
                          oxygen+
                          # OPTIMAL+
                          rand_season+
                          event_season+
                          num_re,
                        data=copd,
                        n.knots = 8,
                        kappa=  10000,
                        cross.validation = T,
                        maxit = 1000)
print(sha.gap)
plot(sha.gap)

sha.death <- frailtyPenal(Surv(time = gaptime,event=death) ~ cluster(ID)+ 
                            # trt +
                            age +
                            male +
                            packyears +
                            # BMI +
                            # fev1+
                            fev1pp+
                            # fev1fvcratio+
                            # fvc+
                            nowsmk+
                            oxygen+
                            OPTIMAL+
                            num_re,
                          data=copd,
                          n.knots = 8,
                          kappa=   1,
                          cross.validation = T,
                          maxit = 1000)
print(sha.death)
plot(sha.death)

joint.frailty <- frailtyPenal(formula = Surv(time=gaptime,event=recurrent_event) ~ cluster(ID)+ 
                                trt +
                                # age +
                                male +
                                # packyears +
                                # BMI +
                                fev1pp+
                                # fev1fvcratio+
                                # nowsmk+
                                oxygen+
                                # OPTIMAL+
                                rand_season+
                                event_season+
                                num_re+
                                terminal(death),
                              formula.terminalEvent = ~
                                # trt +
                                age +
                                male +
                                packyears +
                                # BMI +
                                # fev1+
                                fev1pp+
                                # fev1fvcratio+
                                # fvc+
                                nowsmk+
                                oxygen+
                                OPTIMAL+
                                num_re,
                              data=copd,
                              n.knots = 8,
                              init.B=c(sha.gap$coef,sha.death$coef),
                              # kappa=c(1000,1000),
                              # kappa=c(1,1),
                              kappa=c(sha.gap$kappa,sha.death$kappa),
                              jointGeneral = FALSE,
                              recurrentAG=FALSE, # gap-time
                              maxit=250)
print(joint.frailty)
plot(joint.frailty)

# data(dataNCC)
# head(dataNCC)
# modJoint.ncc <- frailtyPenal(Surv(t.start,t.stop,event)~cluster(id)+cov1
#                              +cov2+terminal(death)+wts(ncc.wts), formula.terminalEvent=~cov1+cov2,
#                              data=dataNCC,n.knots=8,kappa=c(1.6e+10, 5.0e+03),recurrentAG=TRUE, RandDist="LogN") 
# summary(modJoint.ncc)
##################################

# grabage
# try frailtypack on a tody dataset
#######################
#toy example
data("readmission")
readmission %>% glimpse()
# #-- Gap-time
# mod.sha.gap <- frailtyPenal(Surv(time, event) ~ cluster(id) + dukes +
#
#                               charlson + sex + chemo, n.knots = 9, kappa = 1,
#
#                               data = readmission, cross.validation = TRUE)
# mod.sha.gap$kappa/1e10
#
# plot(mod.sha.gap)

# modJoint.gap <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+charlson+
#                                terminal(death),formula.terminalEvent=~sex+dukes+charlson,
#                              data=readmission,n.knots=14,kappa=c(9.55e+9,1.41e+12),
#                              recurrentAG=FALSE,jointGeneral = T)
# 
# modJoint.gap2 <- frailtyPenal(Surv(time,event)~cluster(id)+sex+dukes+charlson+
#                                terminal(death),formula.terminalEvent=~sex+dukes+charlson,
#                              data=readmission,n.knots=10,kappa=c(9.55e+9,1.41e+12),
#                              recurrentAG=FALSE,jointGeneral = T)