# visualization
####################################################
library(tidyverse)
copd_visual <- read_rds('SPPH501_copd.rds')
copd_summary_stats <- copd_visual
copd_summary_stats[,c(22,23,30:32)]  <- lapply(copd_summary_stats[,c(22,23,30:32)],as.factor)
summary(copd_summary_stats[,1:7] %>% distinct())[,-1] # basic
summary(copd_summary_stats[,c(1,8:17)] %>% distinct())[,-1] # lung
summary(copd_summary_stats[,c(1,18:21)] %>% distinct())[,-1] # randomization date
summary(copd_summary_stats[,22:27]) # event
summary(copd_summary_stats[,28:33])

library(reReg)

copd_visual <- copd_visual %>% 
  filter(!( period < max_period & recurrent_event==0))

reObj <- with(copd_visual, Recur(time=event_time, id=ID, event=recurrent_event,death))
plot(reObj)

plot(reObj,cex=1.5,xlab= 'Time in days', ylab='Patients',
     main= 'Event plot',
     terminal.name='Death',
     recurrent.name = 'Hospital readmission')

data(readmission, package = "frailtypack")
readmission <- subset(readmission, !(id %in% c(60, 109, 280)))
tmp2 <- readmission[1:3,] %>% select(id,t.stop,event,death)
attach(tmp2)
reObj2 <- Recur(t.stop, id, event, death)
plot(reObj2)
detach(tmp2)

