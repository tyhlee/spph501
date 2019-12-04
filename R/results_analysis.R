library(tidyverse)
library(frailtypack)
library(parfm)
select <-dplyr::select
copd <- read_rds('copd_modeling_final_corrected.rds')
sha.gap <- frailtyPenal(Surv(gaptime,recurrent_event) ~ cluster(ID)+ 
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
                          num_re+
                          event_season,
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

sha.gap.weibull <- parfm(Surv(gaptime,recurrent_event) ~ 
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
                          num_re+
                          event_season,
                          cluster='ID',
                          data=copd,
                          dist='weibull',
                          frailty='gamma')

sha.death.weibull <- parfm(Surv(gaptime,death) ~ 
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
                           cluster='ID',
                           data=copd,
                           dist='weibull',
                           frailty='gamma')

###########################
# RE
tmp.name <- names(sha.gap$coef) 
tmp.name[match(c('male1','nowsmk1','oxygen1','OPTIMAL1'),tmp.name)] <- c('male','nowsmk','oxygen','OPTIMAL')
tmp1 <- data.frame(Variable=tmp.name,
                   cbind(exp(sha.gap$coef),
                         exp(sha.gap$coef-qnorm(1-0.05/2)*(sqrt(diag(sha.gap$varH)))),
                         exp(sha.gap$coef+qnorm(1-0.05/2)*(sqrt(diag(sha.gap$varH))))),
                   baseline='frailtypack-splines')

tmp2 <- data.frame(Variable=tmp.name,
                   cbind(exp(coef(sha.gap.weibull)),parfm::ci.parfm(sha.gap.weibull)),
                   baseline='parfm-Weibull')

look <- read_csv("DEC/__ADD_PARAM_RE_GAMMA.csv")
look <- look[-c(1,2),]
matching <- sapply(tmp.name,function(x)grep(x ,substring(look$Label,3)))
matching[[1]] <- 1 # manual
matching <- matching %>% unlist()
tmp <- data.frame(Variable = tmp.name,
                  Estimate = 1,
                  Lower = 1,
                  Upper= 1,
                  Baseline = 'NLMIXED-AFT_Weibull')

tmp[-which(is.na(match(tmp.name,names(matching)))),c(-1,-5)] <- (cbind(exp(look$Estimate),
                                                                       exp(look$Estimate-qnorm(1-0.05/2)*look$StandardError),
                                                                       exp(look$Estimate+qnorm(1-0.05/2)*look$StandardError)))[matching,]


names(tmp) <- names(tmp1) <- names(tmp2) <- c("Variable","Estimate","Lower","Upper","Baseline")
tmp3 <- rbind(tmp,tmp1,tmp2)

# names(tmp1) <- names(tmp2) <- c("Variable","Estimate","Lower","Upper","Baseline")
# tmp3 <- rbind(tmp1,tmp2)

p.re = ggplot(tmp3, 
           aes(Estimate,Variable,xmax=Upper,xmin=Lower)) + 
  coord_cartesian(xlim=c(0,15 )) +
  #scale_alpha_discrete(range = c(0.2, 1)) +
  geom_vline(xintercept = 1.0, linetype=2, alpha=0.75) +
  geom_errorbarh(alpha=0.5, color='black') + 
  # geom_point(aes(size = RelConf)) +
  geom_point(aes(Estimate),size=3)+
  # geom_point(data = subset(remresdf, Study=="RE Model"), size=7) +
  scale_size(range = c(2, 5), guide=FALSE) +
  theme_bw() + 
  theme(text = element_text(size=16))+
  facet_grid(~Baseline) +
  ylab("") +
  xlab("HR")

# summary(sha.gap);parfm::ci.parfm(sha.gap.weibull)

frailty.name <- data.frame(Model=c('NLMIXED-AFT_Weibull','frailtypack-splines','parfm-Weibull'))
frailty.re <- cbind(frailty.name,rbind(read_csv("DEC/__PARAM_RE_GAMMA.csv")[22,c(2,3)],
                       c(sha.gap$theta,sqrt(sha.gap$varTheta[1])),
                       sha.gap.weibull[which(rownames(sha.gap.weibull)=='theta'),c(1,2)]))
xtable::xtable(frailty.re)
# saveRDS(frailty.re,'frailty_re.rds')
###############################
# DEATH
tmp.name <- names(sha.death$coef) 
tmp.name[match(c('male1','nowsmk1','oxygen1','OPTIMAL1'),tmp.name)] <- c('male','nowsmk','oxygen','OPTIMAL')
tmp1 <- data.frame(Variable=tmp.name,
                   cbind(exp(sha.death$coef),
      exp(sha.death$coef-qnorm(1-0.05/2)*(sqrt(diag(sha.death$varH)))),
      exp(sha.death$coef+qnorm(1-0.05/2)*(sqrt(diag(sha.death$varH))))),
      baseline='frailtypack-splines')

tmp2 <- data.frame(Variable=tmp.name,
                   cbind(exp(coef(sha.death.weibull)),parfm::ci.parfm(sha.death.weibull)),
           baseline='parfm-Weibull')

look <- read_csv("DEC/__ADD_PARAM_DEATH_GAMMA.csv")
look <- look[-c(1,2,16),]
# look2 <- read_csv("Nov_corrected/__ADD_PARAM_DEATH_GAMMA.csv")
matching <- sapply(tmp.name,function(x)grep(x ,substring(look$Parameter,5)))
matching[[1]] <- 1 # manual
matching <- matching %>% unlist()
tmp <- data.frame(Variable = tmp.name,
                  Estimate = 1,
                  Lower = 1,
                  Upper= 1,
                  Baseline = 'NLMIXED-AFT_Weibull')

tmp[-which(is.na(match(tmp.name,names(matching)))),c(-1,-5)] <- (cbind(exp(look$Estimate),
                                                              exp(look$Estimate-qnorm(1-0.05/2)*look$StandardError),
                                                              exp(look$Estimate+qnorm(1-0.05/2)*look$StandardError)))[matching,]


names(tmp) <- names(tmp1) <- names(tmp2) <- c('Variable',"Estimate","Lower","Upper","Baseline")
tmp3 <- rbind(tmp,tmp1,tmp2)
# names(tmp1) <- names(tmp2) <- c('Variable',"Estimate","Lower","Upper","Baseline")
# tmp3 <- rbind(tmp1,tmp2)

p.death = ggplot(tmp3, 
           aes(Estimate,Variable,xmax=Upper,xmin=Lower)) + 
  coord_cartesian(xlim=c(0,20 )) +
  #scale_alpha_discrete(range = c(0.2, 1)) +
  geom_vline(xintercept = 1.0, linetype=2, alpha=0.75) +
  geom_errorbarh(alpha=0.5, color='black') + 
  # geom_point(aes(size = RelConf)) +
  geom_point(aes(Estimate),size=3)+
  # geom_point(data = subset(remresdf, Study=="RE Model"), size=7) +
  scale_size(range = c(2, 5), guide=FALSE) +
  theme_bw() + 
  theme(text = element_text(size=16))+
  facet_grid(~Baseline) +
  ylab("") +
  xlab("HR")

frailty.name <- data.frame(Model=c('NLMIXED-AFT_Weibull','frailtypack-splines','parfm-Weibull'))

look <- read_csv("DEC/__PARAM_DEATH_GAMMA.csv")

frailty.death <- cbind(frailty.name,rbind(look[which(look$Parameter=='theta'),c(2,3)],
                                       c(sha.death$theta,sqrt(sha.death$varTheta[1])),
                                       sha.death.weibull[which(rownames(sha.death.weibull)=='theta'),c(1,2)]))
xtable::xtable(frailty.death)
# saveRDS(frailty.death,'frailty_death.rds')

######################################
# Joint

look <- read_csv("Nov_corrected/__JOINT_ADDITIONAL.csv")
colnames(look)[1] <- 'Parameter'
# look.death <- look[grep('^d_*',look$Parameter),]
# look.re <- look[grep('^re_*',look$Parameter)[which((grep('b_*',look$Parameter) %in% grep('b_d_*',look$Parameter))==FALSE)],]
look.death <- look[grep('^d_*',look$Parameter),]
look.re <- rbind(look[grep('^re_*',look$Parameter),],
                 look[grep('^CON*',look$Parameter),])
look <- look.re

tmp <- data.frame(Variable = look$Parameter,
                  Estimate = 1,
                  Lower = 1,
                  Upper= 1,
                  Baseline = 'joint-weibull-normal')

tmp[,c(2,3,4)] <- (cbind(exp(look$Estimate),
                         exp(look$Estimate-qnorm(1-0.05/2)*look$StandardError),
                         exp(look$Estimate+qnorm(1-0.05/2)*look$StandardError)))

# tmp$Variable <- substring(tmp$Variable,3)
tmp$Variable <- gsub('re_',"",tmp$Variable)

p.re.joint <- ggplot(tmp, 
                        aes(Estimate,Variable,xmax=Upper,xmin=Lower)) + 
  coord_cartesian(xlim=c(0,10 )) +
  #scale_alpha_discrete(range = c(0.2, 1)) +
  geom_vline(xintercept = 1.0, linetype=2, alpha=0.75) +
  geom_errorbarh(alpha=0.5, color='black') + 
  # geom_point(aes(size = RelConf)) +
  geom_point(aes(Estimate),size=3)+
  # geom_point(data = subset(remresdf, Study=="RE Model"), size=7) +
  scale_size(range = c(2, 5), guide=FALSE) +
  theme_bw() + 
  theme(text = element_text(size=16))+
  ylab("") +
  xlab("HR")

# tmp <- look.re[grep('b_num_re*',look.re$Parameter),]
# tmp.lag <- lag(tmp$Estimate)
# tmp.lag[1] <- 0
# re <-exp(tmp$Estimate- tmp.lag)

look <- look.death

tmp <- data.frame(Variable = look$Parameter,
                  Estimate = 1,
                  Lower = 1,
                  Upper= 1,
                  Baseline = 'joint-weibull-normal')
tmp[,c(2,3,4)] <- (cbind(exp(look$Estimate),
                         exp(look$Estimate-qnorm(1-0.05/2)*look$StandardError),
                         exp(look$Estimate+qnorm(1-0.05/2)*look$StandardError)))

# tmp$Variable <- substring(tmp$Variable,3)
tmp$Variable <- gsub('d_',"",tmp$Variable)

p.death.joint <- ggplot(tmp, 
       aes(Estimate,Variable,xmax=Upper,xmin=Lower)) + 
  coord_cartesian(xlim=c(0,10 )) +
  #scale_alpha_discrete(range = c(0.2, 1)) +
  geom_vline(xintercept = 1.0, linetype=2, alpha=0.75) +
  geom_errorbarh(alpha=0.5, color='black') + 
  # geom_point(aes(size = RelConf)) +
  geom_point(aes(Estimate),size=3)+
  # geom_point(data = subset(remresdf, Study=="RE Model"), size=7) +
  scale_size(range = c(2, 5), guide=FALSE) +
  theme_bw() + 
  theme(text = element_text(size=16))+
  ylab("") +
  xlab("HR")

look <- read_csv("Nov_corrected/__JOINT_ADDITIONAL.csv")

frailty.joint <- look[look$Label %in% c('k','v'),c(1:3)]
xtable::xtable(frailty.joint)
# tmp <- look.death[grep('b_d_num_*',look.death$Parameter),]
# tmp.lag <- lag(tmp$Estimate)
# tmp.lag[1] <- 0
# death <-exp(tmp$Estimate- tmp.lag)
# 
# re.death <- rbind(re,death)
# colnames(re.death) <- c("1st","2nd","3rd")

rm(look,look.death,look.re,tmp,tmp1,tmp2,tmp3)

# random effects
copd <- copd %>% 
  mutate(ID.char = as.character(ID))

copd.ID1 <- copd %>% 
  mutate(num_re = as.numeric(num_re)) %>% 
  filter(num_re > 1) %>% 
  select(ID.char,num_re) %>% 
  select(ID.char) %>% 
  unlist() %>% 
  unique()

IDs <- match(copd.ID1,(copd$ID.char %>% unique()))
# Mon Dec  2 03:58:00 2019 ------------------------------

tmp <- read_csv("DEC/__Z_RE.csv")
gg <- tmp[,'Estimate'] %>% unlist()
# gg <- ifelse(g>0.99999,0.99999,g)
# gg <- pnorm(g)
# gg <- qgamma(p = gg,shape = 1/2.6,scale=2.6)
h <- hist(gg, breaks = 20, density = 10,
          col = "lightgray", xlab = "Accuracy", main = "Overall") 
xfit <- seq(min(gg), max(gg), length = 100) 
yfit <- dnorm(xfit, mean=0,sd=1)
# yfit <- dgamma(xfit, shape=1/2.6,scale=2.6)
yfit <- yfit * diff(h$mids[1:2]) * length(g) 
lines(xfit, yfit, col = "black", lwd = 2)


tmp <- read_csv("Nov_corrected/__JOINT_Z.csv")
tmp <- read_csv("DEC/__Z_RE_GAMMA.csv")

# p=cdf("NORMAL",z); */
#   /* 	if p > .99999 then p=0.99999; */
#   /* 	g2=quantile('GAMMA',p,1/theta); */
#   /* 	g=g2 * theta; */
g <- tmp[IDs,'Estimate'] %>% unlist()
gg <- ifelse(g>0.99999,0.99999,g)
gg <- pnorm(g)
gg <- qgamma(p = gg,shape = 1/2.6,scale=2.6)
h <- hist(gg, breaks = 20, density = 10,
          col = "lightgray", xlab = "Accuracy", main = "Overall") 
# xfit <- seq(0, max(gg), length = 1000) 
xfit <- seq(0, 10, length = 1000) 
# yfit <- dnorm(xfit, mean=0,sd=1)
yfit <- dgamma(xfit, shape=1/2.6,scale=2.6)
yfit <- yfit * diff(h$mids[1:2]) * length(g) 
lines(xfit, yfit, col = "black", lwd = 2)

save.image('result_analysis.RData')
