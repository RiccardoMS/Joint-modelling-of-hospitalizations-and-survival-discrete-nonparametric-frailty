############################################################################################
################################## Classical Cox Analysis ##################################
############################################################################################
## [1] Data preparation
## [2] Example plot
## [3] Summary statistics
## [4] KM curves
## [5] Log rank tests
## [6] Cox PH models

## clear workspace
rm(list = ls())

## load packages
library(data.table)
library(survival)
library(survminer)
library(ggplot2)

## load dataset
load("With_Adherence_Dataset_FFU/ACE_Inhibitors.RData")

## arrange Dataset
# total number of hospitalizations
new[, tot_hosp:=max(hosp,na.rm = TRUE), by=COD_REG]

# number of hospitalizations within 1st year of follow up
new[, tot_hosp_1st:=max(hosp[which(data_prest-data_rif_ev<=365)],na.rm = TRUE), by=COD_REG]

# number of comorbidities at first discharge
new[!is.na(hosp), comorbidity:=rowSums(.SD), .SDcols = 36:55]

# max number of comorbidities observed within 1st year of follow up
new[, max_como_1st:=max(comorbidity[which(data_prest-data_rif_ev<=365)],na.rm = TRUE), by=COD_REG]

# rearrange dataset 
data <- subset(new,hosp==1,select = c(COD_REG,death,timeOUT,SESSO,eta_Min,tot_hosp_1st,tot_hosp,ADERENTE,ADH_LEV,comorbidity, max_como_1st))

# arrange variables
data$COD_REG= factor(data$COD_REG)
data$status=factor(data$death, labels=c('censored','death'))
data$ADERENTE=factor(data$ADERENTE)
data$ADH_LEV= factor(data$ADH_LEV)
data$timeOUT = data$timeOUT - 365

## Example plot
data_head=head(data)
data_head$time_y=data_head$timeOUT/365
ggplot(data=data_head,aes(x=COD_REG,y=as.double(time_y))) +
      geom_bar(stat='identity',width=0.05)+
      geom_point(aes(color=status,shape=status),size=3)+
      xlab("ID")+
      ylab("years")+
      coord_flip()

## Summary statistics for data
summary(data)

## Kaplan_meier curves
KM <- survfit(Surv(timeOUT,death==1)~1,data=data)
summary(KM)
# Survival Probability
ggsurvplot(KM,data=data,risk.table=TRUE,risk.table.col="strata",legend.title="",
           ggtheme = theme_minimal(),break.time.by=360,censor=F)
# Cumulative Incidence
ggsurvplot(KM,data=data,risk.table=TRUE,risk.table.col="strata",
           ggtheme = theme_minimal(),break.time.by=360,fun='event',censor=F)
# Cumulative Hazard
ggsurvplot(KM,data=data,risk.table=TRUE,risk.table.col="strata",
           ggtheme = theme_minimal(),break.time.by=360,fun='cumhaz',censor=F)

# Investigate some factors
# Sex
KM.sex <- survfit(Surv(timeOUT,death==1)~SESSO,data=data)
summary(KM.sex)
ggsurvplot(KM.sex,data=data,conf.int=T,risk.table=TRUE,risk.table.col="strata",
           ggtheme = theme_minimal(),break.time.by=360,legend.labs=c('Female','Male'),legend.title="",
           palette=c('orchid2','dodgerblue2'),title="",pval = T,censor=F)

# Adherence
KM.ADH <- survfit(Surv(timeOUT,death==1)~ADERENTE,data=data)
summary(KM.ADH)
ggsurvplot(KM.ADH,data=data,conf.int=T,risk.table=TRUE,risk.table.col="strata",
           ggtheme = theme_minimal(),break.time.by=360,legend.labs=c('Non Aderente','Aderente'),legend.title="",
           palette=c('red','green'),title="",pval = T,censor=F)

# Adherence Levels
KM.ADH_LEV <- survfit(Surv(timeOUT,death==1)~ADH_LEV,data=data)
summary(KM.ADH)
ggsurvplot(KM.ADH_LEV,data=data,conf.int=T,risk.table=TRUE,risk.table.col="strata",
           ggtheme = theme_minimal(),break.time.by=360,legend.labs=c('[0.00;0.25)','[0.25;0.50)','[0.50;0.75)','[0.75;1.00]'),legend.title="",
           palette='heat',title="",censor=F,risk.table.height=0.35,pval=T)

survdiff(Surv(timeOUT,death==1)~ADH_LEV,data=data)
HR_12=(115/88.2)/(99/93.4)
HR_23=(99/93.4)/(160/122.5)
HR_34=(160/122.5)/(344/413.9)
HR_14=(115/88.2)/(344/413.9)
HR_KM_LEVADH <- c(HR_12,HR_23,HR_34,HR_14)

## LOG-RANK TESTS
# sex 
survdiff(Surv(timeOUT,death==1)~SESSO,data=data)
# Adherence
survdiff(Surv(timeOUT,death==1)~ADERENTE,data=data)
# AdhLev
survdiff(Surv(timeOUT,death==1)~ADH_LEV,data=data)
adh_lev_pairwise = ifelse(data$ADH_LEV==4,1,0)
survdiff(Surv(timeOUT,death==1)~adh_lev_pairwise,data=data)
# Eta_min
eta_Min_fac=cut(data$eta_Min,breaks = c(17,66,80,98))
survdiff(Surv(timeOUT,death==1)~eta_Min_fac,data=data)
# comorbidity
como_fac=cut(data$comorbidity,breaks = c(-1,3,7))
survdiff(Surv(timeOUT,death==1)~como_fac,data=data)
como_fac=cut(data$comorbidity,breaks = c(-1,1,3,7))
survdiff(Surv(timeOUT,death==1)~como_fac,data=data)
# maxcomorbidity
max_como_fac=cut(data$max_como_1st,breaks = c(-1,1,10))
survdiff(Surv(timeOUT,death==1)~max_como_fac,data=data)
max_como_fac=cut(data$max_como_1st,breaks = c(-1,1,3,10))
survdiff(Surv(timeOUT,death==1)~max_como_fac,data=data)
# tot_hosp_1Y
hosp_1Y_fac=cut(data$tot_hosp_1st,breaks = c(0,3,14))
survdiff(Surv(timeOUT,death==1)~hosp_1Y_fac,data=data)

## COX PH: full model and variable selection
# Adherence and maxComo
cox.mod.1 <- coxph(Surv(timeOUT,death)~SESSO + eta_Min + ADERENTE + tot_hosp_1st  + max_como_1st,data=data)
summary(cox.mod.1)
extractAIC(cox.mod.1)
# Adherence and comorbidity
cox.mod.2 <- coxph(Surv(timeOUT,death)~SESSO + eta_Min + ADERENTE + tot_hosp_1st + comorbidity,data=data)
summary(cox.mod.2)
extractAIC(cox.mod.2)


#final
cox.mod <- coxph(Surv(timeOUT,death)~SESSO + eta_Min + ADERENTE + tot_hosp_1st + max_como_1st,data=data)
summary(cox.mod)

#full patient history
cox.mod.tothosp <- coxph(Surv(timeOUT,death)~SESSO + eta_Min + ADERENTE + tot_hosp + max_como_1st,data=data)
summary(cox.mod.tothosp)

# Hazard ratios
x11()
ggforest(cox.mod,data=data)

# Baseline Hazard
x11()
plot(survfit(cox.mod,data=data),col = 'grey',xlab = "time",ylab = "Probability")

# Sex effect
sex_df <-with(data,
              data.frame(SESSO=c('F','M'),
                         eta_Min=rep(median(data$eta_Min,na.rm=T),2),
                         ADERENTE=factor(rep(1,2)),
                         tot_hosp_1st=rep(median(data$tot_hosp_1st, na.rm=T),2),
                         max_como_1st=rep(median(data$max_como_1st,na.rm = T),2)))
cox_sex <- survfit(cox.mod,newdata=sex_df)
ggsurvplot(cox_sex,data=sex_df, conf.int = T,legend.labs=c("Female","Male"),
           ggtheme = theme_minimal(),break.time.by =365, xlab =" Time ( days )", censor =F)

# Adherence Effect
adh_df <-with(data,
              data.frame(SESSO=factor(rep('F',2)),
                         eta_Min=rep(median(data$eta_Min,na.rm=T),2),
                         ADERENTE=factor(c(1,0)),
                         tot_hosp_1st=rep(median(data$tot_hosp_1st, na.rm=T),2),
                         max_como_1st=rep(median(data$max_como_1st,na.rm = T),2)))
cox_adh <- survfit(cox.mod,newdata=adh_df)
ggsurvplot(cox_adh, data=adh_df, conf.int = T,legend.labs=c("Adherent","Non adherent"),
           ggtheme = theme_minimal(),break.time.by =365, xlab =" Time ( days )", censor =F,
           palette = c('green','red'))

# Age_Min effect
age_df <-with(data,
              data.frame(SESSO=factor(rep('F',9)),
                         eta_Min=floor(seq(18,98,length.out=9)),
                         ADERENTE=rep(1,9),
                         tot_hosp_1st=rep(median(data$tot_hosp_1st, na.rm=T),9),
                         max_como_1st=rep(median(data$max_como_1st,na.rm = T),9)))
age_df$ADERENTE <- factor(age_df$ADERENTE)
cox_age<- survfit(cox.mod,newdata=age_df)
ggsurvplot(cox_age, data=age_df, conf.int = F,legend.labs=levels(factor(age_df$eta_Min)),
           ggtheme = theme_minimal(),break.time.by =365, xlab =" Time ( days )", censor =F,
           palette = "Blues")

# Tot_hosp_1st effect
hosp_df <-with(data,
              data.frame(SESSO=factor(rep('F',9)),
                         eta_Min=rep(median(data$eta_Min,na.rm=T),9),
                         ADERENTE=rep(1,9),
                         tot_hosp_1st=floor(seq(1,14,length.out=9)),
                         max_como_1st=rep(median(data$max_como_1st,na.rm = T),9)))
hosp_df$ADERENTE <- factor(hosp_df$ADERENTE)
cox_hosp<- survfit(cox.mod,newdata=hosp_df)
ggsurvplot(cox_hosp, data=hosp_df, conf.int = F,legend.labs=levels(factor(hosp_df$tot_hosp_1st)),
           ggtheme = theme_minimal(),break.time.by =365, xlab =" Time ( days )", censor =F,
           palette = "Reds")

# max_como_1st effect
como_df <-with(data,
               data.frame(SESSO=factor(rep('F',9)),
                          eta_Min=rep(median(data$eta_Min,na.rm=T),9),
                          ADERENTE=rep(1,9),
                          tot_hosp_1st=rep(median(data$tot_hosp_1st, na.rm=T),9),
                          max_como_1st=floor(seq(0,10,length.out=9))))
como_df$ADERENTE <- factor(como_df$ADERENTE)
cox_como<- survfit(cox.mod,newdata=como_df)
ggsurvplot(cox_como, data=como_df, conf.int = F,legend.labs=levels(factor(como_df$max_como_1st)),
           ggtheme = theme_minimal(),break.time.by =365, xlab =" Time ( days )", censor =F,
           palette = "Greens")


# Diagnostics
x11()
ggcoxdiagnostics(cox.mod,type='martingale')
ggcoxdiagnostics(cox.mod, type='deviance')
ggcoxdiagnostics(cox.mod, type='schoenfeld')
