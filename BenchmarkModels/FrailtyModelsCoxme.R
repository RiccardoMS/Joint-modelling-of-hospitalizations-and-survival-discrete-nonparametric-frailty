############################################################################################
################################ Cox Mixed Effect Models ###################################
############################################################################################
## [1] Data preparation
## [2] Cox models
## [3] Coxme models
## [4] Plot 6: HR comparison

## clean workspace
rm(list = ls())

## load packages
library(data.table)
library(survival)
library(survminer)
library(coxme)

## load data
load("dataRec.RData")
load("dataDeath.RData")

## Arrange variables
data$COD_REG= factor(data$COD_REG)
data$SESSO=factor(data$SESSO)
data$ADERENTE=factor(data$ADERENTE)
data$etaEvent=as.double(data$etaEvent)
data$status = factor('na')
for( i in 1:dim(data)[1]){
  if(data[i]$event==1)
    data[i]$status="hospitalization"
  else {
    if(data[i]$cens==0)
      data[i]$status="censored"
    else
      data[i]$status="dead"
  }
}
dataDeath$timeEvent<-dataDeath$timeEvent - 365

## Cox PH - Recurrent Events
Cox.Rec <- coxph(Surv(GapEvent,event)~ SESSO + ADERENTE + etaEvent + comorbidity, data=data)
summary(Cox.Rec)
print(Cox.Rec)

## Cox PH - Terminal Events
Cox.Death <- coxph(Surv(GapEvent,cens)~ SESSO + ADERENTE + etaEvent + comorbidity, data=dataDeath)
summary(Cox.Death)
print(Cox.Death)

## Diagnostics
ggcoxdiagnostics(Cox.Rec,type='martingale')
ggcoxdiagnostics(Cox.Rec,type='schoenfeld')
ggcoxdiagnostics(Cox.Death, type='martingale')
ggcoxdiagnostics(Cox.Death,type='schoenfeld')

## Coxme - Recurrent Events
coxme.Rec <- coxme(Surv(GapEvent,event) ~ SESSO + ADERENTE + etaEvent + comorbidity +(1|COD_REG), data = data)
coxme.Rec

## Coxme - Terminal Events Events
coxme.Death <- coxme(Surv(GapEvent,cens) ~ SESSO + ADERENTE + etaEvent + comorbidity + (1|COD_REG), data = dataDeath)
coxme.Death

## Plot: HR CI comparison Cox vs Coxme
df_cox_rec <- data.frame(
  x   = seq(0,0.5,length.out=length(Cox.Rec$coefficients))-0.05,
  l95 = exp(Cox.Rec$coefficients-1.96*sqrt(diag(Cox.Rec$var))),
  hr  = exp(Cox.Rec$coefficients),
  u95 = exp(Cox.Rec$coefficients+1.96*sqrt(diag(Cox.Rec$var))))

df_cox_death <- data.frame(
  x   = seq(0,0.5,length.out=length(Cox.Death$coefficients))-0.05,
  l95 = exp(Cox.Death$coefficients-1.96*sqrt(diag(Cox.Death$var))),
  hr  = exp(Cox.Death$coefficients),
  u95 = exp(Cox.Death$coefficients+1.96*sqrt(diag(Cox.Death$var))))

df_coxme_rec <- data.frame(
  x   = seq(0,0.5,length.out=length(coxme.Rec$coefficients))+0.05,
  l95 = exp(coxme.Rec$coefficients-1.96*sqrt(diag(vcov(coxme.Rec)))),
  hr  = exp(coxme.Rec$coefficients),
  u95 = exp(coxme.Rec$coefficients+1.96*sqrt(diag(vcov(coxme.Rec)))))

df_coxme_death <- data.frame(
  x   = seq(0,0.5,length.out=length(coxme.Death$coefficients))+0.05,
  l95 = exp(coxme.Death$coefficients-1.96*sqrt(diag(vcov(coxme.Death)))),
  hr  = exp(coxme.Death$coefficients),
  u95 = exp(coxme.Death$coefficients+1.96*sqrt(diag(vcov(coxme.Death)))))

# Rec
library(ggpubr)
g1<-ggplot(data=df_cox_rec[1,]) +
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#F8766D",size=1)+
  geom_point(aes(x=x,y=hr),color="#F8766D",size=3)+
  geom_segment(data=df_coxme_rec[1,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_rec[1,], aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_rec[1,], aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#00BFC4",size=1)+
  geom_point(data=df_coxme_rec[1,], aes(x=x,y=hr),color="#00BFC4",size=3)+
  theme(legend.position="none")+
  scale_x_continuous(name="",breaks = c(0.0),labels=c("Sex [M]"), limits = c(-0.1,0.1))+
  ylab("")
g2<-ggplot(data=df_cox_rec[2,]) +
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#F8766D",size=1)+
  geom_point(aes(x=x,y=hr),color="#F8766D",size=3)+
  geom_segment(data=df_coxme_rec[2,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_rec[2,], aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_rec[2,], aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#00BFC4",size=1)+
  geom_point(data=df_coxme_rec[2,], aes(x=x,y=hr),color="#00BFC4",size=3)+
  theme(legend.position="none")+
  scale_x_continuous(name="",breaks = c(0.1666667),labels=c("Adherence [1]"), limits = c(0.06666667,0.26666667))+
  ylab("")
g3<-ggplot(data=df_cox_rec[3,]) +
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#F8766D",size=1)+
  geom_point(aes(x=x,y=hr),color="#F8766D",size=3)+
  geom_segment(data=df_coxme_rec[3,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_rec[3,], aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_rec[3,], aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#00BFC4",size=1)+
  geom_point(data=df_coxme_rec[3,], aes(x=x,y=hr),color="#00BFC4",size=3)+
  theme(legend.position="none")+
  scale_x_continuous(name="",breaks = c(0.3333),labels=c("AgeEvent"), limits = c(0.2333,0.4333))+
  ylab("")
g4<-ggplot(data=df_cox_rec[4,]) +
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#F8766D",size=1)+
  geom_point(aes(x=x,y=hr),color="#F8766D",size=3)+
  geom_segment(data=df_coxme_rec[4,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_point(data=df_coxme_rec[4,], aes(x=x,y=hr),color="#00BFC4",size=3)+
  geom_segment(data=df_coxme_rec[4,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_rec[4,], aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_rec[4,], aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#00BFC4",size=1)+
  theme(legend.position="none")+
  scale_x_continuous(name="",breaks = c(0.5),labels=c("Comorbidity"), limits = c(0.4,0.6))+
  ylab("")

g5<-ggplot(data=df_cox_death[1,]) +
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#F8766D",size=1)+
  geom_point(aes(x=x,y=hr),color="#F8766D",size=3)+
  geom_segment(data=df_coxme_death[1,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[1,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[1,], aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[1,], aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#00BFC4",size=1)+
  geom_point(data=df_coxme_death[1,], aes(x=x,y=hr),color="#00BFC4",size=3)+
  theme(legend.position="none")+
  scale_x_continuous(name="",breaks = c(0.0),labels=c("Sex [M]"), limits = c(-0.1,0.1))+
  ylab("")
g6<-ggplot(data=df_cox_death[2,]) +
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#F8766D",size=1)+
  geom_point(aes(x=x,y=hr),color="#F8766D",size=3)+
  geom_segment(data=df_coxme_death[2,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[2,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[2,], aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[2,], aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#00BFC4",size=1)+
  geom_point(data=df_coxme_death[2,], aes(x=x,y=hr),color="#00BFC4",size=3)+
  theme(legend.position="none")+
  scale_x_continuous(name="",breaks = c(0.1666667),labels=c("Adherence [1]"), limits = c(0.06666667,0.26666667))+
  ylab("")
g7<-ggplot(data=df_cox_death[3,]) +
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#F8766D",size=1)+
  geom_point(aes(x=x,y=hr),color="#F8766D",size=3)+
  geom_segment(data=df_coxme_death[3,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[3,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[3,], aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[3,], aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#00BFC4",size=1)+
  geom_point(data=df_coxme_death[3,], aes(x=x,y=hr),color="#00BFC4",size=3)+
  theme(legend.position="none")+
  scale_x_continuous(name="",breaks = c(0.3333),labels=c("AgeEvent"), limits = c(0.2333,0.4333))+
  ylab("")
g8<-ggplot(data=df_cox_death[4,]) +
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#F8766D",size=1)+
  geom_segment(aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#F8766D",size=1)+
  geom_point(aes(x=x,y=hr),color="#F8766D",size=3)+
  geom_segment(data=df_coxme_death[4,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[4,], aes(x=x,y=l95,xend=x,yend=u95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[4,], aes(x=x-0.01,y=l95,xend=x+0.01,yend=l95),color="#00BFC4",size=1)+
  geom_segment(data=df_coxme_death[4,], aes(x=x-0.01,y=u95,xend=x+0.01,yend=u95),color="#00BFC4",size=1)+
  geom_point(data=df_coxme_death[4,], aes(x=x,y=hr),color="#00BFC4",size=3)+
  theme(legend.position="none")+
  scale_x_continuous(name="",breaks = c(0.5),labels=c("Comorbidity"), limits = c(0.4,0.6))+
  ylab("")

pl1<-ggarrange(g1,g2,g3,g4,nrow=1,labels = "")
annotate_figure(pl1, top = text_grob("Recurrent Events", 
                                        color = "black", face = "bold", size = 10))
pl2<-ggarrange(g5,g6,g7,g8,nrow=1,labels = "")
annotate_figure(pl2, top = text_grob("Terminal Events", 
                                     color = "black", face = "bold", size = 10))

rm(df_cox_rec,df_cox_death,df_coxme_rec,df_coxme_death,g1,g2,g3,g4,g5,g6,g7,g8,pl1,pl2)


