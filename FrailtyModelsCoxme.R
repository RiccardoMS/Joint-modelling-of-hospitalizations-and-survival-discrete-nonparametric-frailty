############################################################################################
################################ Cox Mixed Effect Models ###################################
############################################################################################
## [1] Data preparation
## [2] Plot 1: recurrent events encoding example
## [3] Plot 2: multiple events & composite endpoint
## [4] Plot 3: Passing to Gap Times may cause loss of information
## [5] Plot 4: Gap times vs Number event
## [6] Plot 5: Gap Times Histogram
## [7] Cox models
## [8] Coxme models
## [9] Plot 6: HR comparison

## clean workspace
rm(list = ls())

## load packages
library(data.table)
library(survival)
library(survminer)
library(coxme)

## load data
load("10_days_of_Hell/dataRecHell.RData")
load("10_days_of_Hell/dataDeathHell.RData")

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

## Plot 1: recurrent events encoding example
require(gridExtra)
df_plot1<-data.frame(
  COD_REG=c("subject 1","subject 1","subject 1","subject 2","subject 3",
            "subject 3","subject 3","subject 3","subject 4"),
  event =c("event","event","censoring","censoring","event","event",
           "event","censoring","censoring"),
  time_y=c(0.2,0.3,0.7,2,0.35,0.7,1.3,2,1),
  left=rep(0.0,9),
  right=c(0.7,0.7,0.7,2.0,2.0,2.0,2.0,2.0,1.0)
)

df_plot2<-data.frame(
  COD_REG=c("subject 1","subject 2","subject 3","subject 4"),
  event =c("event","censoring","event","censoring"),
  time_y=c(0.2,2,0.35,1),
  left=c(0,0,0,0)
)

plot1 <-ggplot(data=df_plot1,aes(x=COD_REG,y=as.double(time_y),group=COD_REG)) +
  geom_segment(aes(x=COD_REG,y=left,xend=COD_REG,yend=right))+
  geom_point(aes(color=event),size=3)+
  scale_color_manual(values=c("white", "black"))+
  xlab("")+
  ylab("years")+
  coord_flip()+
  geom_hline(yintercept = 2,linetype= "dotted")+
  theme(legend.title = element_blank())

plot2 <- 
  ggplot(data=df_plot2,aes(x=COD_REG,y=as.double(time_y),group=COD_REG)) +
  geom_segment(aes(x=COD_REG,y=left,xend=COD_REG,yend=time_y))+
  geom_point(aes(color=event),size=3)+
  scale_color_manual(values=c("white", "black"))+
  xlab("")+
  ylab("years")+
  coord_flip()+
  geom_hline(yintercept = 2,linetype= "dotted")+
  theme(legend.title = element_blank())

grid.arrange(plot1, plot2, ncol=1)

rm(df_plot1,df_plot2,plot1,plot2)


## Plot 2: multiple events & composite endpoint
data_head=data[1:18,]
data_head$time_y=data_head$timeEvent/365
data_head$event=factor(data_head$event, labels=c('terminal','recurrent'))
data_head$left=rep(0.0,18)
ggplot(data=data_head,aes(x=COD_REG,y=as.double(time_y),group=COD_REG)) +
  geom_segment(aes(x=COD_REG,y=left,xend=COD_REG,yend=time_y))+
  geom_point(aes(color=event,shape=status),size=3)+
  xlab("ID")+
  ylab("Time [years]")+
  coord_flip()
rm(data_head)

## Plot 3: Passing to Gap Times may cause loss of information
data_plot3=data[COD_REG %in% c(10007000,10000717),]
data_plot3$time_y=data_plot3$GapEvent
data_plot3$event=factor(data_plot3$event, labels=c('terminal','recurrent'))
data_plot3$left=rep(0.0,11)
data_plot3$index=c(1.5,2.0,2.5,3.0,3.5,4.0,4.5,6.5,7.0,7.5,8.0)
ggplot(data=data_plot3,aes(x=index,y=as.double(time_y))) +
  geom_segment(aes(x=index,y=left,xend=index,yend=time_y))+
  geom_point(aes(color=event,shape=status),size=3)+
  xlab("patient")+
  ylab("Time [days]")+
  scale_x_discrete(breaks=NULL)+
  annotate("text",x=3.0,y=-20,label="10000717", size=3)+
  annotate("text",x=7.25,y=-20,label="10007000",size=3)
rm(data_plot3)

## Plot 4: Gap times vs Number event
library(plyr)
print(table(data$hosp<=9))
hist(data[hosp>9]$hosp)
indiv.stats <- ddply(data[hosp<=9,],c("COD_REG","hosp"),summarize,time2event = GapEvent)
indiv.stats  <- indiv.stats[order(indiv.stats$hosp, -indiv.stats$time2event),]
indiv.stats$order <- c(1:length(indiv.stats$COD_REG))
indiv.stats <- na.omit(indiv.stats)
medianepi <- ddply(indiv.stats,"hosp",summarize, medep = median(order))
medianepi <- na.omit(medianepi)
ggplot(data=indiv.stats, aes(x=time2event, y=order)) +
  geom_rect(mapping = aes(xmin=0,xmax=indiv.stats$time2event,ymin=indiv.stats$order,
                          ymax=indiv.stats$order+1,fill=factor(hosp))) +
  scale_y_continuous(breaks=NULL) +
  xlab("Time [days]") + ylab("") + labs(fill="Event\nNumber") + theme_gray()
rm(indiv.stats,medianepi)

## Plot 5: Gap Times Histogram
p2<-ggplot(dataDeath, aes(x=GapEvent)) +
  geom_histogram(position="identity", binwidth=30)+
  xlab("Terminal Gap time [days]")+
  ylab("count")

p1<-ggplot(data, aes(x=GapEvent)) +
  geom_histogram(position="identity", binwidth=30)+
  xlab("Gap time [days]")+
  ylab("count")
  
ggarrange(p1,p2,labels = c("",""),nrow=1) 


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

## Plot 5: HR CI comparison
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


