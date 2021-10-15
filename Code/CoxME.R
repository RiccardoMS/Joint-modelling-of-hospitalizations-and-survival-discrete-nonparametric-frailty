################################ COX MIXED EFFECT MODEL ###########################################
## Load packages
rm(list = ls())
library(data.table)
library(survival)
library(survminer)
library(coxme)

## Load Dataset
load("With_Adherence_Dataset/ACE_Inhibitors.RData")

## Arrange Dataset
# histogram
per_hist=new[is.na(hosp)==F,max(hosp),by=COD_REG]
hist(per_hist$V1)

# time at hospitalization events
new[!is.na(hosp), timeEvent:= data_prest - data_rif_ev]

# eta at hospitalization events
new[!is.na(hosp), etaEvent:= eta_Min]

# numbers of comorbidity at hospitalization events discharge
new[!is.na(hosp), comorbidity:=rowSums(.SD), .SDcols = 36:55]

# flag for event
new[!is.na(hosp), event:= 1]

# flag for type of censoring
new[!is.na(hosp), cens:=NA]

## build dataset
# key: COD_REG
# flag: event
# time: timeEvent
# patient level: SESSO, ADERENTE
# event level: eta_event, comorbidity
data <- subset(new,hosp>=1,select = c(COD_REG,event,timeEvent,SESSO,ADERENTE,etaEvent,comorbidity,cens))
names(data)
# add censoring event per patient
codici<- unique(data$COD_REG)
for(i in 1:length(codici)){
  paz_corrente <- codici[i]
  temp <- data.frame(paz_corrente,0,unique(new[COD_REG==paz_corrente]$timeOUT),unique(new[COD_REG==paz_corrente]$SESSO), 
                     unique(new[COD_REG==paz_corrente]$ADERENTE),
                     min(new[COD_REG==paz_corrente]$eta_Min) + 
                     as.integer(format(unique(new[COD_REG==paz_corrente]$data_studio_out), format="%Y")) - 
                     as.integer(format(unique(new[COD_REG==paz_corrente]$data_rif_ev), format="%Y")),
                     tail(data[COD_REG==paz_corrente]$comorbidity,n=1),unique(new[COD_REG==paz_corrente]$death))
  names(temp) <- names(data)
  attributes(temp$timeEvent)<-attributes(data$timeEvent)
  data <- rbind(data,temp)
}
# sort data
data <- data[order(COD_REG),]

## Arrange variables
data$COD_REG= factor(data$COD_REG)
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
data$ADERENTE=factor(data$ADERENTE)


## Plot: multiple events & death as non-informative censoring
data_head=data[1:15,]
data_head$time_y=data_head$timeEvent/365
data_head$event=factor(data_head$event, labels=c('terminal','recurrent'))
ggplot(data=data_head,aes(x=COD_REG,y=as.double(time_y),group=COD_REG)) +
  geom_line()+
  geom_point(aes(color=event,shape=status),size=3)+
  xlab("ID")+
  ylab("years")+
  coord_flip()

## Pass to gap times between events
data[,GapEvent:=as.integer(timeEvent)-as.integer(shift(timeEvent)),by=COD_REG]
data<-data[!is.na(GapEvent)]

## COXME model
coxme.mod <- coxme(Surv(GapEvent,event) ~ SESSO + ADERENTE + etaEvent + comorbidity + (1|COD_REG), data = data)
coxme.mod
