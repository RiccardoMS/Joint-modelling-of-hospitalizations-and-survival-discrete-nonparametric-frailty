########################## COX MODELS FOR RECURRENT EVENTS AND DEATH #######################
## Load packages
rm(list = ls())
library(data.table)
library(data.table)
library(survival)
library(survminer)
library(ggplot2)

#Load data
load("With_Adherence_Dataset_FFU/ACE_Inhibitors.RData")

## Arrange Dataset
# time at hospitalization events
new[!is.na(hosp), timeEvent:= data_prest - data_rif_ev]

# eta at hospitalization events
new[!is.na(hosp), etaEvent:= eta_Min]

# numbers of comorbidity at hospitalization events discharge
new[!is.na(hosp), comorbidity:=rowSums(.SD), .SDcols = 36:55]

# flag for event
new[!is.na(hosp), event:= 1]

# flag for type of censoring
new[!is.na(hosp), cens:=0]

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
data$SESSO=factor(data$SESSO)
data$ADERENTE=factor(data$ADERENTE)
data$etaEvent=as.double(data$etaEvent)

# Gap times
data[,check:=as.integer(timeEvent)-as.integer(shift(timeEvent,n=-1)),by=COD_REG]
data<-data[!(event==1 & check==0)]
data[,GapEvent:=as.integer(timeEvent)-as.integer(shift(timeEvent)),by=COD_REG]
data<-data[!is.na(GapEvent) & GapEvent!=0]

## First 1500 patients
#codici<-levels(as.factor(data$COD_REG))
#codici<-codici[1:1500]
#data<-data[COD_REG %in% codici]


# data for terminal events
dataDeath <- data[event==0]


## remove unused structures
rm(new)
rm(temp)
gc()

# Define Data DEath
dataDeath <- data[event==0]

############################# RECURRENT EVENTS #############################################
Cox.Rec <- coxph(Surv(GapEvent,event)~ SESSO + ADERENTE + etaEvent + comorbidity, data=data)
summary(Cox.Rec)
print(Cox.Rec)

Cox.Death <- coxph(Surv(GapEvent,cens)~ SESSO + ADERENTE + etaEvent + comorbidity, data=dataDeath)
summary(Cox.Death)
print(Cox.Death)

x11()
ggforest(Cox.Rec,data=data)

x11()
ggforest(Cox.Death,data=dataDeath)

x11()
ggcoxdiagnostics(Cox.Rec,type='martingale')
ggcoxdiagnostics(Cox.Death, type='martingale')
