######################## JOINT FRAILTY MODEL by Rondeau,2007 ###################################
## Load packages
rm(list = ls())
library(data.table)
library(frailtypack)

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

#data[,check:=as.integer(timeEvent)-as.integer(shift(timeEvent,n=-1)),by=COD_REG]
#data<-data[!(event==1 & check==0)]
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


## Cox Model via frailtypack (penalized MLE) -- Hospitalization
mod.cox.gap <- frailtyPenal(Surv(GapEvent,event)~ SESSO + ADERENTE + 
                            scale(etaEvent) + scale(comorbidity),n.knots=12,kappa=1,data=data,
                            cross.validation = TRUE)
print(mod.cox.gap)
summary(mod.cox.gap)

## Cox Model with random effect -- Hospitalization
mod.cox.gap.shared <- frailtyPenal(Surv(GapEvent,event)~ cluster(COD_REG) + SESSO + ADERENTE + 
                            scale(etaEvent) + scale(comorbidity),n.knots=12,kappa=1,data=data,
                            cross.validation = TRUE)
print(mod.cox.gap.shared)
summary(mod.cox.gap.shared)

## Cox Model with random effect -- Death
mod.cox.death.shared <- frailtyPenal(Surv(GapEvent,cens)~ cluster(COD_REG) + SESSO + ADERENTE + 
                                     scale(etaEvent) + scale(comorbidity),n.knots=12,kappa=1,data=dataDeath,
                                   cross.validation = TRUE)
print(mod.cox.death.shared)
summary(mod.cox.death.shared)

## smoothing parameters
kappa1 <- mod.cox.gap.shared$kappa
kappa2 <- mod.cox.death.shared$kappa

## Joint frailty -- Gamma 
# Not functioning including EtaEvent??

# no etaevent; nknots=8;kappa1*1e15; 
#modJoint.gap <- frailtyPenal(Surv(GapEvent,event)~cluster(COD_REG)+SESSO+ADERENTE+scale(etaEvent)+scale(comorbidity)+
#                             terminal(cens),formula.terminalEvent=~SESSO+ADERENTE++scale(etaEvent)+scale(comorbidity),
#                             data=data,n.knots=8,kappa=c(kappa1,kappa2))
# refine
#modJoint.gap1 <- frailtyPenal(Surv(GapEvent,event)~cluster(COD_REG)+SESSO+ADERENTE+comorbidity+
#                               terminal(cens),formula.terminalEvent=~SESSO+ADERENTE+comorbidity,
#                             data=data,n.knots=20,kappa=c(1,1),nb.gh=32,nb.gl=32,
#                             init.B = modJoint.gap$coef)

#print(modJoint.gap1)
#summary(modJoint.gap1)

## Plots
#plot(modJoint.gap1, event = "Recurrent", type.plot = "Survival", conf.bands
#     = TRUE, pos.legend="bottomright", cex.legend = 0.7, color = 2, median=TRUE,
#     Xlab = "Time", Ylab = "Survival Probability")

## Joint Frailty -- LogNormal
modJoint.logNorm <- frailtyPenal(Surv(GapEvent,event)~cluster(COD_REG)+SESSO+ADERENTE+scale(etaEvent)+ scale(comorbidity)+
                               terminal(cens),formula.terminalEvent=~SESSO+ADERENTE+scale(etaEvent)+scale(comorbidity),
                             data=data,n.knots=10,kappa=c(kappa1,kappa2),RandDist = "LogN")


print(modJoint.logNorm)
summary(modJoint.logNorm)

## Plots
plot(modJoint.logNorm, event = "Both", type.plot = "Survival", conf.bands
     = TRUE, pos.legend="bottomright", cex.legend = 0.7, color = 2, median=TRUE,
     Xlab = "Time", Ylab = "Survival Probability")

