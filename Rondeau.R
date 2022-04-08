######################## JOINT FRAILTY MODEL by Rondeau,2007 ###################################
## Load packages
rm(list = ls())
library(data.table)
library(frailtypack)

## load data
load("10_days_of_Hell/dataRecHell.RData")
load("10_days_of_Hell/dataDeathHell.RData")

## Arrange variables
data$COD_REG= factor(data$COD_REG)
data$SESSO=factor(data$SESSO)
data$ADERENTE=factor(data$ADERENTE)
data$etaEvent=as.double(data$etaEvent)/100
dataDeath$etaEvent=as.double(dataDeath$etaEvent)/100
dataDeath$timeEvent<-dataDeath$timeEvent - 365


## Cox Model via frailtypack (penalized MLE) -- Hospitalization
mod.cox.gap <- frailtyPenal(Surv(GapEvent,event)~ SESSO + ADERENTE + 
                            scale(etaEvent) + scale(comorbidity),n.knots=12,kappa=1,data=data,
                            cross.validation = TRUE)
print(mod.cox.gap)
summary(mod.cox.gap)

## Cox Model with random effect -- Hospitalization
mod.cox.gap.shared <- frailtyPenal(Surv(GapEvent,event)~ cluster(COD_REG) + SESSO + ADERENTE + 
                            etaEvent + comorbidity,n.knots=12,kappa=1,data=data,
                            cross.validation = TRUE)
print(mod.cox.gap.shared)
summary(mod.cox.gap.shared)

## Cox Model with random effect -- Death
mod.cox.death.shared <- frailtyPenal(Surv(GapEvent,cens)~ cluster(COD_REG) + SESSO + ADERENTE + 
                                     etaEvent + comorbidity,n.knots=12,kappa=1,data=dataDeath,
                                   cross.validation = TRUE)
print(mod.cox.death.shared)
summary(mod.cox.death.shared)

## smoothing parameters
kappa1 <- mod.cox.gap.shared$kappa
kappa2 <- mod.cox.death.shared$kappa

## Joint frailty -- Gamma 
# Not functioning including EtaEvent??

# no etaevent; nknots=8;kappa1*1e15; 
#modJoint.gap <- frailtyPenal(Surv(GapEvent,event)~cluster(COD_REG)+SESSO+ADERENTE+scale(etaEvent)+comorbidity+
#                             terminal(cens),formula.terminalEvent=~SESSO+ADERENTE++scale(etaEvent)+comorbidity,
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
modJoint.logNorm <- frailtyPenal(Surv(GapEvent,event)~cluster(COD_REG)+SESSO+ADERENTE+etaEvent+comorbidity+
                               terminal(cens),formula.terminalEvent=~SESSO+ADERENTE+etaEvent+comorbidity,
                             data=data,n.knots=12,kappa=c(kappa1,kappa2),RandDist = "LogN", )


print(modJoint.logNorm)
summary(modJoint.logNorm)

## Plots
plot(modJoint.logNorm, event = "Both", type.plot = "Survival", conf.bands
     = TRUE, pos.legend="bottomright", cex.legend = 0.7, color = 2, median=TRUE,
     Xlab = "Time", Ylab = "Survival Probability")

save(modJoint.logNorm, file="10_days_of_Hell/Rondeau.RData")

## Summary
RondeauSummary<-data.frame(
          Estimate=c(0.044,-0.263,-0.014,0.118,0.055,-0.211,0.041,0.137),
          StdDev=c(0.023,0.022,0.001,0.007,0.081,0.087,0.003,0.023),
          L95=c(0.995,0.736,0.984,1.109,0.901,0.683,1.036,1.096),
          HR=c(1.044,0.769,0.986,1.125,1.056,0.810,1.042,1.147),
          U95=c(1.093,0.802,0.988,1.141,1.238,0.960,1.048,1.199))
save(RondeauSummary, file="10_days_of_Hell/RondeauSummary.RData")
  
  
  
  
  
  
  
  
  
  

