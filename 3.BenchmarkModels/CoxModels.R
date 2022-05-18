###########################################################################################
################################# Cox PH models  ##########################################
###########################################################################################
## [1] Data preparation
## [2] Recurrent Events model
## [3] Terminal Events model
## [4] Post processing

## clean workspace
rm(list = ls())

## load packages
library(data.table)
library(survival)
library(survminer)
library(ggplot2)

## load data
load("dataRec.RData")
load("dataDeath.RData")

## Arrange variables
data$COD_REG= factor(data$COD_REG)
data$SESSO=factor(data$SESSO)
data$ADERENTE=factor(data$ADERENTE)
data$etaEvent=as.double(data$etaEvent)

## Recurrent Events
Cox.Rec <- coxph(Surv(GapEvent,event)~ SESSO + ADERENTE + etaEvent + comorbidity, data=data)
summary(Cox.Rec)
print(Cox.Rec)

## Terminal Events
Cox.Death <- coxph(Surv(GapEvent,cens)~ SESSO + ADERENTE + etaEvent + comorbidity, data=dataDeath)
summary(Cox.Death)
print(Cox.Death)

## Hazard Ratios
x11()
ggforest(Cox.Rec,data=data)
x11()
ggforest(Cox.Death,data=dataDeath)

## Martingale Diagnostics
x11()
ggcoxdiagnostics(Cox.Rec,type='martingale')
ggcoxdiagnostics(Cox.Death, type='martingale')
