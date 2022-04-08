################################################################################################
##################################### HR COMPARISON ############################################
################################################################################################
# Clean workspace
rm(list=ls())

# Load Libraries
library(survival)
library(coxme)
library(frailtypack)
library(ggplot2)
library(dplyr)

# Load models
load("Saved_Extension/NDGaussian.RData")
load("FFU_ACE_runs/Tenth_run_FFU_ACE/FFU_try.RData")
load("Saved_Coxme/FinalCoxmeRec.RData")
load("Saved_Coxme/FinalCoxmeDeath.RData")
load("Saved_Rondeau/RondeauSummary.RData")

# Build DataFrame for visualization
Visual <- data.frame() # names: x,l95,hr,up95, modelName,coeffName, procName
temp <- data.frame(x=rep(3,8),l95=FFU_try$summary$L95[1:8],hr=FFU_try$summary$HR[1:8],
                u95=FFU_try$summary$U95[1:8],modelName=rep("Ng",8),
                coeffName=rep(c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity"),2),
                procName=rep(c("Hosp","Death"),each=4))
Visual<-rbind(Visual,temp)

final<-Saved[[length(Saved)]]
temp1 <- data.frame(x=rep(4,4),l95=summary(final$modelR)$conf.int[,3],
                   hr=summary(final$modelR)$conf.int[,1],
                   u95=summary(final$modelR)$conf.int[,4],modelName=rep("Discrete \nGaussian",4),
                   coeffName=rep(c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity")),
                   procName=rep("Hosp",4))
temp2 <- data.frame(x=rep(4,4),l95=summary(final$modelT)$conf.int[,3],
                    hr=summary(final$modelT)$conf.int[,1],
                    u95=summary(final$modelT)$conf.int[,4],modelName=rep("Discrete \nGaussian",4),
                    coeffName=rep(c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity")),
                    procName=rep("Death",4))
temp<-rbind(temp1,temp2)
Visual<-rbind(Visual,temp)

load("Saved_Extension/NDUniform.RData")
final<-Saved[[length(Saved)]]
temp1 <- data.frame(x=rep(5,4),l95=summary(final$modelR)$conf.int[,3],
                    hr=summary(final$modelR)$conf.int[,1],
                    u95=summary(final$modelR)$conf.int[,4],modelName=rep("Discrete \nUniform",4),
                    coeffName=rep(c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity")),
                    procName=rep("Hosp",4))
temp2 <- data.frame(x=rep(5,4),l95=summary(final$modelT)$conf.int[,3],
                    hr=summary(final$modelT)$conf.int[,1],
                    u95=summary(final$modelT)$conf.int[,4],modelName=rep("Discrete \nUniform",4),
                    coeffName=rep(c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity")),
                    procName=rep("Death",4))
temp<-rbind(temp1,temp2)
Visual<-rbind(Visual,temp)

temp1 <- data.frame(x=rep(1,4),l95=exp(coxme.Rec$coefficients-1.96*sqrt(diag(vcov(coxme.Rec)))),
                    hr=exp(coxme.Rec$coefficients),
                    u95=exp(coxme.Rec$coefficients+1.96*sqrt(diag(vcov(coxme.Rec)))),
                    modelName=rep("Disjoint",4),
                    coeffName=rep(c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity")),
                    procName=rep("Hosp",4))
temp2 <- data.frame(x=rep(1,4),l95=exp(coxme.Death$coefficients-1.96*sqrt(diag(vcov(coxme.Death)))),
                    hr=exp(coxme.Death$coefficients),
                    u95=exp(coxme.Death$coefficients+1.96*sqrt(diag(vcov(coxme.Death)))),
                    modelName=rep("Disjoint",4),
                    coeffName=rep(c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity")),
                    procName=rep("Death",4))
temp<-rbind(temp1,temp2)
Visual<-rbind(Visual,temp)

temp <- data.frame(x=rep(2,8),l95=RondeauSummary$L95,hr=RondeauSummary$HR,
                   u95=RondeauSummary$U95,modelName=rep("Rondeau",8),
                   coeffName=rep(c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity"),2),
                   procName=rep(c("Hosp","Death"),each=4))
Visual<-rbind(Visual,temp)

# Plot 
Visual %>%
 mutate(modelName=factor(modelName, levels = c("Disjoint","Rondeau","Ng","Discrete \nGaussian", "Discrete \nUniform"))) %>%
 mutate(coeffName=factor(coeffName, levels = c("Sex[M]","Adherence[1]","AgeEvent","Comorbidity"))) %>%
 mutate(procName=factor(procName, levels = c("Hosp","Death"))) %>%
  ggplot() +
  theme(axis.ticks.x=element_blank())+
  facet_grid(procName ~ coeffName, scales = "free_y")+
  geom_segment(aes(x=x,y=l95,xend=x,yend=u95, color=modelName),size=1)+
  geom_point(aes(x=x,y=hr,color=modelName),size=1.5)+
  geom_segment(aes(x=x-0.05,y=l95,xend=x+0.05,yend=l95,color=modelName),size=1)+
  geom_segment(aes(x=x-0.05,y=u95,xend=x+0.05,yend=u95,color=modelName),size=1)+
  ylab("")+
  xlab("")




