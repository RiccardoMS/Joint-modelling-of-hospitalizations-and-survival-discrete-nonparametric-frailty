##############################################################################################
############################### NG Postprocessing ############################################
##############################################################################################
## Load packages
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggthemes)
library(latex2exp)

## Results Data frame
loadData=TRUE
if(loadData){
load("Ng_runs/First_run/Omega.RData")
load("Ng_runs/First_run/Phi.RData")
Omega1<-Omega0
Phi1  <-Phi0
load("Ng_runs/Second_run/Omega.RData")
load("Ng_runs/Second_run/Phi.RData")
Omega2<-Omega0
Phi2  <-Phi0
load("Ng_runs/Third_run/Omega.RData")
load("Ng_runs/Third_run/Phi.RData")
Omega3<-Omega0
Phi3  <-Phi0
load("Ng_runs/Fourth_run/Omega.RData")
load("Ng_runs/Fourth_run/Phi.RData")
Omega4<-Omega0
Phi4  <-Phi0
load("Ng_runs/Fifth_run/Omega.RData")
load("Ng_runs/Fifth_run/Phi.RData")
Omega5<-Omega0
Phi5  <-Phi0
load("Ng_runs/Sixth_run/Omega.RData")
load("Ng_runs/Sixth_run/Phi.RData")
Omega6<-Omega0
Phi6  <-Phi0
load("Ng_runs/Seventh_run/Omega.RData")
load("Ng_runs/Seventh_run/Phi.RData")
Omega7<-Omega0
Phi7  <-Phi0
load("Ng_runs/Eigth_run/Omega.RData")
load("Ng_runs/Eigth_run/Phi.RData")
Omega8<-Omega0
Phi8  <-Phi0
load("Ng_runs/Ninth_run/Omega.RData")
load("Ng_runs/Ninth_run/Phi.RData")
Omega9<-Omega0
Phi9  <-Phi0
load("Ng_runs/Tenth_run/Omega.RData")
load("Ng_runs/Tenth_run/Phi.RData")
Omega10<-Omega0
Phi10  <-Phi0
NgResults=data.frame()
for(i in 1:10){
  temp<-data.frame(iteration=i*10,
                   beta1=eval(parse(text = paste("Omega",i,sep="")))[1],
                   beta2=eval(parse(text = paste("Omega",i,sep="")))[2],
                   beta3=eval(parse(text = paste("Omega",i,sep="")))[3],
                   beta4=eval(parse(text = paste("Omega",i,sep="")))[4],
                   gamma1=eval(parse(text = paste("Omega",i,sep="")))[5],
                   gamma2=eval(parse(text = paste("Omega",i,sep="")))[6],
                   gamma3=eval(parse(text = paste("Omega",i,sep="")))[7],
                   gamma4=eval(parse(text = paste("Omega",i,sep="")))[8],
                   thetau=eval(parse(text = paste("Phi",i,sep="")))[1],
                   thetav=eval(parse(text = paste("Phi",i,sep="")))[2],
                   rho=eval(parse(text = paste("Phi",i,sep="")))[3])
  NgResults<-rbind(NgResults,temp)
} 
save(NgResults, file="NgResults.RData")} else {
  load("Saved_Ng/NgResults.RData")
}

## Plot Estimates across iterations
tempB <- NgResults %>%
  select(iteration,beta1,beta2,beta3, beta4) %>%
  gather(key="variable",value="value",-iteration)
g1<-ggplot(data=tempB, aes(x=iteration, y=value))+
  geom_line(aes(color=variable))+
  ylab("")+
  labs(color="Betas")
tempG <- NgResults %>%
  select(iteration,gamma1,gamma2,gamma3, gamma4) %>%
  gather(key="variable",value="value",-iteration)
g2<-ggplot(data=tempG, aes(x=iteration, y=value))+
  geom_line(aes(color=variable))+
  ylab("")+
  labs(color="Gammas")
tempPhi <- NgResults %>%
  select(iteration,thetau,thetav,rho) %>%
  gather(key="variable",value="value",-iteration)

g3<-ggplot(data=tempPhi, aes(x=iteration, y=value))+
  geom_line(aes(color=variable), size=1.5)+
  scale_color_discrete(labels=c(TeX("$\\rho$"),TeX("$\\theta_{u}^{2}$"),TeX("$\\theta_{v}^{2}$")))+
  theme_fivethirtyeight() +
  theme(axis.title = element_text())+
  labs(title = "RE parameters across Iteration",
       x = "Iteration",
       y = "Value",
       caption = "Drug: ACE Inhibitors",
       color="Parameter")
  


## Plot Random Effects estimates
NgRE<-data.frame(u=Omega6[9:3240],v=Omega6[3241:6472],score=Omega6[9:3240]+Omega6[3241:6472])
g4<-ggplot(data=NgRE, aes(x=u,y=v,color=score))+
    geom_point()+
    scale_color_gradient(low = "blue", high = "red" )+
    xlab("u")+
    ylab("v")+
    theme_fivethirtyeight() +
    theme(axis.title = element_text())+
    labs(title = "Random Effects Estimates",
       x = "u",
       y = "v",
       caption = "Drug: ACE Inhibitors",
       color="Frailty\nScore")+
    annotate("rect", xmin = 0.5, xmax = 1.0, ymin = 1.5, ymax = 3,
           alpha = .3)+
    annotate("rect", xmin = -0.25, xmax = 0.25, ymin = -0.5, ymax = 0.5,
           alpha = .2)+
    annotate("rect", xmin = -0.0, xmax = 0.5, ymin = 0.0, ymax = 1.5,
           alpha = .25)+
    annotate("rect", xmin = -0.5, xmax = 0.0, ymin = -1.5, ymax = 0.0,
           alpha = .1)+
    annotate("text",x=-0.25,y=-1.7,label="Protected",color="blue")+
    annotate("text",x=-0.1,y=0.7,label="Neutral", color="purple")+
    annotate("text",x=0.4,y=-0.1,label="At Risk", color="magenta")+
    annotate("text",x=0.75:1.0,y=1.3,label="Highly at Risk",color="red")
  