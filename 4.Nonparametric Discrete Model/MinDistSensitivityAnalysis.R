###############################################################################################
####################### MinDist Sensitivity Analysis ##########################################
###############################################################################################

## Script build several discrete nonparam models with different MinDist values
## to builds a smooth loglikelihood curve 

## clear workspace
rm(list=ls())

## load libraries
library(data.table)
library(survival)
library(gganimate)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(gifski)
library(tidyr)
library(mvtnorm)
library(survminer)

## load data
load("dataRec.RData")
load("dataDeath.RData")

## data preprocessing
# ordering (not necessary?)
data <- data[order(data$GapEvent),]
dataDeath <- dataDeath[order(dataDeath$GapEvent),]

# arrange data
ID_rec          <- factor(data$COD_REG)
ADERENTE_rec    <- factor(data$ADERENTE)
SESSO_rec       <- factor(data$SESSO)
etaEvent_rec    <- as.double(data$etaEvent)
comorbidity_rec <- as.double(data$comorbidity)

ID_term          <- factor(dataDeath$COD_REG)
ADERENTE_term    <- factor(dataDeath$ADERENTE)
SESSO_term       <- factor(dataDeath$SESSO)
etaEvent_term    <- as.double(dataDeath$etaEvent)
comorbidity_term <- as.double(dataDeath$comorbidity)

# Model Matrices
X1 <- model.matrix(~ SESSO_rec + ADERENTE_rec + etaEvent_rec + comorbidity_rec )[,-1]
X2 <- model.matrix(~ SESSO_term + ADERENTE_term + etaEvent_term + comorbidity_term )[,-1]

# Response vectors
time1 <- data$GapEvent
time2 <- dataDeath$GapEvent

# Number of patients
N = length(unique(ID_rec))

# Number of unique times
nt1 = length(unique(data$GapEvent))
nt2 = length(unique(dataDeath$GapEvent))

# Number of observations
nobs1 = nrow(data)
nobs2 = nrow(dataDeath)

# Number of covariates
ncov1 = ncol(X1)
ncov2 = ncol(X2)

# Cumulative Hazards: data frame encoding
cumhaz1 = as.data.frame( cbind( hazard=rep( 0, nt1 ), time = sort( unique( time1 ))))
cumhaz2 = as.data.frame( cbind( hazard=rep( 0, nt2 ), time = sort( unique( time2 ))))

# Instantaneuos hazards
haz1 = as.data.frame( cbind( hazard=rep( 0, nt1 ), time = sort( unique( time1 ))))
haz2 = as.data.frame( cbind( hazard=rep( 0, nt2 ), time = sort( unique( time2 ))))

# Number of events per patient   
D1 <- table( ID_rec[ data$event == 1 ] )
orderD1<-match(unique(ID_rec),unique(data$COD_REG)[order(unique(data$COD_REG))])
D1<- D1[orderD1]
D2 <- table( ID_term[dataDeath$cens == 1])
orderD2<-match(unique(ID_rec),unique(dataDeath$COD_REG)[order(unique(dataDeath$COD_REG))])
D2<-D2[orderD2]

# Risk sets
risk_index1 <- matrix( 0, nrow = nobs1, ncol = nt1 )
risk_index2 <- matrix(0, nrow = nobs2, ncol = nt2)

# Time lists
time_list1 <- sapply( 1:nt1, function(x) !is.na(match(time1, cumhaz1$time[x])))
time_list2 <- sapply( 1:nt2, function(x) !is.na(match(time2, cumhaz2$time[x])))

# Number of ties
m1 <- sapply( 1:dim(time_list1)[2], function(x) sum(data$event[time_list1[,x]]))
m2 <- sapply( 1:dim(time_list2)[2], function(x) sum(dataDeath$cens[time_list2[,x]]))

# fill risk index 
for( l in 1:nt1 )
{
  risk_index1[ which( time1 >= cumhaz1$time[ l ]), l ] <-  1 
}
for( l in 1:nt2 )
{
  risk_index2[ which( time2 >= cumhaz2$time[ l ]), l ] <-  1 
}

# Encode ID as numeric 
groups1 <- match(ID_rec, unique(ID_rec))
groups2 <- match(ID_term, unique(ID_rec))

## Initialization parameters
library(MASS)
Sigma <- matrix(c(0.12,0.0,0.0,1.4),nrow = 2, ncol = 2)
mu=c(0,0)

## Support Functions
# LogLikelihood
LOGL<-function(Z,w,haz1,cumhaz1,beta,haz2,cumhaz2, gamma,P,E_formula1,E_formula2,E_haz1,E_haz2){
  lw<-sum(Z%*%matrix(log(w)))
  lr<- sum(Z*(E_haz1-E_formula1))
  ld<- sum(Z*(E_haz2-E_formula2))
  return (lw+lr+ld)
}

# Initialize saving structure
Saved <- list()

## MinDist Candidates
candidates <- seq(0.1,1.0,by=0.025)
candNumber <- 1

## Initialization setting
Uniform=FALSE

## MinDist wrapper
for(MinDist in candidates){
  
  ## Internal Saving structure
  SavedInt <- list()
  
  ## Set Convergence Indices
  # Count
  count <- 0
  
  # epsilon
  eps_conv = 1e-3
  eps      = 1e5
  
  # Grid Initialization
  set.seed(210197)
  if (!Uniform){
    K<-1000
    P <- mvrnorm(K,mu,Sigma)
    w<-dmvnorm(P,mu,Sigma)
    w<-w/sum(w)
    P_show<-P
    P_show[,1]<-P[,1]-rep(w%*%P[,1],length(P[,1]))
    P_show[,2]<-P[,2]-rep(w%*%P[,2],length(P[,2]))
  } else {
    temp_u<-seq(-3*sqrt(Sigma[1,1]),3*sqrt(Sigma[1,1]), by=MinDist)
    temp_v<-seq(-3*sqrt(Sigma[2,2]),3*sqrt(Sigma[2,2]), by=MinDist)
    P     <- unname(data.matrix(expand.grid(temp_u,temp_v)))
    w     <- rep(1/dim(P)[1],dim(P)[1])
    P_show<-P
    P_show[,1]<-P[,1]-rep(w%*%P[,1],length(P[,1]))
    P_show[,2]<-P[,2]-rep(w%*%P[,2],length(P[,2]))
    K<-dim(P)[1]
  }

  # Grid shrinking
  is_near<-TRUE
  while(is_near){
  D<-dist(P)
  D<-as.matrix(D)
  D[upper.tri(D)]<-10
  diag(D)<-10
  out<-which(D == min(D), arr.ind = TRUE)
  if(D[out][1]<MinDist){
    #merge
    P[out[1,2],]=(P[out[1,2],]+P[out[1,1],])/2
    P<-P[-out[1,1],]
    #update weights
    w[out[1,2]]<-w[out[1,2]]+w[out[1,1]]
    w<-w[-out[1,1]]
    w<-w/sum(w)
    K<-K-1
  }
  else {
    is_near=FALSE
  }
}

  # Assign patient to random frailty, built frailties vectors
  P_index <- sample(1:K,size=N,replace = T, prob = w)
  P_off1  <- P[,1][P_index[groups1]]
  P_off2  <- P[,2][P_index[groups2]]

  # Estimate Initial Recurrent model, Cumulative Hazard and Hazard
  cox1 <- coxph(Surv(time1,event)~SESSO_rec + ADERENTE_rec + etaEvent_rec +
                comorbidity_rec + offset(P_off1),data=data)
  beta<-cox1$coefficients
  s1   <- survfit(cox1,data=data)
  cumhaz1$hazard = s1$cumhaz
  haz1$hazard = diff(c(0,cumhaz1$hazard)) 
  for(j in 1:length(haz1$hazard)){
  if(haz1$hazard[j]==0)
    haz1$hazard[j]<-haz1$hazard[j-1]
}

  # Estimate Initial Recurrent model, Cumulative Hazard and Hazard
  cox2 <- coxph(Surv(time2,cens)~SESSO_term + ADERENTE_term + etaEvent_term +
                comorbidity_term + offset(P_off2),data=dataDeath)
  gamma <-cox2$coefficients
  s2   <- survfit(cox2,data=dataDeath)
  cumhaz2$hazard = s2$cumhaz
  haz2$hazard = diff(c(0,cumhaz2$hazard)) 
  for(j in 1:length(haz2$hazard)){
  if(haz2$hazard[j]==0)
    haz2$hazard[j]<-haz2$hazard[j-1]
}

  # Initialize structures for computations
  numerator <- rep( 0, K )
  Z <- E_formula1 <- E_formula2 <- E_haz1<-E_haz2<- matrix( 0, nrow = N, ncol = K)
  E_part1 <- E_part2<-rep( 0, N)

  # Start loop
  while (eps_conv < eps & count < 100 ){
  
  # Save current w estimates
  w_old <- w
  
  # Save current Z estimates
  Z_old <- Z
  
  # Save current P estimates
  P_old <-P
  
  # Grid Shrinking
  is_near=TRUE
  while(is_near){
    D<-dist(P)
    D<-as.matrix(D)
    D[upper.tri(D)]<-10
    diag(D)<-10
    out<-which(D == min(D), arr.ind = TRUE)
    if(D[out][1]<MinDist){
      #merge
      P[out[1,2],]=(P[out[1,2],]+P[out[1,1],])/2
      P<-P[-out[1,1],]
      #update weights
      w[out[1,2]]<-w[out[1,2]]+w[out[1,1]]
      w<-w[-out[1,1]]
      w<-w/sum(w)
      K<-K-1
    }
    else{
      is_near<-FALSE
    }
  }
  
  # Clean Structures
  Z <- E_formula1 <- E_formula2 <- E_haz1 <- E_haz2<- matrix( 0, nrow = N, ncol = K)
  numerator <- rep(0,K)
  
  # Expectation Step
  for(i in 1:N){
    
    current_patient1 <- groups1==i
    current_patient2 <- groups2==i
    
    ebz1 <- exp( X1[current_patient1,] %*% beta )
    ebz2 <- exp( X2[current_patient2,] %*% gamma)
    
    tRij <- match(time1[current_patient1], cumhaz1$time)
    H01t <- cumhaz1$hazard[tRij]
    lh01t<-  log(haz1$hazard[tRij])
    
    tDi <- match(time2[current_patient2], cumhaz2$time)
    H02t <- cumhaz2$hazard[tDi]
    lh02t<-  log(haz2$hazard[tDi])
    
    E_part1[i] <- ifelse( ncov1 > 0,
                          sum( H01t*ebz1 ),
                          sum( H01t ) )
    E_part2[i] <- ifelse( ncov2 > 0,
                          H02t*ebz2 ,
                          H02t)
    
    for(l in 1:K){
      E_formula1[i,l] <- ifelse( ncov1 > 0,
                                 sum( H01t*ebz1*exp(P[l,1])),
                                 sum( H01t*exp(P[l,1]) ) )
      E_formula2[i,l] <- ifelse( ncov2 > 0,
                                 H02t*ebz2*exp(P[l,2]),
                                 H02t*exp(P[l,2]))
      E_haz1[i,l]     <- sum(data$event[current_patient1]*(lh01t+log(ebz1)+P[l,1]))
      
      E_haz2[i,l]     <- dataDeath$cens[current_patient2]*(lh02t+log(ebz2)+P[l,2])
      
      pivot <- min(as.numeric(D1)[i]*(P[l,1]) - E_formula1[i,l] +
                     + as.numeric(D2)[i]*(P[l,2]) - E_formula2[i,l])
      numerator[l] <- w[l]*exp(as.numeric(D1)[i]*(P[l,1]) - E_formula1[i,l] +
                                 + as.numeric(D2)[i]*(P[l,2]) - E_formula2[i,l])
      
      if(max(numerator)==0)
        numerator <- 1e-16/((as.numeric(D1)[i]*P[,1]-E_formula1[i,]+as.numeric(D2)[i]*P[,2]-E_formula2[i,])/pivot)
    }
    
    Z[i,] <- numerator/sum(numerator)
    
  }
  
  # Maximization Step
  # Latent partition
  belonging <- as.numeric( apply(Z, 1, which.max) )
  
  # Grid Shrinking - Unassigned points
  t<-table(factor(belonging, levels = 1:K))
  to_elim<-which(as.numeric(t)==0)
  if(length(to_elim)>0){
    Z<-Z[,-to_elim]
    E_formula1<-E_formula1[,-to_elim]
    E_formula2<-E_formula2[,-to_elim]
    E_haz1<-E_haz1[,-to_elim]
    E_haz2<-E_haz2[,-to_elim]
    numerator <-numerator[-to_elim]
    P<-P[-to_elim,]
    K<-K-length(to_elim)
  }
  
  # Vector of proportions
  if(K>1)
    w <- (colSums(Z))/ N
  else{
    w <- 1 
    P <- matrix(c(0,0),nrow = 1,ncol = 2)
  }
  
  ## Unconstrained Optimization
  P[,1]   <- log(( as.numeric(D1) %*% Z )/( E_part1 %*% Z))
  P[,2]   <- log(( as.numeric(D2) %*% Z )/( E_part2 %*% Z))
  
  P_show    <-P
  P_show[,1]<-P[,1]-w%*%P[,1]
  P_show[,2]<-P[,2]-w%*%P[,2]
  
  # Frailties update
  P_off1 <-log((Z%*%exp(matrix(P[,1])))[groups1])
  P_off2 <-log((Z%*%exp(matrix(P[,2])))[groups2])
  
  # Estimate Betas
  temp_model1 <- coxph(Surv(time1,event)~ SESSO_rec + ADERENTE_rec + etaEvent_rec +
                         comorbidity_rec + offset(P_off1),data=data,method = "breslow")
  beta       <- temp_model1$coef
  
  # Estimate Gammas
  temp_model2 <- coxph(Surv(time2,cens)~ SESSO_term + ADERENTE_term + etaEvent_term +
                         comorbidity_term + offset(P_off2),data=dataDeath, method = "breslow")
  gamma       <- temp_model2$coef 
  
  # Estimate Cumulative Hazard and Hazard - Recurrent
  s1   <- survfit(temp_model1,data=data)
  cumhaz1$hazard = s1$cumhaz
  haz1$hazard = diff(c(0,cumhaz1$hazard)) 
  for(j in 1:length(haz1$hazard)){
    if(haz1$hazard[j]==0)
      haz1$hazard[j]<-haz1$hazard[j-1]
  }
  
  # Estimate Cumulative Hazard and Hazard - Terminal
  s2   <- survfit(temp_model2,data=dataDeath)
  cumhaz2$hazard = s2$cumhaz
  haz2$hazard = diff(c(0,cumhaz2$hazard)) 
  for(j in 1:length(haz2$hazard)){
    if(haz2$hazard[j]==0)
      haz2$hazard[j]<-haz2$hazard[j-1]
  }
  
  # Convergence
  if(length(w)==length(w_old))
    eps <- max(abs(w - w_old),na.rm=T)
  else if (length(w)==1)
    eps <- 0
  else
    eps<-1
  
  # Count
  count <- count + 1
  
  # Print
  #print(P)
  #print(w)
  #print(table(belonging))
  
  # AIC computation
  LogL= -LOGL(Z=Z,w=w, haz1 = haz1,beta = beta,cumhaz1 = cumhaz1,
                         haz2 = haz2, gamma = gamma, cumhaz2 = cumhaz2,
                         P=P,E_formula1 = E_formula1,E_formula2 = E_formula2,
                         E_haz1 = E_haz1, E_haz2=E_haz2)
  # Save 
  temp_list<-list("modelR"=temp_model1,"modelT"=temp_model2,"w"=w,"P"=P,"P_show"=P_show,"cumhaz1"=cumhaz1,"cumhaz2"=cumhaz2,"LogL"=LogL,"Table"=table(belonging))
  SavedInt[[count]]<-temp_list
}
  
  # Save
  Saved[[candNumber]]<-SavedInt
  
  # update candNumber
  candNumber <- candNumber + 1
}  

############################ PLOTS: Gaussian Initialization #################################
load("C:/Users/aughi/Desktop/TESI/Code/MinDistSensitivityGaussian.RData")
AICbyMinDist<-data.frame()
for( i in 1:(length(Saved)-1)){
  temp <- data.frame(MinDist=candidates[i],
                     aic=(Saved[[i]])[[length(Saved[[i]])]]$LogL)
  AICbyMinDist<-rbind(AICbyMinDist, temp)
}

AICbyMinDist %>%
  ggplot(aes(x=MinDist, y=aic))+
  geom_line(size=1.0)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text())+
  labs(title = "AIC vs MinDist",
       x = "MinDist",
       y = "AIC",
       caption = "Drug: ACE Inhibitors \nInitialization: Gaussian")+
  annotate("rect", xmin = 0.475, xmax = 0.725, ymin = 91000, ymax = 91700,
           alpha = .3)+
  annotate("rect", xmin = 0.750, xmax = 0.825, ymin = 90000, ymax = 90400,
         alpha = .3)+
  annotate("rect", xmin = 0.85, xmax = 0.975, ymin = 90400, ymax = 90800,
           alpha = .3)+
  annotate("rect", xmin = 0.1, xmax = 0.45, ymin = 91000, ymax = 95000,
           alpha = .0, color="red", fill="red")+
  annotate("text",x=0.275,y=90900, label="transition",color="red")+
  annotate("text",x=0.600,y=90900, label="stable [1]",color="blue")+
  annotate("text",x=0.79,y=90700, label="stable [2]",color="blue")+
  annotate("text",x=0.925,y=91000, label="stable [3]",color="blue")


  

FinalResultsStable<-data.frame()
for(number in 16:26){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                          w=S$w,
                          P1=(S$P_show)[,1],
                          P2=(S$P_show)[,2],
                          LogL=rep(S$LogL,length(S$w)),
                          score=(S$P_show)[,1]+(S$P_show)[,2])
      FinalResultsStable<-rbind(FinalResultsStable,tempDF)
}


FinalResultsStable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist))+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Stable [1]")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none")


FinalResultsUnstable<-data.frame()
for(number in 1:15){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist))+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Transition")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none")

# 3
FinalResultsUnstable<-data.frame()
for(number in 27:30){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist))+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Stable [2]")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none")

#4
FinalResultsUnstable<-data.frame()
for(number in 31:37){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist))+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Stable [3]")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none")

#RepresentativePoints
FinalResultsUnstable<-data.frame()
for(number in c(21,29,34)){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist), ncol = 1)+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Representative values of MinDist in the stability Region")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none")

############################ PLOTS: Uniform Initialization #################################
load("C:/Users/aughi/Desktop/TESI/Code/MindistSensitivityUniform.RData")
candidates <- seq(0.1,1.0,by=0.025)
AICbyMinDist<-data.frame()
for( i in 1:(length(Saved)-1)){
  temp <- data.frame(MinDist=candidates[i],
                     aic=(Saved[[i]])[[length(Saved[[i]])]]$LogL)
  AICbyMinDist<-rbind(AICbyMinDist, temp)
}

AICbyMinDist %>%
  ggplot(aes(x=MinDist, y=aic))+
  geom_line(size=1.0)+
  theme_fivethirtyeight()+
  theme(axis.title = element_text())+
  labs(title = "AIC vs MinDist",
       x = "MinDist",
       y = "AIC",
       caption = "Drug: ACE Inhibitors \nInitialization: Uniform")+
  annotate("rect", xmin = 0.48, xmax = 0.52, ymin = 92400, ymax = 92850,
           alpha = .0, color="purple",fill="purple")+
  annotate("rect", xmin = 0.1, xmax = 0.24, ymin = 92100, ymax = 95600,
           alpha = .0, color="red", fill="red")+
  annotate("rect", xmin = 0.25, xmax = 0.474, ymin = 92100, ymax = 93700,
           alpha = .0, color="red", fill="red",linetype="dashed")+
  annotate("rect", xmin = 0.525, xmax = 0.725, ymin = 91200, ymax = 92500,
           alpha = .0, color="red", fill="red")+
  annotate("rect", xmin = 0.726, xmax = 1.0, ymin = 91000, ymax = 92000,
           alpha = .3)+
  annotate("text",x=0.17,y=91950, label="Descent",color="red")+
  annotate("text",x=0.47,y=91400, label="Descent",color="red")+
  annotate("text",x=0.53,y=93050, label="Elbow",color="purple")+
  annotate("text",x=0.85,y=92150, label="Stable",color="blue")
  

# Elbow
FinalResultsUnstable<-data.frame()
for(number in 17){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist), ncol = 1)+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Elbow - MinDist=0.50")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none")  
# Steep Descent
FinalResultsUnstable<-data.frame()
for(number in 1:6){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist))+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Steep descent")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none") 

# Unstable 
FinalResultsUnstable<-data.frame()
for(number in 7:16){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist))+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "High noise")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none") 
# Second Descent Region
FinalResultsUnstable<-data.frame()
for(number in 18:26){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist))+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Descent")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none") 
# Stable
FinalResultsUnstable<-data.frame()
for(number in 27:37){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist))+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Stable")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none") 

# Representative
# RepresentativePoints
FinalResultsUnstable<-data.frame()
for(number in c(17,32)){
  S<-Saved[[number]][[length(Saved[[number]])]]
  tempDF <-data.frame(MinDist=rep(candidates[number], length(S$w)),
                      w=S$w,
                      P1=(S$P_show)[,1],
                      P2=(S$P_show)[,2],
                      LogL=rep(S$LogL,length(S$w)),
                      score=(S$P_show)[,1]+(S$P_show)[,2])
  FinalResultsUnstable<-rbind(FinalResultsUnstable,tempDF)
}


FinalResultsUnstable %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~factor(MinDist, labels = c("Elbow - 0.50", "Minimum - 0.875")), ncol = 1)+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank(),
        title = element_text())+
  labs(title = "Representative values of MinDist")+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                      guide ="none")

