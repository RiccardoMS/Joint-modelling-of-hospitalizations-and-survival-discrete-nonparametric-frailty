########################## Joint Frailty model via GLMM, by Shu kai ###########################
## Load packages
rm(list = ls())
library(data.table)

#Load data
load("With_Adherence_Dataset/ACE_Inhibitors.RData")

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
data$SESSO=factor(data$SESSO)
data$ADERENTE=factor(data$ADERENTE)
data$etaEvent=as.double(data$etaEvent)

## Pass to gap times between events
data[,GapEvent:=as.integer(timeEvent)-as.integer(shift(timeEvent)),by=COD_REG]
data<-data[!is.na(GapEvent) & GapEvent!=0]

# Note: double recordings of same event cause the FrailtyPenal to crash (gapEvent == 0)

## remove unused structures
rm(new)
rm(temp)
gc()


############################ ALGORITHM ####################################################
dataDeath <- data[event==0]

# Variables
ID_rec          <- factor(data$COD_REG)
ADERENTE_rec    <- factor(data$ADERENTE)
SESSO_rec       <- factor(data$SESSO)
etaEvent_rec    <- as.double(data$etaEvent)
comorbidity_rec <- as.numeric(data$comorbidity)

ID_term          <- factor(dataDeath$COD_REG)
ADERENTE_term    <- factor(dataDeath$ADERENTE)
SESSO_term       <- factor(dataDeath$SESSO)
etaEvent_term    <- as.double(dataDeath$etaEvent)
comorbidity_term <- as.numeric(dataDeath$comorbidity)

# Design parameters
n_beta = 4            # recurrent event covariates
n_gamma = 4           # terminal event covariates
n_subjects = 2916     # number of units subject to random effect 

# Initial values
# Omega
beta0  <- rep(0,n_beta)
gamma0 <- rep(0,n_gamma)
u0     <- rep(0,n_subjects)
v0     <- rep(0,n_subjects)
Omega0 <- c(beta0,gamma0,u0,v0)
# Phi
theta_u<- 0.1
theta_v<- 0.1
rho    <- 0.1
Phi0   <- c(theta_u, theta_v, rho)

# General Data Structures: not modified in the looèp, only defined once
# Design Matrices [X1,X2,Z1,Z2]
X1 = model.matrix(~ SESSO_rec + ADERENTE_rec + scale(etaEvent_rec) + comorbidity_rec)[,-1]
X2 = model.matrix(~ SESSO_term + ADERENTE_term + scale(etaEvent_term) + comorbidity_term)[,-1]
Z1 = matrix(0,dim(X1)[1],n_subjects)
Z2 = diag(n_subjects)

codici <- unique(ID_term)
for(i in 1:length(codici)){
    current <- codici[i]
    for(j in 1:dim(X1)[1]){
      if(ID_rec[j]==current)
        Z1[j,i]=1
    }
}

# clean memory
rm(SESSO_rec,ADERENTE_rec,etaEvent_rec,comorbidity_rec,SESSO_term,ADERENTE_term,etaEvent_term,comorbidity_term)
gc()

# Observation Matrix  [OBS]
fill1<-matrix(0,dim(X1)[1],dim(X1)[2])
fill2<-matrix(0,dim(X2)[1],dim(X2)[2])
fill3<-matrix(0,dim(Z1)[1],dim(Z1)[2])
fill4<-matrix(0,dim(Z2)[1],dim(Z2)[2])
OBS1<-cbind(X1,fill1)
OBS2<-cbind(fill2,X2)
OBS3<-cbind(Z1,fill3)
OBS4<-cbind(fill4,Z2)
OBS <-cbind(rbind(OBS1,OBS2),rbind(OBS3,OBS4))
rm(OBS1,OBS2,OBS3,OBS4,fill1,fill2,fill3,fill4)

# Matrices for Phi update
K1<-matrix(0,2*n_subjects,2*n_subjects)
diag(K1[1:n_subjects,1:n_subjects]) <- 1
K2<-matrix(0,2*n_subjects,2*n_subjects)
K2[col(K2)+row(K2)-ncol(K2)==1L] <- 1
K3<-matrix(0,2*n_subjects,2*n_subjects)
diag(K3[(n_subjects+1):nrow(K3),(n_subjects+1):ncol(K3)]) <- 1

# Dummy vector
DELTA_R <- data$event
DELTA_D <- dataDeath$cens

# Convergence 
conv1 <- FALSE      # External loop
maxit1<- 30
maxit2<- 30
it1   <- 0
treshold <- 1e-3

# Start loop
while(!conv1 & it1<maxit1){
  it2 = 0
  conv2 <- FALSE 
  while(!conv2 & it2<maxit2){
    # G first   (depend on beta and u --> in the loop)
    Q1      <- matrix(0,dim(X1)[1],dim(X1)[1])
    diag(Q1)<- exp(X1%*%beta0 + Z1%*%u0)
    E1      <- matrix(0,dim(X1)[1],dim(X1)[1])
    diag(E1)<- DELTA_R/sum(diag(Q1))
    F1      <- matrix(0,dim(X1)[1],dim(X1)[1])
    F1[lower.tri(F1)] <- 1
    S1      <- matrix(0,dim(X1)[1],dim(X1)[1])
    diag(S1)<- cumsum(diag(E1))
    
    Q2      <- matrix(0,dim(X2)[1],dim(X2)[1])
    diag(Q2)<- exp(X2%*%gamma0 + Z2%*%v0)
    E2      <- matrix(0,dim(X2)[1],dim(X2)[1])
    diag(E2)<- DELTA_D/sum(diag(Q2))
    F2      <- matrix(0,dim(X2)[1],dim(X2)[1])
    F2[lower.tri(F2)] <- 1
    S2      <- matrix(0,dim(X2)[1],dim(X2)[1])
    diag(S2)<- cumsum(diag(E2))
    
    D1<- Q1%*%S1 - Q1%*%F1%*%(E1%*%E1)%*%t(F1)%*%t(Q1)
    D2<- Q2%*%S2 - Q2%*%F2%*%(E2%*%E2)%*%t(F2)%*%t(Q2)
    fill1<-matrix(0,dim(X1)[1],dim(X2)[1])
    fill2<-matrix(0,dim(X2)[1],dim(X1)[1])
    
    D<-rbind(cbind(D1,fill1),cbind(fill2,D2))
    rm(D1,D2,fill1,fill2)
    gc()
    
    G_first = t(OBS)%*%D%*%OBS
    
    # G second (depend on theta_u,theta_v,rho --> in the loop)
    G_second = matrix(0, n_beta+n_gamma+2*n_subjects,n_beta+n_gamma+2*n_subjects)
    diag(G_second[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects),n_beta+n_gamma+1:(n_beta+n_gamma+n_subjects)])=theta_v^2
    diag(G_second[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects),(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)])=theta_u^2
    diag(G_second[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects),(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)])=-theta_v*theta_u*rho
    diag(G_second[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects),(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)])=-theta_v*theta_u*rho
    G_second = (1/(theta_u^2*theta_v^2*(1-rho^2)))*G_second
    
    # G
    G = G_first + G_second
    
    # Gradient dL/dOmega [upd_step]
    dL1_dEta <-DELTA_R - Q1%*%F1%*%E1%*%rep(1,dim(X1)[1])
    dL1_dZeta<-DELTA_D - Q2%*%F2%*%E2%*%rep(1,dim(X2)[1])
    dBeta <-t(X1)%*%dL1_dEta 
    dGamma<-t(X2)%*%dL1_dZeta
    dU    <-t(Z1)%*%dL1_dEta  - (theta_v^2*u0 - rho*theta_u*theta_v*v0)/(theta_u^2*theta_v^2*(1-rho)^2)
    dV    <-t(Z2)%*%dL1_dZeta - (theta_v^2*v0 - rho*theta_u*theta_v*u0)/(theta_u^2*theta_v^2*(1-rho)^2)
    upd_step <- c(dBeta,dGamma,dU,dV)
    
    # clean memory
    rm(D,E1,E2,F1,F2,G_first,G_second,Q1,Q2,S1,S2)
    gc()
    
    # update Omega
    invG   <- chol2inv(chol(G))
    Omega0 <- Omega0 + invG%*%upd_step
    
    # update iteration
    it2 <- it2+1
    
    # check convergence
    if(norm(invG%*%upd_step,type="2")<treshold)
      conv2=TRUE
  }
  # update Phi
  Q    <-matrix(Omega0[(n_beta+n_gamma+1):length(Omega0)],2*n_subjects,1)%*%matrix(Omega0[(n_beta+n_gamma+1):length(Omega0)],1,2*n_subjects)
  T1   <- sum(diag(K1%*%(invG[(n_beta+n_gamma+1):nrow(invG),(n_beta+n_gamma+1):ncol(invG)] + Q)))
  T2   <- 0.5*sum(diag(K2%*%(invG[(n_beta+n_gamma+1):nrow(invG),(n_beta+n_gamma+1):ncol(invG)] + Q)))
  T3   <- sum(diag(K3%*%(invG[(n_beta+n_gamma+1):nrow(invG),(n_beta+n_gamma+1):ncol(invG)] + Q)))
  
  rm(Q,T1,T2,T3)
  gc()
  
  converged <- norm(Phi0 - c(T1/n_subjects,T3/n_subjects, T2/sqrt(T1*T3)),type="2")
  Phi0 <- c(T1/n_subjects,T3/n_subjects, T2/sqrt(T1*T3))
  
  # update iteration
  it1 <- it1+1
  
  # check convergence
  if(converged<treshold)
    conv1=TRUE
}

# Std Errors for obtained estimates


# CI for HR, p value, Survival probabilities
