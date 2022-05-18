############################################################################################
##################        Ng, Tawiah, MacLachlan, Gopalan         ##########################
############################################################################################

# Load packages
print("load packages")
library(data.table)
library(Rcpp)
library(parallel)

# Load Data
print("load data")
load("dataRec.RData")
load("dataDeath.RData")

# RCpp function to speed up matrix multiplication (check functions path)
print("load cpp functions")
library(Rcpp)
sourceCpp("MatProd.cpp")
sourceCpp("REMLupdate.cpp")

# Order data according to observed gap times
print("order gap times")
order1 <- order(data$GapEvent)
order2 <- order(dataDeath$GapEvent)
data <- data[order1,]
dataDeath <- dataDeath[order2,]

# Variables
print("define variables")
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

# Design parameters
n_beta = 4                         # recurrent event covariates
n_gamma = 4                        # terminal event covariates
n_subjects = dim(dataDeath)[1]     # number of units subject to random effect 

# Initial values: set to zero and 
print("set initial values")
# Omega
#beta0  <- rep(0,n_beta)
#gamma0 <- rep(0,n_gamma)
#u0     <- rep(0,n_subjects)
#v0     <- rep(0,n_subjects)
#Omega0 <- c(beta0,gamma0,u0,v0)
#rm(beta0,gamma0,u0,v0)
load("Ninth_run/Omega.RData")
# Phi
#theta_u<- 0.7
#theta_v<- 0.7
#rho    <- 0.5
#Phi0   <- c(theta_u, theta_v, rho)
#rm(theta_u,theta_v,rho)
load("Ninth_run/Phi.RData")

# General Data Structures: not modified in the loop, only defined once
# Design Matrices [X1,X2,Z1,Z2]
print("design matrices")
X1 = model.matrix(~ SESSO_rec + ADERENTE_rec + etaEvent_rec + comorbidity_rec)[,-1]
X2 = model.matrix(~ SESSO_term + ADERENTE_term + etaEvent_term + comorbidity_term)[,-1]
Z2 = diag(n_subjects)
Z2 <-Z2[order2,]

loadZ1 <- TRUE
if(loadZ1){
load("RE_design_matrix_FFU.RData")
} else {
print("computing Z1..")
Z1 <- matrix(0,dim(X1)[1],n_subjects)
codici <- levels(ID_rec)
for(i in 1:length(codici)){
    current <- codici[i]
    for(j in 1:dim(X1)[1]){
      if(ID_rec[j]==current)
        Z1[j,i]=1
    }
  }
save(Z1, file="RE_design_matrix_FFU.RData")
}

# Clean memory
rm(SESSO_rec,ADERENTE_rec,etaEvent_rec,comorbidity_rec,SESSO_term,ADERENTE_term,etaEvent_term,comorbidity_term)
gc()

# Censoring vectors
DELTA_R <- data$event
DELTA_D <- dataDeath$cens

# Convergence 
conv1 <- FALSE      # External loop
maxit1<- 100
maxit2<- 10
it1   <- 0
treshold <- 1e-3
treshold1<- 1e-6

# Start loop
print("start loop!")
while(!conv1 & it1<maxit1){
  cat(paste("\nCurrent outer it: ", it1))
  print("Current Estimates:")
  print(paste("Beta: ", Omega0[1:n_beta]))
  print(paste("Gamma: ", Omega0[(n_beta+1):(n_beta+n_gamma)]))
  print(paste("u0 NaN: ", table(is.nan(Omega0[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)]))))
  print(paste("v0 NaN: ", table(is.nan(Omega0[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)]))))
  print(paste("Phi: ", Phi0))
  
  it2 = 0
  conv2 <- FALSE 
  while(!conv2 & it2<maxit2){
    cat(paste("\nCurrent inner it: ", it2))
    print("Current Estimates:")
    print(paste("Beta: ", Omega0[1:n_beta]))
    print(paste("Gamma: ", Omega0[(n_beta+1):(n_beta+n_gamma)]))
    print(paste("u0 NaN: ", table(is.nan(Omega0[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)]))))
    print(paste("v0 NaN: ", table(is.nan(Omega0[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)]))))
    print(paste("Phi: ", Phi0))
    
    # G first   (depend on beta and u --> in the loop)
    Q1      <- matrix(0,dim(X1)[1],dim(X1)[1])
    diag(Q1)<- exp(MatProd(X1,Omega0[1:n_beta]) + MatProd(Z1,Omega0[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)]))
    E1      <- matrix(0,dim(X1)[1],dim(X1)[1])
    diag(E1)<- DELTA_R/rev(cumsum(rev(diag(Q1))))
    F1      <- matrix(0,dim(X1)[1],dim(X1)[1])
    F1[lower.tri(F1)] <- 1
    diag(F1)<-1
    S1      <- matrix(0,dim(X1)[1],dim(X1)[1])
    diag(S1)<- cumsum(diag(E1))
    D1<- MatProd(Q1,S1) - MatProd(Q1,MatProd(F1,MatProd(E1^2,MatProd(t(F1),Q1))))
    dL1_dEta <-DELTA_R - MatProd(Q1,MatProd(F1,MatProd(E1,rep(1,dim(X1)[1]))))
    rm(Q1,E1,F1,S1)
    gc()
    
    Q2      <- matrix(0,dim(X2)[1],dim(X2)[1])
    diag(Q2)<- exp(MatProd(X2,Omega0[(n_beta+1):(n_beta+n_gamma)])+ MatProd(Z2,Omega0[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)]))
    E2      <- matrix(0,dim(X2)[1],dim(X2)[1])
    diag(E2)<- DELTA_D/rev(cumsum(rev(diag(Q2))))
    F2      <- matrix(0,dim(X2)[1],dim(X2)[1])
    F2[lower.tri(F2)] <- 1
    diag(F2)<-1
    S2      <- matrix(0,dim(X2)[1],dim(X2)[1])
    diag(S2)<- cumsum(diag(E2))
    D2<- MatProd(Q2,S2) - MatProd(Q2,MatProd(F2,MatProd(E2^2,MatProd(t(F2),Q2))))
    dL1_dZeta<-DELTA_D - MatProd(Q2,MatProd(F2,MatProd(E2,rep(1,dim(X2)[1]))))
    rm(Q2,E2,F2,S2)
    gc()
    
    # G2
    temp <- matrix(0,2,2)
    temp[1,1]<-1/(Phi0[1]*(1-Phi0[3]^2))
    temp[1,2]<- -Phi0[3]/(sqrt(Phi0[1]*Phi0[2])*(1-Phi0[3]^2))
    temp[2,1]<- -Phi0[3]/(sqrt(Phi0[1]*Phi0[2])*(1-Phi0[3]^2))
    temp[2,2]<-1/(Phi0[2]*(1-Phi0[3]^2))
    G2 <- kronecker(temp,diag(nrow=n_subjects))
    
    # G
    G = matrix(0, n_beta+n_gamma+2*n_subjects,n_beta+n_gamma+2*n_subjects)
    G[1:n_beta,1:n_beta]=MatProd(t(X1),MatProd(D1,X1))
    G[(n_beta+1):(n_beta+n_gamma),(n_beta+1):(n_beta+n_gamma)]=MatProd(t(X2),MatProd(D2,X2))
    G[1:n_beta,(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)]=MatProd(t(X1),MatProd(D1,Z1))
    G[(n_beta+1):(n_beta+n_gamma),(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)]=MatProd(t(X2),MatProd(D2,Z2))
    G[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects),1:n_beta]=MatProd(t(Z1),MatProd(D1,X1))
    G[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects),(n_beta+1):(n_beta+n_gamma)]=MatProd(t(Z2),MatProd(D2,X2))
    G[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects),(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)]=MatProd(t(Z1),MatProd(D1,Z1))
    G[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects),(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)]=MatProd(t(Z2),MatProd(D2,Z2))
    rm(D1,D2)
    gc()
    
    G[(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects),(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects)] = G[(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects),(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects)] + G2
    rm(G2)
    gc()
    
    # Gradient dL/dOmega [upd_step]
    dBeta <-MatProd(t(X1),dL1_dEta)
    dGamma<-MatProd(t(X2),dL1_dZeta)
    dU    <-MatProd(t(Z1),dL1_dEta)  - Omega0[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)]/(Phi0[1]*(1-Phi0[3]^2)) + (Phi0[3]/(sqrt(Phi0[1]*Phi0[2])*(1-Phi0[3]^2)))*Omega0[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)]
    dV    <-MatProd(t(Z2),dL1_dZeta) - Omega0[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)]/(Phi0[2]*(1-Phi0[3]^2)) + (Phi0[3]/(sqrt(Phi0[1]*Phi0[2])*(1-Phi0[3]^2)))*Omega0[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)]
    upd_step <- c(dBeta,dGamma,dU,dV)
    
    # clean memory
    rm(dBeta,dGamma,dU,dV)
    gc()
    
    # update Omega
    x      <- qr.solve(G,upd_step, tol=1e-10)
    Omega0 <- Omega0 + x
    # update iteration
    it2 <- it2+1
    
    # check convergence
    if(norm(x,type="2")<treshold1)
      conv2=TRUE
  }

  ## update Phi
  # invert G
  invG   <- chol2inv(chol(G))
  invG_bb<- invG[(n_beta+n_gamma+1):nrow(invG),(n_beta+n_gamma+1):ncol(invG)]  
  
  # update
  Phi1 <- REMLupdate(Omega0[(n_beta+n_gamma+1):length(Omega0)],invG_bb,n_subjects)
  converged <- norm(Phi0 - Phi1,type="2")
  Phi0 <- Phi1
  gc()
  
  # update iteration
  it1 <- it1+1
  
  # check convergence
  if(converged<treshold)
    conv1=TRUE
}

save(Omega0, file="Omega.RData")
save(Phi0, file="Phi.RData")

# Clean Memory
rm(G)

## Variance-covariance Matrices for obtained estimates
# Beta & Gamma
varBeta <- invG[1:n_beta,1:n_beta]
varGamma<- invG[(n_beta+1):(n_beta+n_gamma),(n_beta+1):(n_beta+n_gamma)]
# Theta_u2,Theta_v2,rho
Epsilon <- kronecker(matrix(c(Phi0[1],Phi0[3]*sqrt(Phi0[2]*Phi0[1]),Phi0[3]*sqrt(Phi0[2]*Phi0[1]),
                              Phi0[2]),2,2),diag(nrow=n_subjects))
dEpsInv_dThetaU2<-kronecker(matrix(c((-1)/(Phi0[1]^2*(1-Phi0[3]^2)),Phi0[3]/(2*sqrt(Phi0[2])*(1-Phi0[3]^2)*Phi0[1]^3),Phi0[3]/(2*sqrt(Phi0[2])*(1-Phi0[3]^2)*Phi0^3),
                                     0),2,2),diag(nrow=n_subjects))
dEpsInv_dThetaV2<-kronecker(matrix(c(0,Phi0[3]/(2*sqrt(Phi0[1])*(1-Phi0[3]^2)*Phi0[2]^3),Phi0[3]/(2*sqrt(Phi0[1])*(1-Phi0[3]^2)*Phi0[2]^3),
                                     (-1)/(Phi0[2]^2*(1-Phi0[3]^2))),2,2),diag(nrow=n_subjects))
dEpsInv_dRho<-kronecker(matrix(c((2*Phi0[3])/(Phi0[1]*(1-Phi0[3]^2)^2),-(Phi0[3]^2+1)/(sqrt(Phi0[1]*Phi0[2])*(1-Phi0[3]^2)^2),-(Phi0[3]^2+1)/(sqrt(Phi0[1]*Phi0[2])*(1-Phi0[3]^2)^2),
                                 (2*Phi0[3])/(Phi0[2]*(1-Phi0[3]^2)^2)),2,2),diag(nrow=n_subjects))

J1    <- MatProd(invG[(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects),(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects)],dEpsInv_dThetaU2)
J2    <- MatProd(Epsilon,dEpsInv_dThetaU2)
J3    <- MatProd(invG[(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects),(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects)],dEpsInv_dThetaV2)
J4    <- MatProd(Epsilon,dEpsInv_dThetaV2)
J5    <- MatProd(invG[(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects),(n_beta+n_gamma+1):(n_beta+n_gamma+2*n_subjects)],dEpsInv_dRho)
J6    <- MatProd(Epsilon,dEpsInv_dRho)

A     <- matrix(0,3,3)
A[1,1]<- sum(diag((J1-J2)))^2
A[1,2]<- sum(diag((MatProd(J1,J3)+MatProd(J2,J4)-2*MatProd(J1,J4))))
A[1,3]<- sum(diag((MatProd(J1,J5)+MatProd(J2,J6)-2*MatProd(J1,J6))))
A[2,2]<- sum(diag((J3-J4)))^4
A[2,3]<- sum(diag((MatProd(J3,J5)+MatProd(J4,J6)-2*MatProd(J3,J6))))
A[3,3]<- sum(diag((J5-J6)))^2
A[3,1]<- A[1,3]
A[2,1]<- A[1,2]
A[3,2]<- A[2,3]
VarPhi<- 2*chol2inv(chol(A))
rm(A,J1,J2,J3,J4,J5,J6,Epsilon,dEpsInv_dRho,dEpsInv_dThetaU2,dEpsInv_dThetaV2)


# CI for HR, Survival probabilities
# HR and 95pc CI
summary <- data.frame(
  Estimate=c(Omega0[1:(n_beta)],Omega0[(n_beta+1):(n_beta+n_gamma)],Phi0),
  StDev   =c(sqrt(diag(varBeta)),sqrt(diag(varGamma)),sqrt(diag(VarPhi))),
  HR      =exp(c(Omega0[1:(n_beta)],Omega0[(n_beta+1):(n_beta+n_gamma)],Phi0)),
  L95     = exp(c(Omega0[1:(n_beta)],Omega0[(n_beta+1):(n_beta+n_gamma)],Phi0) - 1.96*c(sqrt(diag(varBeta)),sqrt(diag(varGamma)),sqrt(diag(VarPhi)))),
  U95     = exp(c(Omega0[1:(n_beta)],Omega0[(n_beta+1):(n_beta+n_gamma)],Phi0) + 1.96*c(sqrt(diag(varBeta)),sqrt(diag(varGamma)),sqrt(diag(VarPhi)))))
row.names(summary) <- c("Beta1","Beta2","Beta3","Beta4","Gamma1","Gamma2","Gamma3","Gamma4","ThetaU2","ThetaV2","rho")
summary

# Save Results
FFU_try <- list(beta=Omega0[1:n_beta],gamma=Omega0[(n_beta+1):(n_beta+n_gamma)],
                 u=Omega0[(n_beta+n_gamma+1):(n_beta+n_gamma+n_subjects)],
                 v=Omega0[(n_beta+n_gamma+n_subjects+1):(n_beta+n_gamma+2*n_subjects)],
                 thetaU2=Phi0[1],thetaV2=Phi0[2],rho=Phi0[3],
                 HR=summary$HR,summary=summary,VarBeta=varBeta,varGamma=varGamma,VarPhi=VarPhi)
save(FFU_try,file='FFU_ACE.RData')
