#####################################################################################################
#########################  Nonparametric Discrete Bivariate Frailty #################################
#####################################################################################################
rm(list=ls())

# load packages
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

## data loading
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

## Set Convergence Indices
# Count
count <- 0

# epsilon
eps_conv = 1e-3
eps      = 1e5


## Set Environmental variables
# Shrinkage: if True, build a grid of K points according to Mu and Sigma, then shrinks it
#            if False, RangeK defines the fixed values of K for which the model will be evaluated
   Shrinkage = TRUE
   K= 1000
   if (!Shrinkage)
     RangeK <- 2:7
# Constrained: if True, applied constrained optimization in P update
   Constrained=FALSE

# Print Ng estimates for initialization of Mu and Sigma
load("FFU_ACE_runs/Tenth_run_FFU_ACE/FFU_try.RData")

# Initialize Mu and Sigma
library(MASS)
Sigma <- matrix(c(0.12,0.0,0.0,1.4),nrow = 2, ncol = 2)
mu=c(0,0)

## Estimation
if(Shrinkage){
  
  # Grid Initialization
  P <- mvrnorm(K,mu,Sigma)
  w<-dmvnorm(P,mu,Sigma)
  w<-w/sum(w)
  P_show<-P
  P_show[,1]<-P[,1]-rep(w%*%P[,1],length(P[,1]))
  P_show[,2]<-P[,2]-rep(w%*%P[,2],length(P[,2]))

  # Grid shrinking
  is_near<-TRUE
    while(is_near){
   D<-dist(P)
   D<-as.matrix(D)
   D[upper.tri(D)]<-10
   diag(D)<-10
   out<-which(D == min(D), arr.ind = TRUE)
   if(D[out][1]<0.25){
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

  # Initialize H0R(t): Breslow estimator
  #YY1     = ((exp( as.numeric(X1 %*% beta) + P_off1)) %*% risk_index1)[1,]
  #haz1    = m1/YY1
  #cumhaz1$hazard = cumsum(haz1)

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

  # Initializa H0D(t): Breslow estimator
  #YY2     = ((exp( as.numeric(X2 %*% gamma) + P_off2)) %*% risk_index2)[1,]
  #haz2    = m2/YY2
  #cumhaz2$hazard = cumsum(haz2)

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

  ## Support Functions
  # Bisection
  Bisection<-function (f, a, b, num = 10, eps = 1e-05, flag=T) 
  { out=f(a)
  h = abs(b - a)/num
  i = 0
  j = 0
  a1 = b1 = 0
  while (i <= num) {
    a1 = a + i * h
    b1 = a1 + h
    if (f(a1) == 0) {
      if(flag){
      print(a1)
      print(f(a1))}
      out=a1
    }
    else if (f(b1) == 0) {
      if(flag){
      print(b1)
      print(f(b1))}
      out=b1
    }
    else if (f(a1) * f(b1) < 0) {
      repeat {
        if (abs(b1 - a1) < eps) 
          break
        x <- (a1 + b1)/2
        if (f(a1) * f(x) < 0) 
          b1 <- x
        else a1 <- x
      }
      if(flag){
      print(j + 1)}
      j = j + 1
      if(flag){
      print((a1 + b1)/2)
      print(f((a1 + b1)/2))}
      out=(a1 + b1)/2
    }
    i = i + 1
  }
  if( flag){
  if (j == 0) 
    print("finding root is fail")
  else print("finding root is successful")
  }
  return (out)
}
  # LogLikelihood
  LOGL<-function(Z,w,haz1,cumhaz1,beta,haz2,cumhaz2, gamma,P,E_formula1,E_formula2,E_haz1,E_haz2){
  lw<-sum(Z%*%matrix(log(w)))
  lr<- sum(Z*(E_haz1-E_formula1))
  ld<- sum(Z*(E_haz2-E_formula2))
  return (lw+lr+ld)
}

  # Initialize saving structure
  Saved <- list()

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
    if(D[out][1]<0.25){
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
  
  # Define auxiliary functions if Constrained P optimization
  if(Constrained & K>1){
  solveP1 <- function(p1,K,Z,D1,E_part1,w){
    out <-rep(0,K)
    out[1] <- p1
    for(l in 2:K){
      out[l] <-log((sum(Z[,l]*D1)+(w[l]/w[1])*(-sum(Z[,1]*D1)+sum(Z[,1]*(E_part1*exp(p1)))))/(sum(Z[,l]*E_part1)))
    }
    return (out)
  }
  
  solveP2 <- function(p2,K,Z,D,E_part2,w){
    out <-rep(0,K)
    out[1] <- p2
    for(l in 2:K){
      out[l] <-log((sum(Z[,l]*D)+(w[l]/w[1])*(-sum(Z[,1]*D)+sum(Z[,1]*(E_part2*exp(p2)))))/(sum(Z[,l]*E_part2)))
    }
    return (out)
  }
  
  to_opt1<-function(p1){
    out = solveP1(p1,K,Z,D1,E_part1,w)
    return (sum(out*w))
  }
  
  to_opt2<-function(p2){
    out = solveP2(p2,K,Z,D2,E_part2,w)
    return (sum(out*w))
  }
  
  CheckLim1<- function(K,Z,D1,E_part1,w){
    out<- (max(log((rep(Z[,1]%*%D1,K-1)-(w[1]/w[2:K])*(D1%*%Z[,2:K]))/rep(Z[,1]%*%E_part1,K-1)),na.rm = T))
    if(is.finite(out))
      return (out+1e-06)
    else
      return (-5)
    }
  
  CheckLim2<- function(K,Z,D2,E_part2,w){
    out<-(max(log((rep(Z[,1]%*%D2,K-1)-(w[1]/w[2:K])*(D2%*%Z[,2:K]))/rep(Z[,1]%*%E_part2,K-1)),na.rm = T))
    if(is.finite(out))
      return (out+1e-06)
    else
      return (-5)
    }
  
  ## Constrained Optimization
  # Check limits
  inf1<-CheckLim1(K,Z,D1,E_part1,w)
  inf2<-CheckLim2(K,Z,D2,E_part2,w)

  # Optimize
  P[1,1] <- Bisection(to_opt1,inf1,10,flag = F)
  P[1,2] <- Bisection(to_opt2,inf2,10,flag = F)
  
  # Support points abscissa
  P[,1] <- solveP1(P[1,1],K,Z,D1,E_part1,w)
 
  # Support points ordinata
  P[,2] <- solveP2(P[1,2],K,Z,D2,E_part2,w)
  }
  
  ## Unconstrained Optimization
  else{
    P[,1]   <- log(( as.numeric(D1) %*% Z )/( E_part1 %*% Z))
    P[,2]   <- log(( as.numeric(D2) %*% Z )/( E_part2 %*% Z))
  }
  P_show    <-P
  P_show[,1]<-P[,1]-w%*%P[,1]
  P_show[,2]<-P[,2]-w%*%P[,2]
  
  
  #P_off1 <- (P[,1][belonging][groups1]) 
  #P_off2 <- (P[,2][belonging][groups2])
  
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
  
  # Estimate hospitalization baseline hazard and cumulative baseline hazard
  #YY1     = as.numeric( (exp( beta %*% t(X1) + P_off1)) %*% risk_index1 )
  #haz1    = m1/YY1
  #cumhaz1$hazard = cumsum(haz1)
  
  # Estimate Cumulative Hazard and Hazard - Recurrent
  s1   <- survfit(temp_model1,data=data)
  cumhaz1$hazard = s1$cumhaz
  haz1$hazard = diff(c(0,cumhaz1$hazard)) 
  for(j in 1:length(haz1$hazard)){
    if(haz1$hazard[j]==0)
      haz1$hazard[j]<-haz1$hazard[j-1]
  }
  
  # Estimate death baseline hazard and cumulative baseline hazard
  #YY2     = as.numeric( (exp( gamma %*% t(X2) + P_off2 )) %*% risk_index2 )
  #haz2    = m2/YY2
  #cumhaz2$hazard = cumsum(haz2)
  
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
  LogL= 3*length(w)-LOGL(Z=Z,w=w, haz1 = haz1,beta = beta,cumhaz1 = cumhaz1,
               haz2 = haz2, gamma = gamma, cumhaz2 = cumhaz2,
               P=P,E_formula1 = E_formula1,E_formula2 = E_formula2,
               E_haz1 = E_haz1, E_haz2=E_haz2)
  # Save 
  temp_list<-list("modelR"=temp_model1,"modelT"=temp_model2,"w"=w,"P"=P,"P_show"=P_show,"cumhaz1"=cumhaz1,"cumhaz2"=cumhaz2,"LogL"=LogL,"Table"=table(belonging))
  Saved[[count]]<-temp_list
  
 }
  
  # No Shrinking
} else {
  
  K_old <- K
  
  for(currentK in RangeK){
  # Grid Initialization
  K <- currentK
  P <- mvrnorm(K,mu,Sigma)
  w<-dmvnorm(P,mu,Sigma)
  w<-w/sum(w)
  
  # To Proper K
  #while(K>currentK){
  #  D<-dist(P)
  #  D<-as.matrix(D)
  #  D[upper.tri(D)]<-10
  #  diag(D)<-10
  #  out<-which(D == min(D), arr.ind = TRUE)
  #  P[out[1,2],]=(P[out[1,2],]+P[out[1,1],])/2
  #  P<-P[-out[1,1],]
  #  w[out[1,2]]<-w[out[1,2]]+w[out[1,1]]
  #  w<-w[-out[1,1]]
  #  w<-w/sum(w)
  #  K<-K-1
  #}
  
  P_show<-P
  P_show[,1]<-P[,1]-rep(w%*%P[,1],length(P[,1]))
  P_show[,2]<-P[,2]-rep(w%*%P[,2],length(P[,2]))
  
  # Assign patient to random frailty, built frailties vectors
  P_index <- sample(1:K,size=N,replace = T, prob = w)
  P_off1  <- P[,1][P_index[groups1]]
  P_off2  <- P[,2][P_index[groups2]]
  
  # Initialize H0R(t): Breslow estimator
  #YY1     = ((exp( as.numeric(X1 %*% beta) + P_off1)) %*% risk_index1)[1,]
  #haz1    = m1/YY1
  #cumhaz1$hazard = cumsum(haz1)
  
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
  
  # Initializa H0D(t): Breslow estimator
  #YY2     = ((exp( as.numeric(X2 %*% gamma) + P_off2)) %*% risk_index2)[1,]
  #haz2    = m2/YY2
  #cumhaz2$hazard = cumsum(haz2)
  
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
  
  ## Support Functions
  # Bisection
  Bisection<-function (f, a, b, num = 10, eps = 1e-05, flag=T) 
  { out=f(a)
  h = abs(b - a)/num
  i = 0
  j = 0
  a1 = b1 = 0
  while (i <= num) {
    a1 = a + i * h
    b1 = a1 + h
    if (f(a1) == 0) {
      if(flag){
        print(a1)
        print(f(a1))}
      out=a1
    }
    else if (f(b1) == 0) {
      if(flag){
        print(b1)
        print(f(b1))}
      out=b1
    }
    else if (f(a1) * f(b1) < 0) {
      repeat {
        if (abs(b1 - a1) < eps) 
          break
        x <- (a1 + b1)/2
        if (f(a1) * f(x) < 0) 
          b1 <- x
        else a1 <- x
      }
      if(flag){
        print(j + 1)}
      j = j + 1
      if(flag){
        print((a1 + b1)/2)
        print(f((a1 + b1)/2))}
      out=(a1 + b1)/2
    }
    i = i + 1
  }
  if( flag){
    if (j == 0) 
      print("finding root is fail")
    else print("finding root is successful")
  }
  return (out)
  }
  # LogLikelihood
  LOGL<-function(Z,w,haz1,cumhaz1,beta,haz2,cumhaz2, gamma,P,E_formula1,E_formula2,E_haz1,E_haz2){
    lw<-sum(Z%*%matrix(log(w)))
    lr<- sum(Z*(E_haz1-E_formula1))
    ld<- sum(Z*(E_haz2-E_formula2))
    return (lw+lr+ld)
  }
  
  # Initialize saving structure
  Saved <- list()
  
  ## Set Convergence Indices
  # Count
  count <- 0
  
  # epsilon
  eps_conv = 1e-3
  eps      = 1e5
  
  # Start loop
  while (eps_conv < eps & count < 100){
    
    # Save current w estimates
    w_old <- w
    
    # Save current Z estimates
    Z_old <- Z
    
    # Save current P estimates
    P_old <-P
    
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
    
    # Vector of proportions
    if(K>1)
      w <- (colSums(Z))/ N
    else{
      w <- 1 
      P <- matrix(c(0,0),nrow = 1,ncol = 2)
    }
    
    # Define auxiliary functions if Constrained P optimization
    if(Constrained & K>1){
      solveP1 <- function(p1,K,Z,D1,E_part1,w){
        out <-rep(0,K)
        out[1] <- p1
        for(l in 2:K){
          out[l] <-log((sum(Z[,l]*D1)+(w[l]/w[1])*(-sum(Z[,1]*D1)+sum(Z[,1]*(E_part1*exp(p1)))))/(sum(Z[,l]*E_part1)))
        }
        return (out)
      }
      
      solveP2 <- function(p2,K,Z,D,E_part2,w){
        out <-rep(0,K)
        out[1] <- p2
        for(l in 2:K){
          out[l] <-log((sum(Z[,l]*D)+(w[l]/w[1])*(-sum(Z[,1]*D)+sum(Z[,1]*(E_part2*exp(p2)))))/(sum(Z[,l]*E_part2)))
        }
        return (out)
      }
      
      to_opt1<-function(p1){
        out = solveP1(p1,K,Z,D1,E_part1,w)
        return (sum(out*w))
      }
      
      to_opt2<-function(p2){
        out = solveP2(p2,K,Z,D2,E_part2,w)
        return (sum(out*w))
      }
      
      CheckLim1<- function(K,Z,D1,E_part1,w){
        out<- (max(log((rep(Z[,1]%*%D1,K-1)-(w[1]/w[2:K])*(D1%*%Z[,2:K]))/rep(Z[,1]%*%E_part1,K-1)),na.rm = T))
        if(is.finite(out))
          return (out+1e-06)
        else
          return (-5)
      }
      
      CheckLim2<- function(K,Z,D2,E_part2,w){
        out<-(max(log((rep(Z[,1]%*%D2,K-1)-(w[1]/w[2:K])*(D2%*%Z[,2:K]))/rep(Z[,1]%*%E_part2,K-1)),na.rm = T))
        if(is.finite(out))
          return (out+1e-06)
        else
          return (-5)
      }
      
      ## Constrained Optimization
      # Check limits
      inf1<-CheckLim1(K,Z,D1,E_part1,w)
      inf2<-CheckLim2(K,Z,D2,E_part2,w)
      
      # Optimize
      P[1,1] <- Bisection(to_opt1,inf1,10,flag = F)
      P[1,2] <- Bisection(to_opt2,inf2,10,flag = F)
      
      # Support points abscissa
      P[,1] <- solveP1(P[1,1],K,Z,D1,E_part1,w)
      
      # Support points ordinata
      P[,2] <- solveP2(P[1,2],K,Z,D2,E_part2,w)
    }
    
    ## Unconstrained Optimization
    else{
      P[,1]   <- log(( as.numeric(D1) %*% Z )/( E_part1 %*% Z))
      P[,2]   <- log(( as.numeric(D2) %*% Z )/( E_part2 %*% Z))
    }
    
    P_Show    <-P
    P_show[,1]<-P[,1]-w%*%P[,1]
    P_show[,2]<-P[,2]-w%*%P[,2]
    
    
    #P_off1 <- (P[,1][belonging][groups1]) 
    #P_off2 <- (P[,2][belonging][groups2])
    
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
    
    # Estimate hospitalization baseline hazard and cumulative baseline hazard
    #YY1     = as.numeric( (exp( beta %*% t(X1) + P_off1)) %*% risk_index1 )
    #haz1    = m1/YY1
    #cumhaz1$hazard = cumsum(haz1)
    
    # Estimate Cumulative Hazard and Hazard - Recurrent
    s1   <- survfit(temp_model1,data=data)
    cumhaz1$hazard = s1$cumhaz
    haz1$hazard = diff(c(0,cumhaz1$hazard)) 
    for(j in 1:length(haz1$hazard)){
      if(haz1$hazard[j]==0)
        haz1$hazard[j]<-haz1$hazard[j-1]
    }
    
    # Estimate death baseline hazard and cumulative baseline hazard
    #YY2     = as.numeric( (exp( gamma %*% t(X2) + P_off2 )) %*% risk_index2 )
    #haz2    = m2/YY2
    #cumhaz2$hazard = cumsum(haz2)
    
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
    LogL= 3*length(w)-LOGL(Z=Z,w=w, haz1 = haz1,beta = beta,cumhaz1 = cumhaz1,
                           haz2 = haz2, gamma = gamma, cumhaz2 = cumhaz2,
                           P=P,E_formula1 = E_formula1,E_formula2 = E_formula2,
                           E_haz1 = E_haz1, E_haz2=E_haz2)
    # Save 
    temp_list<-list("modelR"=temp_model1,"modelT"=temp_model2,"w"=w,"P"=P,"P_show"=P_show,"cumhaz1"=cumhaz1,"cumhaz2"=cumhaz2,"LogL"=LogL,"Table"=table(belonging))
    Saved[[count]]<-temp_list
    
  }
  
  assign(paste("Saved", currentK, sep = ""), Saved)
  }
}

############################################################################################
################################## Posterior inference #####################################
############################################################################################
## Grid inizialization
Pex <- mvrnorm(1000,mu,Sigma)
wex<-dmvnorm(Pex,mu,Sigma)
wex<-wex/sum(wex)
GridInEx <- data.frame("P1"=Pex[,1],"P2"=Pex[,2],"w"=wex, score=Pex[,1]+Pex[,2])
ggplot(data=GridInEx, aes(x=P1,y=P2,color=score))+
  theme_fivethirtyeight() +
  theme(axis.title = element_text())+
  geom_point(aes(size=w))+
  scale_color_gradient(low = "white", high = "red", guide = "none" )+
  scale_size(guide = "none")+
  xlab("u")+
  ylab("v")+
  labs(title = "Discrete Distribution - Initial Grid",
       caption = "Drug: ACE Inhibitors")

# Grid Evolution
ResultsByIteration <- data.frame()
FinalResults       <- data.frame()
for(i in 1:count){
  tempDF<-data.frame(Iteration=rep(i, length(Saved[[i]]$w)),
                     w=Saved[[i]]$w,
                     P1=(Saved[[i]]$P_show)[,1],
                     P2=(Saved[[i]]$P_show)[,2],
                     LogL=rep(Saved[[i]]$LogL,length(Saved[[i]]$w)),
                     score=(Saved[[i]]$P_show)[,1]+(Saved[[i]]$P_show)[,2])
  ResultsByIteration<-rbind(ResultsByIteration,tempDF)
  if(i==count)
    FinalResults<-rbind(FinalResults,tempDF)
}

graph1 = ResultsByIteration %>%
   ggplot(aes(x=P1, y=P2, size=w, color=score)) +
   geom_point(alpha = 0.7, stroke = 0) +
   theme_fivethirtyeight() +
   scale_size(range=c(2,12), guide="none") +
   theme(axis.title = element_text())+
   labs(title = "Discrete Distribution - Evolution",
       x = "u",
       y = "v",
       caption = "Drug: ACE Inhibitors") +
  scale_color_gradient(low = "white", high = "red", guide = "none" ) 

graph1.animation = graph1 +
  transition_time(Iteration) +
  labs(subtitle = "Iteration: {frame_time}")

animate(graph1.animation, height = 500, width = 800, fps = 30, duration = 10,
        end_pause = 60, res = 100)

anim_save("bubbleplotShrink.gif")

# Final Results
FinalResults["Point"]<-paste("P",rownames(FinalResults),sep="")
FinalResults["score"]<-FinalResults$P1+FinalResults$P2
plot1<-FinalResults %>%
  ggplot(aes(x=P1, y=P2)) +
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=factor(score,
                                                    labels = c("0.45","0.24","0.23","0.01","0.08"))),colour="black",pch=21) +
  geom_text(hjust=2,vjust=2.5,aes(size=0.01,label=Point))+
  theme(axis.title = element_text())+
  xlim(-0.5,1.25)+
  ylim(-1.5,4)+
  theme_fivethirtyeight() +
  scale_size(range=c(2,12), guide="none") +
  labs(title = "Discrete Distribution of Random Effects",
       x = "u",
       y = "v",
       caption = "Drug: ACE Inhibitors",
       fill = "Probability") +
  scale_fill_brewer(palette="Reds",
                    guide = guide_legend(override.aes = list(size =4) ))
plot1


## Survival Baselines Per class
library(tidyverse)
detach("package:MASS", unload=TRUE)
# Recurrent events
Ps<-Saved[[length(Saved)]]$P_show[,1]
temp<-data.frame("time"=Saved[[length(Saved)]]$cumhaz1$time,
                 "baseline"=Saved[[length(Saved)]]$cumhaz1$hazard)
for(i in 1:length(Ps)){
 temp[,ncol(temp)+1]<-Saved[[length(Saved)]]$cumhaz1$hazard*exp(Ps[i])
 colnames(temp)[ncol(temp)]<-paste("P",i,sep = "")
}
temp[,ncol(temp)+1]<-s1$std.err
colnames(temp)[ncol(temp)]<-"Err"
temp <- temp %>%
         select(time,baseline,P1,P2,P3,P4,P5) %>%
         gather(key="variable",value="value",-time)

ggplot(data=temp, aes(x=time, y=exp(-value)))+
  geom_line(aes(color=variable), size=1.0)+
  theme_fivethirtyeight() +
  theme(axis.title = element_text())+
  labs(title = "Stratified Baselines for hospitalization",
       x = "Time [days]",
       y = "Survival Probability",
       caption = "Drug: ACE Inhibitors",
       colour="Latent Class")+
  scale_color_manual(values = c("#FFFFFF","#A50F15","#DE2D26","#FEE5D9","#FB6A4A","#FCAE91"))


# Terminal events
Ps<-Saved[[length(Saved)]]$P_show[,2]
temp<-data.frame("time"=Saved[[length(Saved)]]$cumhaz2$time,
                 "baseline"=Saved[[length(Saved)]]$cumhaz2$hazard)
for(i in 1:length(Ps)){
  temp[,ncol(temp)+1]<-Saved[[length(Saved)]]$cumhaz2$hazard*exp(Ps[i])
  colnames(temp)[ncol(temp)]<-paste("P",i,sep = "")
}
temp[,ncol(temp)+1]<-s2$std.err
colnames(temp)[ncol(temp)]<-"Err"
temp <- temp %>%
  select(time,baseline,P1,P2,P3,P4,P5) %>%
  gather(key="variable",value="value",-time)

ggplot(data=temp, aes(x=time, y=exp(-value)))+
  geom_line(aes(color=variable), size=1.0)+
  theme_fivethirtyeight() +
  theme(axis.title = element_text())+
  labs(title = "Stratified Baselines for death",
       x = "Time [days]",
       y = "Survival Probability",
       caption = "Drug: ACE Inhibitors",
       colour="Latent Class")+
  scale_color_manual(values = c("#FFFFFF","#A50F15","#DE2D26","#FEE5D9","#FB6A4A","#FCAE91"))



# AIC
ResultsByIteration %>%
  ggplot(aes(x=Iteration, y=LogL)) +
  geom_line(alpha = 0.7) +
  theme_fivethirtyeight() +
  labs(title = "Minus LogLikelihood Across Iterations",
       x = "Iteration",
       y = "LogL",
       caption = "Drug: ACE Inhibitors") +
  
  scale_color_brewer(palette = "Set2")


##Compute Standard errors
summary(Saved[[length(Saved)]]$modelR)
summary(Saved[[length(Saved)]]$modelT)



######################################### RANGE ###########################################
ResultsByIteration <- data.frame()
FinalResults       <- data.frame()
number<-5
Saved<-eval(parse(text = paste("Saved",number,sep="")))
count<-length(Saved)
for(i in 1:count){
  tempDF<-data.frame(Iteration=rep(i, length(Saved[[i]]$w)),
                     w=Saved[[i]]$w,
                     P1=(Saved[[i]]$P_show)[,1],
                     P2=(Saved[[i]]$P_show)[,2],
                     LogL=rep(Saved[[i]]$LogL,length(Saved[[i]]$w)),
                     score=(Saved[[i]]$P_show)[,1]+(Saved[[i]]$P_show)[,2])
  ResultsByIteration<-rbind(ResultsByIteration,tempDF)
  if(i==count)
    FinalResults<-rbind(FinalResults,tempDF)
}

graph1 = ResultsByIteration %>%
  ggplot(aes(x=P1, y=P2, size=w, color=score)) +
  geom_point(alpha = 0.7, stroke = 0) +
  theme_fivethirtyeight() +
  scale_size(range=c(2,12), guide="none") +
  theme(axis.title = element_text())+
  labs(title = "Discrete Distribution - Evolution",
       x = "u",
       y = "v",
       caption = "Drug: ACE Inhibitors") +
  scale_color_gradient(low = "white", high = "red", guide = "none" ) 

graph1.animation = graph1 +
  transition_time(Iteration) +
  labs(subtitle = "Iteration: {frame_time}")

animate(graph1.animation, height = 500, width = 800, fps = 30, duration = 10,
        end_pause = 60, res = 100)

anim_save("bubbleplotRange5.gif")

## AICc
AICbyIteration<-data.frame()
FinalResults<-data.frame()
for(number in RangeK){
  Saved<-eval(parse(text = paste("Saved",number,sep="")))
  count<-length(Saved)
  for(i in 1:count){
    tempDF<-data.frame(Iteration=i,
                       AICc=Saved[[i]]$LogL,
                       modelK=paste("K=",number))
    AICbyIteration<-rbind(AICbyIteration,tempDF)
    if(i==count){
      tempDF <-data.frame(Iteration=rep(i, length(Saved[[i]]$w)),
                                  w=Saved[[i]]$w,
                                  P1=(Saved[[i]]$P_show)[,1],
                                  P2=(Saved[[i]]$P_show)[,2],
                                  LogL=rep(Saved[[i]]$LogL,length(Saved[[i]]$w)),
                                  score=(Saved[[i]]$P_show)[,1]+(Saved[[i]]$P_show)[,2],
                                  modelK=paste("K=",number, sep = ""))
      FinalResults<-rbind(FinalResults,tempDF)
    }
  }
}

AICbyIteration %>%
  ggplot(aes(x=Iteration, y=AICc, colour=modelK)) +
  geom_line(alpha = 0.7,size=1.0) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text())+
  labs(title = "AICc Across Iterations",
       x = "Iteration",
       y = "AICc",
       caption = "Drug: ACE Inhibitors") +
  scale_color_brewer(palette = "Set2")

ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=1.6, color="white", size=3.5)+
  theme_minimal()

FinalResults %>%
  ggplot(aes(x=modelK, y=LogL)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme(axis.title = element_blank())+
  geom_text(aes(label=round(LogL, digits = 0)), vjust=1.6, color="white", size=3.5)+
  theme_fivethirtyeight() +
  theme(axis.title = element_text())

FinalResults %>%
  ggplot(aes(x=P1, y=P2)) +
  facet_wrap(~modelK)+
  geom_point(alpha = 1, stroke = 0, aes(size=w,
                                        fill=score),colour="black",pch=21) +
  theme(axis.title = element_blank())+
  scale_size(range=c(2,12), guide="none") +
  scale_fill_gradient(low = "blue",high = "red",
                    guide ="none")
