########################## TimeDep Adherence Preprocessing #################################
######################################################################################################
### Preprocessing
## Load packages
rm(list = ls())
library(data.table)

# Load data
load("ACE_Inhibitors.RData")

## Arrange Dataset
# time at hospitalization events
new[!is.na(hosp), timeEvent:= data_prest - data_rif_ev]

# eta at hospitalization events
new[!is.na(hosp), etaEvent:= eta_Min]

# numbers of comorbidity at hospitalization events discharge
new[!is.na(hosp), comorbidity:=rowSums(.SD), .SDcols = 35:54]

# flag for event
new[!is.na(hosp), event:= 1]

# flag for type of censoring
new[!is.na(hosp), cens:=0]

# add censoring event per patient
data <- subset(new,(hosp>=1)||tipo_prest==30,select = c(COD_REG,event,hosp,timeEvent,SESSO,etaEvent,comorbidity,cens,data_rif_ev,data_prest,qt_pharma))
names(data)
codici<- unique(data$COD_REG)
for(i in 1:length(codici)){
  paz_corrente <- codici[i]
  temp <- data.frame(paz_corrente,0,(max(new[(COD_REG==paz_corrente) & (!is.na(hosp))]$hosp)+1),unique(new[COD_REG==paz_corrente]$timeOUT),unique(new[COD_REG==paz_corrente]$SESSO),
                     min(new[COD_REG==paz_corrente]$eta_Min) + 
                       as.integer(format(unique(new[COD_REG==paz_corrente]$data_studio_out), format="%Y")) - 
                       as.integer(format(unique(new[COD_REG==paz_corrente]$data_rif_ev), format="%Y")),
                     tail(data[COD_REG==paz_corrente & !is.na(hosp)]$comorbidity,n=1),unique(new[COD_REG==paz_corrente]$death),
                     unique(new[COD_REG==paz_corrente]$data_rif_ev),unique(new[COD_REG==paz_corrente]$data_studio_out),NA)
  names(temp) <- names(data)
  attributes(temp$timeEvent)<-attributes(data$timeEvent)
  data <- rbind(data,temp)
}
# sort data
data <- data[order(COD_REG),]

# imputation not available qt pharma
data<- data[!(is.na(hosp) & is.na(qt_pharma))]

## Define Timdep Adherence
data[!is.na(event),ADERENZA:=0]
codici <- unique(new$COD_REG)
for (i in 1:length ( codici )){
  paz_corrente = codici[i] 
  temp_hosp <- data[(COD_REG==paz_corrente) & (!is.na(event)),]
  temp_prescr <- data[(COD_REG==paz_corrente) & (is.na(event)),]
  
  for(k in 1:dim(temp_hosp)[1]){
    if(k==1){
      data[(COD_REG==paz_corrente) & (hosp==k)]$ADERENZA=0
    } else{
      #data_rif  =data[(COD_REG==paz_corrente) & (hosp==k)]$data_rif_ev
      data_start=data[(COD_REG==paz_corrente) & (hosp==(k-1))]$data_rif_ev
      data_stop =data[(COD_REG==paz_corrente) & (hosp==k)]$data_prest
      hosp_length= as.integer(data_stop-data_start)
      
      #pdc = (data[(COD_REG==paz_corrente) & (hosp==(k-1))]$ADERENZA)*(as.integer(data_start-data_rif))                                                       # Distinct days covered
      pdc=0
      vet_inizio = temp_prescr$data_prest    
      vet_fine   = temp_prescr$data_prest+temp_prescr$qt_pharma     
  
      for (j in 1:length(vet_inizio)){
        if(vet_inizio[j]<data_start)
           vet_inizio[j] = data_start
        if(vet_fine[j]>data_stop)
           vet_fine[j] = data_stop
        if(vet_fine[j]>vet_inizio[j]){
          pdc = pdc + (vet_fine[j]-vet_inizio[j])
          data_start = vet_fine[j]}
      }
      data[COD_REG==paz_corrente & hosp==k]$ADERENZA =  as.double(pdc/hosp_length)
    }
  }
}
data[,ADERENTE:=ifelse(ADERENZA>0.8,1,0)]
##########################################################################################
data<-data[!(is.na(hosp)),]
## Pass to gap times between events
data[,check:=as.integer(timeEvent)-as.integer(shift(timeEvent,n=-1)),by=COD_REG]
data<-data[!(event==1 & check==0)]
data[,GapEvent:=as.integer(timeEvent)-as.integer(shift(timeEvent)),by=COD_REG]
data<-data[(!is.na(GapEvent) & GapEvent!=0)]
# data for terminal events
dataDeath <- data[event==0]

save(data, file="dataRecTimeDep.RData")
save(dataDeath, file="dataDeathTimeDep.RData")

## remove unused structures
rm(new)
rm(temp)
gc()
