#############################################################################################
############################## Data Preprocessing ###########################################
#############################################################################################
## [1] Cohort Selection
## [2] Additional variables definition

## clean workspace
rm( list = ls ())

## load libraries
library(data.table)

## load data
load("LombardyDataset.RData")

## arrange data
data<-LombardyDataset
dim(data)
names(data)
as.factor(data$COD_REG)

## cohort selection
# hospitalizations and pharmacological prescriptions
selection = data[data$tipo_prest==41 | data$tipo_prest==30]
as.factor(selection$tipo_prest)

# final date
selection = selection[selection$data_rif_ev <= "2011-12-31"]
as.factor(selection$COD_REG)
max(selection$data_rif_ev)

# at least one year survival
selection = selection[which(selection$data_studio_out-selection$data_rif_ev>=365)]
#selection = selection[which(selection$data_prest - selection$data_rif_ev<365)]

# ATC and labels for pharmacological prescriptions
selection$ATC = selection$class_prest
selection$ATC[which(selection$tipo_prest==41)] = NA
selection$classe_pharma = rep("NA", dim(selection)[1])
for (i in 1:dim(selection)[1]){
  if(grepl("C09A",selection[i]$ATC,fixed=TRUE)|grepl("C09B",selection[i]$ATC,fixed=TRUE)|grepl("C09X",selection[i]$ATC,fixed=TRUE))
   {selection[i]$classe_pharma="ACE"}
  else if(grepl("C09C",selection[i]$ATC,fixed=TRUE)|grepl("C09D",selection[i]$ATC,fixed=TRUE))
   {selection[i]$classe_pharma="ARB"}
  else if(grepl("C07",selection[i]$ATC,fixed=TRUE))
   {selection[i]$classe_pharma="BB"}
  else if(grepl("C03D",selection[i]$ATC,fixed=TRUE)|grepl("C03E",selection[i]$ATC,fixed=TRUE))
   {selection[i]$classe_pharma="AA"}
  else if(grepl("C03C",selection[i]$ATC,fixed=TRUE)|grepl("C03BA08",selection[i]$ATC,fixed=TRUE))
   {selection[i]$classe_pharma="DIU"}
  else
   {selection[i]$classe_pharma=NA}
}

# keep only selected drugs
selection = selection[which((!is.na(selection$classe_pharma)&selection$tipo_prest==30)|(selection$tipo_prest==41))]


# hospitalizations and prescriptions counter
selection = selection[,hosp:= sequence(.N), by=list(selection$COD_REG,selection$tipo_prest)][]
selection = selection[,pharma:= sequence(.N), by=list(selection$COD_REG,selection$tipo_prest)][]
for (i in 1:length(selection$pharma)){
  if (selection[i]$tipo_prest==41)
    selection[i]$pharma=NA
  else
    selection[i]$hosp=NA
}

# at least a hospitalization/pharmacological prescription
selection=selection[,.SD[(max(hosp,na.rm=TRUE)>=1&max(pharma,na.rm=TRUE)>=1)],by=COD_REG]
#save(selection,file="SelectedDataFFU.RData")

## Additional variables
# follow up time & censoring dummy
selection[,timeOUT:=data_studio_out - data_rif_ev]
selection[,death:= ifelse(desc_studio_out=="DECEDUTO",1,0)]

#Length of Stay and quantity of prescription
selection[tipo_prest==41,LOS:=qt_prest_Sum]
selection[tipo_prest==30,qt_pharma:=qt_prest_Sum]

# Date of Admission
selection[tipo_prest==41,dataADM:=data_prest-LOS]

#COMBO & DDD
levels(as.factor(selection$ATC))
ddd=c(5,40,15,75,400,NA,NA,NA,NA,15,160,20,160,150,75,10,5,600,37.5,NA,NA,NA,50,10,10,4,
      2.5,15,7.5,2.5,15,2,30,30,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,50,600,80,150,8,40,20,
      NA,NA,NA,NA,NA,NA,NA,NA,150)
for (i in 1:length(ddd)){
  if(!is.na(ddd[i]))
     selection[ATC==levels(as.factor(ATC))[i],DDD:=ddd[i]][ATC==levels(as.factor(ATC))[i],COMBO:=0]
  else
     selection[ATC==levels(as.factor(ATC))[i],DDD:=ddd[i]][ATC==levels(as.factor(ATC))[i],COMBO:=1]
}

# reorder columns
setcolorder(selection,
            c("COD_REG","SESSO",
              "data_rif_ev","data_studio_out", "desc_studio_out",
              "hosp","dataADM","LOS",
              "pharma","ATC","classe_pharma","qt_pharma","DDD","COMBO",
              "death","timeOUT"))
save(selection, file = "SelectedDataFFU.RData")
 




