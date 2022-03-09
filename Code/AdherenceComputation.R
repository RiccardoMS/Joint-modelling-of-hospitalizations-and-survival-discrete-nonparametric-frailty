#############################################################################################
############################# Adherence Computation #########################################
#############################################################################################
## [1] Adherence Computation: Compute Adherence related variables for chosen drug subdataset

## clear workspace
rm( list = ls ())

# load packages
library(data.table)

# load chosen drug subdataset
load("C:/Users/aughi/Desktop/TESI/Code/Split_Dataset/ACE_Inhibitors.RData")

# Compute Adherence
new = subset(new,!(tipo_prest==30 & is.na(qt_pharma)))
new[,ADERENZA:=0]
setcolorder(new,
            c("COD_REG","SESSO",
              "data_rif_ev","data_studio_out", "desc_studio_out",
              "hosp","dataADM","LOS",
              "pharma","ATC","classe_pharma","qt_pharma","DDD","COMBO",
              "death","timeOUT","ADERENZA"))

codici <- unique(new$COD_REG)
for (i in 1:length ( codici )){
  pdc = 0                                                            # Distinct days covered
  paz_corrente = codici[i]                                           # Current patient
  data_rif= new[COD_REG == paz_corrente , unique(data_rif_ev)]       # Reference date
  data_stop=(data_rif+365)                                           # End date
  
  vet_inizio = new[COD_REG==paz_corrente & tipo_prest==30, data_prest]     # Dates of prescription
  vet_fine   = new[COD_REG==paz_corrente & tipo_prest==30, data_prest+qt_pharma]     # End of prescription

  for (j in 1:length(vet_inizio)){
    if(vet_inizio[j]<data_rif)
      vet_inizio[j] = data_rif
    if(vet_fine[j]>data_stop)
      vet_fine[j] = data_stop
    
    if(vet_fine[j]>vet_inizio[j])
       pdc = pdc + (vet_fine[j]-vet_inizio[j])
       data_rif = vet_fine[j]
  }
  
  new[COD_REG==paz_corrente]$ADERENZA =  as.double(pdc/365)
  
}

# binary adherence
new[,ADERENTE:=0]
new[ADERENZA>=0.8]$ADERENTE=1

# adherence classes
new[,ADH_LEV:=0]
new[ADERENZA<0.25]$ADH_LEV=1
new[ADERENZA>=0.25 & ADERENZA<0.5]$ADH_LEV=2
new[ADERENZA>=0.5 & ADERENZA<0.75]$ADH_LEV=3
new[ADERENZA>=0.75]$ADH_LEV=4

# save
save(new, file="With_Adherence_Dataset/ACE_Inhibitors.RData")
