# Clear
rm( list = ls ())
library(data.table)

# Dataset Splitting according to drug taken
load("SelectedDataFFU.RData")

# Chosen Drug
classe = "ACE"

# Patients
patients = unique(selection[classe_pharma == classe]$COD_REG)

# New Dataset
new = selection[(COD_REG %in% patients & tipo_prest==41)|(COD_REG %in% patients & classe_pharma == classe)]
unique(new$classe_pharma)

save(new, file="Split_Dataset_FFU/ACE_Inhibitors.RData")
