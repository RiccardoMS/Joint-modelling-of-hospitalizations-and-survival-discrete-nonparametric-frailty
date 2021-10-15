# Dataset Splitting according to drug taken
load("SelectedData.RData")

# Chosen Drug
classe = "DIU"

# Patients
patients = unique(selection[classe_pharma == classe]$COD_REG)

# New Dataset
new = selection[(COD_REG %in% patients & tipo_prest==41)|(COD_REG %in% patients & classe_pharma == classe)]
unique(new$classe_pharma)

save(new, file="Split_Dataset/Diuretics.RData")
