#############################################################################################
################################ Dataset Splitting ##########################################
#############################################################################################
## [1] Dataset Splitting: create subdatasets with patients undergoing a certain treatment

## clear workspace
rm( list = ls ())

## load libraries
library(data.table)

## load data
load("SelectedDataFFU.RData")

## chosen drug
classe = "ACE"

## patients under chosen treatment
patients = unique(selection[classe_pharma == classe]$COD_REG)

## split dataset
new = selection[(COD_REG %in% patients & tipo_prest==41)|(COD_REG %in% patients & classe_pharma == classe)]
unique(new$classe_pharma)

## save
save(new, file="Split_Dataset_FFU/ACE_Inhibitors.RData")
