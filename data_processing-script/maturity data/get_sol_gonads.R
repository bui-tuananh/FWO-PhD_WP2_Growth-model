# input: getSmartFishData.R
# output: sol_gonads_2004_2023.rds

# ***
# SCRIPT TO SELECT SMARTFISH GONAD DATA
# ***

library(tidyverse)

# Read data ####
datdir <- "./data"
sampleFull <- read_rds(file.path(datdir,paste0("smartFishDataFull_2004_2023.rds")))

sol <- sampleFull[which(sampleFull$SAM.SpeciesFaoCode=="SOL" & !is.na(sampleFull$SPE.WeightGonads)),]
nrow(sol)
table(sol$TRI.Year, sol$HAU.IcesAreaFull)
table(sol$TRI.Year, sol$HAU.IcesAreaFull, sol$SPE.MaturityStage)
sol <- sol[, unlist(lapply(sol, function(x)!all(is.na(x))))] # remove columns with only NAs


# Save gonad data ####
write_rds(sol, file.path(datdir, paste0("sol_gonads_",min(sampleFull$TRI.Year),"_",max(sampleFull$TRI.Year),".rds")))


