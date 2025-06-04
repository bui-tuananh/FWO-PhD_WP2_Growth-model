# BEFORE: Access to SmartFish database
# AFTER: smartFishDataFull_2004_2023.rds

# ***
# SCRIPT TO DOWNLOAD SMARTFISH DATA
# [link](https://gitlab.ilvo.be/d1ict/d1/smartfish/-/blob/master/demoPackage.R)
# ***

# INSTALL Smartfish ####
# follows instruction on the link above
#library(remotes)

# Load libraries ####
library(RODBC)
library(dplyr)
# install.packages("Rcpp") # require to update the Rcpp package
library(Rcpp)
library(odbc)
library(DBI)
library(smartfish)
library(tidyverse)
# sessionInfo()

# Read data using the smartfish package ####
# set var
p <- 'NDGP'
pt <- ''
yFrom <- 1995  #2000
yTo <- 2024
#ds <- ("")
ds <- (list('trip','haul','segment','sample','sampleLengthFrequency','specimen','specimenAging'))

# get smartfish data
smartFishData <- getSmartFishData(p, pt, yFrom, yTo, ds)

# parse into dataframes
trip <- smartFishData[["trip"]]
tripGear <- smartFishData[["tripGear"]]
tripGearSelectiveDevice <- smartFishData[["tripGearSelectiveDevice"]]
tripGearSubType <- smartFishData[["tripGearSubType"]]
haul <- smartFishData[["haul"]]
haulRemark <- smartFishData[["haulRemark"]]
haulExcludeForProjectType <- smartFishData[["haulExcludeForProjectType"]]
haulEventSensorParameter <- smartFishData[["haulEventSensorParameter"]]
segment <- smartFishData[["segment"]]
segmentTripGear <- smartFishData[["segmentTripGear"]]
segmentLitter <- smartFishData[["segmentLitter"]]
segmentQuantityVariant <- smartFishData[["segmentQuantityVariant"]]
sample <- smartFishData[["sample"]]
sampleCountRange <- smartFishData[["sampleCountRange"]]
sampleQuantityVariant <- smartFishData[["sampleQuantityVariant"]]
sampleLengthFrequency <- smartFishData[["sampleLengthFrequency"]]
specimen <- smartFishData[["specimen"]]
specimenAging <- smartFishData[["specimenAging"]]
specimenSurvival <- smartFishData[["specimenSurvival"]]


# Join datasets 
sampleFull <- trip
sampleFull <- inner_join(sampleFull, haul, by="TripID")
sampleFull <- inner_join(sampleFull, segment, by="HaulID")  
sampleFull <- inner_join(sampleFull, sample, by="SegmentID")
sampleFull <- inner_join(sampleFull, specimen, by="SampleID") # "SPE.WeightGonads"
sampleFull <- inner_join(sampleFull, specimenAging, by="SpecimenID") # "SPE.WeightGonads"

# Save SmartFish data ####
datdir <- "./data"
write_rds(sampleFull, file = file.path(datdir,paste0("smartFishDataFull_",min(sampleFull$TRI.Year),"_",max(sampleFull$TRI.Year),".rds")))
