library(NicheMapR)
library(R.matlab)
library(lubridate)
library(dplyr)
library(scales)


###### Regional function for egg data example. Here we will pass a dataframe with lat/longs
##### from each region and species into this function
source("R/Functions/simulations.R")
egg_results <- regional_simulation_egg(longitude = 151.004, #145.003
                                       latitude =  -33.8,#-37.8, 
                                       species = "Lampropholis_delicata",
                                       spatial = "Data/ERA_5/Australia/ERA5", 
                                       region = "AUS")
egg_results  


