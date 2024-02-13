########################### Micro climate parameters
ystart <- 2016 # start date of simulation
yfinish <- 2020 # end date of simulation
scenario <- 0 # climate scenario
nyears <- yfinish - ystart + 4
windfac = 1 # wind factor
RUF = 0.004 #rufness height
REFL <-.15 #  reflectance of soild 
ERR <- 1.5 # model error tolerance
minshade <- 0 # minimum shade
maxshade <- 25 # maximum shade
Usrhyt <- 0.005 # lizard height above ground









########################### Ecotherm parameters
Ww_g <- 1.8 # in between W_p and W_i
shape <- 3 # lizard
CT_min <- 8.8 # critical thermal minimum (deg C) - Assumed
T_RB_min <- 12.6 # minimum emergence (retreat to bask) body temperature (deg C) observed in field - Matthews et al., 2023 lab minimum
T_B_min <-22.6 # minimum basking body temperature (deg C) observed in field - Matthews et al., 2023 field minimum
T_F_min <- 26.3 # minimum feeding/foraging body temperature (deg C) - 
T_pref <- 31 # preferred body temperature (deg C) - Zhang et al., 2023
T_F_max <-36.5 # maximum feeding/foraging temperature (deg C) - 
CT_max <-43 # critical thermal maximum (deg C) - Zhang et al., 2023

# morph, behav and water loss
alpha_max <- 0.80 # Mathews et al., 2016 
alpha_min <- .85 # Mathews et al., 2016 
diurn <- 1
nocturn <- 0
crepus <- 0
shade_seek <- 1
burrow <- 1
climb <- 0
shdburrow <- 2
mindepth <- 1
maxdepth <- 6
pct_wet <- 0.1 
pct_eyes <- 0
crepus <- 0
