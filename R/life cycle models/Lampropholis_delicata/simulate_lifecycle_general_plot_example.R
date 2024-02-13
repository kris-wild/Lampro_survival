library(NicheMapR)
library(R.matlab)
library(lubridate)
library(dplyr)
library(scales)
genus <- "Lampropholis"
species <- "delicata"
species1 <- paste0(genus, "_", species)
species2 <- paste0(genus, " ", species)
species3 <- species2

source(paste0('R/life cycle models/', species1, 
              '/parameters_biophys.R'))
source(paste0('R/life cycle models/', species1, 
              '/run.DEB.sim.R'))

# micro and ecto mods
micro <- readRDS('output/microclimates/micro_Bris_2016_2021.Rda')
ecto <- readRDS('output/life cycle models/ecto_nichemapR_micro_Bris_2016_2021.Rda')


# region info - change 
#region ="Sydney, NSW (Long:151.10 Lat: -33.77)" # 
region = "Brisbane, QLD (Long:153.02 Lat: -27.469)"
#region ="Christchurch, NZ (Long:172.6306, Lat: - 43.5320)"
#region = "Auckland, NZ (Long:174.7645, Lat -36.8509)"
#region = "Honolulu, HI (Long -157.85 lat 21.3"


# micro file and dates
micro_orig <- micro
dates <- micro$dates
dates2 <- micro$dates2
soil <- as.data.frame(micro$soil)
metout <- as.data.frame(micro$metout)
humid <- as.data.frame(micro$humid)
rainfall <- micro$RAINFALL
soilpot <- as.data.frame(cbind(dates,as.data.frame(micro$soilpot)))
metout <- as.data.frame(cbind(dates,as.data.frame(micro$metout))) # above ground microclimatic conditions, min shade
soil <- as.data.frame(cbind(dates,as.data.frame(micro$soil))) # soil temperatures, minimum shade
soilmoist <- as.data.frame(cbind(dates,as.data.frame(micro$soilmoist))) # soil moisture, minimum shade
humid <- as.data.frame(cbind(dates,as.data.frame(micro$humid))) # soil humidity, minimum shade
environ <- as.data.frame(ecto$environ)


# calculating dew for rainfall
dew <- metout$DEW
dew[environ$TC < T_RB_min] <- 0
dew.agg <- aggregate(dew, by = list(format(dates, "%d/%m/%Y")), FUN = 'sum')
dew.agg$Group.1 <- as.Date(dew.agg$Group.1, format = "%d/%m/%Y")
dew.agg <- dew.agg[order(dew.agg$Group.1), ]
dew.agg <- dew.agg$x#[-1]
dew.agg[dew.agg > 0] <- 1
rainfall <- rainfall + dew.agg
micro$RAINFALL <- rainfall


if(exists("E.Hj")==TRUE){rm(E.Hj)}
if(exists("L.j")==TRUE){rm(L.j)}

pars <- readMat(paste0('R/DEB/', species1, '/results_', species1, '.mat'))
par.names <- unlist(labels(pars$par))
for(i in 1:length(par.names)){
  assign(par.names[i], unlist(pars$par[i]))
}
stages <- 4
if(exists("T.L")==FALSE){T_L <- T.L <- 173.15}else{T_L <- T.L}
if(exists("T.H")==FALSE){T_H <- T.H <- 373.15}else{T_H <- T.H}
if(exists("T.AL")==FALSE){T_AL <- T.AL <- 5E04}else{T_AL <- T.AL}
if(exists("T.AH")==FALSE){T_AH <- T.AH <- 9E04}else{T_AH <- T.AH}
if(exists("T.REF")==FALSE){T_REF <- T.REF <- 20 + 273.15}else{T_REF <- T.REF}


######## Grab these from MATLAB OUTPUT
E.0 <- 1145
L.b <- 0.4213

# assign possible missing parameters
if(exists("E.Hj")==FALSE){E.Hj <- E.Hb}
if(exists("E.He")==FALSE){E.He <- E.Hb}
if(exists("L.j")==FALSE){L.j <- L.b}
p.Am <- p.M / kap*z
E.m <- p.Am / v
F.m <- p.Am / kap.X # redefining F.m to max possible value
z.mult <- 1         # DEB body size scaling parameter

# morph, behav and water loss
soilnode <- 7
rundeb <- 1 # run the DEB model?
reset <- 0
metab_mode <- 0

V_i <- z ^ 3
Ww_g <- (E.m * V_i) / mu.E * 23.9 / d.E + V_i
L_i <- z / del.M * 10

# energy content
#E_insect_g <- 23850 # J/g, insects, Buckley 2008 citing Reichle 1971; Andrews and Asato 1977
E_insect_g <- 25600 # J/g, insects, Waldschmidt 1986 cited by Angilletta 2001
V_prey <- 0.082 # cm3 from Pianka 1969, mean of largest 10 items
V_prey <- 0.2 # cm3 from Pianka 1969, largest item

# energy per insect
f_dry.insect <- 0.3 # g dry / g wet
rho.insect <- 1 # g / cm3
W_insect_dry <- V_prey * rho.insect * f_dry.insect # cm3 to g to dry
E_insect <- W_insect_dry * E_insect_g

## Stomach volume
V_stomach <- V_prey # largest prey item
V_stomach <- 0.33E-5 * L_i ^ 2.64 # Avery 1973 Lacerta vivipara 
V_stomach <- 0.048 * Ww_g  ^ 0.99 # Avery 1973 Lacerta vivipara
E_stomach <- E_insect_g * 0.3 * V_stomach
E_sm <- E_stomach / V_i # J / cm3 structure
E_sm <- 2000

d_Egg <- 0.5 # Dry mass fraction of egg (0-1)
pct_H_X <- 70 # Water content of food (%)
pct_H_P <- 30 # Water in faeces (product) (%)
pct_H_R <- 10 # Minimum tolerated dehydration (% of wet mass)
aestivate <- 0
aestdepth <- 5
depress <- 0.7
raindrink <- .1 # mm

K <- 0.05
X_max <- 3
pct_H_X <- 70
X <- 100000
F_m_ref <- p.Am/kap.X * 10000

# egg
V_init <- 3e-9
E_init <- E.0 / V_init
E_H_init <- 0
stage <- 0


# reproduction parameters
viviparous <- 0 # live bearing (1) or egg laying (0)
clutchsize <- 4 # how many eggs per clutch? 
#clutch_ab <- c(0.21354, -4.39046) # regression results from shine paper - in raw_data folder "SVL_Clutchsize.csv" 
photostart <- 3 
photofinish <- 2 
batch <- 1

# overwrite nitrogenous waste indices with those of uric acid (currently ammonia by default)
n.NC <- 1
n.NH <- 4/5
n.NO <- 3/5
n.NN <- 4/5




############################
# Simulation start
startday <- 1
ystart <- 2016
for(ystart in 2016:2016){
datelay <- as.POSIXct(paste0("1/12/", ystart), format = "%d/%m/%Y")
simout <- run.DEB.sim(datelay, 4)
ecto <- simout$ecto
longlat <-micro$longlat 


# retrieve outputs
environ<-as.data.frame(ecto$environ) # behaviour, Tb and environment
metout<-as.data.frame(ecto$metout) # behaviour, Tb and environment
enbal<-as.data.frame(ecto$enbal) # heat balance outputs
masbal<-as.data.frame(ecto$masbal) # mass balance outputs
debout<-as.data.frame(ecto$debout) # DEB model outputs
yearout <- as.data.frame(ecto$yearout) # whole life cycle summary
yearsout <- as.data.frame(ecto$yearsout) # annual summaries
environ <- cbind(environ,metout$SOLR) # add solar radiation for activity window plots
colnames(environ)[ncol(environ)] <- "Solar"
soil <- as.data.frame(ecto$soil)

dates <- micro$dates
dates <- dates[simout$start:simout$finish]
desiccation <- (1 - debout$PCT_DESIC / 100)


###################################
#####  setting up plot 4 panel: mom repo output and thermal biology
par(mfrow = c(4,1))
par(oma = c(2,2,1,1) + 0.1) # margin spacing stuff
par(mar = c(2,4,1,1) + 0.1) # margin spacing stuff

####  Plot 1: DEB PLOT 
plot(dates, debout$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     ylim = c(0, 6))
points(dates, debout$V, type = 'l', col = 'dark green', lwd = 2)
points(dates, debout$WETMASS-debout$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(dates, debout$WETMASS-debout$WETGONAD-debout$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = dates[which(debout$E_H>E.Hb)[1]], lty = 2, col = 'grey')
abline(v = dates[which(debout$E_H>E.Hp)[1]], lty = 2, col = 'grey')
legend("topleft",
        c('repro. buffer', 'food in gut', 'reserve', 'structure'), lty = rep(1, 4),
        col = c("pink", "brown", "grey", "dark green"), bty = 'n')
text(0, max(debout$WETMASS) * 1, labels = "embryo", cex = 0.85)
text((which(debout$E_H > E.Hp)[1] - which(debout$E_H > E.Hp)[1] * .5) / 24 ,
     max(debout$WETMASS) * 1, labels = "immature", cex = 0.85)
text(which(debout$E_H > E.Hp)[1] * 1.2 / 24, max(debout$WETMASS) * 1,
     labels = "adult", cex = 0.85)
title(main = paste0(species2," ", region), line = -1)

####  Plot2: wet mass and reproduction events (observed vs predicted)
plot(dates, debout$WETMASS * desiccation, type = 'l', xlab = 'date',
     ylab = paste0('wet mass, g'), col = 'black', lwd = 1,
     ylim = c(0, 6))
abline(v = (dates)[which(debout$E_H>E.Hb)[1]], lty = 2, col = 'grey')
abline(v = (dates)[which(debout$E_H>E.Hp)[1]], lty = 2, col = 'grey')

# predicted reproduction points
reproevents <- which(debout$E_B[2:length(debout$E_B)] - debout$E_B[1:(length(debout$E_B) - 1)] <= -E.0)
points(dates[reproevents], (debout$WETMASS * desiccation)[reproevents], pch = 17, col = 'blue', cex = 1.2)
text((dates)[which(debout$E_H>E.Hb)[1]], 3, 'hatch', cex = 0.7, srt = 90)
text((dates)[which(debout$E_H>E.Hp)[1]], 3, 'mature', cex = 0.7, srt = 90)


#### plot 3 - SVL and predicted reproduction events
plot(dates, debout$L_W, type = 'l', xlab = 'date',
     ylab = paste0('SVL, mm'), col = 'black', lwd = 1,
     ylim = c(0, 65))
abline(v = (dates)[which(debout$E_H>E.Hb)[1]], lty = 2, col = 'grey')
abline(v = (dates)[which(debout$E_H>E.Hp)[1]], lty = 2, col = 'grey')
text((dates)[which(debout$E_H>E.Hb)[1]], 30, 'hatch', cex = 0.7, srt = 90)
text((dates)[which(debout$E_H>E.Hp)[1]], 30, 'mature', cex = 0.7, srt = 90)
reproevents <- which(debout$E_B[2:length(debout$E_B)] - debout$E_B[1:(length(debout$E_B) - 1)] <= -E.0)
points(dates[reproevents], (debout$L_W)[reproevents], pch = 17, col = 'blue', cex = 1.2)

####  Thermal Plot
transparent_red <- rgb(1, 0, 0, alpha = 0.5)
plot(dates, environ$TSUB, type = 'l', ylab = "Temperature (C)")
points(dates, environ$TC, type = 'l', col = transparent_red)
legend("topleft",
       c('Tsub', 'Tb'), lty = rep(1, 2),
       col = c("black", "red"), bty = 'n')


#############################################
### Nest data and plots
#### reproduction events: grabbing nest date and following day for mass change
# then estimating clutch size basied off the change of E_B of mom and 
# then deviding E_B/E.0
# E_B = energy density at birth & E.0 = energy costs for one egg
repro_dates <- dates[reproevents]
repro_dates <- as.POSIXct(repro_dates, format = "%Y-%m-%d %H:%M:%S")
repro_dates_final <- sort(c(repro_dates, repro_dates + days(1)))
debout$dates <- dates # add dates


# Create a data frame with dates and their corresponding nest identifier
nest_data <- data.frame(date = repro_dates_final,
                        nest = rep(paste0('nest_', 
                                          seq_along(repro_dates)), 
                                   each=2))

# Ensure debout and nest_data have dates in the same format
debout$dates <- as.POSIXct(debout$dates, format = "%Y-%m-%d %H:%M:%S")

# Filter your data based on the dates
clutch_filtered <- debout %>% 
  filter(dates %in% repro_dates_final)

# Join the nest_data with clutch_filtered based on the repro dates
clutch_data <- clutch_filtered %>% 
  left_join(nest_data, by = c("dates" = "date")) %>% 
  dplyr::select(nest, dates, E_B, WETMASS) %>% 
  group_by(nest) %>%
  mutate(E_B_diff = E_B - lag(E_B),
         WETMASS_diff_g = WETMASS - lag(WETMASS)) %>% 
  filter(!is.na(E_B_diff) & !is.na(WETMASS_diff_g)) %>% 
  mutate(clutch_size = round(abs(E_B_diff/E.0), digits = 0),
         mean_egg_wetmass_g = abs(WETMASS_diff_g)/clutch_size)

#### nest temperature plots
# extract soil data
soil$dates <- dates 
soil <- soil %>% dplyr::select(dates, D10cm, D15cm)
soil$dates <- as.POSIXct(soil$dates, format = "%Y-%m-%d %H:%M:%S")

# Initialize an empty data frame to hold the filtered results
nest_data <- data.frame(dates = as.POSIXct(character()),
                            D10cm = numeric(), 
                            D15cm = numeric(),
                            nest = character())

# Loop through each date in repo_dates and keep track of the nest number
nest_number <- 1
for (date in repro_dates) {
  # Calculate the date 35 days later
  end_date <- date + 35*24*60*60  # 35 days in seconds
  # Filter soil for dates between date and end_date
  temp_data <- soil[soil$dates >= date & soil$dates <= end_date, ]
  # Create the nest name for this group with month and year
  nest_name <- paste("nest", nest_number, sep = "_")
  # Add the nest name as a new column to temp_data
  temp_data$nest <- nest_name
  # Combine the filtered data
  nest_data <- rbind(nest_data, temp_data)
  # Calculate the mean temperature for the current nest
  mean_temp <- mean(temp_data$D10cm)
  # Base R plotting
  plot(temp_data$dates, temp_data$D10cm, type = 'l', 
       col = 'red', 
       main = paste("Nest temp", nest_name, region), 
       xlab = "Date", ylab = "Temperature (C) at 10cm", xaxt = "n")
  axis.POSIXct(1, temp_data$dates)  # Add POSIXct dates to the x-axis
  # Add a dashed line at the mean temperature
  abline(h = mean_temp, col = "darkgrey", lty = "dashed")  # red dashed line for the mean
  # Add the mean value as text annotation in the top left corner of the plot
  margin_offset <- par("mar")  # Get the current margin settings
  text(x = min(temp_data$dates), 
       y = max(temp_data$D10cm, na.rm = TRUE) - 
         margin_offset[2] * 0.05, 
       labels = paste("Mean Nest Temperature:", 
                      round(mean_temp, 2)),
       pos = 4, cex = 0.8, col = "Black", font = 2)
  # Increment the nest_number for the next loop iteration
  nest_number <- nest_number + 1}



# deb output
debout$dates <- dates # add dates

#### reproduction events: grabbing nest date and following day for mass change
# then estimating clutch size basied off the change of E_B of mom and 
# then deviding E_B/E.0
# E_B = energy density at birth & E.0 = energy costs for one egg
repro_dates <- dates[reproevents]
repro_dates <- as.POSIXct(repro_dates, format = "%Y-%m-%d %H:%M:%S")
repro_dates_final <- sort(c(repro_dates, repro_dates + days(1)))

# Create a data frame with dates and their corresponding nest identifier
nest_data <- data.frame(date = repro_dates_final,
                        nest = rep(paste0('nest_', 
                                          seq_along(repro_dates)), 
                                   each=2))

# Ensure debout and nest_data have dates in the same format
debout$dates <- as.POSIXct(debout$dates, format = "%Y-%m-%d %H:%M:%S")

# Filter your data based on the dates
clutch_filtered <- debout %>% 
  filter(dates %in% repro_dates_final)

# Join the nest_data with clutch_filtered based on the repro dates
clutch_data <- clutch_filtered %>% 
  left_join(nest_data, by = c("dates" = "date")) %>% 
  dplyr::select(nest, dates, E_B, WETMASS) %>% 
  group_by(nest) %>%
  mutate(E_B_diff = E_B - lag(E_B),
         WETMASS_diff_g = WETMASS - lag(WETMASS)) %>% 
  filter(!is.na(E_B_diff) & !is.na(WETMASS_diff_g)) %>% 
  mutate(clutch_size = round(abs(E_B_diff/E.0), digits = 0),
         mean_egg_wetmass_g = abs(WETMASS_diff_g)/clutch_size)

}
# Simulation end
############################









##############
# environ - make df with dates to work with Tb, activity and shade
environ['dates'] <- dates
environ$year <- lubridate::year(environ$dates)
environ$month <- as.factor(lubridate::month(environ$dates))
environ$months <- as.factor(lubridate::month(environ$dates, label = TRUE))
environ$date <- lubridate::date(environ$dates)
environ$hr <- as.factor(lubridate::hour(environ$dates))
environ <- environ %>% arrange(dates)  

# filter environ df by year
environ_2018 <- environ %>% filter(year == "2018")
environ_2019 <- environ %>% filter(year == "2019")
environ_2020 <- environ %>% filter(year == "2020")
environ_2021 <- environ %>% filter(year == "2021")

##############
# metout - make behaviour df for activity plots
metout$dates <- dates  
metout$year <- lubridate::year(metout$dates)
metout$month <- as.factor(lubridate::month(metout$dates))
metout$months <- as.factor(lubridate::month(metout$dates, label = TRUE))
metout$date <- lubridate::date(metout$dates)
metout$hr <- as.factor(lubridate::hour(metout$dates))
metout <- metout %>% arrange(dates)  
# filter metout df by year
metout2018 <- metout %>% filter(year == "2018")
metout2019 <- metout %>% filter(year == "2019")
metout2020 <- metout %>% filter(year == "2020")
metout2021 <- metout %>% filter(year == "2021")

##############
# Mass balance df
masbal <- cbind(dates, masbal)
masbal$year <-  lubridate::year(dates)
masbal$month <- as.factor(lubridate::month(masbal$dates))
masbal$months <- as.factor(lubridate::month(masbal$dates, label = TRUE))
masbal$date <- lubridate::date(masbal$dates)
masbal$hr <- as.factor(lubridate::hour(masbal$dates))
masbal <- masbal %>% arrange(dates)  
# filter mass balance df by year
masbal_2018 <- masbal %>% filter(year == "2018")
masbal_2019 <- masbal %>% filter(year == "2019")
masbal_2020 <- masbal %>% filter(year == "2020")
masbal_2021 <- masbal %>% filter(year == "2021")

#############
# TB, Activity, Shade, Depth for 2017, 2018, 2019
# setting plots
par(mfrow = c(1,1))

# 2018
with(environ_2018, plot(TC ~ dates, ylab = "", xlab="2018", col = 'black', ylim = c(-50, 40), type = "l", yaxt = 'n', main = paste0('2018 Tb, Activity, Shade, Depth')))
with(environ_2018, points(ACT * 2 + 1 ~ dates, type = "p", pch = 16, col = "orange"))
with(environ_2018, points(SHADE / 9 -16 ~ dates, type = "l", col = "dark green"))
with(environ_2018, points(DEP/10 - 40 ~ dates, type = "l", col = "brown"))
abline(ecto$T_F_min, 0, lty = 2, col = 'blue')
abline(T_pref, 0, lty = 2, col = 'orange')
abline(ecto$T_F_max, 0, lty = 2, col = 'red')
ytick<-seq(10, 40, by=5)
axis(side=2, at=ytick, labels = TRUE)
mtext(text = c('A', 'B', 'I'), side = 2, line = 1, at = c(5, 3, 1))
ytick<-seq(-16, -6, by=2)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = seq(0, 100, 20), side = 2, line = 1, at = seq(-16, -6, 2), las = 2)
ytick<-seq(-45, -25, by=5)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = c(-50, 0, 50, 100, 150), side = 2, line = 1, at = c(-45, -40, -35, -30, -25), las = 2)
abline(h = -40, lty = 2, col = 'grey')
mtext(text = c('body temperature (°C)', 'activity', 'shade (%)', 'depth (cm)'), side = 2, line = 2.5, at = c(30, 5, -7, -35))
text(0.1, c(ecto$T_F_max + 1, ecto$T_F_min + 1), c('T_F_max', 'T_F_min'), col = c('red', 'blue'), cex = 0.75)

# 2018 - Annual activity window with shade options ranging from 0% to 90%
forage <- subset(environ_2018, ACT == 2) # get foraging hours
bask <- subset(environ_2018, ACT == 1) # get basking hours
night <- subset(environ_2018, Solar == 0) # get night hours
with(night, plot(as.numeric(hr) ~ date, ylab = "Hour of Day", xlab = "Day of Year", pch = 15, cex = 2,col = 'dark blue', main = paste0('2018 - Annual activity window with shade options ranging from 0% to 90%'))) # nighttime hours
with(bask, points(as.numeric(hr) ~ date, pch = 15, cex = 2, col = 'light blue')) # basking Tbs
with(forage, points(as.numeric(hr) ~ date, pch = 15, cex = 2, col = 'orange')) # foraging Tbs




# 2019
with(environ_2019, plot(TC ~ dates, ylab = "", xlab="2019", col = 'black', ylim = c(-50, 40), type = "l", yaxt = 'n', main = paste0('2019 Tb, Activity, Shade, Depth')))
with(environ_2019, points(ACT * 2 + 1 ~ dates, type = "p", pch = 16, col = "orange"))
with(environ_2019, points(SHADE / 9 -16 ~ dates, type = "l", col = "dark green"))
with(environ_2019, points(DEP/10 - 40 ~ dates, type = "l", col = "brown"))
abline(ecto$T_F_min, 0, lty = 2, col = 'blue')
abline(T_pref, 0, lty = 2, col = 'orange')
abline(ecto$T_F_max, 0, lty = 2, col = 'red')
ytick<-seq(10, 40, by=5)
axis(side=2, at=ytick, labels = TRUE)
mtext(text = c('A', 'B', 'I'), side = 2, line = 1, at = c(5, 3, 1))
ytick<-seq(-16, -6, by=2)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = seq(0, 100, 20), side = 2, line = 1, at = seq(-16, -6, 2), las = 2)
ytick<-seq(-45, -25, by=5)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = c(-50, 0, 50, 100, 150), side = 2, line = 1, at = c(-45, -40, -35, -30, -25), las = 2)
abline(h = -40, lty = 2, col = 'grey')
mtext(text = c('body temperature (°C)', 'activity', 'shade (%)', 'depth (cm)'), side = 2, line = 2.5, at = c(30, 5, -7, -35))
text(0.1, c(ecto$T_F_max + 1, ecto$T_F_min + 1), c('T_F_max', 'T_F_min'), col = c('red', 'blue'), cex = 0.75)

# 2019 - Annual activity window with shade options ranging from 0% to 90%
forage <- subset(environ_2019, ACT == 2) # get foraging hours
bask <- subset(environ_2019, ACT == 1) # get basking hours
night <- subset(environ_2019, Solar == 0) # get night hours
with(night, plot(as.numeric(hr) ~ date, ylab = "Hour of Day", xlab = "Day of Year", pch = 15, cex = 2,col = 'dark blue', main = paste0('2019 - Annual activity window with shade options ranging from 0% to 90%'))) # nighttime hours
with(bask, points(as.numeric(hr) ~ date, pch = 15, cex = 2, col = 'light blue')) # basking Tbs
with(forage, points(as.numeric(hr) ~ date, pch = 15, cex = 2, col = 'orange')) # foraging Tbs


# 2020
with(environ_2020, plot(TC ~ dates, ylab = "", xlab="2019", col = 'black', ylim = c(-50, 40), type = "l", yaxt = 'n', main = paste0('2020 Tb, Activity, Shade, Depth')))
with(environ_2020, points(ACT * 2 + 1 ~ dates, type = "p", pch = 16, col = "orange"))
with(environ_2020, points(SHADE / 9 -16 ~ dates, type = "l", col = "dark green"))
with(environ_2020, points(DEP/10 - 40 ~ dates, type = "l", col = "brown"))
abline(ecto$T_F_min, 0, lty = 2, col = 'blue')
abline(T_pref, 0, lty = 2, col = 'orange')
abline(ecto$T_F_max, 0, lty = 2, col = 'red')
ytick<-seq(10, 40, by=5)
axis(side=2, at=ytick, labels = TRUE)
mtext(text = c('A', 'B', 'I'), side = 2, line = 1, at = c(5, 3, 1))
ytick<-seq(-16, -6, by=2)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = seq(0, 100, 20), side = 2, line = 1, at = seq(-16, -6, 2), las = 2)
ytick<-seq(-45, -25, by=5)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = c(-50, 0, 50, 100, 150), side = 2, line = 1, at = c(-45, -40, -35, -30, -25), las = 2)
abline(h = -40, lty = 2, col = 'grey')
mtext(text = c('body temperature (°C)', 'activity', 'shade (%)', 'depth (cm)'), side = 2, line = 2.5, at = c(30, 5, -7, -35))
text(0.1, c(ecto$T_F_max + 1, ecto$T_F_min + 1), c('T_F_max', 'T_F_min'), col = c('red', 'blue'), cex = 0.75)

# 2020 - Annual activity window with shade options ranging from 0% to 90%
forage <- subset(environ_2020, ACT == 2) # get foraging hours
bask <- subset(environ_2020, ACT == 1) # get basking hours
night <- subset(environ_2020, Solar == 0) # get night hours
with(night, plot(as.numeric(hr) ~ date, ylab = "Hour of Day", xlab = "Day of Year", pch = 15, cex = 2,col = 'dark blue', main = paste0('2020 - Annual activity window with shade options ranging from 0% to 90%'))) # nighttime hours
with(bask, points(as.numeric(hr) ~ date, pch = 15, cex = 2, col = 'light blue')) # basking Tbs
with(forage, points(as.numeric(hr) ~ date, pch = 15, cex = 2, col = 'orange')) # foraging Tbs


