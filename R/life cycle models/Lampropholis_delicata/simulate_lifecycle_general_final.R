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
micro <- readRDS('output/microclimates/micro_HI_2016_2021.Rda')
ecto <- readRDS('output/life cycle models/ecto_nichemapR_micro_HI_2016_2021.Rda')




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
clutchsize <- 5 # how many eggs per clutch? 
#clutch_ab <- c(0.21354, -4.39046) # regression results from shine paper - in raw_data folder "SVL_Clutchsize.csv" 
photostart <- 3 
photofinish <- 2 
batch <- 1
#clutchsize = 6
#minclutch = 3
#maxclutch = 8

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
longlat <- micro$longlat 


# retrieve outputs
environ<-as.data.frame(ecto$environ) # behaviour, Tb and environment
metout<-as.data.frame(ecto$metout) # behaviour, Tb and environment
enbal<-as.data.frame(ecto$enbal) # heat balance outputs
masbal<-as.data.frame(ecto$masbal) # mass balance outputs
debout<-as.data.frame(ecto$debout) # DEB model outputs
yearout <- as.data.frame(ecto$yearout) # whole life cycle summary
R_0 = yearout$R0
yearsout <- as.data.frame(ecto$yearsout) # annual summaries
environ <- cbind(environ,metout$SOLR) # add solar radiation for activity window plots
colnames(environ)[ncol(environ)] <- "Solar"
soil <- as.data.frame(ecto$soil)
# arranging dates
dates <- micro$dates
dates <- dates[simout$start:simout$finish]
debout$dates <- dates # add dates

##### adding wet mass with dessication
desiccation <- (1 - debout$PCT_DESIC / 100)
debout$Wet_mass_g <- debout$WETMASS*desiccation 
debout$dates <- as.POSIXct(debout$dates, format = "%Y-%m-%d %H:%M:%S")


#############################################
### Nest data and plots
#### reproduction events: grabbing nest date and following day for mass change
# then estimating clutch size basied off the change of E_B of mom and 
# then deviding E_B/E.0
# E_B = energy density at birth & E.0 = energy costs for one egg
reproevents <- which(debout$E_B[2:length(debout$E_B)] - 
                       debout$E_B[1:(length(debout$E_B) - 1)] <= -E.0)
repro_dates <- dates[reproevents]
repro_dates <- as.POSIXct(repro_dates, format = "%Y-%m-%d %H:%M:%S")
repro_dates_final <- sort(c(repro_dates, repro_dates + days(1)))



# Create a data frame with dates and their corresponding nest identifier
nest_data <- data.frame(date = repro_dates_final,
                        nest = rep(paste0('nest_', 
                                          seq_along(repro_dates)), 
                                   each=2))
# Filter your data based on the dates
clutch_filtered <- debout %>% 
  filter(dates %in% repro_dates_final)

# Join the nest_data with clutch_filtered based on the repro dates
clutch_data <- clutch_filtered %>% 
  left_join(nest_data, by = c("dates" = "date")) %>% 
  dplyr::select(nest, dates, E_B, Wet_mass_g) %>% 
  group_by(nest) %>%
  mutate(E_B_diff = E_B - lag(E_B),
         Wet_mass_diff_g = Wet_mass_g - lag(Wet_mass_g)) %>% 
  filter(!is.na(E_B_diff) & !is.na(Wet_mass_diff_g)) %>% 
  mutate(clutch_size = round(abs(E_B_diff/E.0), digits = 0),
         mean_egg_wetmass_g = abs(Wet_mass_diff_g)/clutch_size)

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
mean_temps_df <- data.frame(nest = character(), 
                            mean_temp_10cm = numeric())
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
  # Store the nest name and mean temperature in the summary dataframe
  mean_temps_df <- rbind(mean_temps_df, 
                         data.frame(nest = nest_name, 
                                    mean_temp_10cm = mean_temp))

  # Increment the nest_number for the next loop iteration
  nest_number <- nest_number + 1
  }


# Join clutch data with mean_temps_df with clutch_data_final
clutch_data_final <- clutch_data %>%
  left_join(mean_temps_df, by = "nest") %>% 
  ungroup() %>% 
  mutate(long = longlat[1],
         lat = longlat[2],
         R_0 = R_0) %>% 
  dplyr::select(long, lat, dates, nest, mean_temp_10cm, #ENV variables
          E_B, E_B_diff, Wet_mass_g, Wet_mass_diff_g, R_0, # Mom variables
          clutch_size, mean_egg_wetmass_g) # egg variables


}
# Simulation end
############################

View(clutch_data_final)



