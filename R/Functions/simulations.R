library(NicheMapR)
library(R.matlab)
library(lubridate)
library(dplyr)
library(scales)
######################### DEB in-silico experiment
######## This function is used for the DEB 'in silico experiments'
######## here you can alter the following as 'multipliers'
######## EG - structural costs, 
######## EO - energy content
######## z_mult - zoom scaling coefficent
######## you are able to manipulate these to capture the temperature size rule & 
######## to see how changing the energy content of the egg effects body mass of hatching
######## ** Note zoom factor only changes E.0 in "run.DEB.sim.R" file
####################################
deb_experiment <- function(micro = micro, species = species,
                           E.G_mult = 1, # structural costs multiplier
                           E.0_mult = 1,# energy content of the egg multiplier
                           z_mult = 1) # DEB scaling factor 
  {
  micro <<-micro
  E.G_mult <<- E.G_mult
  E.0_mult <<- E.0_mult
  z_mult <<- z_mult
  
  function_1 <- function(micro = micro, species = species) {
    source(paste0('R/life cycle models/', species, '/parameters_biophys.R'))
   
    
    
    #############  #############  #############
    ############# ECTOTHERM MODEL FOR ADULT ANIMAL
    #############  #############  #############
    # run ecotherm model 
    ecto <- ectotherm(Ww_g = Ww_g,
                      shape = shape,
                      CT_min = CT_min, 
                      T_RB_min = T_RB_min,
                      T_B_min = T_B_min,
                      T_F_min = T_F_min,
                      T_pref = T_pref,
                      T_F_max = T_F_max,
                      CT_max = CT_max,
                      alpha_max = alpha_max,
                      alpha_min = alpha_min,
                      diurn = diurn,
                      nocturn = nocturn,
                      shade_seek = shade_seek,
                      burrow = burrow,
                      climb = climb,
                      shdburrow = shdburrow,
                      mindepth = mindepth,
                      maxdepth = maxdepth,
                      pct_wet = pct_wet,
                      pct_eyes = pct_eyes)
    ecto <<- ecto
    return(list(micro = micro, ecto = ecto))
  }
  
  function_2 <- function(micro, ecto, E.G_mult = E.G_mult, 
                         E.0_mult = E.0_mult, z_mult = z_mult) {
    #############  #############  #############  #############
    ############# DEB MODEL FOR ADULT ANIMAL at gps location
    #############  #############  #############  #############
    # deb simulation function
    source(paste0('R/life cycle models/', species, '/run.DEB.sim.R'))
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
    
    
    ######## Matlab output
    # assign possible missing parameters
    source(paste0('R/life cycle models/', species, '/parameters_matlab_silico_DEB.R'))
    
    
    ######## energy mulipliers for function: 
    E.G <<- E.G*E.G_mult
    E.0 <<- E.0*E.0_mult
    z_mult <<- z_mult
    
    ############################
    # Simulation start
    startday <- 1
    ystart <- 2016
    for(ystart in 2016:2016){
      datelay <- as.POSIXct(paste0("1/12/", ystart), format = "%d/%m/%Y")
      simout <- run.DEB.sim(datelay, 1) # run simulation for how many years?
      
      # retrieve outputs
      ecto <- simout$ecto
      debout<-as.data.frame(ecto$debout) # DEB model outputs
      # arranging dates 
      dates <- micro_orig$dates
      dates <- dates[simout$start:simout$finish]
      debout$dates <- dates # add dates
      
      ##### adding wet mass with dessication
      desiccation <- (1 - debout$PCT_DESIC / 100)
      debout$Wet_mass_g <- debout$WETMASS*desiccation 
      debout$dates <- as.POSIXct(debout$dates, format = "%Y-%m-%d %H:%M:%S")
      day <- debout$DAY
      
      #### Getting out variables we need: incubation time and morphology
      incubation_time <- day[which(debout$E_H>E.Hb)[1]]
      debout_final <- debout %>% filter(DAY <= 80)
      return(debout_final)
    }
  }
  # Calling function_1
  results_1 <- function_1(micro = micro, species = species)
  micro <- results_1$micro
  ecto <- results_1$ecto
  E.G_mult <<- E.G_mult
  E.0_mult <<- E.0_mult
  z_mult <<- z_mult

  # Calling function_2 with the outputs of function_1
  results_2 <<- function_2(micro= micro, ecto = ecto,
                           E.G_mult = E.G_mult, 
                           E.0_mult = E.0_mult,
                           z_mult = z_mult)
  
  return(results_2)
  
}






##############################  MICRO --->  ECTOTHERM ---> LIFE HISTORY MODEL
########  This function is used for the simulation that is looped for each gps point
######## for each region (AUS, NZ, HA). Region will point at where the spatial data (ERA5)
######## is grabbed. The simulation runs in this order for each GPS location given: 
######## microclimate -> adult ecotherm model -> deb simulation
######## the simulation will return nest date, mean nest temperature, E_b, clutch size, and average
######## clutch mass. Read '#' that will tell you notes for each function
####################################
regional_simulation_egg <- function(longitude, latitude, 
                                    species, # delicata or guchi
                                    region, # AUS, HA, NZ
                                    spatial) {
  
  function_1 <- function(longitude, latitude, 
                         species, # delicata or guchi
                         region, # AUS, HA, NZ
                         spatial) {
    source(paste0('R/life cycle models/', species, '/parameters_biophys.R'))
    loc <- c(longitude, latitude) 
    nyears <- yfinish - ystart + 4
    
    # run microclimate model
    micro <- micro_era5(windfac = windfac, 
                        RUF = RUF, 
                        minshade = minshade, 
                        maxshade = maxshade, 
                        ERR = ERR, 
                        loc = loc, 
                        dstart = paste0('01/01/',ystart), 
                        dfinish = paste0('31/12/',yfinish), 
                        Usrhyt = Usrhyt, 
                        REFL = REFL, 
                        spatial = spatial, 
                        scenario = scenario, 
                        runshade = 1, 
                        cap = 1, 
                        soilgrids = 0)
    
    # Set timezone based on region
    if (region == "AUS") {
      if (sign(loc[1]) == -1) {
        gmtzone <- "+"
      } else {
        gmtzone <- ""
      }
      tz <- paste0("Etc/GMT", gmtzone, round(loc[1] / 15 * -1, 0))
    } else if (region == "HI") {
      tz <- "Pacific/Honolulu"
    } else if (region == "NZ") {
      # Add the timezone setting for New Zealand if different
      tz <- "Pacific/Auckland"
    } else {
      # Default timezone or throw an error
      stop("Region not recognized:", region)
    }
    attr(micro$dates, "tzone") <- tz
    micro <<-micro
    
    
    #############  #############  #############
    ############# ECTOTHERM MODEL FOR ADULT ANIMAL
    #############  #############  #############
    # run ecotherm model 
    ecto <- ectotherm(Ww_g = Ww_g,
                      shape = shape,
                      CT_min = CT_min, 
                      T_RB_min = T_RB_min,
                      T_B_min = T_B_min,
                      T_F_min = T_F_min,
                      T_pref = T_pref,
                      T_F_max = T_F_max,
                      CT_max = CT_max,
                      alpha_max = alpha_max,
                      alpha_min = alpha_min,
                      diurn = diurn,
                      nocturn = nocturn,
                      shade_seek = shade_seek,
                      burrow = burrow,
                      climb = climb,
                      shdburrow = shdburrow,
                      mindepth = mindepth,
                      maxdepth = maxdepth,
                      pct_wet = pct_wet,
                      pct_eyes = pct_eyes)
    ecto <<- ecto
    return(list(micro = micro, ecto = ecto))
  }
  
  function_2 <- function(micro, ecto) {
    #############  #############  #############  #############
    ############# DEB MODEL FOR ADULT ANIMAL at gps location
    #############  #############  #############  #############
    # deb simulation function
    source(paste0('R/life cycle models/', species, 
                  '/run.DEB.sim.R'))
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
    
    
    ######## Matlab output
    # assign possible missing parameters
    source(paste0('R/life cycle models/', species, '/parameters_matlab.R'))
    
    
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
                      clutch_size, mean_egg_wetmass_g)
      return(clutch_data_final)
    }
  }
  # Calling function_1
  results_1 <- function_1(longitude, latitude, species, region, spatial)
  micro <- results_1$micro
  ecto <- results_1$ecto
  
  # Calling function_2 with the outputs of function_1
  results_2 <<- function_2(micro= micro, ecto = ecto)
  
  return(results_2)
  
}





