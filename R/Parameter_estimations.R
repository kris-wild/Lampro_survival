pacman::p_load("dplyr","lme4", "lubridate",  
               "emmeans", "ggExtra", "cowplot", 
               "patchwork", "ggplot2", "stringr")


#######################################################
# Summary for matlab zero-variate data across life stages
######################################################

##### ##### ##### 
##### Load in three datasets
## 1) survival data - empirical dataset from our eggs for 
# at birth estimations
survival_data <- read.csv("Data/Lizard_Db_final_survival_est.csv", 
                 head = T) %>% 
  mutate(temp = as.character(temp))

# separate the two species into df
guich_df_survival <-survival_data %>% 
  filter(species == "guichenoti") 
deli_df_survival <- survival_data %>% 
  filter(species == "delicata")


##### ##### ##### 
##### 2) data base life history data - used for DEB parameter estimates 
# across life history (at puberty & ultimate)
life_history_parameter_data <- read.csv("Data/Dalton_growth.csv", 
                                        head = T) %>% 
  filter(sex == "f") # only select females for further measurements




###############################################
##### Analysis: at birth (Lb, Wwb), a puberty (Lp, Wwp), 
#####           maximum (Li, Wwi, Ri)


##### ##### ##### 
##### At birth (ab) estimates with hatchling/survival data from our study
# 1) look at developmental time & morph differences by species and temp
# early puberty is estimated within 6 months 
data_ab <- survival_data %>% 
  group_by(species, temp) %>% 
  summarise(dev_time_mean = mean(dev_time, na.rm = TRUE),
            Lb_SVL = mean(hatch_svl, na.rm = TRUE)*.1, 
            Wwb_mass = mean(hatch_mass, na.rm = TRUE))
data_ab


##########
### 2) Puberty measurements - L. delicata
# filter values at early puberty from SVL from Forsman & Shine 1995:
# https://academic.oup.com/biolinnean/article-abstract/55/4/273/2706115
deli_data_ap_filt <- life_history_parameter_data %>% 
  filter(species == "delicata" ) %>% 
  filter(svl > 31 & svl <33) # Forsman & Shine Table 1
# final dataframe
deli_data_ap <- deli_data_ap_filt %>% 
  summarise(Lp_SVL = round(mean(svl, na.rm = TRUE), 
                             digits = 2)*.1, 
            Wwp_mass = round(mean(mass, na.rm = TRUE), digits = 2))

deli_data_ap


##########
### 3) maximum morph measurements (Li, Wwi) - L. delicata
deli_data_max <- life_history_parameter_data %>% 
  summarise(Li_SVL = max(svl, na.rm = TRUE), 
            Wwi_mass = max(mass, na.rm = TRUE))
deli_data_max



#######################################################
# save data for MATLAB raw data for estimations
#######################################################

###  LW data: length (SVL cm) and wet weight data (g)
LW_data <- read.csv("Data/Dalton_growth.csv", 
                head = T) %>% 
  filter(sex == "f") %>% # only select females for further measurements
  filter(species == "delicata") %>% 
  filter(!is.na(svl) & !is.na(tl) & !is.na(mass)) 
LW_data_final <- LW_data %>% 
  mutate(SVL = round((svl)*.1, digits = 1), 
         mass_rd = round(LW_data$mass, digits = 1)) %>% # change to CM
  dplyr::select(SVL, mass_rd, mass) %>% 
  #group_by(mass_rd) %>% 
  #rename(mass = mass_rd) %>% 
  #summarise(total_length = mean(total_length))%>% # getting mean mass and TL
  dplyr::select(SVL, mass)
plot(LW_data_final$mass, LW_data_final$SVL)
write.table(LW_data_final, "Matlab/Raw_data/LW_data.csv", 
            sep = ",", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)


###  t-L data: length (SVL cm) and time (d)
tL_data <- read.csv("Data/Dalton_growth.csv", 
                    head = T) %>% 
  filter(sex == "f") %>% # only select females for further measurements
  mutate(SVL = (svl)*.1) %>% # change to CM
  filter(species == "delicata") %>% 
  filter(!is.na(SVL) & !is.na(age)) %>% 
  #group_by(age) %>% 
  #summarise(total_length = mean(total_length)) %>% 
  ungroup() %>% 
  dplyr::select(age, SVL) 
plot(tL_data$age, tL_data$SVL)
write.table(tL_data, "Matlab/Raw_data/tL_data.csv", 
            sep = ",", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)


###  t-W data: wet wet weight (g) and time (d)
tW_data <- read.csv("Data/Dalton_growth.csv", 
                    head = T) %>% 
  filter(sex == "f") %>% # only select females for further measurements
  filter(species == "delicata") %>% 
  filter(!is.na(mass) & !is.na(age)) %>% 
  #group_by(age) %>% 
  #summarise(mass = mean(mass)) %>% 
  ungroup() %>% 
  dplyr::select(age, mass) 
plot(tW_data$age, tW_data$mass)
write.table(tW_data, "Matlab/Raw_data/tW_data.csv", 
            sep = ",", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)


### CO2 data Kar et al., 2023_b for raw data used in matlab
kar_CO2_dat <- read.csv(file = "Data/Kar_etal_2023_b_Long_Tinc_MR.csv") %>%
  mutate(mass = exp(lnmass),
         V_CO2_min_g = exp(lnmr),
         temp = as.character(temp))
kar_O2 <- kar_CO2_dat %>% 
  filter(temp == 28) %>% 
  group_by(id) %>% 
  summarise(mean_mass = mean(mass), 
            mean_VCO2_min_g = mean(V_CO2_min_g)) %>% 
  mutate(mean_VO2_min_g = mean_VCO2_min_g/0.77, # RQ value from paper
         mlO2_min = mean_VO2_min_g*mean_mass) %>% 
  filter(!is.na(mean_mass) & !is.na(mlO2_min)) %>% 
  ungroup() %>% 
  dplyr::select(mean_mass, mlO2_min)
kar_O2
plot(kar_O2$mean_mass, kar_O2$mlO2_min)
mod <- lm(mlO2_min~mean_mass, kar_O2)
summary(mod)
abline(mod)
write.table(tL_data, "Matlab/Raw_data/WJO_data.csv", 
            sep = ",", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)




############################
# Regression of Mass and clutch size from:
#Forsman, A. and Shine, R., 1995. The adaptive significance of colour pattern #polymorphism in the Australian scincid lizard Lampropholis delicata. Biological #Journal of the Linnean Society, 55(4), pp.273-291.
# Figure 3 A|B
############################
svl_clutch <- read.csv("Data/SVL_vs_cltuchsize.csv") 
plot(svl_clutch$SVL_mm, svl_clutch$Clutch_size)
lm<- lm(Clutch_size~SVL_mm, data = svl_clutch)
abline(lm)
summary(lm)




#######################################################
# save metabolsim data (CO2) for CTE calculation in R
#######################################################
# summarise by ID and temp
kar_CO2_O2_dat_summary_id <- kar_CO2_dat %>% 
  group_by(id, temp) %>% 
  summarise(mean_mass = mean(mass), 
          mean_VCO2_min_g = mean(V_CO2_min_g)) %>% 
  mutate(mean_VO2_min_g = mean_VCO2_min_g/0.77, # RQ value from paper
         mlO2_min = mean_VO2_min_g*mean_mass) %>% 
  ungroup()
# Summarize by mean VCO2_min_g by temp for CTE calculation
kar_CO2_dat_CTE <- kar_CO2_O2_dat_summary_id %>% 
  group_by(temp) %>% 
  summarise(mean_VCO2_min_g = mean(mean_VCO2_min_g, na.rm = TRUE),
            mean_mass = mean(mean_mass, na.rm = TRUE))
write.csv(kar_CO2_dat_CTE, file = "Data/CTE_VCO2_temp.csv")



#######################################################
# Reflectance and Body temeprature estimates from Mathews et al., 2023
# Matthews, G., Farquhar, J. E., White, C. R., & Chapple, D. G. (2023). Does thermal biology differ between two colour pattern morphs of a widespread Australian lizard?. Journal of Thermal Biology, 114, 103579.
#######################################################
# ABS  data
ABS_data <- read.csv("Data/Mathews_etal_2023_Refl.csv") 
ABS_summary <- ABS_data %>% 
  group_by(MORPH) %>% 
  summarise(min = 100-min(DORS, na.rm = T),
            mean = 100-mean(DORS, na.rm = T),
            max = 100-max(DORS, na.rm = T))
ABS_summary




# Tb gradient data
T_set_data <- read.csv("Data/Mathews_etal_2023_thermal_gradient.csv")
min(T_set_data$TMIN, na.rm = TRUE)
mean(T_set_data$WMEAN, na.rm = TRUE)
max(T_set_data$TMAX, na.rm = TRUE)
par(mfrow = c(3, 1))
hist(T_set_data$TMIN, xlim = c(10, 40), 
     main = "Min Tb lab", xlab = NULL, ylab = "Frequency")
hist(T_set_data$WMEAN, xlim = c(10, 40), 
     main = "Mean Tb lab", xlab = NULL, ylab = "Frequency")
hist(T_set_data$TMAX,  xlim = c(10, 40), 
     main = "Max Tb lab", xlab = "Temperature (c)", ylab = "Frequency")

# Tb field data
Tb_field_data <- read.csv("Data/Mathews_etal_2023_Field_Tb.csv")
min(Tb_field_data$TB, na.rm = TRUE)
mean(Tb_field_data$TB, na.rm = TRUE)
max(Tb_field_data$TB, na.rm = TRUE)
par(mfrow = c(1, 1))
hist(Tb_field_data$TB, xlim = c(10, 40), 
     main = "Field Tb data", xlab = NULL, ylab = "Frequency")
range(Tb_field_data$TB)










###########
##### 1) mass growth rate of individual with most measurements
# grab data from lizard database
liz_kar_2023_dat <- read.csv(file = "Matlab/Raw_data/Kar_etal_2023.csv") %>% 
  filter(egg_incub_temp == 23)

# find a few individuals that have many remeasurments
liz_db_dat_n_id <- liz_kar_2023_dat %>% 
  group_by(bd_liz_id) %>% 
  summarise(n = n()) %>% 
  filter(n >17) # !!! Sex?

ids_kar_2023 <- unique(liz_db_dat_n_id$bd_liz_id)
ids_kar_2023

# filter these individuals with 19 observations: "ld0244"
kar_growth_filtered <- liz_kar_2023_dat %>% 
  filter(bd_liz_id %in% c("ld0250"))

# arranging dates for days since hatch
kar_growth_filtered$liz_hatch_date <- dmy(kar_growth_filtered$liz_hatch_date)
kar_growth_filtered$bd_date <- mdy(kar_growth_filtered$bd_date)

# Calculate date difference
kar_growth_filtered$date_difference <- as.numeric(kar_growth_filtered$bd_date - kar_growth_filtered$liz_hatch_date)

# rename cols
kar_growth_filtered_final <- kar_growth_filtered %>% 
  rename(ID = bd_liz_id, 
         trt_temp = egg_incub_temp,
         hatch_date = liz_hatch_date, 
         measure_date = bd_date, 
         days_since_hatch = date_difference, 
         mass_g = bd_mass,
         svl_mm = bd_svl,
         tl_mm = bd_tail, 
         hll_mm = bd_hlen,
         notes = process_notes) %>% 
  dplyr::select(ID, trt_temp, hatch_date, measure_date, days_since_hatch, 
         mass_g, svl_mm, tl_mm, hll_mm, notes)


# last two gravid???? see figure 3 kar et al., 2023
plot(kar_growth_filtered_final$days_since_hatch,kar_growth_filtered_final$mass_g)
