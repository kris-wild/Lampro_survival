pacman::p_load("dplyr","lme4", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "latex2exp", "tidybayes", "bayesplot", "rstanarm", "plotrix", "emmeans", "ggExtra", "marginaleffects", "gridExtra", "bayestestR","performance", "cowplot", "patchwork", "png", "magick", "ggplot2")

# import data (be sure to change Y/N to 0/1)
data <- read.csv("Data/Lizard_Db_final.csv", head = T) %>% 
  mutate(temp = as.character(temp),
         hatch_total = as.numeric(hatch_total))

# separate the two species into df
guich_df <-data %>% filter(species == "guichenoti") 
deli_df <- data %>% filter(species == "delicata")

# look at developmental time differences by species and temp
data_hatch <- data %>% 
  group_by(species, temp) %>% 
  summarise(dev_time_mean = mean(dev_time, na.rm = TRUE),
            dev_sd = sd(dev_time, na.rm = TRUE),
            ab_SVL = mean(hatch_svl, na.rm = TRUE), 
            ab_total = mean(hatch_total, na.rm = TRUE),
            ab_mass = mean(hatch_mass, na.rm = TRUE))


###########
# Functions
###########
## pMCMC Function
# Calculates the, p-value or pMCMC value for a posterior distribution
# x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
pmcmc <- function(x){
  2*(1 - max(table(x<0) / nrow(x)))
}


deli_df$temp_cont <- as.numeric(deli_df$temp)

######
# Deli Survival Analysis Numeric Temp 
######
# Deli survival model
Deli_m2_brms <- brm(mortality ~ egg_mass_final + temp_cont+ 
                      temp_cont:egg_mass_final   + # egg mass 
                      (1|clutch),  # clutch random factor
                    family = bernoulli(link = "logit"), 
                    data = deli_df, 
                    iter= 5000, warmup = 1000, 
                    thin = 4, cores = 8)


####################
# Model 1 checks: lags, residuals, r2, summary
####################
# checking lags for this model 
variable.names(post_Deli_m2)
draws <- as.array(Deli_m2_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_egg_mass_final", "b_temp_cont", "b_egg_mass_final:temp_cont"), lags =10)
# plots
plot(Deli_m2_brms)
#R2 and summary of full model
bayes_R2(Deli_m2_brms)
summary(Deli_m2_brms)

# stan plot
stanplot(Deli_m2_brms, 
         type = "areas",
         prob = 0.95)

# overall egg mass on survival
plot(conditional_effects(Deli_m2_brms, "egg_mass_final"), 
     ask = FALSE)

# incubation temperature X mass treatment interaction
plot(conditional_effects(Deli_m2_brms, 
                         "egg_mass_final:temp_cont", points = TRUE))








#######################################################
# Adult summary for matlab zero-variate data
######################################################

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
