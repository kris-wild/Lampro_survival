pacman::p_load(NicheMapR, R.matlab, lubridate, dplyr, scales)


######## ######## ########
#### Lab data from incubation experiments
emperical_mass <- read.csv(file = "Data/Lizard_Db_final.csv") %>% 
  group_by(species, temp) %>% 
  summarise(mean_mass = mean(hatch_mass, na.rm = TRUE), 
            mass_sd = sd(hatch_mass, na.rm = TRUE),
            mean_dev_time = mean(dev_time, na.rm = TRUE), 
            dev_time_sd = sd(dev_time, na.rm = TRUE),
            n = n()) %>% 
  mutate(across(c(mean_mass, mass_sd, mean_dev_time, dev_time_sd), ~round(., 3)))
  


######## ######## ########
#### MICRO CLIMATE FOR IN SILICO EXPERIMENTS
# 23C
micro_23 <- readRDS('output/microclimates/micro_Syd_2016_2021.Rda')
columns_to_change <- setdiff(colnames(micro_23$soil), c("DOY", "TIME"))
micro_23$soil[, columns_to_change] <- 23

# 28C
micro_28 <- readRDS('output/microclimates/micro_Syd_2016_2021.Rda')
columns_to_change <- setdiff(colnames(micro_28$soil), c("DOY", "TIME"))
micro_28$soil[, columns_to_change] <- 28

########## In silico DEB experiment
# 1) Run ecotherm model, z constant and same as we have now, at two constant temps, 23 and 28C.
# - What is the body size at "hatch"/maturity?
# - What is the hatching time?
source('R/Functions/simulations.R')
debout_23 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata")
debout_23_hatch_info <- debout_23 %>% filter(DAY == debout_23$DAY[which(debout_23$E_H>E.Hb)[1]])
round(mean(debout_23_hatch_info$Wet_mass_g), digits = 3)


source('R/Functions/simulations.R')
debout_28 <- deb_experiment(micro = micro_28, species = "Lampropholis_delicata")
debout_28_hatch_info <- debout_28 %>% filter(DAY == debout_28$DAY[which(debout_28$E_H>E.Hb)[1]])
round(mean(debout_28_hatch_info$Wet_mass_g), digits = 3)


###################################
#####  QUICK PLOTS 
par(mfrow = c(2,1))
par(oma = c(2,2,1,1) + 0.1) # margin spacing stuff
par(mar = c(2,4,1,1) + 0.1) # margin spacing stuff

####  Plot 1: DEB PLOT @ 23 C
plot(debout_23$DAY, debout_23$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     xlim = c(0, 50), ylim = c(0, .2))
points(debout_23$DAY, debout_23$V, type = 'l', col = 'dark green', lwd = 2)
points(debout_23$DAY, debout_23$WETMASS-debout_23$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(debout_23$DAY, debout_23$WETMASS-debout_23$WETGONAD-debout_23$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = debout_23$DAY[which(debout_23$E_H>E.Hb)[1]], lty = 2, col = 'grey')
text((debout_23$DAY)[which(debout_23$E_H>E.Hb)[1]], .15, 'hatch', cex = 0.85, srt = 90)
legend("topleft",
       c('repro. buffer', 'food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("pink", "brown", "grey", "dark green"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 23C", line = -1) #(E_G -20%)

####  Plot 2: DEB PLOT @ 23 C
plot(debout_28$DAY, debout_28$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     xlim = c(0, 50), ylim = c(0, .2))
points(debout_28$DAY, debout_28$V, type = 'l', col = 'dark green', lwd = 2)
points(debout_28$DAY, debout_28$WETMASS-debout_28$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(debout_28$DAY, debout_28$WETMASS-debout_28$WETGONAD-debout_28$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = debout_28$DAY[which(debout_28$E_H>E.Hb)[1]], lty = 2, col = 'grey')
text((debout_28$DAY)[which(debout_28$E_H>E.Hb)[1]], .15, 'hatch', cex = 0.85, srt = 90)
legend("topleft",
       c('repro. buffer', 'food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("pink", "brown", "grey", "dark green"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 28 C", line = -1) #(E_G +30%)
