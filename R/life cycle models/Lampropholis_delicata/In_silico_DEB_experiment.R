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

########## In silico DEB experiment - NULL
# 1) Null model at two constant temps, 23 and 28C.
# - What is the body size at "hatch"/maturity?
# - What is the hatching time?.
source('R/Functions/simulations.R')
# NULL: 23C
debout_23 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata")
debout_23_hatch_info <- debout_23 %>% filter(DAY == debout_23$DAY[which(debout_23$E_H>E.Hb)[1]])
null_23_mass <- round(mean(debout_23_hatch_info$Wet_mass_g), digits = 3)
null_23_dev <- unique(debout_23_hatch_info$DAY)

# NULL: 28C
debout_28 <- deb_experiment(micro = micro_28, species = "Lampropholis_delicata")
debout_28_hatch_info <- debout_28 %>% filter(DAY == debout_28$DAY[which(debout_28$E_H>E.Hb)[1]])
null_28_mass <- round(mean(debout_28_hatch_info$Wet_mass_g), digits = 3)
null_28_dev <- unique(debout_28_hatch_info$DAY)

# data
null_23_mass 
null_23_dev
null_28_mass
null_28_dev

# EO 
E.0
E.0_sim  <-  data.frame(E.0_sim =round(rnorm(100, mean=E.0, sd=15)))


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
       c('food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("brown", "grey", "darkgreen"), bty = 'n', cex = 0.5)
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
       c('food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("brown", "grey", "darkgreen"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 28 C", line = -1) 



########## In silico DEB experiment - Structure costs (EG) by temp
# 2) Structure costs (EG) model under varying temperature conditions: 
# Does decreasing EG for cold incubated individuals and 
# increasing EG hot incubated individuals accurately predict 
# - What is the body mass at "hatch"
# - What is the hatching time?
source('R/Functions/simulations.R')

########################### 23C
######## 23C EG 10% decrease
# change EG*0.9
debout_23_EG_10 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",
                                  E.G_mult = 0.9)
debout_23_hatch_info_EG_10 <- debout_23_EG_10 %>% 
  filter(DAY == debout_23_EG_10$DAY[which(debout_23_EG_10$E_H>E.Hb)[1]])
EG_10_dev_23C <- mean(debout_23_hatch_info_EG_10$DAY)
EG_10_mass_23C <- round(mean(debout_23_hatch_info_EG_10$Wet_mass_g), digits = 3)
EG_10_dev_23C
EG_10_mass_23C

######## 23C EG 20% decrease
# change EG*0.8
debout_23_EG_20 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",
                                  E.G_mult = 0.8)
debout_23_hatch_info_EG_20 <- debout_23_EG_20 %>% filter(DAY == debout_23_EG_20$DAY[which(debout_23_EG_20$E_H>E.Hb)[1]])
EG_20_dev_23C <- mean(debout_23_hatch_info_EG_20$DAY)
EG_20_mass_23C <- round(mean(debout_23_hatch_info_EG_20$Wet_mass_g), digits = 3)
EG_20_dev_23C
EG_20_mass_23C

######## 23C EG 30% decrease
# change EG*0.7
debout_23_EG_30 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",
                                  E.G_mult = 1)
debout_23_hatch_info_EG_30 <- debout_23_EG_30 %>% filter(DAY == debout_23_EG_30$DAY[which(debout_23_EG_30$E_H>E.Hb)[1]])
EG_30_dev_23C <- mean(debout_23_hatch_info_EG_30$DAY)
EG_30_mass_23C <- round(mean(debout_23_hatch_info_EG_30$Wet_mass_g), digits = 3)
EG_30_dev_23C
EG_30_mass_23C

########################### 28C 
######## 28C EG 10% increase
# change EG*1.1
debout_28_EG_10 <- deb_experiment(micro = micro_28, species = "Lampropholis_delicata",
                                  E.G_mult = 1.1)
debout_28_hatch_info_EG_10 <- debout_28_EG_10 %>% filter(DAY == debout_28_EG_10$DAY[which(debout_28_EG_10$E_H>E.Hb)[1]])
EG_10_dev_28C <- mean(debout_28_hatch_info_EG_10$DAY)
EG_10_mass_28C <- round(mean(debout_28_hatch_info_EG_10$Wet_mass_g), digits = 3)
EG_10_dev_28C
EG_10_mass_28C


########################### Plots
####  Plot 1A:  23C EO 10%
plot(debout_23_EG_10$DAY, debout_23_EG_10$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     xlim = c(0, 50), ylim = c(0, .2))
points(debout_23_EG_10$DAY, debout_23_EG_10$V, type = 'l', col = 'dark green', lwd = 2)
points(debout_23_EG_10$DAY, debout_23_EG_10$WETMASS-debout_23_EG_10$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(debout_23_EG_10$DAY, debout_23_EG_10$WETMASS-debout_23_EG_10$WETGONAD-debout_23_EG_10$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = debout_23_EG_10$DAY[which(debout_23_EG_10$E_H>E.Hb)[1]], lty = 2, col = 'grey')
text((debout_23_EG_10$DAY)[which(debout_23_EG_10$E_H>E.Hb)[1]], .15, 'hatch', cex = 0.85, srt = 90)
legend("topleft",
       c('food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("brown", "grey", "darkgreen"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 23C (E_G -10%)", line = -1)

####  Plot 1B:  28C EO 10%
plot(debout_28_EG_10$DAY, debout_28_EG_10$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     xlim = c(0, 50), ylim = c(0, .2))
points(debout_28_EG_10$DAY, debout_28_EG_10$V, type = 'l', col = 'dark green', lwd = 2)
points(debout_28_EG_10$DAY, debout_28_EG_10$WETMASS-debout_28_EG_10$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(debout_28_EG_10$DAY, debout_28_EG_10$WETMASS-debout_28_EG_10$WETGONAD-debout_28_EG_10$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = debout_28_EG_10$DAY[which(debout_28_EG_10$E_H>E.Hb)[1]], lty = 2, col = 'grey')
text((debout_28_EG_10$DAY)[which(debout_28_EG_10$E_H>E.Hb)[1]], .15, 'hatch', cex = 0.85, srt = 90)
legend("topleft",
       c('food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("brown", "grey", "darkgreen"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 23C (E_G +10%)", line = -1)


########################### Plots
####  Plot 2A:  23C EO 20%
plot(debout_23_EG_20$DAY, debout_23_EG_20$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     xlim = c(0, 50), ylim = c(0, .2))
points(debout_23_EG_20$DAY, debout_23_EG_20$V, type = 'l', col = 'dark green', lwd = 2)
points(debout_23_EG_20$DAY, debout_23_EG_20$WETMASS-debout_23_EG_20$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(debout_23_EG_20$DAY, debout_23_EG_20$WETMASS-debout_23_EG_20$WETGONAD-debout_23_EG_20$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = debout_23_EG_20$DAY[which(debout_23_EG_20$E_H>E.Hb)[1]], lty = 2, col = 'grey')
text((debout_23_EG_20$DAY)[which(debout_23_EG_20$E_H>E.Hb)[1]], .15, 'hatch', cex = 0.85, srt = 90)
legend("topleft",
       c('food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("brown", "grey", "darkgreen"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 23C (E_G -20%)", line = -1)

####  Plot 2B:  28C EO 20%
plot(debout_28_EG_20$DAY, debout_28_EG_20$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     xlim = c(0, 50), ylim = c(0, .2))
points(debout_28_EG_20$DAY, debout_28_EG_20$V, type = 'l', col = 'dark green', lwd = 2)
points(debout_28_EG_20$DAY, debout_28_EG_20$WETMASS-debout_28_EG_20$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(debout_28_EG_20$DAY, debout_28_EG_20$WETMASS-debout_28_EG_20$WETGONAD-debout_28_EG_20$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = debout_28_EG_20$DAY[which(debout_28_EG_20$E_H>E.Hb)[1]], lty = 2, col = 'grey')
text((debout_28_EG_20$DAY)[which(debout_28_EG_20$E_H>E.Hb)[1]], .15, 'hatch', cex = 0.85, srt = 90)
legend("topleft",
       c('food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("brown", "grey", "darkgreen"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 23C (E_G +20%)", line = -1)


########################### Plots
####  Plot 3A:  23C EO 30%
plot(debout_23_EG_30$DAY, debout_23_EG_30$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     xlim = c(0, 50), ylim = c(0, .2))
points(debout_23_EG_30$DAY, debout_23_EG_30$V, type = 'l', col = 'dark green', lwd = 2)
points(debout_23_EG_30$DAY, debout_23_EG_30$WETMASS-debout_23_EG_30$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(debout_23_EG_30$DAY, debout_23_EG_30$WETMASS-debout_23_EG_30$WETGONAD-debout_23_EG_30$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = debout_23_EG_30$DAY[which(debout_23_EG_30$E_H>E.Hb)[1]], lty = 2, col = 'grey')
text((debout_23_EG_30$DAY)[which(debout_23_EG_30$E_H>E.Hb)[1]], .15, 'hatch', cex = 0.85, srt = 90)
legend("topleft",
       c('food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("brown", "grey", "darkgreen"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 23C (E_G -30%)", line = -1)

####  Plot 3B:  28C EO 30%
plot(debout_28_EG_30$DAY, debout_28_EG_30$WETMASS, type = 'l', xlab = 'date',
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2,
     xlim = c(0, 50), ylim = c(0, .2))
points(debout_28_EG_30$DAY, debout_28_EG_30$V, type = 'l', col = 'dark green', lwd = 2)
points(debout_28_EG_30$DAY, debout_28_EG_30$WETMASS-debout_28_EG_30$WETGONAD, type = 'l',
       lwd = 2, col = 'brown')
points(debout_28_EG_30$DAY, debout_28_EG_30$WETMASS-debout_28_EG_30$WETGONAD-debout_28_EG_30$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = debout_28_EG_30$DAY[which(debout_28_EG_30$E_H>E.Hb)[1]], lty = 2, col = 'grey')
text((debout_28_EG_30$DAY)[which(debout_28_EG_30$E_H>E.Hb)[1]], .15, 'hatch', cex = 0.85, srt = 90)
legend("topleft",
       c('food in gut', 'reserve', 'structure'), lty = rep(1, 4),
       col = c("brown", "grey", "darkgreen"), bty = 'n', cex = 0.5)
text(10, .12, labels = "embryo", cex = 0.85)
title(main = "L. delicata development 23C (E_G +30%)", line = -1)


############## Final data
# mass data
null_23_mass 
null_28_mass
EG_10_mass_23C
EG_10_mass_28C
EG_20_mass_23C
EG_20_mass_28C
EG_30_mass_23C
EG_30_mass_28C


# development data
null_23_dev
null_28_dev
EG_10_dev_23C
EG_10_dev_28C



########## In silico DEB experiment - EO
# 3) Holding EG @20% (+/-) incubation temperature
# how does decreasing energy content of the egg (EO) effect 
# - What is the body size at "hatch"/maturity?
# - What is the hatching time?
#################### 23C EG constant 0.8
######## 23C EO 10% decrease
source('R/Functions/simulations.R')
debout_23_E0_10 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",
                                  E.G_mult = 0.8, # constant multiplier for cool temps
                                  E.0_mult = 0.9)
debout_23_hatch_info_E0_10 <- debout_23_E0_10 %>% filter(DAY == debout_23_E0_10$DAY[which(debout_23_E0_10$E_H>E.Hb)[1]])
E0_10_time_23 <- mean(debout_23_hatch_info_E0_10$DAY)
E0_10_mass_23 <- round(mean(debout_23_hatch_info_E0_10$Wet_mass_g), digits = 3)
E0_10_mass_23
E0_10_time_23

######## 23C EO 20% decrease
debout_23_E0_20 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",
                                  E.G_mult = 0.8, # constant multiplier for cool temps
                                  E.0_mult = 0.8)
debout_23_hatch_info_20 <- debout_23_E0_20 %>% filter(DAY == debout_23_E0_20$DAY[which(debout_23_E0_20$E_H>E.Hb)[1]])
E0_20_time_23 <- mean(debout_23_hatch_info_20$DAY)
E0_20_mass_23 <- round(mean(debout_23_hatch_info_20$Wet_mass_g), digits = 3)
E0_20_mass_23
E0_20_time_23

######## 23C EO 30% decrease
debout_23_E0_30 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",
                                  E.G_mult = 0.8, # constant multiplier for cool temps
                                  E.0_mult = 0.7)
debout_23_hatch_info_30 <- debout_23_E0_30 %>% filter(DAY == debout_23_E0_30$DAY[which(debout_23_E0_30$E_H>E.Hb)[1]])
E0_30_time_23 <- mean(debout_23_hatch_info_30$DAY)
E0_30_mass_23 <- round(mean(debout_23_hatch_info_30$Wet_mass_g), digits = 3)
E0_30_mass_23
E0_30_time_23

######## 23C EO 40% decrease
debout_23_E0_40 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",
                                  E.G_mult = 0.8, # constant multiplier for cool temps
                                  E.0_mult = 0.6)
debout_23_hatch_info_40 <- debout_23_E0_40 %>% filter(DAY == debout_23_E0_40$DAY[which(debout_23_E0_40$E_H>E.Hb)[1]])
E0_40_time_23 <- mean(debout_23_hatch_info_40$DAY)
E0_40_mass_23 <- round(mean(debout_23_hatch_info_40$Wet_mass_g), digits = 3)
E0_40_mass_23
E0_40_time_23

######## 23C EO 50% decrease
debout_23_E0_50 <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",
                                  E.G_mult = 0.8, # constant multiplier for cool temps
                                  E.0_mult = 0.5)
debout_23_hatch_info_50 <- debout_23_E0_50 %>% filter(DAY == debout_23_E0_50$DAY[which(debout_23_E0_50$E_H>E.Hb)[1]])
E0_50_time_23 <- mean(debout_23_hatch_info_50$DAY)
E0_50_mass_23 <- round(mean(debout_23_hatch_info_50$Wet_mass_g), digits = 3)
E0_50_mass_23
E0_50_time_23

# mass change
E0_10_mass_23
E0_20_mass_23
E0_30_mass_23
E0_40_mass_23
E0_50_mass_23
# development time change
E0_10_time_23
E0_20_time_23
E0_30_time_23
E0_40_time_23
E0_50_time_23


#################### 28C EG constant 1.2
######## 28C EO 10% decrease
debout_28_E0_10 <- deb_experiment(micro = micro_28, species = "Lampropholis_delicata",
                                  E.G_mult = 1.2, # constant multiplier for cool temps
                                  E.0_mult = .9)
debout_28_hatch_info_E0_10 <- debout_28_E0_10 %>% filter(DAY == debout_28_E0_10$DAY[which(debout_28_E0_10$E_H>E.Hb)[1]])
E0_10_time_28 <- mean(debout_28_hatch_info_E0_10$DAY)
E0_10_mass_28 <- round(mean(debout_28_hatch_info_E0_10$Wet_mass_g), digits = 3)
E0_10_mass_28
E0_10_time_28

######## 28C EO 20% decrease
debout_28_E0_20 <- deb_experiment(micro = micro_28, species = "Lampropholis_delicata",
                                  E.G_mult = 1.2, # constant multiplier for cool temps
                                  E.0_mult = .8)
debout_28_hatch_info_20 <- debout_28_E0_20 %>% filter(DAY == debout_28_E0_20$DAY[which(debout_28_E0_20$E_H>E.Hb)[1]])
E0_20_time_28 <- mean(debout_28_hatch_info_20$DAY)
E0_20_mass_28 <- round(mean(debout_28_hatch_info_20$Wet_mass_g), digits = 3)
E0_20_mass_28
E0_20_time_28

######## 28C EO 30% decrease
debout_28_E0_30 <- deb_experiment(micro = micro_28, species = "Lampropholis_delicata",
                                  E.G_mult = 1.2, # constant multiplier for cool temps
                                  E.0_mult = .7)
debout_28_hatch_info_30 <- debout_28_E0_30 %>% filter(DAY == debout_28_E0_30$DAY[which(debout_28_E0_30$E_H>E.Hb)[1]])
E0_30_time_28 <- mean(debout_28_hatch_info_30$DAY)
E0_30_mass_28 <- round(mean(debout_28_hatch_info_30$Wet_mass_g), digits = 3)
E0_30_mass_28
E0_30_time_28

######## 28C EO 40% decrease
debout_28_E0_40 <- deb_experiment(micro = micro_28, species = "Lampropholis_delicata",
                                  E.G_mult = 1.2, # constant multiplier for cool temps
                                  E.0_mult = .6)
debout_28_hatch_info_40 <- debout_28_E0_40 %>% filter(DAY == debout_28_E0_40$DAY[which(debout_28_E0_40$E_H>E.Hb)[1]])
E0_40_time_28 <- mean(debout_28_hatch_info_40$DAY)
E0_40_mass_28 <- round(mean(debout_28_hatch_info_40$Wet_mass_g), digits = 3)
E0_40_mass_28
E0_40_time_28

######## 28C EO 50% decrease
debout_28_E0_50 <- deb_experiment(micro = micro_28, species = "Lampropholis_delicata",
                                  E.G_mult = 1.2, # constant multiplier for cool temps
                                  E.0_mult = .5)
debout_28_hatch_info_50 <- debout_28_E0_50 %>% filter(DAY == debout_28_E0_50$DAY[which(debout_28_E0_50$E_H>E.Hb)[1]])
E0_50_time_28 <- mean(debout_28_hatch_info_50$DAY)
E0_50_mass_28 <- round(mean(debout_28_hatch_info_50$Wet_mass_g), digits = 3)
E0_50_mass_28
E0_50_time_28

# mass change
E0_10_mass_28
E0_20_mass_28
E0_30_mass_28
E0_40_mass_28
E0_50_mass_28
# development time change
E0_10_time_28
E0_20_time_28
E0_30_time_28
E0_40_time_28
E0_50_time_28



########################### 23C
######## 23C Zoom mult .9
# change EG*0.9
debout_23_zoom_90pct <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata",z_mult = .9)
debout_23_hatch_info_zoom_90pct <- debout_23_zoom_90pct %>% 
  filter(DAY == debout_23_zoom_90pct$DAY[which(debout_23_zoom_90pct$E_H>E.Hb)[1]])
zoom_90pct_dev_23C <- mean(debout_23_hatch_info_zoom_90pct$DAY)
zoom_90pct_mass_23C <- round(mean(debout_23_hatch_info_zoom_90pct$Wet_mass_g), digits = 3)
zoom_90pct_dev_23C
zoom_90pct_mass_23C


########################### 23C
######## 23C Zoom mult .8
debout_23_zoom_80pct <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata", z_mult = .8)
debout_23_hatch_info_zoom_80pct <- debout_23_zoom_80pct %>% 
  filter(DAY == debout_23_zoom_80pct$DAY[which(debout_23_zoom_80pct$E_H>E.Hb)[1]])
zoom_80pct_dev_23C <- mean(debout_23_hatch_info_zoom_80pct$DAY)
zoom_80pct_mass_23C <- round(mean(debout_23_hatch_info_zoom_80pct$Wet_mass_g), digits = 3)
zoom_80pct_dev_23C
zoom_80pct_mass_23C


########################### 23C
######## 23C Zoom mult .7
debout_23_zoom_70pct <- deb_experiment(micro = micro_23, species = "Lampropholis_delicata", z_mult = .7)
debout_23_hatch_info_zoom_70pct <- debout_23_zoom_70pct %>% 
  filter(DAY == debout_23_zoom_70pct$DAY[which(debout_23_zoom_70pct$E_H>E.Hb)[1]])
zoom_70pct_dev_23C <- mean(debout_23_hatch_info_zoom_70pct$DAY)
zoom_70pct_mass_23C <- round(mean(debout_23_hatch_info_zoom_70pct$Wet_mass_g), digits = 3)
zoom_70pct_dev_23C
zoom_70pct_mass_23C
