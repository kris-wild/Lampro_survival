library(NicheMapR)
genus <- "Lampropholis"
species <- "delicata"
species1 <- paste0(genus, "_", species)
species2 <- paste0(genus, " ", species) 
source(paste0('R/life cycle models/', species1, '/parameters_biophys.R'))
micro <- readRDS('output/microclimates/micro_HI_2016_2021.Rda')



# ecto_nichemapR - ecotherm mod with predicted Tb's
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

# save predicted Tb 
saveRDS(ecto, file = 'output/life cycle models/ecto_nichemapR_micro_HI_2016_2021.Rda')



# predicted Tb (black) vs TE/substrate (red) 
environ <- as.data.frame(ecto$environ)
plot(micro$dates, environ$TSUB, type = 'l', col = 2)
points(micro$dates, environ$TC, type = 'l')












##################### ###################### ###################### 
# set up to check with distrbution of field data from Mathews 2023
environ$dates <- micro$dates
environ$month <- month(micro$dates)
environ$year <- year(micro$dates)
filt <- environ %>% filter(dates >= "2019-09-01 00:00:00" & 
                             dates <= "2019-10-31 00:00:00") 

# bring in Tb data from hanging rock
Tb_dat <- read.csv("Data/Mathews_etal_2023_Field_Tb.csv")

# quick plot for comparison
par(mfrow = c(2,1))
hist(filt$TC)
hist(Tb_dat$TB)
