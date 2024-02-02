library(NicheMapR)
library(mcera5)
loc <- c(-157.85, 21.30)
# SYD c(151.10, -33.77)
# Brisb c(174.7645, -36.8509) 
# NZ c(172.6306, - 43.5320)
# HI c(-157.85, 21.30)
longlat <- loc
ystart <- 2016 # start date needed
yfinish <- 2020
scenario <- 0
nyears <- yfinish - ystart + 4
REFL <-.15 #  0.3 # 
ERR <- 1.5 # model error tolerance
minshade <- 0 # minimum shade
maxshade <-25 # maximum shade
Usrhyt <- 0.005 # lizard height above ground
#spatial = "/Volumes/The Brain/ERA5/Australia_ERA5/"
spatial = 'Data/ERA_5/Hawaii/ERA5'

micro <- micro_era5(windfac = 1, RUF = 0.004, 
                    minshade = minshade, maxshade = maxshade, 
                    ERR = ERR, loc = loc, 
                    dstart = paste0('01/01/', ystart), 
                    dfinish = paste0('31/12/', yfinish), 
                    Usrhyt = Usrhyt, REFL = REFL, 
                    spatial = spatial, 
                    scenario = scenario, runshade = 1, 
                    cap = 1, soilgrids = 0)


# need to change the dates from UTC
if(sign(loc[1]) == -1){
  gmtzone <- "+"
}else{
  gmtzone <- ""
}
tz <- paste0("Etc/GMT", gmtzone, round(loc[1]/15*-1, 0))
attr(micro$dates, "tzone") <- tz


#tz <- "Pacific/Honolulu"
#attr(micro$dates, "tzone") <- tz
saveRDS(micro, file = 'output/microclimates/micro_HI_2016_2021.Rda')



