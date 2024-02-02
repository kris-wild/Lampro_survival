library(ecmwfr)
library(mcera5)
library(lubridate)
library(dplyr)
library(tidync)
setwd("~/Desktop/")
uid <- "233381"
cds_api_key <- "b92ba7e9-f5db-48ce-bce7-62100e55b4a4"
ecmwfr::wf_set_key(user = uid, key = cds_api_key, service = "cds")

# bounding coordinates (in WGS84)

### NZ Whole island 2021
#-45.88673°N 166.54706°E
#-37.71943°N 178.41597°E
#-34.41862°N 172.74527°E
#-46.64581°N 168.83293°E
#xmn <- 166.54706
#xmx <- 178.41597
#ymn <- -46.64581
#ymx <- -34.41862

# ALL OF HAWAII 
#19.48649°N -154.79670°E
#21.79830°N -160.25860°E
#22.23602°N -159.52229°E
#18.91388°N -155.69869°E
xmn <- -160
xmx <- -154
ymn <- 23
ymx <- 18

# temporal extent
st_time <- lubridate::ymd("2016:01:01")
en_time <- lubridate::ymd("2016:12:31")


file_prefix <- "ERA5"
op <- "Era5/"

# build a request (covering multiple years)
req <- build_era5_request(xmin = xmn, xmax = xmx,
                          ymin = ymn, ymax = ymx,
                          start_time = st_time,
                          end_time = en_time,
                          outfile_name = file_prefix)
str(req)
request_era5(request = req, uid = uid, out_path = op)
