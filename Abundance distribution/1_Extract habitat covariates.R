#### getting habitat covariates for SWG ####

library(sf)
library(raster)
library(MODIS) 
library(exactextractr)
library(viridis)
library(tidyverse)

# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection
MODISoptions(localArcPath = "~/.modis")
MODISoptions(outDirPath = "D:/work/SOIB project/Abundance/data/modis")
MODIS::EarthdataLogin(usr = "faizeeali", pwd = "NASAwaleno_1")
MODISoptions(check_earthdata_login = TRUE)
MODISoptions(MODISserverOrder = "LAADS", quiet = FALSE)
# load ebird data
ebird <- read_tsv("ebd_sample_focused.txt")

# get list of tiles required to cover WG
WG <- read_sf("../shapefiles/S_WG.shp")

tiles <- getTile(WG)
tiles@tile
#> [1] "h25v08" "h25v07"

# earliest year of ebird data
begin_year <- format(min(ebird$`OBSERVATION DATE`), "%Y.01.01")
# end date for ebird data
end_year <- format(max(ebird$observation_date), "%Y.12.31")
# download tiles and combine into a single raster for each year
WG_buffered <- WG %>% st_buffer(dist = 10000) %>% st_make_valid()

tifs <- runGdal(product = "MCD12Q1", collection = "061", SDSstring = "01", 
                extent = WG_buffered, 
                begin = "2014.01.01", end = "2023.12.31", 
                outDirPath = "data", job = "modis",
                MODISserverOrder = c("LPDAAC"),
                verbose = TRUE) %>% 
    pluck("MCD12Q1.061") %>% 
    unlist()

# rename tifs to have more descriptive names
new_names <- format(as.Date(names(tifs)), "%Y") %>% 
    sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
    file.path(dirname(tifs), .)
file.rename(tifs, new_names)