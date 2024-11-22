library(tidyverse)
library(raster)
# load the landcover data
landcover <- list.files("tifs/", pattern = ".tif", 
                        full.names = TRUE) 
landcover <- stack(landcover)
# label layers with year
name <- c("y2014", "y2015", "y2016", "y2017", "y2018", "y2019", "y2020", "y2021", "y2022", "y2023")

names(landcover) <- name
landcover

### landcover data upto 2023

library(rgdal)
library(sf)
library(readr)
library(dplyr)

### filter the data only within WG
library(auk)

###zero filling####
data = read.delim("cuce_ebird_w_zf.csv", sep = ",")
WG <- read_sf("../shapefiles/S_WG.shp")
ebd_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)



WG_sf <- st_as_sf(WG)
WG_sf <- st_transform(WG_sf, crs = st_crs(ebd_sf))


points_in_WG <- st_within(ebd_sf, WG_sf, sparse = FALSE)
# subset data frame
ebd_in_WG <- data[points_in_WG, ]
#saving the WG subsetted sampling file  

write_tsv(
    ebd_in_WG,
    "ebd_w_WG_zf.txt")

ebird<- read_tsv("ebd_w_WG_zf.txt")
rm(data, ebd_sf, points_in_WG, WG_sf)
#### binning the landcover data ####
ebird$year <- year(ebird$observation_date)
max_lc_year <- "2023"
a <- c("2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")

ebird <- ebird[ebird$year %in% a,]
ebird_buff <- ebird %>% 
    distinct(year,
             locality_id, latitude, longitude) %>% 
    # for 2019 use 2018 landcover data
    mutate(year_lc = paste0("y", year)) %>% 
    # convert to spatial features
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    # transform to modis projection
    st_transform(crs = projection(landcover)) %>% 
    # buffer to create neighborhood around each point
    st_buffer(dist = neighborhood_radius) %>% 
    # nest by year
    nest(data = c(year, locality_id, geometry))
library(exactextractr)
library(purrr)

#### debug function ####

calculate_pland <- function(yr, regions, lc) {
    locs <- st_set_geometry(regions, NULL)
    exact_extract(lc[[yr]], regions, progress = FALSE) %>% 
        purrr::map(~ count(., landcover = value)) %>% 
        tibble(locs, data = .) %>% 
        unnest(data)
}
# iterate over all years extracting landcover for all checklists in each
#### remove years which do not have LC ####

lc_extract <- ebird_buff %>% 
    mutate(pland = map2(year_lc, data, calculate_pland, lc = landcover)) %>% 
    select(pland) %>% 
    unnest(cols = pland)
###set of land cover values within a neighborhood around each checklist location.

pland <- lc_extract %>% 
    # calculate proporiton
    group_by(locality_id, year) %>% 
    mutate(pland = n / sum(n)) %>% 
    ungroup() %>% 
    select(-n) %>% 
    # remove NAs after tallying so pland is relative to total number of cells
    filter(!is.na(landcover))


###convert the numeric landcover codes to more descriptive names

lc_names <- tibble(landcover = 0:15,
                   lc_name = c("pland_00_water", 
                               "pland_01_evergreen_needleleaf", 
                               "pland_02_evergreen_broadleaf", 
                               "pland_03_deciduous_needleleaf", 
                               "pland_04_deciduous_broadleaf", 
                               "pland_05_mixed_forest",
                               "pland_06_closed_shrubland", 
                               "pland_07_open_shrubland", 
                               "pland_08_woody_savanna", 
                               "pland_09_savanna", 
                               "pland_10_grassland", 
                               "pland_11_wetland", 
                               "pland_12_cropland", 
                               "pland_13_urban", 
                               "pland_14_mosiac", 
                               "pland_15_barren"))
pland <- pland %>% 
    inner_join(lc_names, by = "landcover") %>% 
    arrange(landcover) %>% 
    select(-landcover)

# tranform to wide format, filling in implicit missing values with 0s%>% 
pland <- pland %>% 
    pivot_wider(names_from = lc_name, 
                values_from = pland, 
                values_fill = list(pland = 0))

# save
write_csv(pland, "data/modis_pland_location-year.csv")



