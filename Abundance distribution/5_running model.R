#### Modeling Relative Abundance ####

#Data preparation

library(lubridate)
library(sf)
library(raster)
library(dggridR)
library(pdp)
library(mgcv)
library(fitdistrplus)
library(viridis)
library(fields)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# set random number seed to insure fully repeatable results
set.seed(1)
# setup output directory for saved results
if (!dir.exists("output")) {
    dir.create("output")
}

# ebird data
ebird <- read_csv("cuce_ebird_s_zf.csv") %>% 
    mutate(protocol_type = factor(protocol_type, 
                                  levels = c("Stationary" , "Traveling"))) %>%
    # remove observations with no count
    filter(!is.na(observation_count))
ebird$year <- format(ebird$observation_date, "%Y")
ebird <- ebird[ebird$year %in% c("2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"), ]
ebird$year <- as.numeric(ebird$year)

# modis habitat covariates
habitat <- read_csv("data/pland-elev_location-year.csv") %>% 
    mutate(year = as.integer(year))

# combine ebird and habitat data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# prediction surface
pred_surface <- read_csv("data/pland-elev_prediction-surface.csv")
# latest year of landcover data
max_lc_year <- pred_surface$year[1]
r <- raster("data/prediction-surface.tif")

# load gis data for making maps
map_proj <- st_crs(4326)
ne_land <- read_sf("../shapefiles/S_WG.shp") %>% 
    st_transform(crs = map_proj) %>% 
    st_geometry()
ne_country_lines  <- ne_countries(country = "India", returnclass = "sf") %>% 
    st_geometry()
ne_state_lines <- ne_states(country = "India") %>% 
    st_transform(crs = map_proj) %>% 
    st_geometry()

#### Spatiotemporal subsampling has alrady been done ####

## Test-train split ##

#split the data into 80% of checklists for training and 20% for testing

#At this stage, we’ll also retain only the variables that we’ll use as covariates in the models.

hab_covs <- c("pland_04_deciduous_broadleaf", 
              "pland_05_mixed_forest",
              "pland_12_cropland",
              "pland_13_urban",
              "pland_02_evergreen_broadleaf",
              "pland_10_grassland",
              "elevation_median")

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
    x <- hms(x, quiet = TRUE)
    hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
ebird_habitat <- ebird_habitat %>% 
    mutate(
        observation_count = as.integer(observation_count),
        # effort_distance_km to 0 for non-travelling counts
        effort_distance_km = if_else(protocol_type != "Traveling", 
                                     0, effort_distance_km),
        # convert time to decimal hours since midnight
        time_observations_started = time_to_decimal(time_observations_started),
        # split date into year and day of year
        year = year(observation_date),
        day_of_year = yday(observation_date)
    )

ebird_split <- ebird_habitat %>% 
    # select only the columns to be used in the model
    select(observation_count,
           # effort covariates
           day_of_year, time_observations_started, duration_minutes,
           effort_distance_km, number_observers, protocol_type,
           # habitat covariates
           hab_covs) 
# split 80/20
ebird_split <- ebird_split %>% 
    split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
map_int(ebird_split, nrow)
#>  test train 1731  6633 


#### Exploratory data analysis ####
p <- par(mfrow = c(1, 2))
# counts with zeros
hist(ebird_habitat$observation_count, main = "Histogram of counts", 
     xlab = "Observed count")
# counts without zeros
pos_counts <- keep(ebird_habitat$observation_count, ~ . > 0)
hist(pos_counts, main = "Histogram of counts > 0", 
     xlab = "Observed non-zero count")
par(p)

#### Abundance models ####

### Using GAMs ###
## preparation ##
#joins the ends of the relationship, i.e. 0 and 24 are both midnight

# gam parameters
# degrees of freedom for smoothing
k <- 11 # as we have more habitat covariates in our data as compared to the book example

# degrees of freedom for cyclic time of day smooth
k_time <- 7 

# continuous predictors
# hold out time to treat seperately since it's cyclic
continuous_covs <- ebird_split$train %>% 
    select(-observation_count, -protocol_type, -time_observations_started) %>% 
    names()

# create model formula for predictors

gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                            var = continuous_covs, k = k) %>% 
    str_flatten(collapse = " + ") %>% 
    str_glue(" ~ ", .,
             " + protocol_type + ",
             "s(time_observations_started, bs = \"cc\", k = {k})", 
             k = k_time) %>% 
    as.formula()
# model formula including response
gam_formula <- update.formula(observation_count ~ ., gam_formula_rhs)
gam_formula

#observation_count ~ s(day_of_year, k = 7) + s(duration_minutes, 
#                                              k = 7) + s(effort_distance_km, k = 7) + s(number_observers, 
#                                                                                        k = 7) + s(pland_04_deciduous_broadleaf, k = 7) + s(pland_05_mixed_forest, 
#                                                                                                                                            k = 7) + s(pland_12_cropland, k = 7) + s(pland_13_urban, 
#                                                                                                                                                                                     k = 7) + s(pland_02_evergreen_broadleaf, k = 7) + s(pland_10_grassland, 
#                                                                                                                                                                                                                                         k = 7) + s(elevation_median, k = 7) + protocol_type + s(time_observations_started, 
#                                                                                                                                                                                                                                                                                                 bs = "cc", k = 7)                                                                                    k = 7) + s(pland_12_cropland, k = 7) + s(pland_13_urban, 
#                                                                                                                                       k = 7) + s(pland_02_evergreen_broadleaf, k = 7) + s(pland_10_grassland, 
#                                                                                                                                                                                           k = 7) + s(elevation_median, k = 7) + protocol_type + s(time_observations_started, 
#                                                                                                                                                                                                                                                   bs = "cc", k = 7)

#### Now we’ll use this formula to fit GAM models, testing the following three count response distributions: ####
# Zero-inflated Poisson:
# Negative binomial:
# Tweedie distribution:

# explicitly specify where the knots should occur for time_observations_started
# this ensures that the cyclic spline joins the variable at midnight
# this won't happen by default if there are no data near midnight
time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))

# zero-inflated poisson
m_ziplss <- gam(list(gam_formula,      # count model
                     gam_formula[-2]), # presence model
                data = ebird_split$train, 
                family = "ziplss", 
                knots = time_knots,
                method = "REML")
## model did not merge####
m_tw <- gam(gam_formula,
            data = ebird_split$train, 
            family = "tw",
            knots = time_knots)
# negative binomial
m_nb <- gam(gam_formula,
            data = ebird_split$train, 
            family = "nb",
            knots = time_knots)
save(m_nb, file = "mn_b.Rdata")
## this model merged and it was the best predicting model for the book example as well 

#### use this model to map GHCF relative abundance IN SWG ####

### first we need to bring effort variables into this prediction surface
## to the pediction surface we created  - 
## standard eBird checklist: a 1 km, 1 hour traveling count at the peak time of day for detecting this species##

