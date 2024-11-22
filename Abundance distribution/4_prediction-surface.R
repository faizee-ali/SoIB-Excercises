#### make the prediction surface ####
## to gerenate a regular grid of habitat covariates over which to make predictions

## step 1 generate an empty raster template 
WG <- read_sf("../shapefiles/S_WG.shp")

agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
    aggregate(agg_factor) 
r <- WG %>% 
    st_transform(crs = projection(r)) %>% 
    rasterize(r, field = 1) %>% 
    # remove any empty cells at edges
    trim()
r <- writeRaster(r, filename = "data/prediction-surface.tif", overwrite = TRUE)


### step 2 ; for each cell of this raster, we’ll calculate the PLAND metrics using the same approach

# get cell centers and create neighborhoods
r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
    st_as_sf() %>% 
    transmute(id = row_number())
r_cells <- st_buffer(r_centers, dist = neighborhood_radius)

# extract landcover values within neighborhoods, only needed most recent year
lc_extract_pred <- landcover[[paste0("y", max_lc_year)]] %>% 
    exact_extract(r_cells, progress = FALSE) %>% 
    purrr::map(~ count(., landcover = value)) %>% 
    tibble(id = r_cells$id, data = .) %>% 
    unnest(data)

# calculate the percent for each landcover class
pland_pred <- lc_extract_pred %>% 
    
    group_by(id) %>% 
    mutate(pland = n / sum(n)) %>% 
    ungroup() %>% 
    select(-n) %>% 
    # remove NAs after tallying so pland is relative to total number of cells
    filter(!is.na(landcover))

# convert names to be more descriptive
pland_pred <- pland_pred %>% 
    inner_join(lc_names, by = "landcover") %>% 
    arrange(landcover) %>% 
    select(-landcover)

# tranform to wide format, filling in implicit missing values with 0s
pland_pred <- pland_pred %>% 
    pivot_wider(names_from = lc_name, 
                values_from = pland, 
                values_fill = list(pland = 0)) %>% 
    mutate(year = max_lc_year) %>% 
    select(id, year, everything())

# join in coordinates
pland_coords <- st_transform(r_centers, crs = 4326) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    cbind(id = r_centers$id, .) %>% 
    rename(longitude = X, latitude = Y) %>% 
    inner_join(pland_pred, by = "id")

### import elevation data and crop 

elev <- raster("data/cut_n00e060.tif")
# crop, buffer WG by 10 km to provide a little wiggly room
elev <- WG %>% 
    st_buffer(dist = 10000) %>% 
    st_transform(crs = projection(elev)) %>% 
    crop(elev, .) %>% 
    projectRaster(crs = projection(landcover))

###Now we extract the elevation values within the neighborhood of each checklist location

# buffer each checklist location
ebird_buff_noyear <- ebird %>% 
    distinct(locality_id, latitude, longitude) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(elev)) %>% 
    st_buffer(dist = neighborhood_radius)

# extract elevation values and calculate median and sd
locs <- st_set_geometry(ebird_buff_noyear, NULL) %>% 
    mutate(id = row_number())
elev_checklists <- exact_extract(elev, ebird_buff_noyear, progress = FALSE) %>% 
    map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                     elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get locality_id
    bind_cols(locs, .)


###calculate the elevation covariates for the prediction surface.

# extract and calculate median and sd
elev_pred <- exact_extract(elev, r_cells, progress = FALSE) %>% 
    map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                     elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get locality_id
    bind_cols(st_drop_geometry(r_cells), .)

#### finally combine these elevation covariates with the land cover covariates.
# checklist covariates
pland_elev_checklist <- inner_join(pland, elev_checklists, by = "locality_id")
write_csv(pland_elev_checklist, "data/pland-elev_location-year.csv")

# prediction surface covariates
pland_elev_pred <- inner_join(pland_coords, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/pland-elev_prediction-surface.csv")
glimpse(pland_elev_pred)
