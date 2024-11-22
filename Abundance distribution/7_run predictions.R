#### runing the predictions ####
pred_model <- m_nb
# assuming the peak time of day for detecting is around 7 AM.

# create a dataframe of covariates with a range of start times
seq_tod <- seq(0, 24, length.out = 300)
tod_df <- ebird_split$train %>% 
    # find average pland habitat covariates
    select(starts_with("pland"), elevation_median) %>% 
    summarize_all(mean, na.rm = TRUE) %>% 
    ungroup() %>% 
    # use standard checklist
    mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
           duration_minutes = 60,
           effort_distance_km = 1,
           number_observers = 1,
           protocol_type = "Traveling",
           ) %>% 
    cbind(time_observations_started = seq_tod)

# predict at different start times
pred_tod <- predict(pred_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
    as_tibble() %>% 
    # calculate backtransformed confidence limits
    transmute(time_observations_started = seq_tod,
              pred = pred_model$family$linkinv(fit),
              pred_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
              pred_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit))

# find optimal time of day
t_peak <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# plot the partial dependence plot
ggplot(pred_tod) +
    aes(x = time_observations_started, y = pred,
        ymin = pred_lcl, ymax = pred_ucl) +
    geom_ribbon(fill = "grey80", alpha = 0.5) +
    geom_line() +
    geom_vline(xintercept = t_peak, color = "blue", linetype = "dashed") +
    labs(x = "Hours since midnight",
         y = "Predicted relative abundance",
         title = "Effect of observation start time on Wood Thrush reporting",
         subtitle = "Peak detectability shown as dashed blue line")



t_peak <- 7

# add effort covariates to prediction surface
pred_surface_eff <- pred_surface %>% 
    mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
           time_observations_started = t_peak,
           duration_minutes = 60,
           effort_distance_km = 1,
           number_observers = 1,
           protocol_type = "Traveling")

# predict
pred <- predict(pred_model, newdata = pred_surface_eff, 
                type = "link", 
                se.fit = TRUE) %>% 
    as_tibble() %>% 
    # calculate confidence limits and back transform
    transmute(abd = pred_model$family$linkinv(fit),
              abd_se = pred_model$family$linkinv(se.fit),
              abd_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
              abd_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit)) %>%
    # add to prediction surface
    bind_cols(pred_surface_eff, .) %>% 
    select(latitude, longitude, abd, abd_se, abd_lcl, abd_ucl)
###convert this data frame to spatial features using sf
r <- raster("data/prediction-surface.tif")
r_pred <- pred %>% 
    # convert to spatial features
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    select(abd, abd_se) %>% 
    st_transform(crs = projection(r)) %>% 
    # rasterize
    rasterize(r)
r_pred <- r_pred[[-1]]

# save the rasters
tif_dir <- "output"
if (!dir.exists(tif_dir)) {
    dir.create(tif_dir)
}
writeRaster(r_pred[["abd"]], 
            filename = file.path(tif_dir, "abundance-model_abd_cuce.tif"),
            overwrite = TRUE)
writeRaster(r_pred[["abd_se"]], 
            filename = file.path(tif_dir, "abundance-model_se_cuce.tif"), 
            overwrite = TRUE)

####Finally, let’s make a map!####
# latest year of landcover data
max_lc_year <- pred_surface$year[1]

# load gis data for making maps
map_proj <- projection(r)
WG <- read_sf("../shapefiles/S_WG.shp") %>% 
    st_transform(crs = map_proj) %>% 
    st_geometry()
ne_country_lines  <- ne_countries(country = "India", returnclass = "sf") %>% 
    st_geometry()
ne_state_lines <- ne_states(country = "India") %>% 
    st_transform(crs = map_proj) %>% 
    st_geometry()

# any expected abundances below this threshold are set to zero
zero_threshold <- 0.05

# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj, method = "ngb")



par(mfrow = c(2, 1))
for (nm in names(r_pred)) {
    r_plot <- r_pred_proj[[nm]]
    
    par(mar = c(3.5, 0.25, 0.25, 0.25))
    # set up plot area
    plot(WG, col = NA, border = NA)
    plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
    
    # modified plasma palette
    plasma_rev <- rev(plasma(25, end = 0.9))
    gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
    pal <- c(gray_int(4)[2], plasma_rev)
    
    # abundance vs. se
    if (nm == "abd") {
        title <- "Gray-headed Canary Flycatcher Relative Abundance"
        # set very low values to zero
        r_plot[r_plot <= zero_threshold] <- NA
        # log transform
        r_plot <- log10(r_plot)
        # breaks and legend
        mx <- ceiling(100 * cellStats(r_plot, max)) / 100
        mn <- floor(100 * cellStats(r_plot, min)) / 100
        brks <- seq(mn, mx, length.out = length(pal) + 1)
        lbl_brks <- sort(c(-2:2, mn, mx))
        lbls <- round(10^lbl_brks, 2)
    } else {
        title <- "Gray-headed Canary Flycatcher Abundance Uncertainty (SE)"
        # breaks and legend
        mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
        mn <- floor(1000 * cellStats(r_plot, min)) / 1000
        brks <- seq(mn, mx, length.out = length(pal) + 1)
        lbl_brks <- seq(mn, mx, length.out = 5)
        lbls <- round(lbl_brks, 2)
    }
    
    # abundance
    plot(r_plot, 
         col = pal, breaks = brks, 
         maxpixels = ncell(r_plot),
         legend = FALSE, add = TRUE)
    
    # legend
    par(new = TRUE, mar = c(0, 0, 0, 0))
    image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
               smallplot = c(0.25, 0.75, 0.06, 0.09),
               horizontal = TRUE,
               axis.args = list(at = lbl_brks, 
                                labels = lbls,
                                fg = "black", col.axis = "black",
                                cex.axis = 0.75, lwd.ticks = 0.5,
                                padj = -1.5),
               legend.args = list(text = title,
                                  side = 3, col = "black",
                                  cex = 1, line = 0))
}

nm <- "abd_se"
r_plot <- r_pred_proj[[nm]]

par(mar = c(0.25, 0.25, 0.25, 0.25))# set up plot area
plot(WG)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))

gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

# abundance vs. se
if (nm == "abd") {
    title <- "Gray-headed Canary Flycatcher Relative Abundance"
    # set very low values to zero
    r_plot[r_plot <= zero_threshold] <- NA
    # log transform
    r_plot <- log10(r_plot)
    # breaks and legend
    mx <- ceiling(100 * cellStats(r_plot, max)) / 100
    mn <- floor(100 * cellStats(r_plot, min)) / 100
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- sort(c(-2:2, mn, mx))
    lbls <- round(10^lbl_brks, 2)
} else {
    title <- "Gray-headed Canary Flycatcher Abundance Uncertainty (SE)"
    # breaks and legend
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    mn <- floor(1000 * cellStats(r_plot, min)) / 1000
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- seq(mn, mx, length.out = 5)
    lbls <- round(lbl_brks, 2)
}
title <- "Wood Thrush Abundance Uncertainty (SE)"
# breaks and legend
mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
mn <- floor(1000 * cellStats(r_plot, min)) / 1000

# set very low values to zero
r_plot[r_plot <= zero_threshold] <- NA
# log transform
r_plot <- log10(r_plot)
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- seq(mn, mx, length.out = 5)
lbls <- round(lbl_brks, 2)

# abundance
plot(r_plot, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_plot),
     legend = FALSE, add = TRUE)

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
           smallplot = c(0.25, 0.75, 0.06, 0.09),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, 
                            labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))
