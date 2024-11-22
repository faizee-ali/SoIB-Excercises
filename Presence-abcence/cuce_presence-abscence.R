#### Spatial thinning
#library(auk)
library(sf)
library(dggridR)
library(lubridate)
library(tidyverse)
library(rnaturalearth)

# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)

#> Resolution: 13, Area (km^2): 31.9926151554038, Spacing (km): 5.58632116604266, CLS (km): 6.38233997895802
# read in data
ebird <- read_csv("cuce/summer/cuce_ebird_s_zf.csv") %>% 
    # get hexagonal cell id and week number for each checklist
    mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)

# sample a single checklist from each grid cell

ebird_ss <- ebird %>% 
    group_by(species_observed, cell) %>% 
    sample_n(size = 1) %>% 
    ungroup()

write_tsv(ebird_ss, "cuce-PA-summer.txt")

cuce_s_df <- read_tsv("summer/cuce-PA-summer.txt")
cuce_s_df <- st_as_sf(cuce_s_df, coords = c("longitude", "latitude"), crs = 4326)

#### make a presence-absence map ####

asia <- ne_countries(continent = "Asia", returnclass = "sf") %>% 
    st_geometry()

india <- ne_countries(country = "India", returnclass = "sf") %>% 
    st_geometry()

p <- ggplot() +
    # Background: Asia map
    geom_sf(data = asia, fill = "grey80", color = "white") +
    # Highlight: India map
    geom_sf(data = india, fill = "grey70", color = "white") +
    # Points: Presence/Absence
    geom_sf(data = cuce_s_df, aes(geometry = geometry, color = species_observed), alpha = 0.3, size = 0.3) +
    # Title and theme adjustments
    labs(
        title = "Gray-headed Canary Flycatcher Summer Observations",
        color = "Observation Status"  # Legend title
    ) +
    
    theme(
        legend.position = c(0.3, 0.05),  # Place legend inside (x = 30%, y = 05%)
        legend.justification = c("center", "bottom"),  # Align legend inside top-right corner
        legend.background = element_rect(fill = alpha("white", 0.7), color = "black"),  # Add transparent background
        legend.key = element_rect(fill = "transparent", color = NA),  # Transparent legend keys
        plot.title = element_text(
            hjust = 0.5,         
            size = 20,           
            face = "bold",       
            margin = margin(b = 10)  # Add some space below the title
        )) +
    # Add coordinate limits
    coord_sf(
        xlim = c(68, 123),  # Restrict longitude range
        ylim = c(-11, 40)
    ) +
    # Define colors for presence/absence in the legend
    scale_color_manual(
        values = c("TRUE" = "orange", "FALSE" = "grey20"),
        labels = c("Present", "Not reported")
    )
ggsave(p, filename = "../figures/cuce-s-PA.png",
       width =8, height = 7, units = "in", dpi = 900)


#### winter ####

# read in data
ebird <- read_csv("cuce/winter/cuce_ebird_w_zf.csv") %>% 
    # get hexagonal cell id and week number for each checklist
    mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)

# sample a single checklist from each grid cell

ebird_ss <- ebird %>% 
    group_by(species_observed, cell) %>% 
    sample_n(size = 1) %>% 
    ungroup()

write_tsv(ebird_ss, "cuce-PA-winter.txt")

#read in data
cuce_w_df <- read_tsv("winter/cuce-PA-winter.txt")
cuce_w_df <- st_as_sf(cuce_w_df, coords = c("longitude", "latitude"), crs = 4326)


#### make a presence-absence map ####


#plot
p <- ggplot() +
    # Background: Asia map
    geom_sf(data = asia, fill = "grey80", color = "white") +
    # Highlight: India map
    geom_sf(data = india, fill = "grey70", color = "white") +
    # Points: Presence/Absence
    geom_sf(data = cuce_w_df, aes(geometry = geometry, color = species_observed),alpha = 0.3, size = 0.3) +
    # Title and theme adjustments
    labs(
        title = "Gray-headed Canary Flycatcher Winter Observations",
        color = "Observation Status"  # Legend title
    ) +
    
    theme(
        legend.position = c(0.3, 0.05),  # Place legend inside (x = 30%, y = 05%)
        legend.justification = c("center", "bottom"),  # Align legend inside top-right corner
        legend.background = element_rect(fill = alpha("white", 0.7), color = "black"),  # Add transparent background
        legend.key = element_rect(fill = "transparent", color = NA),  # Transparent legend keys
        plot.title = element_text(
            hjust = 0.5,         
            size = 20,           
            face = "bold",       
            margin = margin(b = 10)  # Add some space below the title
        )) +
    # Add coordinate limits
    coord_sf(
        xlim = c(68, 123),  # Restrict longitude range
        ylim = c(-11, 40)
    ) +
    # Define colors for presence/absence in the legend
    scale_color_manual(
        values = c("TRUE" = "violet", "FALSE" = "grey20"),
        labels = c("Present", "Not reported")
    )
ggsave(p, filename = "../figures/cuce-w-PA.png",
       width =8, height = 7, units = "in", dpi = 900)
