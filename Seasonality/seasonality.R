library(auk)
library(sf)
library(tidyverse)
library(dplyr)
library(lubridate)


library(readr)
#library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
#library(forcats)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(gghalves)
library(stringr)

###zero filling####
data = read.delim("ebd_zf_sampl.txt", sep = "\t")
WG <- read_sf("../shapefiles/S_WG.shp")
ebd_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

points_in_WG <- st_within(ebd_sf, WG, sparse = FALSE)
# subset data frame
ebd_in_WG <- data[points_in_WG, ]


## tidy-up

time_to_decimal <- function(x) {
    x <- hms(x, quiet = TRUE)
    hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
ebd_in_WG <- ebd_in_WG %>% 
    mutate(
        # convert X to NA
        observation_count = if_else(observation_count == "X", 
                                    NA_character_, observation_count),
        observation_count = as.integer(observation_count),
        # effort_distance_km to 0 for non-travelling counts
        effort_distance_km = if_else(protocol_type != "Traveling", 
                                     0, effort_distance_km),
        # convert time to decimal hours since midnight
        time_observations_started = time_to_decimal(time_observations_started),
        # split date into year and day of year
        year = year(observation_date),
        day_of_year = yday(observation_date),
        month = month(observation_date)
    )
rm(data, ebd_sf, points_in_WG)
#saving the WG subsetted sampling file  

write_tsv(
    ebd_in_WG,
    "ebd_cuce_WG_zf.txt")


### generating detection frequencies ###

# Group by month and species, calculate detections and total checklists
monthly_freq <- ebd_in_WG %>%
    group_by(month) %>%
    summarize(
        total_detections = sum(species_observed),
        total_checklists = n(),
        frequency = total_detections / total_checklists
    ) %>%
    ungroup()
monthly_freq$perct <- monthly_freq$frequency*100
monthly_freq <- na.omit(monthly_freq)    
    
# Plot the frequency of detections for each species
p<- ggplot(monthly_freq, aes(x = month, y = perct)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = seq(1, 12, by = 1), limits = c(1, 12)) +
    labs(
        title = "Frequency of Detections in %",
        x = "Month",
        y = "Frequency of Detections",
        color = "Species"
    )
plot(p)
ggsave(p, filename = paste0("cuce-detection", ".png"),
       width = 8, height = 5, units = "in", dpi = 300)

#### violin density plots for cuce along elevation according to seasons

f_data <- read.delim("all-data-merged.txt", sep = "\t")
#remove other species
cuce_data <- f_data %>% filter(COMMON.NAME == "Gray-headed Canary-Flycatcher")



#### melting the above data into a format useful for making the desired voilin plot####
# I want a 2 column dataframe with one column  - season and the other
# with the elevation where the species was found, so that the plot is
# weighted according to the number individuals detected at those elevaions 

# to do this we will repeat the entries for each species at a particular location
# according to the number of individuals detected



cuce_data <- cuce_data %>%
    mutate(season = case_when(
        month %in% 3:6 ~ "summer",
        month %in% 7:10 ~ "monsoon",
        month %in% c(11, 12, 1, 2) ~ "winter",
        TRUE ~ NA_character_  # Default case, though it shouldn't occur
    ))
cuce_data <- cuce_data %>%
    mutate(OBSERVATION.COUNT = case_when(
        OBSERVATION.COUNT == "X" ~ "1",  # Replace 'X' with 1
        TRUE ~ OBSERVATION.COUNT          # Keep other values as they are
    )) 
cuce_data$OBSERVATION.COUNT <- as.numeric(cuce_data$OBSERVATION.COUNT)

plot_data <- data.frame(
    season=character(),
    elevation=numeric()
)
for (i in 1:nrow(cuce_data)) {
    print(i)
    for (j in 1:cuce_data$OBSERVATION.COUNT[i]) {
        print(j)
        row <- data.frame(season = cuce_data$season[i], elevation = cuce_data$elevation[i])
        plot_data <- rbind(plot_data, row)
    }
}

save(plot_data, file = "focus data for stats.Rdata")
load("focus data for stats.Rdata")

sample_size = plot_data %>% group_by(season) %>% summarize(num=n())


# Plot

p <- plot_data %>%
    left_join(sample_size) %>%
    mutate(myaxis = paste0(season, "\n", "n=", num)) %>%
    ggplot( aes(x=myaxis, y=elevation, fill=season)) +
    geom_violin(width=1) +
    geom_violin(width=1, fill="transparent", draw_quantiles = c(0.05, 0.95), color = "magenta", linewidth = 0.5)+
    geom_violin(width=1, fill="transparent", draw_quantiles = c(0.5), color = "grey", linewidth = 0.5)+
    scale_fill_viridis(discrete = TRUE) +
    theme_ipsum() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 14),
          axis.text.y = element_text(angle = 90, hjust = 0.5, size = 18),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          panel.grid.major.x = element_line(color = "gray")) +
    xlab(element_text("Season", size = 20)) +
    ylab(element_text("Elevation (m)"))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10))


ggsave(p, filename = "seasonal_cuce_violin.jpeg",
       width = 10, height = 8, units = "in", dpi = 300)


