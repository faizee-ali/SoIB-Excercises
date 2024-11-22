library(auk)
library(tidyverse)

#testing
cuce <- read.delim("data/ebd_gyhcaf/ebd_gyhcaf1_smp_relSep-2024.txt")
ebd_top<- read_tsv("ebd_IN_unv_smp_relApr-2024_sampling.txt", n_max = 50)

#### start ####
auk_set_ebd_path("data/ebd_yewlap/ebd_yewlap2_smp_relSep-2024.txt")
auk_ebd("data/ebd_gyhcaf/ebd_gyhcaf1_smp_relSep-2024.txt")

## filters ##

#filtering ebd and sampling file together
ebd_filtered <- auk_ebd("data/ebd_yewlap/ebd_yewlap2_smp_relSep-2024.txt", file_sampling = "data/ebd_sampling_relSep-2024/ebd_sampling_relSep-2024.txt") %>% 
    auk_protocol(c("Traveling", "Stationary")) %>% 
    auk_complete() %>% 
    auk_duration(duration = c(1, 120)) %>%
    auk_distance(distance = c(0, 2.5)) %>%
    auk_filter(file = "yama/f_ebd.txt", file_sampling = "yama/f_sed.txt")

####keeping only the relevant sampling data####
cuce <- read_tsv("f_ebd.txt")

#find the extent
bounding_box <- cuce %>%
    summarise(
        xmin = min(LONGITUDE, na.rm = TRUE),
        xmax = max(LONGITUDE, na.rm = TRUE),
        ymin = min(LATITUDE, na.rm = TRUE),
        ymax = max(LATITUDE, na.rm = TRUE)
    )

# Apply the bounding box filter

sampling <- auk_sampling(file = "f_sed.txt")
small_sampling <- auk_bbox(sampling, bbox = c(bounding_box$xmin, bounding_box$ymin, bounding_box$xmax, bounding_box$ymax)) %>%
    auk_filter(file = "ebd_sample_focused.txt", overwrite = TRUE)

ebd <- auk_ebd(file = "f_ebd.txt")
small_ebd <- auk_bbox(ebd, bbox = c(bounding_box$xmin, bounding_box$ymin, bounding_box$xmax, bounding_box$ymax)) %>%
    auk_filter(file = "ebd_verify_focused.txt")


#### zero filling ####

ebd_only <- read_ebd("ebd_verify_focused.txt")
sed_only <- read_sampling("ebd_sample_focused.txt")

ebd_zf <- auk_zerofill(ebd_only, sampling_events = sed_only)
ebd_zf

ebd_zf_df <- collapse_zerofill(ebd_zf)

write_tsv(ebd_zf_df, file = "ebd_zf_sampl.txt")

#winter data
ebd_zf_df <- auk_ebd("ebd_verify_focused.txt", file_sampling ="ebd_sample_focused.txt") %>% 
    auk_date(c("*-11-01", "*-02-28")) %>%
    auk_filter(file = "winter/ebd_cuce_winter.txt", file_sampling = "winter/sampling_cuce_winter.txt")%>% 
    auk_zerofill(collapse = TRUE)

#tidy-up
arrange(ebd_zf_df, desc(observation_count)) %>% 
    select(checklist_id, observation_count) %>% 
    head(10)

zf_count <- ebd_zf_df %>% 
    mutate(observation_count = if_else(observation_count == "X", 
                                       NA_character_, observation_count),
           observation_count = as.integer(observation_count),
           effort_distance_km = if_else(protocol_type == "Stationary", 
                                        0, effort_distance_km))
write_csv(zf_count, "winter/ebird_zf.csv")
rm(ebd_zf_df, zf_count)
###summer data
ebd_zf_df <- auk_ebd("ebd_verify_focused.txt", file_sampling ="ebd_sample_focused.txt") %>% 
    auk_date(c("*-04-01", "*-07-28")) %>%
    auk_filter(file = "summer/ebd_cuce_winter.txt", file_sampling = "summer/sampling_cuce_winter.txt")%>% 
    auk_zerofill(collapse = TRUE)

#tidy-up
arrange(ebd_zf_df, desc(observation_count)) %>% 
    select(checklist_id, observation_count) %>% 
    head(10)

zf_count <- ebd_zf_df %>% 
    mutate(observation_count = if_else(observation_count == "X", 
                                       NA_character_, observation_count),
           observation_count = as.integer(observation_count),
           effort_distance_km = if_else(protocol_type == "Stationary", 
                                        0, effort_distance_km))
write_csv(zf_count, "summer/ebird_zf.csv")
