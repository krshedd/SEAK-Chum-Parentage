# Join OceanAK and Finsight Stream-Specimens Data
# Jodi wants to know sample date and spawning state data to help troubleshoot genotyping failures
# Kyle Shedd
# 2024-05-09

library(tidyverse)
oceanak <- readr::read_csv(file = "../OceanAK/AHRP Salmon Biological Data 20240105_122404.368459.csv") %>% 
  janitor::clean_names()
finsight <- readr::read_csv(file = "../Stream Specimens/finsight_stream_specimens_SEAK_chum_fitness_2013-2023.csv")

finsight %>% 
  dplyr::rename(dna_tray_code = sample_tray_id,
                dna_tray_well_code = sample_cell) %>% 
  dplyr::mutate(dna_tray_well_code = as.numeric(dna_tray_well_code)) %>% 
  dplyr::select(dna_tray_code, dna_tray_well_code, spawning_state) %>% 
  dplyr::right_join(oceanak, by = c("dna_tray_code", "dna_tray_well_code")) %>% 
  dplyr::mutate(sample_date = lubridate::as_date(sample_date)) %>% 
  dplyr::select(collection_id:fish_id, spawning_state, sample_date:tissue_type, dna_tray_code, dna_tray_well_code, location_code:sw_age) %>% 
  readr::write_csv(file = "../OceanAK/AHRP Salmon Biological Data 20240105_122404.368459_with_spawning_state.csv")

# End