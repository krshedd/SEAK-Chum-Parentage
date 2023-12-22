# Add existing extractions from 2013 to gentoyping list
# Kyle Shedd
# 2023-12-22

# eCM001 only had 2014-2022, no including 2013 (already extracted) or 2023 (no otolith reads yet)

library(tidyverse)
.username = readLines("~/R/usr_pw.txt", n = 1)
.password = readLines("~/R/usr_pw.txt" , n = 2)[[2]]

(
  extractions_2013 <- GCLr::get_extraction_info(
    sillyvec = c(
      "CMADMCR13",
      "CMPROSCR13",
      "CMSAWCR13",
      "CMFISHCR13",
      "CMFISHCRT13"
    ),
    username = .username,
    password = .password,
    file = "../OceanAK/extraction_wells_for_genotyping_2013.csv"
  ) %>% 
    janitor::clean_names() %>% 
    dplyr::rename(plate_id = fk_plate_id,
                  tissue_type = tissuetype)
)

# see what we have
extractions_2013 %>% 
  dplyr::group_by(plate_id, silly_code) %>% 
  dplyr::summarise(min_fish_no = min(fish_no),
                   max_fish_no = max(fish_no)) %>% 
  dplyr::arrange(plate_id, dplyr::desc(max_fish_no)) %>% 
  print(n = Inf)

# drop extra, cherry-picked plates
extractions_2013 %>% 
  dplyr::filter(plate_id <= 56364) %>% 
  dplyr::distinct(plate_id, silly_code) %>% 
  tidyr::nest(silly_code = silly_code) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(silly_code = toString(silly_code)) %>% 
  dplyr::mutate(silly_code = stringr::str_remove_all(string = silly_code, pattern = "\"")) %>% 
  dplyr::mutate(silly_code = stringr::str_remove_all(string = silly_code, pattern = "\\)")) %>% 
  dplyr::mutate(silly_code = stringr::str_remove_all(string = silly_code, pattern = "\\(")) %>% 
  dplyr::mutate(silly_code = stringr::str_remove_all(string = silly_code, pattern = "c")) %>% 
  dplyr::mutate(silly_code = stringr::str_replace_all(string = silly_code, pattern = ", ", replacement = "/")) %>% 
  readr::write_csv(file = "../OceanAK/extraction_plates_for_genotyping_2013.csv")

extractions_2013 %>% 
  dplyr::filter(plate_id == 56364)

# end