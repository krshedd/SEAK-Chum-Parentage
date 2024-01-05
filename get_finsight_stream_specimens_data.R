# get Finsight StreamSpecimens data from HWI_STREAM_RAW, Tim Frawley's copy since HatcheryWild.org is down
# Kyle Shedd
# 2024-01-05

library(tidyverse)
library(lubridate)
library(GCLr)

rm(list = ls())

# Begin user input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.username = readLines("~/R/usr_pw.txt", n = 1)  # LOKI username
.password = readLines("~/R/usr_pw.txt" , n = 2)[[2]]  # LOKI password

# Currently set up to read all pink data from all years
# region = "SE" # SE or PWS

# End user input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get from Database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start.time <- Sys.time() 

options(java.parameters = "-Xmx10g")

url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud

drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")

drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")

con <- RJDBC::dbConnect(drv, url = url, user = .username, password = .password)

data_qry <- paste0("SELECT * FROM dwasl.HWI_STREAM_RAW@dwprod_readonly")  # get everything!

# data_qry <- paste0("SELECT * FROM dwasl.HWI_STREAM_LOAD2023@dwprod_readonly")  # get from a specific year

dataAll0 <- RJDBC::dbGetQuery(con, data_qry)   

discon <- RJDBC::dbDisconnect(con)  # disconnect from OceanAK

stop.time <- Sys.time()

fulltime <- stop.time - start.time

print(fulltime)

# Inspect Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

glimpse(dataAll0)   # what does out data look like?

# quick count of samples from 2023, should just be Fish Creek
dataAll0 %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(sample_date = lubridate::as_date(survey_date),
                sample_year = lubridate::year(sample_date)) %>% 
  dplyr::filter(sample_year == 2023) %>% 
  dplyr::count(stream_name)

# filter for all SEAK Chum fitness stream data
dataAll0 %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    sample_date = lubridate::as_date(survey_date),
    sample_year = lubridate::year(sample_date)
  ) %>%
  dplyr::filter(
    region == "SE",
    species == 450,
    grepl(x = stream_name, pattern = "Fish") |
      grepl(x = stream_name, pattern = "Admiralty") |
      grepl(x = stream_name, pattern = "Sawmill") |
      grepl(x = stream_name, pattern = "Prospect")
  ) %>%
  dplyr::count(stream_name, sample_year) %>%
  tidyr::pivot_wider(names_from = sample_year, values_from = n)

# save all SEAK Chum fitness stream data with `janitor::clean_names()`
dataAll0 %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    sample_date = lubridate::as_date(survey_date),
    sample_year = lubridate::year(sample_date)
  ) %>%
  dplyr::filter(
    region == "SE",
    species == 450,
    grepl(x = stream_name, pattern = "Fish") |
      grepl(x = stream_name, pattern = "Admiralty") |
      grepl(x = stream_name, pattern = "Sawmill") |
      grepl(x = stream_name, pattern = "Prospect")
  ) %>% 
  readr::write_csv(file = "../Stream Specimens/finsight_stream_specimens_SEAK_chum_fitness_2013-2023.csv")
# end