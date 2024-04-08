# get Raw Warehouse Data from dwasl.SALMON_BIO_FACT (i.e., MTAL data, including scale ages!!!)
# however, note that Eric L. added `fw_age` and `sw_age` to AKFINADM.V_SALMON_BIO_FACT_GEN_TISS as of 2024-01-05
# Kyle Shedd
# 2024-01-05; updated 2024-04-08

library(tidyverse)
library(lubridate)
library(GCLr)

rm(list = ls())

# Begin user input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.username = readLines("~/R/usr_pw.txt", n = 1)  # LOKI username
.password = readLines("~/R/usr_pw.txt" , n = 2)[[2]]  # LOKI password

# End user input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get from Database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start.time <- Sys.time() 

options(java.parameters = "-Xmx10g")

url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud

drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")

drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")

con <- RJDBC::dbConnect(drv, url = url, user = .username, password = .password)

# data_qry <- paste0("SELECT * FROM dwasl.SALMON_BIO_FACT@dwprod_readonly T WHERE T.BATCH_NUMBER LIKE 'HWI%' AND T.CARD_NUMBER IS NOT NULL")  # get all HWI data, but only for fish with scale ages, need to modify `T.CARD_NUMBER IS NOT NULL` to get everything

data_qry <- paste0("SELECT * FROM dwasl.SALMON_BIO_FACT@dwprod_readonly T WHERE T.STREAM IN '111-50-10690' AND T.SAMPLE_YEAR IN '2023'")  # get all HWI data, but only for fish with scale ages, need to modify `T.CARD_NUMBER IS NOT NULL` to get everything

dataAll0 <- RJDBC::dbGetQuery(con, data_qry)   

discon <- RJDBC::dbDisconnect(con)  # disconnect from OceanAK

stop.time <- Sys.time()

fulltime <- stop.time - start.time

print(fulltime)

# Inspect Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

glimpse(dataAll0)   # what does out data look like?

# is 2023 data in there?
dataAll0 %>% 
  dplyr::count(STREAM, SAMPLE_YEAR) %>% 
  tidyr::pivot_wider(names_from = STREAM, values_from = n)

# quick count of samples from 2023, should just be Fish Creek
dataAll0 %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(sample_date = lubridate::as_date(sample_date)) %>% 
  dplyr::count(sample_year, species, gcl_location_code) %>% 
  tidyr::pivot_wider(names_from = sample_year, values_from = n)

# which streams have age data?
dataAll0 %>% 
  janitor::clean_names() %>% 
  dplyr::filter(!is.na(sw_age)) %>% 
  dplyr::count(sample_year, gcl_location_code, sw_age)

# ditch NA columns
data_noNA <- dataAll0 %>% 
  dplyr::select_if(~any(!is.na(.)))

glimpse(data_noNA)

# fix NA locations, write out SEAK fitness stream scale age + length data for Julia McMahon (UAF)
data_noNA %>% 
  dplyr::mutate(GCL_LOCATION_CODE = dplyr::case_when(is.na(GCL_LOCATION_CODE) & STREAM == "111-41-10050" ~ "Admiralty Creek",
                                                     is.na(GCL_LOCATION_CODE) & STREAM == "111-50-10690" ~ "Fish Creek - Douglas Island",
                                                     is.na(GCL_LOCATION_CODE) & STREAM == "111-33-10100" ~ "Prospect Creek",
                                                     is.na(GCL_LOCATION_CODE) & STREAM == "115-20-10520" ~ "Sawmill Creek",
                                                     TRUE ~ GCL_LOCATION_CODE)) %>% 
  readr::write_csv("../Stream Specimens/StreamSpecimens_SEAK_2013-2022_for_JuliaMcMahon.csv")

# end