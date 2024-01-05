# get "OceanAK" data from Eric Lardizabal's genetics view of the Salmon Biological Fact Warehouse, includes scale ages
# Kyle Shedd
# updated 2024-01-05

library(tidyverse)
library(lubridate)

rm(list = ls())

# Begin user input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.username = readLines("~/R/usr_pw.txt", n = 1)  # LOKI username
.password = readLines("~/R/usr_pw.txt" , n = 2)[[2]]  # LOKI password

# Currently set up to read all pink data from all years
species = "CM"  # P = pink, CM = chum

streams <- c("ADMCR", "PROSCR", "SAWCR", "FISHCR")  # these are the back half of the silly codes

yrs <- 13:23  # 2 digit years

sillyvec <- paste0(species, rep(streams, each = length(yrs)), yrs)  # put it all together to get all possible silly codes

(sillyvec <-
    c(sillyvec, paste0(
      species, rep(streams[-1], each = length(21:23)), "T", 21:23
    ), "CMFISHCRT13")  # add in the Tagged collections for FISH, PROSCR and SAWCR 2021 and 2022 + Fish Creek 2013 (but doesn't appear to be included here)
)

not_tissues = "Otolith"  # OceanAK has 1 row of data per tissue, including both otolith and genetic tissues, we do not want the otolith stuff (this is different than otolith reads)

# End user input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start.time <- Sys.time() 

options(java.parameters = "-Xmx10g")

url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud

drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")

drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")

con <- RJDBC::dbConnect(drv, url = url, user = .username, password = .password)

data_qry <-
  paste(
    "SELECT * FROM AKFINADM.V_SALMON_BIO_FACT_GEN_TISS WHERE SILLY_CODE IN (",
    paste0("'", sillyvec, "'", collapse = ","),
    ") AND TISSUE_TYPE NOT IN (",
    paste0("'", not_tissues, "'", collapse = ","),
    ")",
    sep = ""
  )

dataAll0 <- RJDBC::dbGetQuery(con, data_qry)   

discon <- RJDBC::dbDisconnect(con)  # disconnect from OceanAK

stop.time <- Sys.time()

fulltime <- stop.time - start.time

print(fulltime)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dplyr::glimpse(dataAll0)   # what does out data look like?

dataAll0 %>% 
  dplyr::count(SILLY_CODE, TISSUE_TYPE)  # sample size by silly code and tissue type, should only be hearts and occassionally other random tissues

# sample size by stream by year
dataAll0 %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(sample_date = lubridate::as_date(sample_date),
                sample_year = lubridate::year(sample_date)) %>% 
  dplyr::count(sample_year, location_code) %>% 
  tidyr::pivot_wider(names_from = sample_year, values_from = n)


# Write out the data to "V:\Analysis\5_Coastwide\Multispecies\Alaska Hatchery Research Program\SEAK Chum\OceanAK" with timestamp
readr::write_csv(
  x = dataAll0,
  file = paste0(
    "../OceanAK/AHRP Salmon Biological Data ",
    Sys.time() %>% stringr::str_remove_all(pattern = "-") %>% stringr::str_remove_all(pattern = ":") %>% stringr::str_replace(pattern = " ", replacement = "_"),
    ".csv"
  )
)

# End