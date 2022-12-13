library(tidyverse)
library(lubridate)

rm(list = ls())

# Begin user input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

username = "krshedd"  # LOKI username

password = ""  # LOKI passowrd

# Currently set up to read all pink data from all years
species = "CM"  # P = pink, CM = chum

streams <- c("ADMCR", "PROSCR", "SAWCR", "FISHCR")  # these are the back half of the silly codes

yrs <- 13:22  # 2 digit years

sillyvec <- paste0(species, rep(streams, each = length(yrs)), yrs)  # put it all together to get all possible silly codes

not_tissues = "Otolith"  # OceanAK has 1 row of data per tissue, including both otolith and genetic tissues, we do not want the otolith stuff (this is different than otolith reads)

source("~/../R/Functions.GCL.R")  # set to your own path as necessary

# End user input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start_time <- proc.time()

while(!require(RJDBC)){install.packages("RJDBC")}

if(!file.exists(path.expand("~/R"))){
  
  dir<-path.expand("~/R")
  
  dir.create(dir)
  
  bool <- file.copy(from="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to=path.expand("~/R/ojdbc8.jar"))
  
} else {
  
  if(!file.exists(path.expand("~/R/ojdbc8.jar"))){
    
    bool <- file.copy(from="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to=path.expand("~/R/ojdbc8.jar"))
    
  }
  
}

start.time <- Sys.time() 

options(java.parameters = "-Xmx10g")

if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
  
  drv <- JDBC("oracle.jdbc.OracleDriver",classPath="C:/Program Files/R/RequiredLibraries/ojdbc8.jar"," ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
  
} else {
  
  drv <- JDBC("oracle.jdbc.OracleDriver",classPath=path.expand("~/R/ojdbc8.jar")," ")
  
}

url <-LOKI_URL.GCL()

con <- dbConnect(drv,url=url,user=username,password=password)

data_qry <-
  paste(
    "SELECT * FROM AKFINADM.V_SALMON_BIO_FACT_GEN_TISS WHERE SILLY_CODE IN (",
    paste0("'", sillyvec, "'", collapse = ","),
    ") AND TISSUE_TYPE NOT IN (",
    paste0("'", not_tissues, "'", collapse = ","),
    ")",
    sep = ""
  )

dataAll0 <- dbGetQuery(con, data_qry)   

discon <- dbDisconnect(con)  # disconnect from OceanAK

proc.time() - start_time  # how long did it take?

glimpse(dataAll0)   # what does out data look like?

dataAll0 %>% 
  count(SILLY_CODE, TISSUE_TYPE)  # sample size by silly code and tissue type, should only be hearts and occassionally other random tissues


# Write out the data to "V:\Analysis\5_Coastwide\Multispecies\Alaska Hatchery Research Program\SEAK Chum\OceanAK" with timestamp
write_csv(
  x = dataAll0,
  file = paste0(
    "../OceanAK/AHRP Salmon Biological Data ",
    Sys.time() %>% str_remove_all(pattern = "-") %>% str_remove_all(pattern = ":") %>% str_replace(pattern = " ", replacement = "_"),
    ".csv"
  )
)
