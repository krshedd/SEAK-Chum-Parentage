# OceanAK data pull
# Kyle Shedd
# Tue Mar 05 15:08:00 2019

rm(list = ls())
date()
library(tidyverse)
library(lubridate)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### OceanAK ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# List of silly's for OceanAK filter
streams <- c("ADMCR", "PROSCR", "SAWCR", "FISHCR")
yrs <- c(13:14, 17:18)
writeClipboard(paste(paste0("CM", rep(streams, each = 5), yrs), collapse = ";"))

oceanak <- read_csv(file = "../OceanAK/PedigreeData_AHRP - Salmon Biological Data 2_SEAK_2013-2018_no_otoliths.csv") %>% 
  unite(SillySource, `Silly Code`, `Fish ID`, sep = "_", remove = FALSE) %>% 
  unite(TrayCodeID, `DNA Tray Code`, `DNA Tray Well Code`, sep = "_", remove = FALSE)

# dups <- oceanak %>% 
#   group_by(SillySource) %>% 
#   summarise(n = n()) %>% 
#   filter(n > 1) %>% 
#   arrange(desc(n)) %>% 
#   left_join(oceanak)
# nrow(dups)
# table(dups$`Location Code`, dups$`Sample Year`)
# View(dups)
# length(unique(dups$`Sample ID`))
# 
# dups_tray <- oceanak %>% 
#   group_by(TrayCodeID) %>% 
#   summarise(n = n()) %>% 
#   filter(n > 1) %>% 
#   arrange(desc(n)) %>% 
#   left_join(oceanak)
# nrow(dups_tray)
# table(dups_tray$`Location Code`, dups_tray$`Sample Year`)
# View(dups_tray)
# 
# write_csv(x = dups, path = "OceanAK/AHRP - Salmon Biological Data 2_PWS_2013-2017_duplicates.csv")

# samples per stream per year
addmargins(table(oceanak$`Location Code`, oceanak$`Sample Year`))

# otolith code per stream per year
table(oceanak$`Location Code`, oceanak$`Otolith Mark Present`, oceanak$`Sample Year`, useNA = "always")
table(oceanak$`Location Code`, oceanak$`Otolith Mark Status Code`, oceanak$`Sample Year`, useNA = "always")

# add otolith_read logical, stream as factor, and year
oceanak_mod <- oceanak %>% 
  mutate(otolith_read = !is.na(`Otolith Mark Status Code`)) %>% 
  mutate(stream = factor(x = `Location Code`, levels = c("Admiralty Creek", "Fish Creek - Douglas Island", "Prospect Creek", "Sawmill Creek"))) %>% 
  rename(year = `Sample Year`) %>% 
  mutate(origin = case_when(`Otolith Mark Present` == "NO" ~ "natural",
                            `Otolith Mark Present` == "YES" ~ "hatchery")) %>% 
  mutate(origin = factor(origin, levels = c("natural", "hatchery")))

# table of stream, year, and otolith_read
table(oceanak_mod$stream, oceanak_mod$otolith_read, oceanak_mod$year)

# table of samples per stream per year
addmargins(table(oceanak_mod$stream, oceanak_mod$year))

# table of otoliths read by stream and year
oceanak_mod %>%
  filter(otolith_read == TRUE) %>% 
  group_by(stream, year) %>% 
  summarise(freq = n()) %>% 
  spread(year, freq)

# table of otoliths read by stream and year and sex
oceanak_mod %>%
  filter(otolith_read == TRUE & Sex != "U" & !is.na(origin)) %>% 
  mutate(sex = case_when(Sex == "F" ~ "female",
                         Sex == "M" ~ "male")) %>% 
  unite(origin_sex, sex, origin) %>% 
  mutate(origin_sex = factor(x = origin_sex, levels = c("female_natural", "female_hatchery", "male_natural", "male_hatchery"))) %>% 
  group_by(stream, year, origin_sex) %>% 
  summarise(freq = n()) %>% 
  spread(origin_sex, freq) %>% 
  arrange(year)

# table of otolith mark read by stream and year
oceanak_mod %>%
  filter(`Otolith Mark Present` == "NO") %>%  # natural-origin fish only
  group_by(stream, year) %>% 
  summarise(freq = n()) %>% 
  spread(year, freq)

# min and max date within a year
oceanak_mod %>% 
  group_by(year) %>% 
  summarise(begin_date = min(`Sample Date`), end_date = max(`Sample Date`))

# histogram of samples per date per year
oceanak_mod %>% 
  mutate(julian_date = yday(`Sample Date`)) %>% 
  ggplot(aes(x = julian_date)) +
  geom_histogram() +
  facet_grid(year ~ .)

# histogram of samples per date per year per stream
oceanak_mod %>% 
  mutate(julian_date = yday(`Sample Date`)) %>% 
  mutate(stream = case_when(stream == "Fish Creek - Douglas Island" ~ "Fish Creek",
                            TRUE ~ as.character(stream))) %>% 
  filter(!is.na(origin) & Sex != "U") %>% 
  group_by(year, stream, origin, julian_date) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = julian_date, y = n, fill = origin)) +
  geom_col() +
  ylim(0, 345) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(year ~ stream) +
  labs(fill = "Origin") +
  ylab("Number of Samples") +
  xlab("Day of Year") +
  theme(text = element_text(size = 20)) 
  # ggtitle("AHRP SEAK Chum - number of samples") 

# blank
oceanak_mod %>% 
  mutate(julian_date = yday(`Sample Date`)) %>% 
  mutate(stream = case_when(stream == "Fish Creek - Douglas Island" ~ "Fish Creek",
                            TRUE ~ as.character(stream))) %>% 
  filter(!is.na(origin) & Sex != "U") %>% 
  group_by(year, stream, origin, julian_date) %>% 
  summarise(n = n()) %>% 
  mutate(n = case_when(stream == "Fish Creek" & year == 2013 ~ as.double(n),
                       TRUE ~ 0)) %>% 
  ggplot(aes(x = julian_date, y = n, fill = origin)) +
  geom_col() +
  # geom_histogram(binwidth = 1) +
  ylim(0, 345) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(year ~ stream) +
  labs(fill = "Origin") +
  ylab("Number of Samples") +
  xlab("Day of Year") +
  theme(text = element_text(size = 20)) 


# just fish creek 2013
oceanak_mod %>% 
  mutate(julian_date = yday(`Sample Date`)) %>% 
  mutate(stream = case_when(stream == "Fish Creek - Douglas Island" ~ "Fish Creek",
                            TRUE ~ as.character(stream))) %>% 
  filter(!is.na(origin) & Sex != "U") %>% 
  filter(stream == "Fish Creek" & year == 2013) %>% 
  ggplot(aes(x = julian_date, fill = origin)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylim(0, 345) +
  xlim(201, 242) +
  facet_grid(year ~ stream) +
  labs(fill = "Origin") +
  ylab("Number of Samples") +
  xlab("Day of Year") +
  theme(text = element_text(size = 20)) #+
ggtitle("AHRP SEAK Chum - number of samples") 


yday(Sys.Date())  # today's Julian date

# get ggplot rgb
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col2rgb(gg_color_hue(2))
