# script to pull SST data from NOAA OISST netcdf files

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data_dir <- "../../data/sst/"

library(tidyverse)
source("NOAA_OISST_ncdf4.R")

# define years to pull data for
years <- 2021:2024

# function to pull data from one year
pull_one_year <- function(year){
  # make name for data file
  filename <- paste0("sst.day.mean.", year, ".nc")
  this_file <- extractOISSTdaily(fname = paste0(data_dir, filename),
                                 lsmask = paste0(data_dir, "lsmask.oisst.nc"),
                                 lonW = (-69.324207 + 360) %% 360,
                                 lonE = (-69.304127+ 360) %% 360,
                                 latS = 43.859421,
                                 latN = 43.873913,
                                 date1 = paste0(year,"-1-01"),
                                 date2 = paste0(year, "-12-31"),
                                 varname = "sst")[1, 1, ]
  # make data into tidy df
  this_file_df <- data.frame(date = names(this_file), sst = this_file) %>%
    mutate(date = as.Date(date),
           sst = as.numeric(sst))
  rownames(this_file_df) <- NULL
  
  return(this_file_df)
}

all_years_data <- bind_rows(lapply(years, pull_one_year))

# add columns for day, month, year
all_years_data_clean <- all_years_data %>%
  mutate(day = day(date),
         month = month(date),
         year = year(date))

# write temperature data to CSV
write_csv(all_years_data_clean, paste0(data_dir, "sst.csv"))

# visualize data
ggplot(data = all_years_data, aes(x = date, y = anom)) +
  geom_point() +
  geom_line()







