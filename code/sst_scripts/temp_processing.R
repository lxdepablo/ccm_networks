# temp_processing.R
# Henry Li
# last edit: April 10, 2023
# adapted from Dana Morton's code
# used to generate daily mean SST data for two Sorte study sites, to generate temperature anomaly figure for GRFP proposal

###############################################################################
### script for extracting daily mean SST from NOAA daily mean SST files. 
### uses functions from Luke Miller's R script : NOAA_OISST_ncdf4.R 
### see for info and to download required functions: https://lukemiller.org/index.php/2014/11/extracting-noaa-sea-surface-temperatures-with-ncdf4/

# load packages
library(ncdf4)
library(dplyr)
library(heatwaveR)
library(ggplot2)

# set working directory to data location and source Luke Miller's R functions
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./NOAA_OISST_ncdf4.R")
data_dir <- "../../GOM community data local/temperature_data/NOAA_data/"
setwd(data_dir)

# load csv of coordinates to query
coords <- read.csv("../coordinates.csv", header=T)

# load all temperature data files ending in .nc
fileNames <- Sys.glob("./sst_day_mean_data/*.nc") # CHANGE FOR DIFFERENT SITE FOLDERS

#convert longitudes to 360 format: long360 <- (lon180 + 360) %% 360 
coords$lon.360 <- (coords$OceanLon + 360) %% 360 

# create date vector to use in final output table
dates <- seq(as.Date("2010-01-01"), as.Date("2017-12-31"), by = "days") # CHANGE FOR DIFFERENT SITE FOLDERS
dates <- dates[!(format(dates, "%m") == "02" & format(dates, "%d") == "29")]

###############################################################################
# Canoe Beach
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[1,6],lonE= coords[1,6], latS=coords[1,4], latN=coords[1,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[1,6],lonE= coords[1,6], latS=coords[1,4], latN=coords[1,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
canoe_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(canoe_df) <- c("t", "temp")
# write data to .csv file
write.table(canoe_df, "../SST_data/Canoe.Beach.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

canoe_df$t <- as.Date(canoe_df$t)
# process marine heat wave event data
canoe_hot_ts <- ts2clm(canoe_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
canoe_mhw <- detect_event(canoe_hot_ts)
# process marine cold spell event data
canoe_cold_ts <- ts2clm(canoe_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
canoe_mcs <- detect_event(canoe_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(canoe_mhw$event, "../event_data/Canoe.Beach.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(canoe_mcs$event, "../event_data/Canoe.Beach.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Cape Newagen
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[2,6],lonE= coords[2,6], latS=coords[2,4], latN=coords[2,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[2,6],lonE= coords[2,6], latS=coords[2,4], latN=coords[2,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
cape_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(cape_df) <- c("t", "temp")
# write data to .csv file
write.table(cape_df, "../SST_data/Cape.Newagen.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

cape_df$t <- as.Date(cape_df$t)
# process marine heat wave event data
cape_hot_ts <- ts2clm(cape_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
cape_mhw <- detect_event(cape_hot_ts)
# process marine cold spell event data
cape_cold_ts <- ts2clm(cape_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
cape_mcs <- detect_event(cape_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(cape_mhw$event, "../event_data/Cape.Newagen.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(cape_mcs$event, "../event_data/Cape.Newagen.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Chamberlain
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[3,6],lonE= coords[3,6], latS=coords[3,4], latN=coords[3,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[3,6],lonE= coords[3,6], latS=coords[3,4], latN=coords[3,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
chamberlain_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(chamberlain_df) <- c("t", "temp")
# write data to .csv file
write.table(chamberlain_df, "../SST_data/Chamberlain.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

chamberlain_df$t <- as.Date(chamberlain_df$t)
# process marine heat wave event data
chamberlain_hot_ts <- ts2clm(chamberlain_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
chamberlain_mhw <- detect_event(chamberlain_hot_ts)
# process marine cold spell event data
chamberlain_cold_ts <- ts2clm(chamberlain_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
chamberlain_mcs <- detect_event(chamberlain_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(chamberlain_mhw$event, "../event_data/Chamberlain.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(chamberlain_mcs$event, "../event_data/Chamberlain.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Cunner Ledge
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[4,6],lonE= coords[4,6], latS=coords[4,4], latN=coords[4,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[4,6],lonE= coords[4,6], latS=coords[4,4], latN=coords[4,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
cunner_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(cunner_df) <- c("t", "temp")
# write data to .csv file
write.table(cunner_df, "../SST_data/Cunner.Ledge.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

cunner_df$t <- as.Date(cunner_df$t)
# process marine heat wave event data
cunner_hot_ts <- ts2clm(cunner_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
cunner_mhw <- detect_event(cunner_hot_ts)
# process marine cold spell event data
cunner_cold_ts <- ts2clm(cunner_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
cunner_mcs <- detect_event(cunner_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(cunner_mhw$event, "../event_data/Cunner.Ledge.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(cunner_mcs$event, "../event_data/Cunner.Ledge.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Cutler
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[5,6],lonE= coords[5,6], latS=coords[5,4], latN=coords[5,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[5,6],lonE= coords[5,6], latS=coords[5,4], latN=coords[5,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
cutler_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(cutler_df) <- c("t", "temp")
# write data to .csv file
write.table(cutler_df, "../SST_data/Cutler.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

cutler_df$t <- as.Date(cutler_df$t)
# process marine heat wave event data
cutler_hot_ts <- ts2clm(cutler_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
cutler_mhw <- detect_event(cutler_hot_ts)
# process marine cold spell event data
cutler_cold_ts <- ts2clm(cutler_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
cutler_mcs <- detect_event(cutler_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(cutler_mhw$event, "../event_data/Cutler.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(cutler_mcs$event, "../event_data/Cutler.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Dyer Cove
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[6,6],lonE= coords[6,6], latS=coords[6,4], latN=coords[6,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[6,6],lonE= coords[6,6], latS=coords[6,4], latN=coords[6,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
dyer_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(dyer_df) <- c("t", "temp")
# write data to .csv file
write.table(dyer_df, "../SST_data/Dyer.Cove.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

dyer_df$t <- as.Date(dyer_df$t)
# process marine heat wave event data
dyer_hot_ts <- ts2clm(dyer_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
dyer_mhw <- detect_event(dyer_hot_ts)
# process marine cold spell event data
dyer_cold_ts <- ts2clm(dyer_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
dyer_mcs <- detect_event(dyer_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(dyer_mhw$event, "../event_data/Dyer.Cove.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(dyer_mcs$event, "../event_data/Dyer.Cove.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Folly's Point
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[7,6],lonE= coords[7,6], latS=coords[7,4], latN=coords[7,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[7,6],lonE= coords[7,6], latS=coords[7,4], latN=coords[7,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
follys_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(follys_df) <- c("t", "temp")
# write data to .csv file
write.table(follys_df, "../SST_data/Follys.Point.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

follys_df$t <- as.Date(follys_df$t)
# process marine heat wave event data
follys_hot_ts <- ts2clm(follys_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
follys_mhw <- detect_event(follys_hot_ts)
# process marine cold spell event data
follys_cold_ts <- ts2clm(follys_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
follys_mcs <- detect_event(follys_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(follys_mhw$event, "../event_data/Follys.Point.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(follys_mcs$event, "../event_data/Follys.Point.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Grindstone Neck
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[8,6],lonE= coords[8,6], latS=coords[8,4], latN=coords[8,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[8,6],lonE= coords[8,6], latS=coords[8,4], latN=coords[8,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
grindstone_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(grindstone_df) <- c("t", "temp")
# write data to .csv file
write.table(grindstone_df, "../SST_data/Grindstone.Neck.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

grindstone_df$t <- as.Date(grindstone_df$t)
# process marine heat wave event data
grindstone_hot_ts <- ts2clm(grindstone_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
grindstone_mhw <- detect_event(grindstone_hot_ts)
# process marine cold spell event data
grindstone_cold_ts <- ts2clm(grindstone_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
grindstone_mcs <- detect_event(grindstone_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(grindstone_mhw$event, "../event_data/Grindstone.Neck.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(grindstone_mcs$event, "../event_data/Grindstone.Neck.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Hamilton Cove
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[9,6],lonE= coords[9,6], latS=coords[9,4], latN=coords[9,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[9,6],lonE= coords[9,6], latS=coords[9,4], latN=coords[9,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
hamilton_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(hamilton_df) <- c("t", "temp")
# write data to .csv file
write.table(hamilton_df, "../SST_data/Hamilton.Cove.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

hamilton_df$t <- as.Date(hamilton_df$t)
# process marine heat wave event data
hamilton_hot_ts <- ts2clm(hamilton_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
hamilton_mhw <- detect_event(hamilton_hot_ts)
# process marine cold spell event data
hamilton_cold_ts <- ts2clm(hamilton_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
hamilton_mcs <- detect_event(hamilton_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(hamilton_mhw$event, "../event_data/Hamilton.Cove.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(hamilton_mcs$event, "../event_data/Hamilton.Cove.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Horizon Beach
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[10,6],lonE= coords[10,6], latS=coords[10,4], latN=coords[10,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[10,6],lonE= coords[10,6], latS=coords[10,4], latN=coords[10,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
horizon_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(horizon_df) <- c("t", "temp")
# write data to .csv file
write.table(horizon_df, "../SST_data/Horizon.Beach.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

horizon_df$t <- as.Date(horizon_df$t)
# process marine heat wave event data
horizon_hot_ts <- ts2clm(horizon_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
horizon_mhw <- detect_event(horizon_hot_ts)
# process marine cold spell event data
horizon_cold_ts <- ts2clm(horizon_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
horizon_mcs <- detect_event(horizon_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(horizon_mhw$event, "../event_data/Horizon.Beach.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(horizon_mcs$event, "../event_data/Horizon.Beach.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Loblolly Point
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[11,6],lonE= coords[11,6], latS=coords[11,4], latN=coords[11,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[11,6],lonE= coords[11,6], latS=coords[11,4], latN=coords[11,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
loblolly_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(loblolly_df) <- c("t", "temp")
# write data to .csv file
write.table(loblolly_df, "../SST_data/Loblolly.Point.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

loblolly_df$t <- as.Date(loblolly_df$t)
# process marine heat wave event data
loblolly_hot_ts <- ts2clm(loblolly_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
loblolly_mhw <- detect_event(loblolly_hot_ts)
# process marine cold spell event data
loblolly_cold_ts <- ts2clm(loblolly_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
loblolly_mcs <- detect_event(loblolly_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(loblolly_mhw$event, "../event_data/Loblolly.Point.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(loblolly_mcs$event, "../event_data/Loblolly.Point.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Marshall Point
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[12,6],lonE= coords[12,6], latS=coords[12,4], latN=coords[12,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[12,6],lonE= coords[12,6], latS=coords[12,4], latN=coords[12,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
marshall_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(marshall_df) <- c("t", "temp")
# write data to .csv file
write.table(marshall_df, "../SST_data/Marshall.Point.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

marshall_df$t <- as.Date(marshall_df$t)
# process marine heat wave event data
marshall_hot_ts <- ts2clm(marshall_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
marshall_mhw <- detect_event(marshall_hot_ts)
# process marine cold spell event data
marshall_cold_ts <- ts2clm(marshall_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
marshall_mcs <- detect_event(marshall_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(marshall_mhw$event, "../event_data/Marshall.Point.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(marshall_mcs$event, "../event_data/Marshall.Point.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Pemaquid Point
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[13,6],lonE= coords[13,6], latS=coords[13,4], latN=coords[13,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[13,6],lonE= coords[13,6], latS=coords[13,4], latN=coords[13,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
pemaquid_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(pemaquid_df) <- c("t", "temp")
# write data to .csv file
write.table(pemaquid_df, "../SST_data/Pemaquid.Point.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

pemaquid_df$t <- as.Date(pemaquid_df$t)
# process marine heat wave event data
pemaquid_hot_ts <- ts2clm(pemaquid_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
pemaquid_mhw <- detect_event(pemaquid_hot_ts)
# process marine cold spell event data
pemaquid_cold_ts <- ts2clm(pemaquid_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
pemaquid_mcs <- detect_event(pemaquid_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(pemaquid_mhw$event, "../event_data/Pemaquid.Point.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(pemaquid_mcs$event, "../event_data/Pemaquid.Point.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Quoddy Head
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[14,6],lonE= coords[14,6], latS=coords[14,4], latN=coords[14,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[14,6],lonE= coords[14,6], latS=coords[14,4], latN=coords[14,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
quoddy_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(quoddy_df) <- c("t", "temp")
# write data to .csv file
write.table(quoddy_df, "../SST_data/Quoddy.Head.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

quoddy_df$t <- as.Date(quoddy_df$t)
# process marine heat wave event data
quoddy_hot_ts <- ts2clm(quoddy_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
quoddy_mhw <- detect_event(quoddy_hot_ts)
# process marine cold spell event data
quoddy_cold_ts <- ts2clm(quoddy_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
quoddy_mcs <- detect_event(quoddy_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(quoddy_mhw$event, "../event_data/Quoddy.Head.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(quoddy_mcs$event, "../event_data/Quoddy.Head.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Rock Harbor
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[15,6],lonE= coords[15,6], latS=coords[15,4], latN=coords[15,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[15,6],lonE= coords[15,6], latS=coords[15,4], latN=coords[15,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
rock_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(rock_df) <- c("t", "temp")
# write data to .csv file
write.table(rock_df, "../SST_data/Rock.Harbor.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

rock_df$t <- as.Date(rock_df$t)
# process marine heat wave event data
rock_hot_ts <- ts2clm(rock_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
rock_mhw <- detect_event(rock_hot_ts)
# process marine cold spell event data
rock_cold_ts <- ts2clm(rock_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
rock_mcs <- detect_event(rock_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(rock_mhw$event, "../event_data/Rock.Harbor.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(rock_mcs$event, "../event_data/Rock.Harbor.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Rye Beach
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[16,6],lonE= coords[16,6], latS=coords[16,4], latN=coords[16,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[16,6],lonE= coords[16,6], latS=coords[16,4], latN=coords[16,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
rye_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(rye_df) <- c("t", "temp")
# write data to .csv file
write.table(rye_df, "../SST_data/Rye.Beach.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

rye_df$t <- as.Date(rye_df$t)
# process marine heat wave event data
rye_hot_ts <- ts2clm(rye_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
rye_mhw <- detect_event(rye_hot_ts)
# process marine cold spell event data
rye_cold_ts <- ts2clm(rye_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
rye_mcs <- detect_event(rye_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(rye_mhw$event, "../event_data/Rye.Beach.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(rye_mcs$event, "../event_data/Rye.Beach.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Mt Desert Island
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[17,6],lonE= coords[17,6], latS=coords[17,4], latN=coords[17,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[17,6],lonE= coords[17,6], latS=coords[17,4], latN=coords[17,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
mdi_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(mdi_df) <- c("t", "temp")
# write data to .csv file
write.table(mdi_df, "../SST_data/Mt.Desert.Island.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

mdi_df$t <- as.Date(mdi_df$t)
# process marine heat wave event data
mdi_hot_ts <- ts2clm(mdi_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
mdi_mhw <- detect_event(mdi_hot_ts)
# process marine cold spell event data
mdi_cold_ts <- ts2clm(mdi_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
mdi_mcs <- detect_event(mdi_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(mdi_mhw$event, "../event_data/Mt.Desert.Island.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(mdi_mcs$event, "../event_data/Mt.Desert.Island.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# South Addison
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[18,6],lonE= coords[18,6], latS=coords[18,4], latN=coords[18,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[18,6],lonE= coords[18,6], latS=coords[18,4], latN=coords[18,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
south_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(south_df) <- c("t", "temp")
# write data to .csv file
write.table(south_df, "../SST_data/South.Addison.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

south_df$t <- as.Date(south_df$t)
# process marine heat wave event data
south_hot_ts <- ts2clm(south_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
south_mhw <- detect_event(south_hot_ts)
# process marine cold spell event data
south_cold_ts <- ts2clm(south_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
south_mcs <- detect_event(south_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(south_mhw$event, "../event_data/South.Addison.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(south_mcs$event, "../event_data/South.Addison.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# The Glades
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[19,6],lonE= coords[19,6], latS=coords[19,4], latN=coords[19,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[19,6],lonE= coords[19,6], latS=coords[19,4], latN=coords[19,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
glades_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(glades_df) <- c("t", "temp")
# write data to .csv file
write.table(glades_df, "../SST_data/The.Glades.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

glades_df$t <- as.Date(glades_df$t)
# process marine heat wave event data
glades_hot_ts <- ts2clm(glades_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
glades_mhw <- detect_event(glades_hot_ts)
# process marine cold spell event data
glades_cold_ts <- ts2clm(glades_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
glades_mcs <- detect_event(glades_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(glades_mhw$event, "../event_data/The.Glades.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(glades_mcs$event, "../event_data/The.Glades.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Vinalhaven Island
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 34, 37))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[20,6],lonE= coords[20,6], latS=coords[20,4], latN=coords[20,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[20,6],lonE= coords[20,6], latS=coords[20,4], latN=coords[20,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
vinalhaven_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(vinalhaven_df) <- c("t", "temp")
# write data to .csv file
write.table(vinalhaven_df, "../SST_data/Vinalhaven.Island.2010-2017.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

vinalhaven_df$t <- as.Date(vinalhaven_df$t)
# process marine heat wave event data
vinalhaven_hot_ts <- ts2clm(vinalhaven_df, climatologyPeriod = c("2010-01-01", "2017-12-31"))
vinalhaven_mhw <- detect_event(vinalhaven_hot_ts)
# process marine cold spell event data
vinalhaven_cold_ts <- ts2clm(vinalhaven_df, climatologyPeriod = c("2010-01-01", "2017-12-31"), pctile = 10)
vinalhaven_mcs <- detect_event(vinalhaven_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(vinalhaven_mhw$event, "../event_data/Vinalhaven.Island.2010-2017.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(vinalhaven_mcs$event, "../event_data/Vinalhaven.Island.2010-2017.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Appledore Island
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 38, 41))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[21,6],lonE= coords[21,6], latS=coords[21,4], latN=coords[21,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[21,6],lonE= coords[21,6], latS=coords[21,4], latN=coords[21,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
appledore_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(appledore_df) <- c("t", "temp")
# write data to .csv file
write.table(appledore_df, "../SST_data/Appledore.Island.1982-2022.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

appledore_df$t <- as.Date(appledore_df$t)
# process marine heat wave event data
appledore_hot_ts <- ts2clm(appledore_df, climatologyPeriod = c("1982-01-01", "2022-12-31"))
appledore_mhw <- detect_event(appledore_hot_ts)
# process marine cold spell event data
appledore_cold_ts <- ts2clm(appledore_df, climatologyPeriod = c("1982-01-01", "2022-12-31"), pctile = 10)
appledore_mcs <- detect_event(appledore_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(appledore_mhw$event, "../event_data/Appledore.Island.1982-2022.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(appledore_mcs$event, "../event_data/Appledore.Island.1982-2022.cold.events.csv", sep=",", row.names = F, col.names = T)

###############################################################################
# Allen Island
temps1 <- matrix(,nrow = length(fileNames), ncol=59)  #Jan - Feb 28 to account for leap years 
temps2 <- matrix(,nrow = length(fileNames), ncol=306) #March 1 - Dec 31

R= 1

for (i in 1:nrow(temps1)) {
  name.str <- toString(fileNames[i])
  year = as.numeric(substr(name.str, 40, 43))
  temps1[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[22,6],lonE= coords[22,6], latS=coords[22,4], latN=coords[22,4],
                                  date1=paste(year,"-1-01", sep=""),date2= paste(year,'-02-28',sep=""))
  temps2[R,] <- extractOISSTdaily(fileNames[i],"./lsmask.oisst.v2.nc",lonW=coords[22,6],lonE= coords[22,6], latS=coords[22,4], latN=coords[22,4],
                                  date1=paste(year,"-3-01", sep=""),date2= paste(year,'-12-31',sep=""))
  
  R=R+1}  

all.temps <- cbind(temps1,temps2) # combine two matrices into a single matrix
temps <- as.vector(t(all.temps)) # convert table of temperature data into a single vector
allen_df <- cbind.data.frame(dates, temps) # column bind the previously defined date vector and the temperature vector
colnames(allen_df) <- c("t", "temp")
# write data to .csv file
write.table(allen_df, "../SST_data/Allen.Island.2015-2022.daily.mean.SST.csv", sep=",", row.names = F, col.names = T)

allen_df$t <- as.Date(allen_df$t)
# process marine heat wave event data
allen_hot_ts <- ts2clm(allen_df, climatologyPeriod = c("2015-01-01", "2022-12-31"))
allen_mhw <- detect_event(allen_hot_ts)
# process marine cold spell event data
allen_cold_ts <- ts2clm(allen_df, climatologyPeriod = c("2015-01-01", "2022-12-31"), pctile = 10)
allen_mcs <- detect_event(allen_cold_ts, coldSpells = TRUE)
# write data to .csv files
write.table(allen_mhw$event, "../event_data/Allen.Island.2015-2022.hot.events.csv", sep=",", row.names = F, col.names = T)
write.table(allen_mcs$event, "../event_data/Allen.Island.2015-2022.cold.events.csv", sep=",", row.names = F, col.names = T)