###
# MedAtlasModelling
# Christine Spencer, April 2023
#
# R Workflow designed to extract data from bivariate frequency tables in the MedAtlas (2004) database, specifically orientation data 
# (wind direction or wave direction), that requires straightforward trigonometric functions to produce averages of the observations of the measurement
# stations. 
# The data is exported as CSVs from the database, from which point this workflow can be used, by setting the working directory to the folder in which
# the exported CSVs are located. Please note that if multiple seasons are being processed with the goal of analysing seasonal patterns, the CSVs
# of each season/time chunk must be kept in separate folders as this method will apply the functions to every CSV in the given working directory.

# The result will be a CSV file where each measuremenet station from the MedAtlas is summarised with the following attributes:
# latitude, longitude (coordinates for direct input into GIS), U and V components (commonly the format of directional data now available online, for
# easy integration with other datasets), average speed (in both m/s the native unit and converted into kph), the average direction (radians and 
# degrees), and the corrected average angle in degrees, which is the final average direction that should be used.

# This workflow includes toy datasets using the same CSV formatting as the MedAtlasDatabase. If using, set working directory to the folder where the
# data files are stored.
# This workflow also includes an optional workflow for a spatial interpolation of the data using IDW methods.

library(spatstat)
library(raster)
library(rgdal)
library(terra)
library(maptools)
library(rgeos)

#### load data
# toy data
#### this code is all completed in nested loops, for every file in specified directory folder

path="D://Phd/OtherProjects/Paula-CAA2023/MedAtlasData/freqtables/Winter/" # set directory to where files are stored (read intro for specifics)
list_csv_files <- list.files(path=path, pattern="*.csv")

# create vec which turns categorical data from frequency tables to a numeric value for each reading by the measurement station
speed_vec <- (c(seq(0.5,9.5,1),seq(11,19,2)))

# beginning of outer loop
for(i in 1:length(list_csv_files)){
  #for(i in 1:5){
  tmp1 <- read.csv(file=paste(path,list_csv_files[i],sep=""), sep=";", stringsAsFactors = FALSE)
  tmp <- tmp1[,1:17] # extract frequencies from table
  
  # inner loop which converts frequency of observations by speed and direction to individual readings for the given station in the outer loop
  for(j in 1:24){
    if(tmp[j,17]!=0) {
      test <- data.frame(reading=paste(list_csv_files[i],j,sep=""),
                         speed_ms=rep(speed_vec,tmp[j,2:16]),dir_deg=rep(as.numeric(tmp[j,1]),tmp[j,17]))
      if(j==1){ dat <- test }
      else{ dat <- rbind(dat,test) } }
  }
  
  # for each reading of the station given in the outer loop, the angles are transformed into their U and V components
  dat$EW_V_DEG <- dat$speed_ms*(sin(dat$dir_deg))
  dat$NS_U_DEG <- dat$speed_ms*(cos(dat$dir_deg))
  dat$EW_V_RAD <- dat$speed_ms*(sin((dat$dir_deg*pi)/180))
  dat$NS_U_RAD <- dat$speed_ms*(cos((dat$dir_deg*pi)/180))
  
  
  # creation of data.frame which will be exported 
  loc_avg <- data.frame(long=sub("E.*","",sub(".*, ","",tmp1[27,1])),
                        lat=substr(sub("N.*","",tmp1[27,1]),2,nchar(sub("N.*","",tmp1[27,1])))) # extract long and lat (x and y from csv)
  loc_avg$EW_V_RAD <-(sum(dat$EW_V_RAD)/nrow(dat))*-1
  loc_avg$NS_U_RAD <- (sum(dat$NS_U_RAD)/nrow(dat))*-1
  loc_avg$speed_ms <- sqrt(loc_avg$NS_U_RAD^2+loc_avg$EW_V_RAD^2)
  loc_avg$speed_kph <- loc_avg$speed_ms*3.6 # convert to kph
  loc_avg$dir_rad <- atan2(loc_avg$NS_U_RAD,loc_avg$EW_V_RAD)
  loc_avg$dir_deg <- (loc_avg$dir_rad*180)/pi
  if(loc_avg$dir_deg<180){ loc_avg$dir_deg_corr <- loc_avg$dir_deg+180 } 
  if (loc_avg$dir_deg>180) {loc_avg$dir_deg_corr <- loc_avg$dir_deg-180}
  if(i==1){ offshore_data <- loc_avg }
  else { 
    offshore_data <- rbind(offshore_data,loc_avg)
  }
}

# produces a data.frame where each row is the summary of each station CSV exported from the MedAtlas database
offshore_data

# write csv
write.csv(offshore_data,file="offshore_data_mean_summer.csv")
