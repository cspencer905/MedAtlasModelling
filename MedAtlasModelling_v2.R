###
# MedAtlasModelling
# Christine Spencer, April 2023, updated September 2023
#
# R Workflow designed to extract data from bivariate frequency tables in the MedAtlas (2004) database, specifically orientation data 
# (wind direction or wave direction), that requires straightforward trigonometric functions to produce averages of the observations of the measurement
# stations. 
# The data is exported as CSVs from the database, from which point this workflow can be used, by setting the working directory to the folder in which
# the exported CSVs are located. Please note that if multiple seasons are being processed with the goal of analysing seasonal patterns, the CSVs
# of each season/time chunk must be kept in separate folders as this method will apply the functions to every CSV in the given working directory.
#
# The result will be a CSV file where each measurement station from the MedAtlas is summarised with the following attributes:
# latitude, longitude (coordinates for direct input into GIS), U and V components (commonly the format of directional data now available online, for
# easy integration with other datasets), average speed (in both m/s the native unit and converted into kph), the average direction (radians and 
# degrees), and the corrected average angle in degrees, which is the final average direction that should be used.
#
# This workflow includes toy datasets using the same CSV formatting as the MedAtlasDatabase. If using, set working directory to the folder where the
# data files are stored.
#
# This workflow also includes an optional workflow for a spatial interpolation of the data using thin plate spline methods. It presupposes 
# the summary table produced by the main part of the code was joined to a point shapefile of the offshore data stations, as the data in this
# application was used in a local coordinate system for the Eastern Mediterranean (ESPG: 2100) and used in conjunction with coastal weather 
# station data. If using directly from the summary table of offshore data, ensure that the amalgamated file is called the same object (stations)
# and its coordinates can be accessed by calling @coords.


library(spatstat)
library(raster)
library(rgdal)
library(terra)
library(rgeos)

########## MedAtlas data Extraction #############
## load data
# set path and create list of all tables exported from the MedAtlas database
path="D://MAM_toy"                 # set directory to where files are stored (read intro for specifics)
list_csv_files <- list.files(path=path, pattern="*.csv")

# create vec which turns categorical data from frequency tables to a numeric value for each reading by the measurement station
speed_vec <- (c(seq(0.5,9.5,1),seq(11,19,2)))

#### this code is all completed in nested loops, for every file in specified directory folder
# beginning of outer loop
for(i in 1:length(list_csv_files)){
  tmp1 <- read.csv(file=paste(path,list_csv_files[i],sep=""), sep=";", stringsAsFactors = FALSE)
  tmp <- tmp1[,1:17]        # extract frequencies from table
  
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
  # convert into sin and cosine components - from degree, resultant components are in radians 
  dat$EW_V_RAD <- dat$speed_ms*(sin((dat$dir_deg*pi)/180))
  dat$NS_U_RAD <- dat$speed_ms*(cos((dat$dir_deg*pi)/180))
  
  
  # creation of data.frame which will be exported 
  loc_avg <- data.frame(long=sub("E.*","",sub(".*, ","",tmp1[27,1])),
                        lat=substr(sub("N.*","",tmp1[27,1]),2,nchar(sub("N.*","",tmp1[27,1]))))     # extract long and lat (x and y from csv)
  loc_avg$EW_V_RAD <-(sum(dat$EW_V_RAD)/nrow(dat))*-1
  loc_avg$NS_U_RAD <- (sum(dat$NS_U_RAD)/nrow(dat))*-1                                              # U and V components
  loc_avg$speed_ms <- sqrt(loc_avg$NS_U_RAD^2+loc_avg$EW_V_RAD^2)                                   # average speed in m/s
  loc_avg$speed_kph <- loc_avg$speed_ms*3.6                                                         # convert to kph  
  loc_avg$dir_rad <- atan2(loc_avg$NS_U_RAD,loc_avg$EW_V_RAD)
  loc_avg$dir_deg <- (loc_avg$dir_rad*180)/pi
  if(loc_avg$dir_deg<180){ loc_avg$dir_deg_corr <- loc_avg$dir_deg+180 }                            # corrected average direction
  if (loc_avg$dir_deg>180) {loc_avg$dir_deg_corr <- loc_avg$dir_deg-180}
  if(i==1){ offshore_data <- loc_avg }
  else { 
    offshore_data <- rbind(offshore_data,loc_avg)
  }
}

# produces a data.frame where each row is the summary of each station CSV exported from the MedAtlas database - nrow=n database tables
print(offshore_data)

# write csv
write.csv(offshore_data,file="offshore_data_mean_summer.csv")                 # name file based on directory input
write.csv(offshore_data,file="offshore_data_mean_winter.csv")                 # name file based on directory input

## END OF DATABASE EXTRACTION ##



######## INTERPOLATION #########
# point data interpolation to raster surface
library(gstat)
library(raster)
library(rgdal)
library(fields)
library(sp)
library(mapview)

# load raster file of area of interest
sea <- raster("raster/sea.tif")
sea[sea <0] <- NA                     # if not masked to water (i.e. landmass not NA) turn land into NA
# should have a raster where sea is 0 and landmass is NA

# create point file from station summary table
dat <- read.csv("offshore_data_mean_winter.csv")    # choose seasonal file
dat <- read.csv("offshore_data_mean_summer.csv")    # choose seasonal file

sp.dat <- dat
coordinates(sp.dat) <- cbind(dat$long, dat$lat)
proj4string(sp.dat) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# plot points on map to ensure no error
mapview(sp.dat)

stations <- sp.dat

# add coastal data if available, otherwise save stations object as using writeOGR().

### thin plate spline model - from U and V components
# create model for interpolation to be called in function for U and V components (in radians) of angle
xy <- data.frame(stations@coords)

# U component
v <- as.numeric(stations$NS_U_RAD)
tps <- Tps(xy, v)
int_U_RAD <- interpolate(sea, tps)

# V component
v <- as.numeric(stations$EW_V_RAD)
tps <- Tps(xy, v)
int_V_RAD <- interpolate(sea, tps)

writeRaster(int_U_RAD,"winter_int_U_RAD_spline.tif")    # save seasonal file
writeRaster(int_V_RAD,"winter_int_V_RAD_spline.tif")    # save seasonal file

#writeRaster(int_U_RAD,"summer_int_U_RAD_spline.tif")    # save seasonal file
#writeRaster(int_V_RAD,"summer_int_V_RAD_spline.tif")    # save seasonal file

### Direction raster surface
# convert U and V components back to angle with atan function
int_dir_deg <- (atan2(int_U_RAD,int_V_RAD))*180/pi
int_dir_deg2 <- int_dir_deg
# add 180 to return to 360 degree angle units
int_dir_deg2 <- int_dir_deg2+180
# mask to sea
int_dir_deg_msk <- mask(int_dir_deg2, sea)
# plot result
plot(int_dir_deg_msk)

writeRaster(int_dir_deg2,"winter_uv_dir_deg_spl.tif")            # save seasonal file
writeRaster(int_dir_deg_msk,"winter_uv_dir_deg_spl_mask.tif")    # save seasonal file

#writeRaster(int_dir_deg2,"summer_uv_dir_deg_spl.tif")           # save seasonal file
#writeRaster(int_dir_deg_msk,"summer_uv_dir_deg_spl_mask.tif")   # save seasonal file

### Speed raster surface 
# function with U and V components
int_uv_speed <- (sqrt(int_U_RAD^2+int_V_RAD^2)) *3.6 #3.6 converts speed unit to kph
# mask to sea
int_uv_speed_msk <- mask(int_uv_speed, sea)
# plot result
plot(int_uv_speed_msk)

writeRaster(int_uv_speed,"winter_uv_speed_spline.tif")             # save seasonal file
writeRaster(int_uv_speed_msk,"winter_uv_speed_spline_mask.tif")    # save seasonal file

#writeRaster(int_uv_speed,"summer_uv_speed_spline.tif")            # save seasonal file
#writeRaster(int_uv_speed_msk,"summer_uv_speed_spline_mask.tif")   # save seasonal file

## END OF INTERPOLATION ##

##### END #######