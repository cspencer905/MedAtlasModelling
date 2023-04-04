# MedAtlasModelling
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
