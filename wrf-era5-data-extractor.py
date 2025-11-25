'''
This script will extract data from WRF output files produced from ERA5 datasets
and calculate the extracted variable over a specific location by identifying
the nearest grid point to the given latitude and longitude coordinates. Outputs
will be exported as CSV files for perusal.
'''

# Importing necessary libraries
from netCDF4 import Dataset
from wrf import (getvar, ALL_TIMES, latlon_coords)
import pandas as pd

# Reading the wrfout file
wrfout_file = Dataset(r'C:\Users\pbuli\Desktop\WRF-things\wrfout_d01_2025-11-03_00_00_00')

# Extracting meteorological variable of interest at all time steps. Change variable as needed
prec_acc_c = getvar(wrfout_file, 'PREC_ACC_C', timeidx=ALL_TIMES) 
prec_acc_nc = getvar(wrfout_file, 'PREC_ACC_NC', timeidx=ALL_TIMES) 
total_prec = prec_acc_c + prec_acc_nc # Comment out this line if only one variable is needed

# Taking the lat and lon coordinates
lats, lons = latlon_coords(total_prec)

# Specifying the target latitude and longitude coordinates. Change as needed.
target_lat = 10.2926 # Cebu City
target_lon = 123.9022

# Finding the nearest grid point to the specified coordinates
lat_idx = (abs(lats[:,0].values - target_lat)).argmin()
lon_idx = (abs(lons[0,:].values - target_lon)).argmin()

# Extracting precipitation data for the nearest grid point at all time steps
prec_at_point = total_prec[:, lat_idx, lon_idx].values

# Extracting the time values
time_values = total_prec['Time'].values

# Creating a DataFrame to hold the extracted data
prec_data = {'DateTime': time_values, '1-hr Rainfall Accumulation (mm)': prec_at_point}
prec_df = pd.DataFrame(prec_data)

# Converting the DateTime values to a readable format
prec_df['DateTime'] = pd.to_datetime(prec_df['DateTime'])

# Saving the DataFrame to a CSV file
prec_df.to_csv(r'C:\Users\pbuli\Desktop\WRF-things\CSV\precipitation_at_point.csv', index=False)



# End of program