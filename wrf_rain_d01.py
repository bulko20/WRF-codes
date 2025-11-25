import pytz
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import cartopy.crs as crs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors

from wrf import (to_np, getvar, ALL_TIMES, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

###################
#### FUNCTIONS ####
###################

def read_rgb_file(rgb_file_path):
    with open(rgb_file_path, 'r') as open_file:
        rgb_lines = open_file.readlines()

    converted_colors = []
    for rgb_line in rgb_lines:
        # Split the line into components
        rgb_components = rgb_line.split()

        # Convert the RGB values to the range 0-1 and add to the list
        r, g, b = [int(c) / 255 for c in rgb_components[1:4]]
        converted_colors.append((r, g, b))

    return converted_colors

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Preparing the colormap
color_file_path1 = 'C:/Users/Brian/Desktop/WRF visualization/cloudcover.rgb'
cloud_colors = read_rgb_file(color_file_path1)
cloud_colormap = mcolors.ListedColormap(cloud_colors)

color_file_path2 = 'C:/Users/Brian/Desktop/WRF visualization/PhilRainRate_extended.rgb'
rain_colors = read_rgb_file(color_file_path2)
rain_colors[0] = (rain_colors[0][0], rain_colors[0][1], rain_colors[0][2], 0.0)  # Set alpha to 0.0 for full transparency of the first color in the colormap (white)
rain_colormap = mcolors.ListedColormap(rain_colors)

# Preparing file path to shapefile of Philippine coastline
PH_shapefile_province = 'C:/Users/Brian/Desktop/WRF visualization/gadm41_PHL_1.shp'

# Importing the wrfout file
wrfout_file = Dataset('C:/Users/Brian/Desktop/DMet WRF Outputs/wrfout_d01_2024-10-21_12_00_00')

# Get the rainfall datasets
prec_data_c = getvar(wrfout_file, 'PREC_ACC_C', timeidx = ALL_TIMES)
prec_data_nc = getvar(wrfout_file, 'PREC_ACC_NC', timeidx = ALL_TIMES)

# Get the cloud fraction dataset
cloudfrac = getvar(wrfout_file, 'CLDFRA', timeidx = ALL_TIMES)

# Getting the time values
times = prec_data_c['Time'].values

# Get the cartopy mapping object
cart_proj = get_cartopy(prec_data_c)

for time in times:
    # Extracting precipitation intensity and total cloud cover data
    # for a specific date and time
    prec = prec_data_c.sel(Time = time) + prec_data_nc.sel(Time = time)
    clouds = cloudfrac.sel(Time = time).max(dim='bottom_top') * 100

    # Get the latitude and longitude points
    lats, lons = latlon_coords(prec_data_c)
    
    # Create a figure
    fig = plt.figure(figsize=(25., 25.), dpi=250)
    
    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)
    
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(prec_data_c))
    ax.set_ylim(cartopy_ylim(prec_data_c))
    
    # Add the gridlines
    gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 2),
                             ylocs=np.arange(-90, 91, 2), color='gray',
                             linestyle='--')
    gridlines.top_labels = False
    gridlines.right_labels = False
    gridlines.xlabel_style = {'size': 20}
    gridlines.ylabel_style = {'size': 20}
    
    # Add the coastlines
    ax.coastlines('10m', edgecolor='black', linewidth=3)
    shape_feature_province = ShapelyFeature(Reader(PH_shapefile_province).geometries(),
                                            crs.PlateCarree(), edgecolor='black', 
                                            facecolor='none', linewidth=0.5)
    ax.add_feature(shape_feature_province)
    
    # Make the filled contours for cloud cover and rainfall
    isonephs = plt.contourf(to_np(lons), to_np(lats), to_np(clouds), 
                            levels=np.arange(0, 100, 1),
                            transform=crs.PlateCarree(),
                            cmap=cloud_colormap)

    isohyets = plt.contourf(to_np(lons), to_np(lats), to_np(prec), 
                            levels=np.arange(0, 100, 0.1),
                            transform=crs.PlateCarree(),
                            cmap=rain_colormap, extend='max')
    
    # Add the colorbars
    cax1 = fig.add_axes([0.125, 0.070, 0.77, 0.01])
    norm1 = Normalize(vmin=0, vmax=100)
    colorbar1 = ColorbarBase(cax1, cmap=cloud_colormap, norm=norm1,
                             orientation='horizontal', ticks=np.arange(0, 101, 5))
    colorbar1.set_label('Total cloud cover (%)', fontsize=25)
    colorbar1.ax.tick_params(labelsize=15)

    cax2 = fig.add_axes([0.125, 0.025, 0.77, 0.01])
    norm2 = Normalize(vmin=0, vmax=100)
    colorbar2 = ColorbarBase(cax2, cmap=rain_colormap, norm=norm2,
                             extend='max', orientation='horizontal')
    colorbar2.set_ticks(ticks=[0, 2.5, 7.5, 15.0, 30.0, 60.0, 100.0],
                       labels=['0', '2.5', '7.5', '15.0', '30.0', '60.0', '100.0'])
    colorbar2.set_label('Rainfall intensity (mm/h)', fontsize=25)
    colorbar2.ax.tick_params(labelsize=15)

    # Formatting time before printing on the graphics
    time_dt = pd.to_datetime(time).to_pydatetime()
    utc = pytz.utc
    phst = pytz.timezone('Asia/Manila')
    time_dt_utc = utc.localize(time_dt)
    time_dt_phst = time_dt_utc.astimezone(phst)
    formatted_time = time_dt_phst.strftime('%I:%M %p PhST %d %B %Y')

    # Adding plot title and captions
    title = 'Rainfall intensity and cloud cover'
    caption1 = f'Valid for {formatted_time}'
    caption2 = 'Weather Research & Forecasting (WRF) Model, ver. 4.6.0\nDriving model: GFS 2024-10-21 12:00 UTC | Initialized 2024-10-21 12:00 UTC'
    plt.suptitle(title, x=0.125, y=0.977, fontsize=60, ha='left', va='top')
    plt.figtext(0.125, 0.941, caption1, fontsize=40, ha='left', va='top')
    plt.figtext(0.125, 0.915, caption2, fontsize=25, ha='left', va='top')
    
    plt.show()
    
    break

# End of program