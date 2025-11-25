'''
This script read WRF output files and generates plots for rainfall.
It attempts to use xarray, netcdf4, matplotlib, and cartopy libraries for data handling and visualization
which is contrary to wrf-python that is based on the netcdf4-python library.
'''

# Importing necessary libraries
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import pandas as pd

# Reading the dataset
wrfout = Dataset(r'C:\Users\pbuli\Desktop\WRF-things\wrfouts\wrfout_d01_2025-11-03_00_00_00')

# Extracting the variables of interest using xarray
prec_acc_c = xr.DataArray(wrfout.variables['PREC_ACC_C'][:])
prec_acc_nc = xr.DataArray(wrfout.variables['PREC_ACC_NC'][:])
total_prec = prec_acc_c + prec_acc_nc

# Parse the WRF Times variable (if present) to get datetimes for captions
time_index = None
if 'Times' in wrfout.variables:
    times_raw = wrfout.variables['Times'][:]
    times = []
    for t in times_raw:
        # each t is often an array of chars; attempt safe decoding to a string
        try:
            # handle array of bytes or strings
            s = ''.join([c.decode('utf-8') if isinstance(c, bytes) else str(c) for c in t]).strip()
        except Exception:
            try:
                s = ''.join([chr(c) for c in t]).strip()
            except Exception:
                s = str(t)
        times.append(s)
    # Let pandas parse the time strings (flexible on formats)
    try:
        time_index = pd.to_datetime(times)
    except Exception:
        # fallback: keep raw strings
        time_index = times
else:
    # If Times not available, try common alternatives
    if 'XTIME' in wrfout.variables:
        # XTIME often stores minutes since start; convert to pandas timedelta
        xt = wrfout.variables['XTIME'][:]
        try:
            # assume first value is minutes since start and convert
            base = pd.Timestamp(str(wrfout.creation_date)) if hasattr(wrfout, 'creation_date') else None
            time_index = pd.to_timedelta(xt, unit='m')
        except Exception:
            time_index = list(xt)
    else:
        time_index = None

# Extracting the lat and lon coordinates
lat = xr.DataArray(wrfout.variables['XLAT'][0, :, :])
lon = xr.DataArray(wrfout.variables['XLONG'][0, :, :])

# Compute global min and max across all time steps so color scale is consistent
arr = total_prec.values 
vmin = float(np.nanmin(arr))
vmax = float(np.ceil(np.nanmax(arr)))

# For loop that creates plots for each time step and saves them as PNG files
for time_idx in range(total_prec.shape[0]):
    plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    
    # Adding gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5,
                      linestyle='--', x_inline=False, y_inline=False,
                      xlabel_style={'size': 12}, ylabel_style={'size': 12})
    gl.top_labels = False
    gl.right_labels = False

    # Prepare data and consistent contour levels for color scale
    data = np.asarray(total_prec[time_idx, :, :])
    levels = np.linspace(vmin, vmax, 21)

    # Use contourf for filled contours (with transform for cartopy)
    cf = ax.contourf(lon, lat, data, levels=levels, cmap='YlGnBu', extend='both', transform=ccrs.PlateCarree())

    # Adding a colorbar
    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', shrink=0.7, aspect=20)
    cbar.set_label('Rainfall Intensity (mm/h)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Adding titles and captions
    title = 'Rainfall Intesity During the Passage of Typhoon Tino'
    # Format caption using parsed datetime if available
    if time_index is not None:
        try:
            formatted_time = pd.to_datetime(time_index[time_idx]).strftime('%Y-%m-%d %H:%M UTC')
        except Exception:
            formatted_time = str(time_index[time_idx])
        caption1 = f'Valid for {formatted_time}'
    else:
        caption1 = f'Valid for Time Index {time_idx}'
    caption2 = 'Weather Research & Forecasting (WRF) Model, ver. 4.6.1\nDriving model: ERA5 2025-11-03 00:00 UTC | Initialized 2025-11-03 00:00 UTC'
    ax.set_title(title, x=0.00, y=1.095, fontsize=18, ha='left', va='bottom',)
    plt.figtext(0.125, 0.855, caption1, fontsize=14, ha='left', va='top')
    plt.figtext(0.125, 0.83, caption2, fontsize=12, ha='left', va='top')

    # Saving the plot with higher resolution
    plt.savefig(rf'C:\Users\pbuli\Desktop\WRF-things\Plots\total_precip_time_{time_idx}.png')
    plt.close()

# End of program