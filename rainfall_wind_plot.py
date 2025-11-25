'''
This script reads WRF output files and generates plots showing rainfall (shaded contours)
and 10-m wind as streamlines on the same map for each time step.
'''

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import pandas as pd
from matplotlib.colors import Normalize

# Open WRF file
wrfout = Dataset(r'C:\Users\pbuli\Desktop\WRF-things\wrfouts\wrfout_d01_2025-11-03_00_00_00')

# Read precipitation accumulations and wind components
prec_acc_c = xr.DataArray(wrfout.variables['PREC_ACC_C'][:])
prec_acc_nc = xr.DataArray(wrfout.variables['PREC_ACC_NC'][:])
total_prec = prec_acc_c + prec_acc_nc

u10 = xr.DataArray(wrfout.variables['U10'][:])
v10 = xr.DataArray(wrfout.variables['V10'][:])

# Wind speed
wind_speed = np.sqrt(u10.values**2 + v10.values**2)

# Parse Times (if present) for captions
time_index = None
if 'Times' in wrfout.variables:
    times_raw = wrfout.variables['Times'][:]
    times = []
    for t in times_raw:
        try:
            s = ''.join([c.decode('utf-8') if isinstance(c, bytes) else str(c) for c in t]).strip()
        except Exception:
            try:
                s = ''.join([chr(c) for c in t]).strip()
            except Exception:
                s = str(t)
        times.append(s)
    try:
        time_index = pd.to_datetime(times)
    except Exception:
        time_index = times

# Lat/Lon from dataset
lat = xr.DataArray(wrfout.variables['XLAT'][0, :, :])
lon = xr.DataArray(wrfout.variables['XLONG'][0, :, :])

# Consistent color scale for precipitation across times
arr = total_prec.values
vmin_rf = float(np.nanmin(arr))
vmax_rf = float(np.ceil(np.nanmax(arr)))

# Loop over time steps and create a single plot per time containing both rainfall and wind streamlines
for time_idx in range(total_prec.shape[0]):
    fig = plt.figure(figsize=(12, 9))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.2)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5,
                      linestyle='--', x_inline=False, y_inline=False,
                      xlabel_style={'size': 12}, ylabel_style={'size': 12})
    gl.top_labels = False
    gl.right_labels = False

    # Prepare precipitation for this time
    data_prec = np.asarray(total_prec[time_idx, :, :])
    levels = np.linspace(vmin_rf, vmax_rf, 21)
    cf = ax.contourf(lon, lat, data_prec, levels=levels, cmap='YlGnBu', extend='both', transform=ccrs.PlateCarree())

    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', shrink=0.7, aspect=20, pad=0.02)
    cbar.set_label('Rainfall Intensity (mm/h)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Prepare wind for this time (decimate grid for streamlines if needed)
    u = np.asarray(u10[time_idx, :, :])
    v = np.asarray(v10[time_idx, :, :])
    speed = np.sqrt(u**2 + v**2)

    ny, nx = lon.shape
    # choose a decimation step so streamplot is not overloaded
    step = max(1, int(min(ny, nx) / 40))
    lon_dec = lon.values[::step, ::step]
    lat_dec = lat.values[::step, ::step]
    u_dec = u[::step, ::step]
    v_dec = v[::step, ::step]
    speed_dec = speed[::step, ::step]

    # Normalize wind speed color scale
    vmin_ws = float(np.nanmin(wind_speed))
    vmax_ws = float(np.nanmax(wind_speed))
    norm_ws = Normalize(vmin=vmin_ws, vmax=vmax_ws)

    # Plot streamlines colored by wind speed with normalization
    strm = ax.streamplot(lon_dec, lat_dec, u_dec, v_dec,
                         density=1.2, color=speed_dec, cmap='plasma', norm=norm_ws,
                         linewidth=1.0, arrowsize=1.0)

    # Colorbar for wind speed
    cbar_w = plt.colorbar(strm.lines, ax=ax, orientation='vertical', shrink=0.6, aspect=30, pad=0.01)
    cbar_w.set_label('Wind speed (m/s)', fontsize=9)
    cbar_w.ax.tick_params(labelsize=8)
    '''
    Disclaimer: Try block checking if set_clim is supported. Help was obtained from ChatGPT.
    '''
    # Ensure the colorbar reflects the normalization used for the streamlines by setting
    # the limits on the underlying mappable and updating the colorbar.
    try:
        cbar_w.mappable.set_clim(vmin_ws, vmax_ws)
        cbar_w.update_normal(cbar_w.mappable)
    except Exception:
        # If the underlying object doesn't support set_clim, skip silently.
        pass

    # Titles and captions
    title = 'Rainfall Intensity (shaded) and Wind (streamlines)'
    if time_index is not None:
        try:
            formatted_time = pd.to_datetime(time_index[time_idx]).strftime('%Y-%m-%d %H:%M UTC')
        except Exception:
            formatted_time = str(time_index[time_idx])
        caption1 = f'Valid for {formatted_time}'
    else:
        caption1 = f'Valid for Time Index {time_idx}'

    caption2 = 'Weather Research & Forecasting (WRF) Model, ver. 4.6.1\nDriving model: ERA5 2025-11-03 00:00 UTC | Initialized 2025-11-03 00:00 UTC'
    ax.set_title(title, x=0.00, y=1.098, fontsize=18, ha='left', va='bottom')
    plt.figtext(0.125, 0.83, caption1, fontsize=14, ha='left', va='top')
    plt.figtext(0.125, 0.81, caption2, fontsize=12, ha='left', va='top')

    # Save figure
    outpath = rf'C:\Users\pbuli\Desktop\WRF-things\Plots\precip_wind_time_{time_idx}.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close(fig)