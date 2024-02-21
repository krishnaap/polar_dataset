#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 17:00:06 2024

@author: krishna
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import matplotlib.ticker as mticker
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.util import add_cyclic_point

ds = xr.open_dataset('/media/krishna/FAF039FFF039C325/temp/MODEL.SST.HAD187001-198110.OI198111-202206.nc')
min_val = -5
max_val = 25
norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
for month in range(1, 13):  # This will save 12 months in 1970
    sst_month = ds['SST'].sel(time=slice('1970-{:02d}-01'.format(month), '1970-{:02d}-28'.format(month))).mean(dim='time')

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=1))


    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='black', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True  # Turn on to see labels on the left
    gl.bottom_labels = True  # Turn on to see labels on the bottom
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 30))  # Longitude lines every 30 degrees
    gl.ylocator = mticker.FixedLocator(range(60, 91, 10))  # Latitude lines every 10 degrees from 60N to 90N
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
    
    ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())
    ax.coastlines()

    sst_cyclic, lon_cyclic = add_cyclic_point(sst_month, coord=ds['lon'])
    # Plot the data using pcolormesh
    mesh = ax.pcolormesh(ds['lon'], ds['lat'], sst_month,
                         transform=ccrs.PlateCarree(), cmap='RdYlBu_r',norm=norm)
    land = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                    edgecolor='face',
                                    facecolor='lightgray')
    ax.add_feature(land, zorder=1) # land masking.
    

    cbar = plt.colorbar(mesh, orientation='vertical', pad=0.05, aspect=30, shrink=0.8)     # Add a colorbar with bold label text
    cbar.set_label('Sea Surface Temperature (°C)', weight='bold')

    # Add title with the current month and year
    plt.title('SST for 1970-{:02d}'.format(month))

    # Show the plot
    plt.title('SST for 1970-{:02d}'.format(month))

    # Saving into PNG file
    plt.savefig('/home/krishna/files/out/SST_1970_{:02d}.png'.format(month), dpi=300, bbox_inches='tight') #year and months in the file name
    # plt.savefig('/media/krishna/FAF039FFF039C325/temp/SST_1970_{:02d}.png'.format(month), dpi=300, bbox_inches='tight')
    plt.show()
    
    
    
    
    plt.clf() # Clear the current figure's memory to make room for the next plot
