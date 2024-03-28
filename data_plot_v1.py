#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 02:38:59 2024

@author: krishna
"""



import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import rasterio
import numpy as np

# Path to the GeoTIFF file
geotiff_path = '/media/krishna/FAF039FFF039C325/temp/Sea_Ice_extend_NOAA_data/geotiff_files/1978/12_Dec/N_19781211_concentration_v3.0.tif'

# Read the GeoTIFF file using rasterio
with rasterio.open(geotiff_path) as dataset:
    # Read the data
    data = dataset.read(1)
    
    # Define the Cartopy projection from the dataset's CRS
    crs = ccrs.NorthPolarStereo(central_longitude=0)
    data_crs = ccrs.epsg(dataset.crs.to_epsg())

    # Mask out the invalid data (assumed here to be -9999 or similar)
    data = np.ma.masked_where(data < 0, data)

    # Create a figure with an axes object with a polar stereographic projection
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=crs)
    ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())
    
    # Set the map background color to a blue shade for ocean
    ocean_color = '#2B3A67'
    ax.set_facecolor(ocean_color)

    # Add coastlines, land feature with grey color, and gridlines
    ax.coastlines()
    land_color = 'lightgray'
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor=land_color)
    
    # Set up gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='black', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 30))
    gl.ylocator = mticker.FixedLocator(range(60, 90, 10))
    gl.xlabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
    
    # Plot the data
    img_extent = (dataset.bounds.left, dataset.bounds.right, dataset.bounds.bottom, dataset.bounds.top)
    norm = mcolors.Normalize(vmin=0, vmax=100)  # Normalize data range to 0-100%
    
    # Create a custom colormap for ice data from blue to red
    cmap = mcolors.LinearSegmentedColormap.from_list("", ["#2B3A67", "white"])

    # Apply the normalization and colormap during the plot
    img = ax.imshow(data, origin='upper', extent=img_extent, transform=data_crs, cmap=cmap, interpolation='nearest', norm=norm)
    land = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgray')
    ax.add_feature(land, zorder=1)
    # Add colorbar with the same normalization
    cbar = plt.colorbar(img, orientation='vertical', pad=0.05, aspect=50, shrink=0.8, norm=norm)
    cbar.set_label('Sea Ice Concentration (%)', weight='bold')  # Adjust label as needed

    plt.show()
    # plt.savefig('/path/where/you/want/to/save/plot.png')
# =============================================================================
# Deprecated
# =============================================================================
#%%  Working one with messy colormap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import rasterio
import numpy as np

# Path to the GeoTIFF file
geotiff_path = '/media/krishna/FAF039FFF039C325/temp/Sea_Ice_extend_NOAA_data/geotiff_files/1978/12_Dec/N_19781211_concentration_v3.0.tif'

# Read the GeoTIFF file using rasterio
with rasterio.open(geotiff_path) as dataset:
    
    # Read the data
    data = dataset.read(1)
    
    # Get the metadata
    metadata = dataset.meta
    
    # Define the Cartopy projection from the dataset's CRS
    crs = ccrs.NorthPolarStereo(central_longitude=0)
    data_crs = ccrs.epsg(dataset.crs.to_epsg())

    # Mask out the invalid data (assumed here to be -9999 or similar)
    data = np.ma.masked_where(data < 0, data)

    # Create a figure with an axes object with a polar stereographic projection
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=crs)
    ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())
    
    # Add coastlines, land feature, and gridlines
    ax.coastlines()
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
    
    # Set up gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='black', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.xlocator = mticker.FixedLocator(range(-180, 181, 30))
    gl.ylocator = mticker.FixedLocator(range(60, 90, 10))
    gl.xlabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
    
    # Plot the data
    img_extent = (dataset.bounds.left, dataset.bounds.right, dataset.bounds.bottom, dataset.bounds.top)
    norm = mcolors.Normalize(vmin=0, vmax=100)  # Adjust the vmax to the max data value if it's different

    # When plotting the data, apply the normalization
    ax.imshow(data, origin='upper', extent=img_extent, transform=data_crs, cmap='RdYlBu_r', interpolation='nearest', norm=norm)
    land = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgray')
    ax.add_feature(land, zorder=1)
    
    # When creating the colorbar, pass the same normalization
    cbar = plt.colorbar(ax.images[0], orientation='vertical', pad=0.05, aspect=50, shrink=0.8, norm=norm)
    cbar.set_label('Sea Ice Concentration (%)', weight='bold')  # Adjust label as needed

    # Mask the land

    # Add a colorbar
    # cbar = plt.colorbar(ax.images[0], orientation='vertical', pad=0.05, aspect=50, shrink=0.8)
    # cbar.set_label('Data Value', weight='bold')

    # Save the figure
    # plt.savefig('/path/where/you/want/to/save/plot.png')

    plt.show()

#%%import matplotlib.pyplot as plt

#%%
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import rasterio
from rasterio.plot import show

# Path to your GeoTIFF file
geotiff_path = '/media/krishna/FAF039FFF039C325/temp/Sea_Ice_extend_NOAA_data/geotiff_files/1978/12_Dec/N_19781201_concentration_v3.0.tif'

# Open the GeoTIFF file
with rasterio.open(geotiff_path) as src:
    # Read the dataset's data
    data = src.read(1)
    # Get bounds and coordinate reference system from the GeoTIFF file
    bounds = src.bounds
    crs = src.crs

# Create a figure with Polar Stereographic projection
fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')

# Add gridlines
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = True
gl.bottom_labels = True
gl.xlocator = mticker.FixedLocator(range(-180, 181, 30))
gl.ylocator = mticker.FixedLocator(range(60, 90, 10))
gl.xlabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
gl.ylabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}

# Plot the data
extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
mappable = ax.imshow(data, origin='upper', extent=extent, transform=ccrs.PlateCarree(), cmap='RdYlBu_r')

# Add a colorbar
cbar = plt.colorbar(mappable, orientation='vertical', pad=0.05, aspect=50, shrink=0.8)
cbar.set_label('Data Value', weight='bold')

plt.show()
