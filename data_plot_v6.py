#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 02:38:59 2024

@author: krishna
version: 5
Additions : 3 sections.
1. plots for all months within a time period, as in years
2. finding mean and plotting for a particular month for an entire period years.
3. Finding meand and plotting for selected years
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import rasterio
import numpy as np
import os
import calendar

# Base path where GeoTIFF files are stored
base_path = '/media/krishna/FAF039FFF039C325/temp/Sea_Ice_extend_NOAA_data/geotiff_files/'

#%%============================================================================
# Section 1. Plotting monthly mean for the entire time period (1970-2023) or 
# for a selected years.
# =============================================================================
#_______comment or uncomment as per your requirement
# Set years as a range from 1970 to 2023
years = [str(year) for year in range(1979, 1981)]
# years = ['1978', '1982', '1990', '2001', '2020'] # Or, set years to specific years of interest

# Define a function to process files and calculate the mean
def process_files_for_month(year, month):
    # Initialize an empty list to store data arrays
    monthly_data = []
    
    # Construct the path to the month's directory
    month_folder_name = f"{month:02d}_{calendar.month_abbr[month]}"
    month_path = os.path.join(base_path, year, month_folder_name)
    
    # Check if the month path exists
    if not os.path.exists(month_path):
        print(f"Path does not exist: {month_path}")
        return None
    
    # List all GeoTIFF files in the month's directory
    for file_name in os.listdir(month_path):
        if file_name.endswith('.tif'):
            # Construct the file path
            file_path = os.path.join(month_path, file_name)

            # Read the GeoTIFF file
            with rasterio.open(file_path) as dataset:
                # Read and mask the data
                data = dataset.read(1)
                data = np.ma.masked_where(data < 0, data)
                monthly_data.append(data)
    
    
    # Calculate the mean of the data for the month
    if monthly_data:
        mean_data = np.ma.mean(np.ma.stack(monthly_data), axis=0)
        return mean_data, dataset.meta, dataset.bounds, dataset.crs
    else:
        return None

# Loop through each year and month and calculate the mean
for year in years:
    for month in range(1, 13):  # Months from 1 to 12
        result = process_files_for_month(year, month)
        if result:
            mean_data, metadata, bounds, crs = result
            
            # Plotting
            fig = plt.figure(figsize=(10, 10))
            ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
            ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())

            # Set the map background color to a blue shade for ocean
            ocean_color = '#2B3A67'
            ax.set_facecolor(ocean_color)

            # Add coastlines, land feature with grey color, and gridlines

            land_feature = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='lightgray')
            ax.add_feature(land_feature, zorder=2)
            ax.coastlines(zorder=3)            
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
            img_extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)
            norm = mcolors.Normalize(vmin=0, vmax=100)
            cmap = mcolors.LinearSegmentedColormap.from_list("", ["#2B3A67", "white"])
            img = ax.imshow(mean_data, origin='upper', extent=img_extent, transform=ccrs.epsg(crs.to_epsg()), cmap=cmap, interpolation='nearest', norm=norm)
            cbar = plt.colorbar(img, orientation='vertical', pad=0.05, aspect=50, shrink=0.8, norm=norm)
            cbar.set_label('Sea Ice Concentration (%)', weight='bold')
            
            plt.show()

            output_path = '/media/krishna/FAF039FFF039C325/temp/Sea_Ice_extend_NOAA_data/output'
            output_filename = f'{year}_{month:02d}_mean_sea_ice_concentration.png'
            output_filepath = os.path.join(output_path, output_filename)
            fig.savefig(output_filepath)
            plt.close(fig)
            print(f"Saved: {output_filepath}")
        else:
            print(f"No data to process for {year}-{month:02d}")
#%%============================================================================
# Section 2. Calculating the monthly mean for the entire time period (1970-2023)
# Save mean into variable name: mean_data
# Inputs are at the end of this section.
# =============================================================================

# Function to calculate the monthly mean, same type of function we will use to calculate 


def calculate_monthly_mean(years_list, month):
    global last_bounds, last_crs  # Declare global variables to store bounds and CRS
    data_list = []
    last_bounds = last_crs = None  # Reset these each time function is called

    for year in years_list:
        month_folder_name = f"{month:02d}_{calendar.month_abbr[month]}"
        month_path = os.path.join(base_path, year, month_folder_name)

        if not os.path.isdir(month_path):
            print(f"Directory does not exist, skipping: {month_path}")
            continue

        for file_name in filter(lambda f: f.endswith('.tif'), os.listdir(month_path)):
            file_path = os.path.join(month_path, file_name)

            with rasterio.open(file_path) as src:
                data = src.read(1)
                masked_data = np.ma.masked_less_equal(data, 0)
                data_list.append(masked_data)
                last_bounds = src.bounds  # Store the last bounds
                last_crs = src.crs  # Store the last CRS


    if data_list:
        # Calculate the mean
        stacked_data = np.ma.stack(data_list)
        mean_data = np.ma.mean(stacked_data, axis=0)
        ##
        fig = plt.figure(figsize=(10, 10))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
        ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())

        # Set the map background color to a blue shade for ocean
        ocean_color = '#2B3A67'
        ax.set_facecolor(ocean_color)

        # Add coastlines, land feature with grey color, and gridlines

        land_feature = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='lightgray')
        ax.add_feature(land_feature, zorder=2)
        ax.coastlines(zorder=3)            
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
        img_extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)
        norm = mcolors.Normalize(vmin=0, vmax=100)
        cmap = mcolors.LinearSegmentedColormap.from_list("", ["#2B3A67", "white"])
        img = ax.imshow(mean_data, origin='upper', extent=img_extent, transform=ccrs.epsg(crs.to_epsg()), cmap=cmap, interpolation='nearest', norm=norm)
        cbar = plt.colorbar(img, orientation='vertical', pad=0.05, aspect=50, shrink=0.8, norm=norm)
        cbar.set_label('Sea Ice Concentration (%)', weight='bold')
        
        plt.show()

        output_path = '/media/krishna/FAF039FFF039C325/temp/Sea_Ice_extend_NOAA_data/output'
        output_filename = f'Whole_{month:02d}_mean_sea_ice_concentration.png'  #uncomment during whole time period           
        output_filepath = os.path.join(output_path, output_filename)
        fig.savefig(output_filepath)
        plt.close(fig)
        print(f"Saved: {output_filepath}")
   
        ##
        return np.ma.mean(np.ma.stack(data_list), axis=0)
    else:
        print("No valid data to process.")
        return None

# Set the desired whole time period in this following range.
# Variables to store the last processed bounds and CRS for plotting use
last_bounds = None
last_crs = None

years_range = [str(year) for year in range(1979, 2023)]
month_to_process = 1
mean_data = calculate_monthly_mean(years_range, month_to_process)


#%%============================================================================
# Section 3. Calculating monthly mean for selected years (e.g., 1983, 1990, 2022)
# Save this into variable selected_mean
# Inputs are at the end of this section.
# =============================================================================

# Function to calculate the monthly mean, same function we will use to calculate 
# whole data mean and mean for selected years as well.
def calculate_monthly_mean(years_list, month):
    data_list = []
    last_bounds = last_crs = None  # These variables will store the last read dataset bounds and CRS

    for year in years_list:
        month_folder_name = f"{month:02d}_{calendar.month_abbr[month]}"
        month_path = os.path.join(base_path, year, month_folder_name)

        if not os.path.isdir(month_path):
            print(f"Directory does not exist, skipping: {month_path}")
            continue

        for file_name in filter(lambda f: f.endswith('.tif'), os.listdir(month_path)):
            file_path = os.path.join(month_path, file_name)

            with rasterio.open(file_path) as src:
                data = src.read(1)
                masked_data = np.ma.masked_less_equal(data, 0)
                data_list.append(masked_data)
                last_bounds = src.bounds
                last_crs = src.crs

    if data_list:
        # Calculate the mean
        stacked_data = np.ma.stack(data_list)
        mean_data = np.ma.mean(stacked_data, axis=0)
        
        ##
        fig = plt.figure(figsize=(10, 10))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
        ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())

        # Set the map background color to a blue shade for ocean
        ocean_color = '#2B3A67'
        ax.set_facecolor(ocean_color)

        # Add coastlines, land feature with grey color, and gridlines

        land_feature = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='lightgray')
        ax.add_feature(land_feature, zorder=2)
        ax.coastlines(zorder=3)            
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
        img_extent = (bounds.left, bounds.right, bounds.bottom, bounds.top)
        norm = mcolors.Normalize(vmin=0, vmax=100)
        cmap = mcolors.LinearSegmentedColormap.from_list("", ["#2B3A67", "white"])
        img = ax.imshow(mean_data, origin='upper', extent=img_extent, transform=ccrs.epsg(crs.to_epsg()), cmap=cmap, interpolation='nearest', norm=norm)
        cbar = plt.colorbar(img, orientation='vertical', pad=0.05, aspect=50, shrink=0.8, norm=norm)
        cbar.set_label('Sea Ice Concentration (%)', weight='bold')
        
        plt.show()

        output_path = '/media/krishna/FAF039FFF039C325/temp/Sea_Ice_extend_NOAA_data/output'
        output_filename = f'Selected_{month:02d}_mean_sea_ice_concentration.png' #uncomment during selected time period.
        output_filepath = os.path.join(output_path, output_filename)
        fig.savefig(output_filepath)
        plt.close(fig)
        print(f"Saved: {output_filepath}")
   
        ##
        return mean_data, last_bounds, last_crs  # Return the mean data along with the last bounds and CRS
    else:
        print("No valid data to process.")
        return None, None, None  # Return None for all values if there's nothing to process

# Use the function for selected years and store the additional information
selected_years = ['1983', '1990', '2022']
selected_mean, selected_bounds, selected_crs = calculate_monthly_mean(selected_years, month_to_process)


#%%============================================================================
# Section 4. Anomaly
# Calculate the anomaly
# and plot it
# =============================================================================
anomaly = selected_mean - mean_data

# Plotting the anomaly
fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
ax.set_extent([-180, 180, 55, 90], crs=ccrs.PlateCarree())

# Add land, coastlines, and gridlines
land_feature = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='lightgray')
ax.add_feature(land_feature, zorder=2)
ax.coastlines(zorder=3)
gl = ax.gridlines(draw_labels=True, linewidth=1, color='black', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = True
gl.bottom_labels = True
gl.xlocator = mticker.FixedLocator(range(-180, 181, 30))
gl.ylocator = mticker.FixedLocator(range(60, 90, 10))
gl.xlabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
gl.ylabel_style = {'size': 12, 'color': 'black', 'weight': 'bold'}

# Calculate min and max values for the normalization, ensure they are symmetrical around zero for anomalies
# anomaly_min = np.min(anomaly)
# anomaly_max = np.max(anomaly)
# vmax = max(abs(anomaly_min), abs(anomaly_max))
# vmin = -vmax
vmax = 10
vmin = -10


# Set up the normalization and the colormap
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
cmap = 'RdBu'  # This colormap is good for differences or anomalies

# Plot the data with the specified norm and colormap
img = ax.imshow(anomaly, origin='upper', extent=(selected_bounds.left, selected_bounds.right, selected_bounds.bottom, selected_bounds.top),
                transform=ccrs.epsg(selected_crs.to_epsg()), cmap=cmap, norm=norm, interpolation='nearest')

# Add a colorbar
cbar = plt.colorbar(img, orientation='vertical', pad=0.05, aspect=50, shrink=0.8)
cbar.set_label('Sea Ice Concentration Anomaly (%)', weight='bold')

# Title and display
plt.title('Sea Ice Concentration Anomaly')
plt.show()

# Save the plot
output_filename = 'sea_ice_concentration_anomaly.png'
output_filepath = os.path.join(output_path, output_filename)
fig.savefig(output_filepath, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {output_filepath}")
