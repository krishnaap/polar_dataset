#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 02:19:20 2024

@author: krishna
"""

import requests
import os

# Define the base URL and the range of years to cover
base_url = 'https://noaadata.apps.nsidc.org/NOAA/G02135/north/daily/geotiff/'
start_year = 1990
end_year = 2023
months = ['01_Jan', '02_Feb', '03_Mar', '04_Apr', '05_May', '06_Jun',
          '07_Jul', '08_Aug', '09_Sep', '10_Oct', '11_Nov', '12_Dec']

def download_file(url, path):
    """Download a single file and save it to the specified path."""
    local_filename = os.path.join(path, url.split('/')[-1])
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        os.makedirs(os.path.dirname(local_filename), exist_ok=True)
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print(f"Downloaded: {local_filename}")

for year in range(start_year, end_year + 1):
    for month in months:
        # Construct the URL for each month of each year
        month_url = f"{base_url}{year}/{month}/"
        try:
            response = requests.get(month_url)
            response.raise_for_status()  # Check if the URL exists
            # If the request is successful, parse for TIFF links
            for line in response.iter_lines():
                if line:
                    decoded_line = line.decode('utf-8')
                    if '.tif"' in decoded_line:
                        # Extract the .tif file name
                        start_pos = decoded_line.find('href="') + 6
                        end_pos = decoded_line.find('.tif"', start_pos) + 4
                        tif_file = decoded_line[start_pos:end_pos]
                        # Download the .tif file
                        download_file(f"{month_url}{tif_file}", f"geotiff_files/{year}/{month}")
        except requests.RequestException as e:
            print(f"Failed to access {month_url}: {e}")

#%%

