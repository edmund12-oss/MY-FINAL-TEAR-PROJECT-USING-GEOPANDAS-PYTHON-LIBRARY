#Import necessary python libraries
import geopandas as gpd
import xarray as xr
import numpy as np
import rasterio
from rasterio.enums import Resampling
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import matplotlib.patches as mpatches

# Load your rainfall data from a NetCDF file
rainfall_data_path = r'E:\project execution\geopandas\New thematic maps\rainfall_map\new_rf_data.nc'
rainfall_data = xr.open_dataset(rainfall_data_path)

# Assuming your rainfall data has lat and lon dimensions, and a rainfall variable
# Adjust 'rainfall' below to match the variable name in your dataset
rainfall_2003_2023 = rainfall_data['rainfall']  # Replace 'rainfall' with your variable name
lat = rainfall_data['lat']
lon = rainfall_data['lon']

#Extract rainfall values as numpy arrays to be used in creating a histograph
values = rainfall_2003_2023.values

# Load your shapefile
study_area = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp')

#1735.82299805, 1379.10754395
boundaries = [1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750]  # Adjust these boundaries as needed
cmap = plt.get_cmap('Blues', len(boundaries) - 1)

# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, rainfall_2003_2023, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.9)

# Label each polygon with the respective label from the metadata
for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
    label = row['UGA_adm4_9']
    # Get the centroid of the polygon to place the label
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, bbox=dict(facecolor='none', alpha=0.7, edgecolor='none'))

# Add colorbar
cbar = plt.colorbar(contour, ax=ax, label='Rainfall Values (mm)', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Rainfall Thematic Map') #Set title for your map
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4) #Set limit of longitude values to be displayed
plt.ylim(0, 0.6) #Set limit of latitude values to be displayed
plt.xlabel('Longitude') #Label x axis
plt.ylabel('Latitude') #Label y axis
plt.show()

# Displaying Histograph
plt.hist(values, bins=20, alpha=0.7)
plt.title('Distribution of Rainfall Values (mm)')
plt.xlabel('Rainfall Values (mm)')
plt.ylabel('Frequency')
plt.show()












 # Replace 'your_dem_file.tif' with the path to your DEM file
#dem_path = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil aluminium map(0_20)cm\alu.tif'
#NB: DEM has to be in a projected CS for the slope calculations

#with rasterio.open(dem_path) as dem:
    # Read the DEM data, resampling as needed
    #elevation = dem.read(1, resampling=Resampling.bilinear)

    # Get metadata for later use
    #meta = dem.meta

# Write the slope data to a new file
#slope_path = r'E:\project execution\geopandas\New thematic maps\rainfall_map\rfl.tif'
#with rasterio.open(slope_path, 'w', **meta) as dst:
    #dst.write(rainfall_2003_2023.astype(rasterio.float32), 1)


