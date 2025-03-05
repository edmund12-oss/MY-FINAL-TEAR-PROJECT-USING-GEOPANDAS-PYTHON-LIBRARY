import geopandas as gpd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import matplotlib.patches as mpatches

# Load the NetCDF file
nc_file = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil bulk density(0_20)cm\masked_bulk_density.nc'
ds = xr.open_dataset(nc_file) #Read the NetCDF file using xarray

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
bulk_density = ds['mean_0_20']  # Read the variable of interest within the NetCDF file
lat = ds['lat'] #Read the latitude values within the NetCDF file
lon = ds['lon'] #Read the longitude values within the NetCDF file

#Extract the values for bulk density into a numpy array to be used for displaying the histograph
values = bulk_density.values

#Load study area shapefile
study_area = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp') 

boundaries = [1, 1.1, 1.2, 1.3, 1.4, 1.5]  # Adjust these boundaries as needed
labels = ['#FDE725', '#22A884', '#30678D', '#440154']  # Custom labels for each range
cmap = plt.get_cmap('viridis', len(boundaries) - 1)

# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, bulk_density, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.6)

# Label each polygon with the respective label from the metadata
for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
    label = row['UGA_adm4_9']
    # Get the centroid of the polygon to place the label
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, bbox=dict(facecolor='none', alpha=0.7, edgecolor='none'))

# Add colorbar
cbar = plt.colorbar(contour, ax=ax, label='Soil Bulk Density Mean_0_20 Values (g/cm3)', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Soil Bulk Density (0-20)cm')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

# Displaying Histograph
plt.hist(values, bins=20, alpha=0.7)
plt.title('Distribution of Soil Bulk Density Mean_0_20 Values (g/cm3)')
plt.xlabel('Soil Bulk Density Mean_0_20 Values (g/cm3)')
plt.ylabel('Frequency')
plt.show()