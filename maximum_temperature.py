import xarray as xr
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.enums import Resampling
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import matplotlib.patches as mpatches

# Load the NetCDF file
nc_file = r'C:\Users\hp\Downloads\mg_wind_speed.nc'
ds = xr.open_dataset(nc_file)


# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
temperature_2m = ds['__xarray_dataarray_variable__']# You might need to adjust variable names based on your NetCDF file structure
temperature_2m_2d = temperature_2m.isel(time=0)
lat = ds['lat']
lon = ds['lon']

min_val = temperature_2m_2d.min()
print(min_val)


values = temperature_2m_2d.values

#3.2126546, 1.6348957

#Load study area shapefile
study_area = gpd.read_file(r'E:\Surveying\MARY GIFT\Nabajuzi Shapefile\nabajuzi_shp.shp') 

boundaries = [1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]  # Adjust these boundaries as needed
cmap = plt.get_cmap('plasma', len(boundaries) - 1)

#boundaries = [0, 0.25, 0.5, 0.75, 1]  # Adjust these boundaries as needed
#labels = ['#FDE725', '#22A884', '#30678D', '#440154']  # Custom labels for each range
#categories = ['Highly Suitable (0.75 - 1)', 'Moderately Suitable (0.5 - 0.75)', 'Marginally Suitable (0.25 - 0.5)', 'Unsuitable (0 - 0.25)']
#cmap = plt.get_cmap('coolwarm', len(boundaries) - 1)

#Create a list  of patches for the legend
#patches = [mpatches.Patch(color=labels[i], label=categories[i]) for i in range(len(categories))]


# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots( figsize =(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, temperature_2m_2d, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1.0)

# Label each polygon with the respective label from the metadata
#for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
#   label = row['UGA_adm4_9']
   # Get the centroid of the polygon to place the label
#   centroid = row.geometry.centroid
#   ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, bbox=dict(facecolor='none', alpha=0.7, edgecolor='none'))

# Add a custom legend to the plot
#legend = ax.legend(handles=patches, loc='lower right', title="Legend")

# Add colorbar
cbar = plt.colorbar(contour, ax=ax, label='Maximum Daily Temperature Values (degrees celcius)', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Maximum Daily Temperature Thematic Map')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
#plt.xlim(32.7, 33.4)
#plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

# Displaying Histograph
#plt.hist(values, bins=20, alpha=0.7)
#plt.title('Distribution of Maximum Daily Temperature Values (degrees celcius)')
#plt.xlabel('Maximum Daily Temperature Values (degrees celcius)')
#plt.ylabel('Frequency')
#plt.show()

 # Replace 'your_dem_file.tif' with the path to your DEM file
dem_path = r'E:\Surveying\MARY GIFT\mg_dem\dem_projected.tif'
#NB: DEM has to be in a projected CS for the slope calculations

with rasterio.open(dem_path) as dem:
    # Read the DEM data, resampling as needed
    #elevation = dem.read(1, resampling=Resampling.bilinear)

    # Get metadata for later use
    meta = dem.meta

# Write the slope data to a new file
slope_path = r'E:\Surveying\MARY GIFT\WEATHER\windspeed\mg_windspeed.tif'
with rasterio.open(slope_path, 'w', **meta) as dst:
    dst.write(temperature_2m_2d.astype(rasterio.float32), 1)