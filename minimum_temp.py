import xarray as xr
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.enums import Resampling
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import matplotlib.patches as mpatches

# Load the NetCDF file
nc_file = r'E:\project execution\geopandas\New thematic maps\minimum temperature map\min_t.nc'
ds = xr.open_dataset(nc_file)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
min_temperature_2m = ds['air_temperature_at_2_metres']# You might need to adjust variable names based on your NetCDF file structure
#min_temperature_2m_2d = min_temperature_2m.isel(time=0)
lat = ds['lat']
lon = ds['lon']

#18.52807617, 22.30535889

#Reclassifying according to the suitability levels
#highly_suitable_threshold_mintemp_1, highly_suitable_threshold_mintemp_2 = 29, 32
#moderately_suitable_threshold_mintemp = 21
#marginally_suitable_threshold_mintemp = 18


#highly_suitable_mintemp = np.where((min_temperature_2m_2d >= highly_suitable_threshold_mintemp_1) & (min_temperature_2m_2d <= highly_suitable_threshold_mintemp_2), 1, 0)
#moderately_suitable_mintemp = np.where((min_temperature_2m_2d >= moderately_suitable_threshold_mintemp) & (min_temperature_2m_2d < highly_suitable_threshold_mintemp_1), 0.67, 0)
#marginally_suitable_mintemp = np.where((min_temperature_2m_2d < moderately_suitable_threshold_mintemp) & (min_temperature_2m_2d >= marginally_suitable_threshold_mintemp), 0.33, 0)
#unsuitable_mintemp = np.where(min_temperature_2m_2d < marginally_suitable_threshold_mintemp, 0, 0)

#combined_classes_mintemp = highly_suitable_mintemp + moderately_suitable_mintemp + marginally_suitable_mintemp + unsuitable_mintemp


values = min_temperature_2m.values

#Load study area shapefile
study_area = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp') 

boundaries = [18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23]  # Adjust these boundaries as needed
#cmap = plt.get_cmap('plasma', len(boundaries) - 1)

#boundaries = [0, 0.25, 0.5, 0.75, 1]  # Adjust these boundaries as needed
labels = ['#FDE725', '#22A884', '#30678D', '#440154']  # Custom labels for each range
categories = ['Highly Suitable (0.75 - 1)', 'Moderately Suitable (0.5 - 0.75)', 'Marginally Suitable (0.25 - 0.5)', 'Unsuitable (0 - 0.25)']
cmap = plt.get_cmap('coolwarm', len(boundaries) - 1)

#Create a list  of patches for the legend
patches = [mpatches.Patch(color=labels[i], label=categories[i]) for i in range(len(categories))]


# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, min_temperature_2m, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1.0)

# Label each polygon with the respective label from the metadata
for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
    label = row['UGA_adm4_9']
    # Get the centroid of the polygon to place the label
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, bbox=dict(facecolor='none', alpha=0.7, edgecolor='none'))

# Add a custom legend to the plot
#legend = ax.legend(handles=patches, loc='lower right', title="Legend")

# Add colorbar
cbar = plt.colorbar(contour, ax=ax, label='Minimum Daily Temperature Values (degrees celcius)', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Minimum Daily Temperature Thematic Map')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

# Displaying Histograph
plt.hist(values, bins=20, alpha=0.7)
plt.title('Distribution of Minimum Daily Temperature Values (degrees celcius)')
plt.xlabel('Minimum Daily Temperature Values (degrees celcius)')
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
#slope_path = r'E:\project execution\geopandas\New thematic maps\minimum temperature map\mint.tif'
#with rasterio.open(slope_path, 'w', **meta) as dst:
    #dst.write(min_temperature_2m_2d.astype(rasterio.float32), 1)















