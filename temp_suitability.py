import xarray as xr
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
from rasterio.enums import Resampling

# Load the NetCDF file
nc_file = r'E:\project execution\geopandas\New thematic maps\final_suitability_map\f_suitability.nc'
ds = xr.open_dataset(nc_file)

final_suitability = ds['final suitability']  # Replace 'rainfall' with your variable name
lat = ds['lat']
lon = ds['lon']


# Load your shapefile
study_area = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp')


boundaries = [0, 0.5, 1.5, 2.5, 3.5] 
#labels = ['#FDE725', '#22A884', '#30678D', '#440154']  # Custom labels for each range
categories = ['Highly Suitable', 'Moderately Suitable', 'Marginally Suitable', 'Unsuitable']
#cmap = plt.get_cmap('viridis', len(boundaries) - 1)

labels_2 = ['#FDE725', '#22A884', '#000000', '#440154']

custom_cmap = ['#440154', '#000000', '#22A884', '#FDE725']

#Create a list  of patches for the legend
patches = [mpatches.Patch(color=labels_2[i], label=categories[i]) for i in range(len(categories))]


# Create a BoundaryNorm instance
#norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, final_suitability, colors=custom_cmap, levels=boundaries)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.9)

# Label each polygon with the respective label from the metadata
for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
    label = row['UGA_adm4_9']
    # Get the centroid of the polygon to place the label
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, fontweight='bold', bbox=dict(facecolor='none', alpha=0.7, edgecolor='none'))

# Add a custom legend to the plot
legend = ax.legend(handles=patches, loc='lower right', title="Legend")

# Add colorbar
#cbar = plt.colorbar(contour, ax=ax, label='Suitability Values', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('MAP OF BUIKWE DISTRICT SHOWING DIFFERENT CLASSES OF \n SUITABILITY FOR THE GROWTH OF COCOA')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

# Replace 'your_dem_file.tif' with the path to your DEM file
#dem_path = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil aluminium map(0_20)cm\alu.tif'
#NB: DEM has to be in a projected CS for the slope calculations

#with rasterio.open(dem_path) as dem:
#    # Read the DEM data, resampling as needed
    #elevation = dem.read(1, resampling=Resampling.bilinear)

    # Get metadata for later use
    #meta = dem.meta

# Write the slope data to a new file
#slope_path = r'E:\project execution\geopandas\New thematic maps\final_suitability_map\final_suitability_2.tif'
#with rasterio.open(slope_path, 'w', **meta) as dst:
    #dst.write(combined_classes.astype(rasterio.float32), 1)























