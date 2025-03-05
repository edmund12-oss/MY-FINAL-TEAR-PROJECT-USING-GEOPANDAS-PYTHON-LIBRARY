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
nc_file = r'E:\project execution\geopandas\New thematic maps\LULC\lulc_reclass.nc'
ds = xr.open_dataset(nc_file)

lulc_suitability = ds['lulc']  # Replace 'rainfall' with your variable name
lat = ds['lat']
lon = ds['lon']

threshold = 0.5



highly_suitable = np.where(lulc_suitability < threshold, 1, 0)
unsuitable = np.where(lulc_suitability > threshold, 0, 0)


combined_classes = highly_suitable + unsuitable



print(combined_classes)
# Load your shapefile
study_area = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp')

boundaries = [0, 0.5, 1] 
cmap = plt.get_cmap('viridis', len(boundaries) - 1)

# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, combined_classes, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.9)

# Label each polygon with the respective label from the metadata
for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
    label = row['UGA_adm4_9']
    # Get the centroid of the polygon to place the label
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, bbox=dict(facecolor='none', alpha=0.7, edgecolor='none'))

# Add colorbar
cbar = plt.colorbar(contour, ax=ax, label='Suitability Values', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Suitability Map')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()
