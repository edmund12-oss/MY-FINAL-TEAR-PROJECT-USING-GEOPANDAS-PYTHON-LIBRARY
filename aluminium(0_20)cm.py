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
nc_file = r'C:\Users\hp\Desktop\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil aluminium map(0_20)cm\aluminium_extractable_mean_0_20 (1).nc'
ds = xr.open_dataset(nc_file)


# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
#Aluminium, extractable, predicted mean at 0-20 cm depth
aluminium_extractable = ds['Aluminium, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds['y']
lon = ds['x']

#converting from mg/kg to g/kg
aluminium_extractable_final = aluminium_extractable / 1000


#Reclassifying according to the suitability levels
highly_suitable_threshold_al = 0.5
moderately_suitable_threshold_al = 1
marginally_suitable_threshold_al = 2


highly_suitable_al = np.where(aluminium_extractable_final < highly_suitable_threshold_al, 1, 0)
moderately_suitable_al = np.where((aluminium_extractable_final >= highly_suitable_threshold_al) & (aluminium_extractable_final <= moderately_suitable_threshold_al), 0.67, 0)
marginally_suitable_al = np.where((aluminium_extractable_final > moderately_suitable_threshold_al) & (aluminium_extractable_final <= marginally_suitable_threshold_al), 0.33, 0)
unsuitable_al = np.where(aluminium_extractable_final > marginally_suitable_threshold_al, 0, 0)

combined_classes_al = highly_suitable_al + moderately_suitable_al + marginally_suitable_al + unsuitable_al


#values = aluminium_extractable_final.values

#Load study area shapefile
study_area = gpd.read_file(r'C:\Users\hp\Desktop\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp') 

boundaries = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
#boundaries = [0, 0.25, 0.5, 0.75, 1]  # Adjust these boundaries as needed
labels = ['#FDE725', '#22A884', '#30678D', '#440154']  # Custom labels for each range
categories = ['Highly Suitable (0.75 - 1)', 'Moderately Suitable (0.5 - 0.75)', 'Marginally Suitable (0.25 - 0.5)', 'Unsuitable (0 - 0.25)']
cmap = plt.get_cmap('viridis', len(boundaries) - 1)

# Create a list of patches for the legend
patches = [mpatches.Patch(color=labels[i], label=categories[i]) for i in range(len(categories))]


# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)


# Calculate midpoints of boundaries for tick placement
#ticks = [(boundaries[i] + boundaries[i+1])/2 for i in range(len(boundaries)-1)]



# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, aluminium_extractable_final, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.6)


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
cbar = plt.colorbar(contour, ax=ax, label='Aluminium Extractable Mean_0_20 Values (g/kg)', norm=norm, boundaries=boundaries, ticks=boundaries)
#cbar.set_ticklabels(labels)  # Apply the custom labels to the colorbar
ax.set_title('Soil Aluminium (0-20)cm')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()



#Displaying Histogram
#plt.hist(values, bins=20, alpha=0.7)
#plt.title('Distribution of Aluminium Extractable Mean_0_20 Values (g/kg)')
#plt.xlabel('Aluminium Extractable Mean_0_20 Values (g/kg)')
#plt.ylabel('Frequency')
#plt.show()

# Replace 'your_dem_file.tif' with the path to your DEM file
#dem_path = r'C:\Users\hp\Desktop\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil aluminium map(0_20)cm\alu.tif'
#NB: DEM has to be in a projected CS for the slope calculations

#with rasterio.open(dem_path) as dem:
#    # Read the DEM data, resampling as needed
#    elevation = dem.read(1, resampling=Resampling.bilinear)

    # Get metadata for later use
#    meta = dem.meta

# Write the slope data to a new file
#slope_path = r'C:\Users\hp\Desktop\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil aluminium map(0_20)cm\al.tif'
#with rasterio.open(slope_path, 'w', **meta) as dst:
#    dst.write(aluminium_extractable_final.astype(rasterio.float32), 1)




    