import xarray as xr
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import matplotlib.patches as mpatches

# Load the NetCDF file
nc_file = r'C:\Users\hp\Desktop\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil phosphorous extractable(20_50)cm\masked_phosphorus.nc'
ds = xr.open_dataset(nc_file)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
phosphorus_extractable_20_50 = ds['Phosphorus, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds['lat']
lon = ds['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_phosphorus_1_20_50, highly_suitable_threshold_phosphorus_2_20_50 = 20, 50
moderately_suitable_threshold_phosphorus_20_50 = 10
marginally_suitable_threshold_phosphorus_20_50 = 5


highly_suitable_phosphorus_20_50 = np.where((phosphorus_extractable_20_50 >= highly_suitable_threshold_phosphorus_1_20_50) & (phosphorus_extractable_20_50 <= highly_suitable_threshold_phosphorus_2_20_50), 1, 0)
moderately_suitable_phosphorus_20_50 = np.where((phosphorus_extractable_20_50 >= moderately_suitable_threshold_phosphorus_20_50) & (phosphorus_extractable_20_50 < highly_suitable_threshold_phosphorus_1_20_50), 0.67, 0)
marginally_suitable_phosphorus_20_50 = np.where((phosphorus_extractable_20_50 < moderately_suitable_threshold_phosphorus_20_50) & (phosphorus_extractable_20_50 >= marginally_suitable_threshold_phosphorus_20_50), 0.33, 0)
unsuitable_phosphorus_20_50 = np.where(phosphorus_extractable_20_50 < marginally_suitable_threshold_phosphorus_20_50, 0, 0)

combined_classes_phosphorus_20_50 = highly_suitable_phosphorus_20_50 + moderately_suitable_phosphorus_20_50 + marginally_suitable_phosphorus_20_50 + unsuitable_phosphorus_20_50


#3.48168898, 17.17414665
values = ds['Phosphorus, extractable, predicted mean at 20-50 cm depth'].values

#Load study area shapefile
study_area = gpd.read_file(r'C:\Users\hp\Desktop\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp') 

boundaries = [0, 0.25, 0.5, 0.75, 1]  # Adjust these boundaries as needed
labels = ['#FDE725', '#22A884', '#30678D', '#440154']  # Custom labels for each range
categories = ['Highly Suitable (0.75 - 1)', 'Moderately Suitable (0.5 - 0.75)', 'Marginally Suitable (0.25 - 0.5)', 'Unsuitable (0 - 0.25)']
cmap = plt.get_cmap('viridis', len(boundaries) - 1)

#Create a list  of patches for the legend
patches = [mpatches.Patch(color=labels[i], label=categories[i]) for i in range(len(categories))]

# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, combined_classes_phosphorus_20_50, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.6)

# Label each polygon with the respective label from the metadata
for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
    label = row['UGA_adm4_9']
    # Get the centroid of the polygon to place the label
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, bbox=dict(facecolor='none', alpha=0.7, edgecolor='none'))

# Add a custom legend to the plot
legend = ax.legend(handles=patches, loc='lower right', title="Legend")

# Add colorbar
#cbar = plt.colorbar(contour, ax=ax, label='Soil Phosphorus Mean_20_50 Values (mg/kg)', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Soil Phosphorus (20-50)cm')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

# Displaying Histogram
#plt.hist(values, bins=20, alpha=0.7)
#plt.title('Distribution of Soil Phosphorus Mean_20_50 Values (mg/kg)')
#plt.xlabel('Soil Phosphorus Mean_20_50 Values (mg/kg)')
#plt.ylabel('Frequency')
#plt.show()