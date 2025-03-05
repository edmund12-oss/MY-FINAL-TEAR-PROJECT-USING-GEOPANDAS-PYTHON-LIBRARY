import geopandas as gpd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import matplotlib.patches as mpatches

# Load the NetCDF file
nc_file = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil bedrock depth(0_20)cm\masked_depth.nc'
ds = xr.open_dataset(nc_file)



# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
bedrock_depth = ds['mean_0_200']  # You might need to adjust variable names based on your NetCDF file structure
#bedrock_depth_2d =bedrock_depth.isel(time=0)
lat = ds['lat']
lon = ds['lon']


#Reclassifying according to the suitability levels
highly_suitable_threshold_bd = 120
moderately_suitable_threshold_bd = 90
marginally_suitable_threshold_bd = 50


highly_suitable_bd = np.where(bedrock_depth > highly_suitable_threshold_bd, 1, 0)
moderately_suitable_bd = np.where((bedrock_depth <= highly_suitable_threshold_bd) & (bedrock_depth > moderately_suitable_threshold_bd), 0.67, 0)
marginally_suitable_bd = np.where((bedrock_depth <= moderately_suitable_threshold_bd) & (bedrock_depth >= marginally_suitable_threshold_bd), 0.33, 0)
unsuitable_bd = np.where(bedrock_depth < marginally_suitable_threshold_bd, 0, 0)

combined_classes_bd = highly_suitable_bd + moderately_suitable_bd + marginally_suitable_bd + unsuitable_bd


#180, 200
values = bedrock_depth.values

#Load study area shapefile
study_area = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp') 

boundaries = [0, 0.25, 0.5, 0.75, 1]  # Adjust these boundaries as needed
labels = ['#FDE725', '#22A884', '#30678D', '#440154']  # Custom labels for each range
categories = ['Highly Suitable (0.75 - 1)', 'Moderately Suitable (0.5 - 0.75)', 'Marginally Suitable (0.25 - 0.5)', 'Unsuitable (0 - 0.25)']
cmap = plt.get_cmap('viridis', len(boundaries) - 1)

# Create a list of patches for the legend
patches = [mpatches.Patch(color=labels[i], label=categories[i]) for i in range(len(categories))]


# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, combined_classes_bd, cmap=cmap, levels=boundaries, norm=norm)
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
#cbar = plt.colorbar(contour, ax=ax, label='Soil Bedrock Depth Mean_0_20 Values (cm)', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Soil Bedrock Depth')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

# Displaying Histograph
plt.hist(values, bins=20, alpha=0.7)
plt.title('Distribution of Soil Bedrock Depth (cm)')
plt.xlabel('Soil Bedrock Depth Values (cm)')
plt.ylabel('Frequency')
plt.show()
