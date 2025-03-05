import geopandas as gpd
import rasterio
from rasterio.enums import Resampling
import numpy as np
import matplotlib.pyplot as plt

dem_path = 'dem.tif'

with rasterio.open(dem_path) as dem:
#Reading the dem data and resampling
elevation_data = dem.read(1, resampling=Resampling.bilinear)

#Getting metadata for later use
meta = dem.meta

#Calculate gradients
grad_x, grad_y = np.gradient(elevation_data, edge_order=2)

#Calculate Slope in radians
slope_radians = np.arctan(np.sqrt(grad_x**2 + grad_y**2))

#Convert to degrees
slope_degrees = np.degrees(slope_radians)

plt.imshow(slope_degress, cmap='terrain')
plt.colorbar(label='Slope (degrees)')
plt.title('Slope Map')
plt.show()