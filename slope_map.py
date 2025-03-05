import rasterio
import numpy as np
import matplotlib.pyplot as plt
from rasterio.enums import Resampling

# Replace 'your_dem_file.tif' with the path to your DEM file
dem_path = r'C:\Users\USER\Desktop\project execution\geopandas\Thematic maps\slope map\buikwe dem files\bwe_dem1.tif'
#NB: DEM has to be in a projected CS for the slope calculations

with rasterio.open(dem_path) as dem:
    # Read the DEM data, resampling as needed
    elevation = dem.read(1, resampling=Resampling.bilinear)

    # Get metadata for later use
    meta = dem.meta

# Calculate gradients
grad_x, grad_y = np.gradient(elevation, edge_order=2)

# Calculate slope in radians
slope_radians = np.arctan(np.sqrt(grad_x**2 + grad_y**2))

# Convert to degrees
slope_degrees = np.degrees(slope_radians)

# Update the metadata for the slope dataset
meta.update(dtype=rasterio.float32, count=1)

# Write the slope data to a new file
slope_path = 'slope_map.tif'
with rasterio.open(slope_path, 'w', **meta) as dst:
    dst.write(slope_degrees.astype(rasterio.float32), 1)

#print(slope_degrees)

plt.imshow(slope_degrees, cmap='terrain')
plt.colorbar(label='Slope (degrees)')
plt.title('Slope Map')
plt.show()
