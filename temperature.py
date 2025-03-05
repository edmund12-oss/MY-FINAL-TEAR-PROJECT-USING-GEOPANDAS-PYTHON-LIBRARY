import xarray as xr
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata

# Load the NetCDF file
nc_file = r'C:\Users\USER\Desktop\project execution\geopandas\New thematic maps\temperature map\2001_2021_temp_max_DEA_2.nc'
ds = xr.open_dataset(nc_file)

#print(ds)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
temperature_2m = ds['air_temperature_at_2_metres']  # You might need to adjust variable names based on your NetCDF file structure
temperature_2m_2d = temperature_2m.isel(time=0)
lat = ds['lat']
lon = ds['lon']

#18.61502075, 32.98782349

values = temperature_2m_2d.values

#Load study area shapefile
study_area = gpd.read_file(r'C:\Users\USER\Desktop\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_admin.shp') 

boundaries = [23, 23.2, 23.4, 23.6, 23.8, 24, 24.2, 24.4, 24.6, 24.8, 25]  # Adjust these boundaries as needed
cmap = plt.get_cmap('viridis', len(boundaries) - 1)

# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, temperature_2m_2d, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.6)

# Add colorbar
cbar = plt.colorbar(contour, ax=ax, label='Temperature Values', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Temperature Thematic Map')
plt.show()

# Displaying Histograph
plt.hist(values, bins=20, alpha=0.7)
plt.title('Distribution of temperature Values')
plt.xlabel('Temperature Values')
plt.ylabel('Frequency')
plt.show()

