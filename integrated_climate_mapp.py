import xarray as xr
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata

#RAINFALL MAP

# Load your rainfall data from a NetCDF file
rainfall_data_path = r'C:\Users\hp\Desktop\project execution\geopandas\New thematic maps\rainfall_map\2003_2023_rf_data.nc'
rainfall_data = xr.open_dataset(rainfall_data_path)

# Assuming your rainfall data has lat and lon dimensions, and a rainfall variable
# Adjust 'rainfall' below to match the variable name in your dataset
rainfall_values_2d = rainfall_data['rainfall']  # Replace 'rainfall' with your variable name
lat = rainfall_data['y']
lon = rainfall_data['x']

# If you need to aggregate the rainfall data (e.g., mean over time), do so here
# For example, to get the mean across the time dimension if it exists:
rainfall_sum = rainfall_values_2d.sum(dim='time', skipna=True)
rainfall_2003_2023 = rainfall_sum / 20

highly_suitable_threshold_rf = 1500
moderately_suitable_threshold_rf = 1250
marginally_suitable_threshold_rf = 1000


highly_suitable_rf = np.where(rainfall_2003_2023 >= highly_suitable_threshold_rf, 1, 0)
moderately_suitable_rf = np.where((rainfall_2003_2023 < highly_suitable_threshold_rf) & (rainfall_2003_2023 >= moderately_suitable_threshold_rf), 0.67, 0)
marginally_suitable_rf = np.where((rainfall_2003_2023 < moderately_suitable_threshold_rf) & (rainfall_2003_2023 >= marginally_suitable_threshold_rf), 0.33, 0)
unsuitable_rf = np.where(rainfall_2003_2023 < marginally_suitable_threshold_rf, 0, 0)

combined_classes_rf = highly_suitable_rf + moderately_suitable_rf + marginally_suitable_rf + unsuitable_rf

#MAXIMUM TEMPERATURE

# Load the NetCDF file
nc_file_mxtemp = r'C:\Users\hp\Desktop\project execution\geopandas\New thematic maps\maximum temperature map\final_max_temperature.nc'
ds_mxtemp = xr.open_dataset(nc_file_mxtemp)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
temperature_2m = ds_mxtemp['air_temperature_at_2_metres']# You might need to adjust variable names based on your NetCDF file structure
temperature_2m_2d = temperature_2m.isel(time=0)
#lat = ds_mxtemp['lat']
#lon = ds_mxtemp['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_mxtemp_1, highly_suitable_threshold_mxtemp_2 = 29, 32
moderately_suitable_threshold_mxtemp = 21
marginally_suitable_threshold_mxtemp = 18


highly_suitable_mxtemp = np.where((temperature_2m_2d >= highly_suitable_threshold_mxtemp_1) & (temperature_2m_2d <= highly_suitable_threshold_mxtemp_2), 1, 0)
moderately_suitable_mxtemp = np.where((temperature_2m_2d >= moderately_suitable_threshold_mxtemp) & (temperature_2m_2d < highly_suitable_threshold_mxtemp_1), 0.67, 0)
marginally_suitable_mxtemp = np.where((temperature_2m_2d < moderately_suitable_threshold_mxtemp) & (temperature_2m_2d >= marginally_suitable_threshold_mxtemp), 0.33, 0)
unsuitable_mxtemp = np.where(temperature_2m_2d < marginally_suitable_threshold_mxtemp, 0, 0)

combined_classes_mxtemp = highly_suitable_mxtemp + moderately_suitable_mxtemp + marginally_suitable_mxtemp + unsuitable_mxtemp

#MINIMUM TEMPERATURE

# Load the NetCDF file
nc_file_mintemp = r'C:\Users\hp\Desktop\project execution\geopandas\New thematic maps\minimum temperature map\final_min_temperature.nc'
ds_mintemp = xr.open_dataset(nc_file_mintemp)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
min_temperature_2m = ds_mintemp['air_temperature_at_2_metres']# You might need to adjust variable names based on your NetCDF file structure
min_temperature_2m_2d = min_temperature_2m.isel(time=0)
#lat = ds_mintemp['lat']
#lon = ds_mintemp['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_mintemp_1, highly_suitable_threshold_mintemp_2 = 29, 32
moderately_suitable_threshold_mintemp = 21
marginally_suitable_threshold_mintemp = 18


highly_suitable_mintemp = np.where((min_temperature_2m_2d >= highly_suitable_threshold_mintemp_1) & (min_temperature_2m_2d <= highly_suitable_threshold_mintemp_2), 1, 0)
moderately_suitable_mintemp = np.where((min_temperature_2m_2d >= moderately_suitable_threshold_mintemp) & (min_temperature_2m_2d < highly_suitable_threshold_mintemp_1), 0.67, 0)
marginally_suitable_mintemp = np.where((min_temperature_2m_2d < moderately_suitable_threshold_mintemp) & (min_temperature_2m_2d >= marginally_suitable_threshold_mintemp), 0.33, 0)
unsuitable_mintemp = np.where(min_temperature_2m_2d < marginally_suitable_threshold_mintemp, 0, 0)

combined_classes_mintemp = highly_suitable_mintemp + moderately_suitable_mintemp + marginally_suitable_mintemp + unsuitable_mintemp

#WIND SPEED

# Load the NetCDF file
nc_file_windspeed = r'C:\Users\hp\Desktop\project execution\geopandas\New thematic maps\windspeed map\wind_speed_new.nc'
ds_windspeed = xr.open_dataset(nc_file_windspeed)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
wind_speed = ds_windspeed['__xarray_dataarray_variable__']# You might need to adjust variable names based on your NetCDF file structure
wind_speed_2d = wind_speed.isel(time=0)
#lat = ds_windspeed['lat']
#lon = ds_windspeed['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_windspeed_1, highly_suitable_threshold_windspeed_2 = 0, 2
moderately_suitable_threshold_windspeed = 4
marginally_suitable_threshold_windspeed = 6


highly_suitable_windspeed = np.where((wind_speed_2d >= highly_suitable_threshold_mintemp_1) & (wind_speed_2d <= highly_suitable_threshold_mintemp_2), 1, 0)
moderately_suitable_windspeed = np.where((wind_speed_2d <= moderately_suitable_threshold_mintemp) & (wind_speed_2d > highly_suitable_threshold_mintemp_2), 0.67, 0)
marginally_suitable_windspeed = np.where((wind_speed_2d > moderately_suitable_threshold_mintemp) & (wind_speed_2d <= marginally_suitable_threshold_mintemp), 0.33, 0)
unsuitable_windspeed = np.where(wind_speed_2d > marginally_suitable_threshold_mintemp, 0, 0)

combined_classes_windspeed = highly_suitable_windspeed + moderately_suitable_windspeed + marginally_suitable_windspeed + unsuitable_windspeed


#tentative weight till I use AHP to determine actual weights
weight_rf = 0.35445804
weight_mxtemp = 0.58391608
weight_mintemp = 0.58391608
weight_windspeed = 0.06162587


climate_map_1 = (weight_rf * combined_classes_rf) 
climate_map_2 = (weight_mxtemp * combined_classes_mxtemp) + (weight_mintemp * combined_classes_mintemp) + (weight_windspeed * combined_classes_windspeed)

#changing the shape of an array
new_array = np.zeros_like(climate_map_1)

new_array[:climate_map_2.shape[0], :climate_map_2.shape[1]] = climate_map_2

final_climate_data = climate_map_1 + new_array



#Load study area shapefile
study_area = gpd.read_file(r'C:\Users\hp\Desktop\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_admin.shp') 


boundaries = [0, 0.25, 0.5, 0.75, 1]  # Adjust these boundaries as needed
cmap = plt.get_cmap('viridis', len(boundaries) - 1)

# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, final_climate_data, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.6)

# Add colorbar
cbar = plt.colorbar(contour, ax=ax, label='Climate Values', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Climate')
plt.show()
