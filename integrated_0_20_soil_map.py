#SOIL PH

import geopandas as gpd
import rasterio
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject
import pyproj
from affine import Affine
import matplotlib.patches as mpatches
from rasterio.transform import from_origin

# Load the NetCDF file
nc_file_ph = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil ph(0_20)cm\masked_ph.nc'
ds_ph = xr.open_dataset(nc_file_ph)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
ph = ds_ph['pH, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_ph['lat']
lon = ds_ph['lon']

#4.69999981, 6.9000001
#Reclassifying according to the suitability levels
highly_suitable_threshold_ph = 6
moderately_suitable_threshold_ph = 5
#marginally_suitable_threshold_ph = 0


highly_suitable_ph = np.where(ph >= highly_suitable_threshold_ph, 1, 0)
moderately_suitable_ph = np.where((ph < highly_suitable_threshold_ph) & (ph >= moderately_suitable_threshold_ph), 0.67, 0)
#marginally_suitable_ph = np.where((grid_temperature < moderately_suitable_threshold) & (grid_temperature >= marginally_suitable_threshold), 0.33, 0)
unsuitable_ph = np.where(ph < moderately_suitable_threshold_ph, 0, 0)

combined_classes_ph = highly_suitable_ph + moderately_suitable_ph + unsuitable_ph


#SOIL ALUMINIUM EXTRACTABLE*

# Load the NetCDF file
nc_file_al = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil aluminium map(0_20)cm\masked_al.nc'
ds_al = xr.open_dataset(nc_file_al)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
aluminium_extractable = ds_al['masked_al']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_al['lat']
lon = ds_al['lon']

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

#SOIL DEPTH TO BEDROCK

# Load the NetCDF file
nc_file_dtb = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil bedrock depth(0_20)cm\masked_depth.nc'
ds_dtb = xr.open_dataset(nc_file_dtb)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
bedrock_depth_2d = ds_dtb['mean_0_200']  # You might need to adjust variable names based on your NetCDF file structure
#bedrock_depth_2d =bedrock_depth.isel(time=0)
#lat = ds_dtb['lat']
#lon = ds_dtb['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_bd = 120
moderately_suitable_threshold_bd = 90
marginally_suitable_threshold_bd = 50


highly_suitable_bd = np.where(bedrock_depth_2d > highly_suitable_threshold_bd, 1, 0)
moderately_suitable_bd = np.where((bedrock_depth_2d <= highly_suitable_threshold_bd) & (bedrock_depth_2d > moderately_suitable_threshold_bd), 0.67, 0)
marginally_suitable_bd = np.where((bedrock_depth_2d <= moderately_suitable_threshold_bd) & (bedrock_depth_2d >= marginally_suitable_threshold_bd), 0.33, 0)
unsuitable_bd = np.where(bedrock_depth_2d < marginally_suitable_threshold_bd, 0, 0)

combined_classes_bd = highly_suitable_bd + moderately_suitable_bd + marginally_suitable_bd + unsuitable_bd


#SOIL BULK DENSITY

# Load the NetCDF file
nc_file_bdensity = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil bulk density(0_20)cm\masked_bulk_density.nc'
ds_bdensity = xr.open_dataset(nc_file_bdensity)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
bulk_density_2d = ds_bdensity['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
#bulk_density_2d =bulk_density.isel(time=0)
#lat = ds_bdensity['lat']
#lon = ds_bdensity['lon']

#1, 1.5

#Reclassifying according to the suitability levels
highly_suitable_threshold_bdensity = 1.2
moderately_suitable_threshold_bdensity = 1.4
marginally_suitable_threshold_bdensity = 1.5


highly_suitable_bdensity = np.where(bulk_density_2d <= highly_suitable_threshold_bdensity, 1, 0)
moderately_suitable_bdensity = np.where((bulk_density_2d > highly_suitable_threshold_bdensity) & (bulk_density_2d <= moderately_suitable_threshold_bdensity), 0.67, 0)
marginally_suitable_bdensity = np.where((bulk_density_2d > moderately_suitable_threshold_bdensity) & (bulk_density_2d <= marginally_suitable_threshold_bdensity), 0.33, 0)
unsuitable_bdensity = np.where(bulk_density_2d > marginally_suitable_threshold_bdensity, 0, 0)

combined_classes_bdensity = highly_suitable_bdensity + moderately_suitable_bdensity + marginally_suitable_bdensity + unsuitable_bdensity

#SOIL NITROGEN

# Load the NetCDF file
nc_file_nitro = r'geopandas/New thematic maps/soil_maps(0_20)cm/soil nitrogen total(0_20)cm/masked_nitrogen_total.nc'
ds_nitro = xr.open_dataset(nc_file_nitro)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
nitrogen_total = ds_nitro['Nitrogen, total, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_nitro['lat']
lon = ds_nitro['lon']

nitrogen_total_final = nitrogen_total / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_nitrogen_1, highly_suitable_threshold_nitrogen_2 = 1.5, 2
moderately_suitable_threshold_nitrogen = 1
marginally_suitable_threshold_nitrogen = 0.5


highly_suitable_nitrogen = np.where((nitrogen_total_final >= highly_suitable_threshold_nitrogen_1) & (nitrogen_total_final <= highly_suitable_threshold_nitrogen_2), 1, 0)
moderately_suitable_nitrogen = np.where((nitrogen_total_final >= moderately_suitable_threshold_nitrogen) & (nitrogen_total_final < highly_suitable_threshold_nitrogen_1), 0.67, 0)
marginally_suitable_nitrogen = np.where((nitrogen_total_final < moderately_suitable_threshold_nitrogen) & (nitrogen_total_final >= marginally_suitable_threshold_nitrogen), 0.33, 0)
unsuitable_nitrogen = np.where(nitrogen_total_final < marginally_suitable_threshold_nitrogen, 0, 0)

combined_classes_nitrogen = highly_suitable_nitrogen + moderately_suitable_nitrogen + marginally_suitable_nitrogen + unsuitable_nitrogen

#SOIL CALCIUM

# Load the NetCDF file
nc_file_ca = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil calcium extractable(0_20)cm\masked_calcium.nc'
ds_ca = xr.open_dataset(nc_file_ca)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
calcium_extractable = ds_ca['Calcium, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_ca['lat']
lon = ds_ca['lon']

#163.02189636, 2207.34765625
#converting to g/kg
calcium_extractable_final = calcium_extractable / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_calcium_1, highly_suitable_threshold_calcium_2 = 0.6, 1.2
moderately_suitable_threshold_calcium_1, moderately_suitable_threshold_calcium_2 = 0.3, 1.5
marginally_suitable_threshold_calcium_1, marginally_suitable_threshold_calcium_2 = 0.15, 2


highly_suitable_calcium = np.where((calcium_extractable_final >= highly_suitable_threshold_calcium_1) & (calcium_extractable_final <= highly_suitable_threshold_calcium_2), 1, 0)
moderately_suitable_calcium = np.where((calcium_extractable_final >= moderately_suitable_threshold_calcium_1) & (calcium_extractable_final < highly_suitable_threshold_calcium_1) & (calcium_extractable_final > highly_suitable_threshold_calcium_2) & (calcium_extractable_final <= moderately_suitable_threshold_calcium_2), 0.67, 0)
marginally_suitable_calcium = np.where((calcium_extractable_final < moderately_suitable_threshold_calcium_1) & (calcium_extractable_final >= marginally_suitable_threshold_calcium_1) & (calcium_extractable_final > moderately_suitable_threshold_calcium_2) & (calcium_extractable_final <= marginally_suitable_threshold_calcium_2), 0.33, 0)
unsuitable_calcium = np.where((calcium_extractable_final < marginally_suitable_threshold_calcium_1) & (calcium_extractable_final > marginally_suitable_threshold_calcium_2), 0, 0)

combined_classes_calcium = highly_suitable_calcium + moderately_suitable_calcium + marginally_suitable_calcium + unsuitable_calcium

#SOIL CARBON ORGANIC

# Load the NetCDF file
nc_file_co = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil carbon organic(0_20)cm\masked_carbon_organic.nc'
ds_co = xr.open_dataset(nc_file_co)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
carbon_organic = ds_co['Carbon, organic, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_co['lat']
lon = ds_co['lon']

#max_value = carbon_organic.max()
#print(max_value)

#8.02501392, 23.53253174

#Converting from g/kg to a percentage

carbon_organic_percentage = (carbon_organic / 1000) * 100

#Reclassifying according to the suitability levels
highly_suitable_threshold_co = 3
moderately_suitable_threshold_co = 2
marginally_suitable_threshold_co = 1


highly_suitable_co = np.where(carbon_organic_percentage > highly_suitable_threshold_co, 1, 0)
moderately_suitable_co = np.where((carbon_organic_percentage <= highly_suitable_threshold_co) & (carbon_organic_percentage >= moderately_suitable_threshold_co), 0.67, 0)
marginally_suitable_co = np.where((carbon_organic_percentage >= marginally_suitable_threshold_co) & (carbon_organic_percentage < moderately_suitable_threshold_co), 0.33, 0)
unsuitable_co = np.where(carbon_organic_percentage < marginally_suitable_threshold_co, 0, 0)

combined_classes_co = highly_suitable_co + moderately_suitable_co + marginally_suitable_co + unsuitable_co

#SOIL CARBON TOTAL

# Load the NetCDF file
nc_file_ct = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil carbon total(0_20)cm\masked_carbon_total.nc'
ds_ct = xr.open_dataset(nc_file_ct)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
carbon_total_2d = ds_ct['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
#carbon_total_2d = carbon_total.isel(time=0)
#lat = ds_ct['lat']
#lon = ds_ct['lon']

#converting from g/kg to percentage 
carbon_total_percentage = carbon_total_2d / 10

#Reclassifying according to the suitability levels
highly_suitable_threshold_ct = 2
moderately_suitable_threshold_ct = 1.2
marginally_suitable_threshold_ct = 0.6


highly_suitable_ct = np.where(carbon_total_percentage > highly_suitable_threshold_ct, 1, 0)
moderately_suitable_ct = np.where((carbon_total_percentage <= highly_suitable_threshold_ct) & (carbon_total_percentage >= moderately_suitable_threshold_ct), 0.67, 0)
marginally_suitable_ct = np.where((carbon_total_percentage >= marginally_suitable_threshold_ct) & (carbon_total_percentage < moderately_suitable_threshold_ct), 0.33, 0)
unsuitable_ct = np.where(carbon_total_percentage < marginally_suitable_threshold_ct, 0, 0)

combined_classes_ct = highly_suitable_ct + moderately_suitable_ct + marginally_suitable_ct + unsuitable_ct


#SOIL CATION EXCHANGE CAPACITY

nc_file_cec = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil cation exchange capacity(0_20)cm\masked_cec.nc'
ds_cec = xr.open_dataset(nc_file_cec)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
cec = ds_cec['Effective Cation Exchange Capacity, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_cec['lat']
lon = ds_cec['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_cec = 20
moderately_suitable_threshold_cec = 10
marginally_suitable_threshold_cec = 5


highly_suitable_cec = np.where(cec > highly_suitable_threshold_cec, 1, 0)
moderately_suitable_cec = np.where((cec <= highly_suitable_threshold_cec) & (cec >= moderately_suitable_threshold_cec), 0.67, 0)
marginally_suitable_cec = np.where((cec >= marginally_suitable_threshold_cec) & (cec < moderately_suitable_threshold_cec), 0.33, 0)
unsuitable_cec = np.where(cec < marginally_suitable_threshold_cec, 0, 0)

combined_classes_cec = highly_suitable_cec + moderately_suitable_cec + marginally_suitable_cec + unsuitable_cec

#SOIL CLAY CONTENT

# Load the NetCDF file
nc_file_clay = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil clay content(0_20)cm\masked_clay_content.nc'
ds_clay = xr.open_dataset(nc_file_clay)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
clay_content_2d = ds_clay['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
#clay_content_2d =clay_content.isel(time=0)
#lat = ds_clay['lat']
#lon = ds_clay['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_clay_1, highly_suitable_threshold_clay_2 = 20, 35
moderately_suitable_threshold_clay = 45
marginally_suitable_threshold_clay_1, marginally_suitable_threshold_clay_2 = 15, 60


highly_suitable_clay = np.where((clay_content_2d >= highly_suitable_threshold_clay_1) & (clay_content_2d <= highly_suitable_threshold_clay_2), 1, 0)
moderately_suitable_clay = np.where((clay_content_2d > highly_suitable_threshold_clay_2) & (clay_content_2d <= moderately_suitable_threshold_clay), 0.67, 0)
marginally_suitable_clay = np.where((clay_content_2d >= marginally_suitable_threshold_clay_1) & (clay_content_2d < highly_suitable_threshold_clay_1) & (clay_content_2d > moderately_suitable_threshold_clay) & (clay_content_2d <= marginally_suitable_threshold_clay_2), 0.33, 0)
unsuitable_clay = np.where((clay_content_2d < marginally_suitable_threshold_clay_1) & (clay_content_2d > marginally_suitable_threshold_clay_2), 0, 0)

combined_classes_clay = highly_suitable_clay + moderately_suitable_clay + marginally_suitable_clay + unsuitable_clay

#SOIL IRON EXTRACTABLE

# Load the NetCDF file
nc_file_fe = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil iron extractable(0_20)cm\masked_iron.nc'
ds_fe = xr.open_dataset(nc_file_fe)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
iron_extractable = ds_fe['Iron, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_fe['lat']
lon = ds_fe['lon']

#converting from mg/kg to g/kg
iron_extractable_final = iron_extractable / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_iron_1, highly_suitable_threshold_iron_2 = 0.1, 0.3
moderately_suitable_threshold_iron_1, moderately_suitable_threshold_iron_2 = 0.05, 0.5
marginally_suitable_threshold_iron_1, marginally_suitable_threshold_iron_2 = 0.02, 0.8


highly_suitable_iron = np.where((iron_extractable_final >= highly_suitable_threshold_iron_1) & (iron_extractable_final <= highly_suitable_threshold_iron_2), 1, 0)
moderately_suitable_iron = np.where((iron_extractable_final >= moderately_suitable_threshold_iron_1) & (iron_extractable_final < highly_suitable_threshold_iron_1) & (iron_extractable_final > highly_suitable_threshold_iron_2) & (iron_extractable_final <= moderately_suitable_threshold_iron_2), 0.67, 0)
marginally_suitable_iron = np.where((iron_extractable_final < moderately_suitable_threshold_iron_1) & (iron_extractable_final >= marginally_suitable_threshold_iron_1) & (iron_extractable_final > moderately_suitable_threshold_iron_2) & (iron_extractable_final <= marginally_suitable_threshold_iron_2), 0.33, 0)
unsuitable_iron = np.where((iron_extractable_final < marginally_suitable_threshold_iron_1) & (iron_extractable_final > marginally_suitable_threshold_iron_2), 0, 0)

combined_classes_iron = highly_suitable_iron + moderately_suitable_iron + marginally_suitable_iron + unsuitable_iron

# SOIL MAGNESIUM

# Load the NetCDF file
nc_file_mg = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil magnesium extractable(0_20)cm\masked_magnesium.nc'
ds_mg = xr.open_dataset(nc_file_mg)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
magnesium_extractable = ds_mg['Magnesium, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_mg['lat']
lon = ds_mg['lon']

#Converting from mg/kg to g/kg
magnesium_extractable_final = magnesium_extractable / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_mg_1, highly_suitable_threshold_mg_2 = 0.2, 0.4
moderately_suitable_threshold_mg_1, moderately_suitable_threshold_mg_2 = 0.1, 0.6
marginally_suitable_threshold_mg_1, marginally_suitable_threshold_mg_2 = 0.05, 0.8


highly_suitable_mg = np.where((magnesium_extractable_final >= highly_suitable_threshold_mg_1) & (magnesium_extractable_final <= highly_suitable_threshold_mg_2), 1, 0)
moderately_suitable_mg = np.where((magnesium_extractable_final >= moderately_suitable_threshold_mg_1) & (magnesium_extractable_final < highly_suitable_threshold_mg_1) & (magnesium_extractable_final > highly_suitable_threshold_mg_2) & (magnesium_extractable_final <= moderately_suitable_threshold_mg_2), 0.67, 0)
marginally_suitable_mg = np.where((magnesium_extractable_final < moderately_suitable_threshold_mg_1) & (magnesium_extractable_final >= marginally_suitable_threshold_mg_1) & (magnesium_extractable_final > moderately_suitable_threshold_mg_2) & (magnesium_extractable_final <= marginally_suitable_threshold_mg_2), 0.33, 0)
unsuitable_mg = np.where((magnesium_extractable_final < marginally_suitable_threshold_mg_1) & (magnesium_extractable_final > marginally_suitable_threshold_mg_2), 0, 0)

combined_classes_mg = highly_suitable_mg + moderately_suitable_mg + marginally_suitable_mg + unsuitable_mg

#SOIL PHOSPHOROUS

# Load the NetCDF file
nc_file_phosphorus = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil phosphorous extractable(0_20)cm\masked_phosphorus.nc'
ds_phosphorus = xr.open_dataset(nc_file_phosphorus)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
phosphorus_extractable = ds_phosphorus['Phosphorus, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_phosphorus['lat']
lon = ds_phosphorus['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_phosphorus_1, highly_suitable_threshold_phosphorus_2 = 20, 50
moderately_suitable_threshold_phosphorus = 10
marginally_suitable_threshold_phosphorus = 5


highly_suitable_phosphorus = np.where((phosphorus_extractable >= highly_suitable_threshold_phosphorus_1) & (phosphorus_extractable <= highly_suitable_threshold_phosphorus_2), 1, 0)
moderately_suitable_phosphorus = np.where((phosphorus_extractable >= moderately_suitable_threshold_phosphorus) & (phosphorus_extractable < highly_suitable_threshold_phosphorus_1), 0.67, 0)
marginally_suitable_phosphorus = np.where((phosphorus_extractable < moderately_suitable_threshold_phosphorus) & (phosphorus_extractable >= marginally_suitable_threshold_phosphorus), 0.33, 0)
unsuitable_phosphorus = np.where(phosphorus_extractable < marginally_suitable_threshold_phosphorus, 0, 0)

combined_classes_phosphorus = highly_suitable_phosphorus + moderately_suitable_phosphorus + marginally_suitable_phosphorus + unsuitable_phosphorus

#SOIL POTASSIUM

# Load the NetCDF file
nc_file_k = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil potassium extractable(0_20)cm\masked_potassium.nc'
ds_k = xr.open_dataset(nc_file_k)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
potassium_extractable = ds_k['Potassium, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_k['lat']
lon = ds_k['lon']

#Converting from mg/kg to g/kg
potassium_extractable_final = potassium_extractable / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_k_1, highly_suitable_threshold_k_2 = 0.15, 0.25
moderately_suitable_threshold_k = 0.10
marginally_suitable_threshold_k = 0.05


highly_suitable_k = np.where((potassium_extractable_final >= highly_suitable_threshold_k_1) & (potassium_extractable_final <= highly_suitable_threshold_k_2), 1, 0)
moderately_suitable_k = np.where((potassium_extractable_final >= moderately_suitable_threshold_k) & (potassium_extractable_final < highly_suitable_threshold_k_1), 0.67, 0)
marginally_suitable_k = np.where((potassium_extractable_final < moderately_suitable_threshold_k) & (potassium_extractable_final >= marginally_suitable_threshold_k), 0.33, 0)
unsuitable_k = np.where(potassium_extractable_final < marginally_suitable_threshold_k, 0, 0)

combined_classes_k = highly_suitable_k + moderately_suitable_k + marginally_suitable_k + unsuitable_k

#SOIL SAND CONTENT

# Load the NetCDF file
nc_file_sand = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil sand content(0_20)cm\masked_sand.nc'
ds_sand = xr.open_dataset(nc_file_sand)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
sand_content_2d = ds_sand['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
#sand_content_2d = sand_content.isel(time=0)
#lat = ds_sand['lat']
#lon = ds_sand['lon']


#Reclassifying according to the suitability levels
highly_suitable_threshold_sand_1, highly_suitable_threshold_sand_2 = 20, 40
moderately_suitable_threshold_sand  = 60
marginally_suitable_threshold_sand_1, marginally_suitable_threshold_sand_2 = 15, 70


highly_suitable_sand = np.where((sand_content_2d >= highly_suitable_threshold_sand_1) & (sand_content_2d <= highly_suitable_threshold_sand_2), 1, 0)
moderately_suitable_sand  = np.where((sand_content_2d > highly_suitable_threshold_sand_2) & (sand_content_2d <= moderately_suitable_threshold_sand), 0.67, 0)
marginally_suitable_sand  = np.where((sand_content_2d >= marginally_suitable_threshold_sand_1) & (sand_content_2d < highly_suitable_threshold_sand_1) & (sand_content_2d > moderately_suitable_threshold_sand) & (sand_content_2d <= marginally_suitable_threshold_sand_2), 0.33, 0)
unsuitable_sand  = np.where((sand_content_2d < marginally_suitable_threshold_sand_1) & (sand_content_2d > marginally_suitable_threshold_sand_2), 0, 0)

combined_classes_sand = highly_suitable_sand + moderately_suitable_sand + marginally_suitable_sand + unsuitable_sand

#SOIL SILT CONTENT 

# Load the NetCDF file
nc_file_silt = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil silt content(0_20)cm\masked_silt.nc'
ds_silt = xr.open_dataset(nc_file_silt)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variables
silt_content_2d = ds_silt['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
#silt_content_2d = silt_content.isel(time=0)
#lat = ds_silt['lat']
#lon = ds_silt['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_silt_1, highly_suitable_threshold_silt_2 = 20, 40
moderately_suitable_threshold_silt = 15
marginally_suitable_threshold_silt = 10


highly_suitable_silt = np.where((silt_content_2d >= highly_suitable_threshold_silt_1) & (silt_content_2d <= highly_suitable_threshold_silt_2), 1, 0)
moderately_suitable_silt = np.where((silt_content_2d >= moderately_suitable_threshold_silt) & (silt_content_2d < highly_suitable_threshold_silt_1), 0.67, 0)
marginally_suitable_silt = np.where((silt_content_2d < moderately_suitable_threshold_silt) & (silt_content_2d >= marginally_suitable_threshold_silt), 0.33, 0)
unsuitable_silt = np.where(silt_content_2d < marginally_suitable_threshold_silt, 0, 0)

combined_classes_silt = highly_suitable_silt + moderately_suitable_silt + marginally_suitable_silt + unsuitable_silt

#SOIL STONE CONTENT

# Load the NetCDF file
nc_file_stone = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil stone content(0_20)cm\masked_stone_content.nc'
ds_stone = xr.open_dataset(nc_file_stone)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
stone_content = ds_stone['Stone content, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_stone['lat']
lon = ds_stone['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_stone = 10
moderately_suitable_threshold_stone = 20
marginally_suitable_threshold_stone = 35


highly_suitable_stone = np.where(stone_content < highly_suitable_threshold_stone, 1, 0)
moderately_suitable_stone = np.where((stone_content >= highly_suitable_threshold_stone) & (stone_content <= moderately_suitable_threshold_stone), 0.67, 0)
marginally_suitable_stone = np.where((stone_content <= marginally_suitable_threshold_stone) & (stone_content > moderately_suitable_threshold_stone), 0.33, 0)
unsuitable_stone = np.where(stone_content > marginally_suitable_threshold_stone, 0, 0)

combined_classes_stone = highly_suitable_stone + moderately_suitable_stone + marginally_suitable_stone + unsuitable_stone

#SOIL SULPHUR EXTRACTABLE

# Load the NetCDF file
nc_file_sulphur = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil sulphur extractable(0_20)cm\masked_sulphur.nc'
ds_sulphur = xr.open_dataset(nc_file_sulphur)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
sulphur_extractable = ds_sulphur['Sulphur, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_sulphur['lat']
lon = ds_sulphur['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_sulphur_1, highly_suitable_threshold_sulphur_2 = 10, 30
moderately_suitable_threshold_sulphur = 5
marginally_suitable_threshold_sulphur = 3


highly_suitable_sulphur = np.where((sulphur_extractable >= highly_suitable_threshold_sulphur_1) & (sulphur_extractable <= highly_suitable_threshold_sulphur_2), 1, 0)
moderately_suitable_sulphur = np.where((sulphur_extractable >= moderately_suitable_threshold_sulphur) & (sulphur_extractable < highly_suitable_threshold_sulphur_1), 0.67, 0)
marginally_suitable_sulphur = np.where((sulphur_extractable < moderately_suitable_threshold_sulphur) & (sulphur_extractable >= marginally_suitable_threshold_sulphur), 0.33, 0)
unsuitable_sulphur = np.where(sulphur_extractable < marginally_suitable_threshold_sulphur, 0, 0)

combined_classes_sulphur = highly_suitable_sulphur + moderately_suitable_sulphur + marginally_suitable_sulphur + unsuitable_sulphur

#SOIL ZINC EXTRACTABLE

# Load the NetCDF file
nc_file_zinc = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil zinc extractable(0_20)cm\masked_zinc.nc'
ds_zinc = xr.open_dataset(nc_file_zinc)

# Assuming the dataset contains 'bulk_density', 'lat', and 'lon' variable
zinc_extractable = ds_zinc['Zinc, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_zinc['lat']
lon = ds_zinc['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_zinc_1, highly_suitable_threshold_zinc_2 = 1.5, 10
moderately_suitable_threshold_zinc = 1
marginally_suitable_threshold_zinc = 0.5


highly_suitable_zinc = np.where((zinc_extractable >= highly_suitable_threshold_zinc_1) & (zinc_extractable <= highly_suitable_threshold_zinc_2), 1, 0)
moderately_suitable_zinc = np.where((zinc_extractable >= moderately_suitable_threshold_zinc) & (zinc_extractable < highly_suitable_threshold_zinc_1), 0.67, 0)
marginally_suitable_zinc = np.where((zinc_extractable < moderately_suitable_threshold_zinc) & (zinc_extractable >= marginally_suitable_threshold_zinc), 0.33, 0)
unsuitable_zinc = np.where(zinc_extractable < marginally_suitable_threshold_zinc, 0, 0)

combined_classes_zinc = highly_suitable_zinc + moderately_suitable_zinc + marginally_suitable_zinc + unsuitable_zinc


#WEIGHTED OVERLAY OF SOIL MAPS(0_20)

#tentative weight till I use AHP to determine actual weights
weight = 0.0526315789473684


#changing shape of climate data to (5479, 4574)
combined_classes_silt_2 = np.zeros_like(combined_classes_zinc)

combined_classes_silt_2[:combined_classes_silt.shape[0], :combined_classes_silt.shape[1]] = combined_classes_silt 


#changing shape of climate data to (5479, 4574)
combined_classes_sand_2 = np.zeros_like(combined_classes_zinc)

combined_classes_sand_2[:combined_classes_sand.shape[0], :combined_classes_sand.shape[1]] = combined_classes_sand 


#changing shape of climate data to (5479, 4574)
combined_classes_ct_2 = np.zeros_like(combined_classes_zinc)

combined_classes_ct_2[:combined_classes_ct.shape[0], :combined_classes_ct.shape[1]] = combined_classes_ct 


#changing shape of climate data to (5479, 4574)
combined_classes_bd_2 = np.zeros_like(combined_classes_zinc)

combined_classes_bd_2[:combined_classes_bd.shape[0], :combined_classes_bd.shape[1]] = combined_classes_bd 



#changing shape of climate data to (5479, 4574)
combined_classes_bdensity_2 = np.zeros_like(combined_classes_zinc)

combined_classes_bdensity_2[:combined_classes_bdensity.shape[0], :combined_classes_bdensity.shape[1]] = combined_classes_bdensity 





#changing shape of climate data to (5479, 4574)
combined_classes_clay_2 = np.zeros_like(combined_classes_zinc)

combined_classes_clay_2[:combined_classes_clay.shape[0], :combined_classes_clay.shape[1]] = combined_classes_clay 






soil_map_0_20 = (weight * combined_classes_ph) + (weight * combined_classes_al) + (weight * combined_classes_nitrogen) + (weight * combined_classes_calcium) + (weight * combined_classes_co) + (weight * combined_classes_cec) + (weight * combined_classes_iron) + (weight * combined_classes_mg) + (weight * combined_classes_phosphorus) + (weight * combined_classes_k) + (weight * combined_classes_stone) + (weight * combined_classes_sulphur) + (weight * combined_classes_zinc) + (weight * combined_classes_silt_2) + (weight * combined_classes_sand_2) + (weight * combined_classes_zinc) + (weight * combined_classes_clay_2) + (weight * combined_classes_ct_2) + (weight * combined_classes_bdensity_2)+ ( + (weight * combined_classes_bd_2))

#Load study area shapefile
study_area = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp') 

#in mg/kg
boundaries = [0, 0.25, 0.5, 0.75, 1]  # Adjust these boundaries as needed
labels = ['#FDE725', '#22A884', '#30678D', '#440154']  # Custom labels for each range
categories = ['Highly Suitable', 'Moderately Suitable', 'Marginally Suitable', 'Unsuitable']
cmap = plt.get_cmap('viridis', len(boundaries) - 1)


#Create a list  of patches for the legend
patches = [mpatches.Patch(color=labels[i], label=categories[i]) for i in range(len(categories))]

# Create a BoundaryNorm instance
norm = mcolors.BoundaryNorm(boundaries, cmap.N, clip=True)

# Calculate midpoints of boundaries for tick placement
ticks = [(boundaries[i] + boundaries[i+1])/2 for i in range(len(boundaries)-1)]
# Plotting the results using matplotlib
fig, ax = plt.subplots(figsize=(12, 8), alpha=0.2)
contour = ax.contourf(lon, lat, soil_map_0_20, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.6)

# Label each polygon with the respective label from the metadata
for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
    label = row['UGA_adm4_9']
    # Get the centroid of the polygon to place the label
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

# Add a custom legend to the plot
legend = ax.legend(handles=patches, loc='upper left', title="Legend")

# Add colorbar
cbar = plt.colorbar(contour, ax=ax, label='Soil Map Values', norm=norm, boundaries=boundaries, ticks=boundaries)
ax.set_title('Soil Map(0_20)')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()