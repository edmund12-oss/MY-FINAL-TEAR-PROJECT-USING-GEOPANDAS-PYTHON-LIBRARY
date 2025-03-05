#IMporting the necessary python libraries
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

#SOIL MAPS 0-20cm DEPTH
#SOIL PH
# Load the NetCDF file
nc_file_ph = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil ph(0_20)cm\ph.nc'
ds_ph = xr.open_dataset(nc_file_ph)

# Assuming the dataset contains 'pH, predicted mean at 0-20 cm depth', 'y', and 'x' variables
ph = ds_ph['pH, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_ph['y']
lon = ds_ph['x']

#4.69999981, 6.9000001
#Reclassifying according to the suitability levels
highly_suitable_threshold_ph_1, highly_suitable_threshold_ph_2 = 6, 7.5
moderately_suitable_threshold_ph = 5
marginally_suitable_threshold_ph = 4.5

highly_suitable_ph = np.where((ph >= highly_suitable_threshold_ph_1) & (ph <= highly_suitable_threshold_ph_2) , 1, 0)
moderately_suitable_ph = np.where((ph < highly_suitable_threshold_ph) & (ph >= moderately_suitable_threshold_ph), 0.67, 0)
marginally_suitable_ph = np.where((ph < moderately_suitable_threshold) & (grid_temperature >= marginally_suitable_threshold), 0.33, 0)
unsuitable_ph = np.where(ph < marginally_suitable_threshold_ph, 0, 0)

combined_classes_ph = highly_suitable_ph + moderately_suitable_ph + marginally_suitable_threshold_ph + unsuitable_ph

#SOIL ALUMINIUM EXTRACTABLE*
# Load the NetCDF file
nc_file_al = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil aluminium map(0_20)cm\aluminium_extractable_mean_0_20 (1).nc'
ds_al = xr.open_dataset(nc_file_al)

# Assuming the dataset contains 'Aluminium, extractable, predicted mean at 0-20 cm depth', 'y', and 'x' variables
aluminium_extractable = ds_al['Aluminium, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_al['y']
lon = ds_al['x']

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
nc_file_dtb = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil bedrock depth(0_20)cm\isda_soil_bedrock_depth.nc'
ds_dtb = xr.open_dataset(nc_file_dtb)

# Assuming the dataset contains 'mean_0_200', 'y', and 'x' variables
bedrock_depth = ds_dtb['mean_0_200']  # You might need to adjust variable names based on your NetCDF file structure
bedrock_depth_2d =bedrock_depth.isel(time=0)
lat = ds_dtb['y']
lon = ds_dtb['x']

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
nc_file_bdensity = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil bulk density(0_20)cm\isda_soil_bulk_density.nc'
ds_bdensity = xr.open_dataset(nc_file_bdensity)

# Assuming the dataset contains 'mean_0_20', 'y', and 'x' variables
bulk_density = ds_bdensity['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
bulk_density_2d =bulk_density.isel(time=0)
lat = ds_bdensity['y']
lon = ds_bdensity['x']

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
nc_file_nitro = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil nitrogen total(0_20)cm\nitrogen_total.nc'
ds_nitro = xr.open_dataset(nc_file_nitro)

# Assuming the dataset contains 'Nitrogen, total, predicted mean at 0-20 cm depth', 'y', and 'x' variables
nitrogen_total = ds_nitro['Nitrogen, total, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_nitro['y']
lon = ds_nitro['x']

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
nc_file_ca = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil calcium extractable(0_20)cm\calcium_extractable.nc'
ds_ca = xr.open_dataset(nc_file_ca)

# Assuming the dataset contains 'Calcium, extractable, predicted mean at 0-20 cm depth', 'y', and 'x' variables
calcium_extractable = ds_ca['Calcium, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_ca['y']
lon = ds_ca['x']

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
nc_file_co = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil carbon organic(0_20)cm\carbon_organic.nc'
ds_co = xr.open_dataset(nc_file_co)

# Assuming the dataset contains 'Carbon, organic, predicted mean at 0-20 cm depth', 'y', and 'x' variable
carbon_organic = ds_co['Carbon, organic, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_co['y']
lon = ds_co['x']

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
nc_file_ct = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil carbon total(0_20)cm\isda_soil_carbon_total.nc'
ds_ct = xr.open_dataset(nc_file_ct)

# Assuming the dataset contains 'mean_0_20', 'y', and 'x' variable
carbon_total = ds_ct['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
carbon_total_2d = carbon_total.isel(time=0)
lat = ds_ct['y']
lon = ds_ct['x']

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

nc_file_cec = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil cation exchange capacity(0_20)cm\cation_exchange_capacity.nc'
ds_cec = xr.open_dataset(nc_file_cec)

# Assuming the dataset contains 'Effective Cation Exchange Capacity, predicted mean at 0-20 cm depth', 'y', and 'x' variable
cec = ds_cec['Effective Cation Exchange Capacity, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_cec['y']
lon = ds_cec['x']

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
nc_file_clay = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil clay content(0_20)cm\isda_soil_clay_content.nc'
ds_clay = xr.open_dataset(nc_file_clay)

# Assuming the dataset contains 'mean_0_20', 'y', and 'x' variable
clay_content = ds_clay['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
clay_content_2d =clay_content.isel(time=0)
lat = ds_clay['y']
lon = ds_clay['x']

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
nc_file_fe = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil iron extractable(0_20)cm\iron_extractable.nc'
ds_fe = xr.open_dataset(nc_file_fe)

# Assuming the dataset contains 'Iron, extractable, predicted mean at 0-20 cm depth', 'y', and 'x' variables
iron_extractable = ds_fe['Iron, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_fe['y']
lon = ds_fe['x']

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
nc_file_mg = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil magnesium extractable(0_20)cm\magnesium_extractable.nc'
ds_mg = xr.open_dataset(nc_file_mg)

# Assuming the dataset contains 'Magnesium, extractable, predicted mean at 0-20 cm depth', 'y', and 'x' variables
magnesium_extractable = ds_mg['Magnesium, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_mg['y']
lon = ds_mg['x']

#Converting from mg/kg to g/kg
magnesium_extractable_final = magnesium_extractable / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_mg_1, highly_suitable_threshold_mg_2 = 0.2, 0.4
moderately_suitable_threshold_mg_1 = 0.1
marginally_suitable_threshold_mg_1 = 0.05

highly_suitable_mg = np.where((magnesium_extractable_final >= highly_suitable_threshold_mg_1) & (magnesium_extractable_final <= highly_suitable_threshold_mg_2), 1, 0)
moderately_suitable_mg = np.where((magnesium_extractable_final >= moderately_suitable_threshold_mg_1) & (magnesium_extractable_final < highly_suitable_threshold_mg_1), 0.67, 0)  #& (magnesium_extractable_final > highly_suitable_threshold_mg_2) & (magnesium_extractable_final <= moderately_suitable_threshold_mg_2), 0.67, 0)
marginally_suitable_mg = np.where((magnesium_extractable_final < moderately_suitable_threshold_mg_1) & (magnesium_extractable_final >= marginally_suitable_threshold_mg_1), 0.33, 0)
unsuitable_mg = np.where((magnesium_extractable_final < marginally_suitable_threshold_mg_1),0, 0)

combined_classes_mg = highly_suitable_mg + moderately_suitable_mg + marginally_suitable_mg + unsuitable_mg


#SOIL PHOSPHOROUS

# Load the NetCDF file
nc_file_phosphorus = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil phosphorous extractable(0_20)cm\phosphorous_extractable.nc'
ds_phosphorus = xr.open_dataset(nc_file_phosphorus)

# Assuming the dataset contains 'Phosphorus, extractable, predicted mean at 0-20 cm depth', 'y', and 'x' variables
phosphorus_extractable = ds_phosphorus['Phosphorus, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_phosphorus['y']
lon = ds_phosphorus['x']

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
nc_file_k = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil potassium extractable(0_20)cm\potassium_extractable.nc'
ds_k = xr.open_dataset(nc_file_k)

# Assuming the dataset contains 'Potassium, extractable, predicted mean at 0-20 cm depth', 'lat', and 'lon' variables
potassium_extractable = ds_k['Potassium, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_k['y']
lon = ds_k['x']

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
nc_file_sand = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil sand content(0_20)cm\isda_soil_sand_content.nc'
ds_sand = xr.open_dataset(nc_file_sand)

# Assuming the dataset contains 'mean_0_20', 'y', and 'x' variables
sand_content = ds_sand['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
sand_content_2d = sand_content.isel(time=0)
lat = ds_sand['y']
lon = ds_sand['x']

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
nc_file_silt = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil silt content(0_20)cm\isda_soil_silt_content.nc'
ds_silt = xr.open_dataset(nc_file_silt)

# Assuming the dataset contains 'mean_0_20', 'y', and 'x' variables
silt_content = ds_silt['mean_0_20']  # You might need to adjust variable names based on your NetCDF file structure
silt_content_2d = silt_content.isel(time=0)
lat = ds_silt['y']
lon = ds_silt['x']

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
nc_file_stone = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil stone content(0_20)cm\stone_content.nc'
ds_stone = xr.open_dataset(nc_file_stone)

# Assuming the dataset contains 'Stone content, predicted mean at 0-20 cm depth', 'y', and 'x' variables
stone_content = ds_stone['Stone content, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_stone['y']
lon = ds_stone['x']

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
nc_file_sulphur = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil sulphur extractable(0_20)cm\sulphur_extractable.nc'
ds_sulphur = xr.open_dataset(nc_file_sulphur)

# Assuming the dataset contains 'Sulphur, extractable, predicted mean at 0-20 cm depth', 'y', and 'x' variables
sulphur_extractable = ds_sulphur['Sulphur, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_sulphur['y']
lon = ds_sulphur['x']

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
nc_file_zinc = r'E:\project execution\geopandas\New thematic maps\soil_maps(0_20)cm\soil zinc extractable(0_20)cm\zinc_extractable.nc'
ds_zinc = xr.open_dataset(nc_file_zinc)

# Assuming the dataset contains 'Zinc, extractable, predicted mean at 0-20 cm depth', 'y', and 'x' variables
zinc_extractable = ds_zinc['Zinc, extractable, predicted mean at 0-20 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_zinc['y']
lon = ds_zinc['x']

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

#weights determined using AHP
weight_al = 0.01320508
weight_bdepth = 0.03284426
weight_bdensity = 0.01716747
weight_calcium = 0.03951289
weight_carbon_organic = 0.04899452
weight_ct = 0.04350074
weight_cec = 0.07660854
weight_clay = 0.03786267
weight_iron = 0.03866181
weight_mg = 0.08365584
weight_nitrogen = 0.09870465
weight_ph = 0.13146055
weight_phosphorus = 0.08984892
weight_potassium = 0.10243938
weight_sand = 0.02174041
weight_silt = 0.02745628
weight_stone = 0.01843849
weight_sulphur = 0.03428968
weight_zinc = 0.04360782




soil_map_0_20 = (weight_ph * combined_classes_ph) + (weight_al * combined_classes_al) + (weight_bdepth * combined_classes_bd) + (weight_bdensity * combined_classes_bdensity) + (weight_nitrogen * combined_classes_nitrogen) + (weight_calcium * combined_classes_calcium) + (weight_carbon_organic * combined_classes_co) + (weight_ct * combined_classes_ct) + (weight_cec * combined_classes_cec) + (weight_clay * combined_classes_clay) + (weight_iron * combined_classes_iron) + (weight_mg * combined_classes_mg) + (weight_phosphorus * combined_classes_phosphorus) + (weight_potassium * combined_classes_k) + (weight_sand * combined_classes_sand) + (weight_silt * combined_classes_silt) + (weight_stone * combined_classes_stone) + (weight_sulphur * combined_classes_sulphur) + (weight_zinc * combined_classes_zinc)




#SOIL MAPS 20-50cm DEPTH





#SOIL PH 20-50

# Load the NetCDF file
nc_file_ph_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil ph(20_50)cm\ph.nc'
ds_ph_20_50 = xr.open_dataset(nc_file_ph_20_50)

# Assuming the dataset contains 'pH, predicted mean at 20-50 cm depth', 'y', and 'x' variables
ph_20_50 = ds_ph_20_50['pH, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_ph_20_50['y']
lon = ds_ph_20_50['x']

#4.69999981, 6.9000001
#Reclassifying according to the suitability levels
highly_suitable_threshold_ph_20_50 = 6
moderately_suitable_threshold_ph_20_50 = 5
#marginally_suitable_threshold_ph = 0


highly_suitable_ph_20_50 = np.where(ph_20_50 >= highly_suitable_threshold_ph_20_50, 1, 0)
moderately_suitable_ph_20_50 = np.where((ph_20_50 < highly_suitable_threshold_ph_20_50) & (ph_20_50 >= moderately_suitable_threshold_ph_20_50), 0.67, 0)
#marginally_suitable_ph = np.where((grid_temperature < moderately_suitable_threshold) & (grid_temperature >= marginally_suitable_threshold), 0.33, 0)
unsuitable_ph_20_50 = np.where(ph_20_50 < moderately_suitable_threshold_ph_20_50, 0, 0)

combined_classes_ph_20_50 = highly_suitable_ph_20_50 + moderately_suitable_ph_20_50 + unsuitable_ph_20_50


#SOIL ALUMINIUM EXTRACTABLE 20-50

# Load the NetCDF file
nc_file_al_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil aluminium map(20_50)cm\aluminium_extractable_mean_0_20 (1).nc'
ds_al_20_50 = xr.open_dataset(nc_file_al_20_50)

# Assuming the dataset contains 'Aluminium, extractable, predicted mean at 20-50 cm depth', 'y', and 'x' variables
aluminium_extractable_20_50 = ds_al_20_50['Aluminium, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_al_20_50['y']
lon = ds_al_20_50['x']

#converting from mg/kg to g/kg
aluminium_extractable_final_20_50 = aluminium_extractable_20_50 / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_al_20_50 = 0.5
moderately_suitable_threshold_al_20_50 = 1
marginally_suitable_threshold_al_20_50 = 2


highly_suitable_al_20_50 = np.where(aluminium_extractable_final_20_50 < highly_suitable_threshold_al_20_50, 1, 0)
moderately_suitable_al_20_50 = np.where((aluminium_extractable_final_20_50 >= highly_suitable_threshold_al_20_50) & (aluminium_extractable_final_20_50 <= moderately_suitable_threshold_al_20_50), 0.67, 0)
marginally_suitable_al_20_50 = np.where((aluminium_extractable_final_20_50 > moderately_suitable_threshold_al_20_50) & (aluminium_extractable_final_20_50 <= marginally_suitable_threshold_al_20_50), 0.33, 0)
unsuitable_al_20_50 = np.where(aluminium_extractable_final_20_50 > marginally_suitable_threshold_al_20_50, 0, 0)

combined_classes_al_20_50 = highly_suitable_al_20_50 + moderately_suitable_al_20_50 + marginally_suitable_al_20_50 + unsuitable_al_20_50

#SOIL DEPTH TO BEDROCK 20-50

# Load the NetCDF file
nc_file_dtb_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil bedrock depth(20_50)cm\isda_soil_bedrock_depth.nc'
ds_dtb_20_50 = xr.open_dataset(nc_file_dtb_20_50)

# Assuming the dataset contains 'mean_0_200', 'y', and 'x' variables
bedrock_depth_20_50 = ds_dtb_20_50['mean_0_200']  # You might need to adjust variable names based on your NetCDF file structure
bedrock_depth_2d_20_50 = bedrock_depth_20_50.isel(time=0)
lat = ds_dtb_20_50['y']
lon = ds_dtb_20_50['x']

#Reclassifying according to the suitability levels
highly_suitable_threshold_bd_20_50 = 120
moderately_suitable_threshold_bd_20_50 = 90
marginally_suitable_threshold_bd_20_50 = 50


highly_suitable_bd_20_50 = np.where(bedrock_depth_2d_20_50 > highly_suitable_threshold_bd_20_50, 1, 0)
moderately_suitable_bd_20_50 = np.where((bedrock_depth_2d_20_50 <= highly_suitable_threshold_bd_20_50) & (bedrock_depth_2d_20_50 > moderately_suitable_threshold_bd_20_50), 0.67, 0)
marginally_suitable_bd_20_50 = np.where((bedrock_depth_2d_20_50 <= moderately_suitable_threshold_bd_20_50) & (bedrock_depth_2d_20_50 >= marginally_suitable_threshold_bd_20_50), 0.33, 0)
unsuitable_bd_20_50 = np.where(bedrock_depth_2d_20_50 < marginally_suitable_threshold_bd_20_50, 0, 0)

combined_classes_bd_20_50 = highly_suitable_bd_20_50 + moderately_suitable_bd_20_50 + marginally_suitable_bd_20_50 + unsuitable_bd_20_50


#SOIL BULK DENSITY 20-50

# Load the NetCDF file
nc_file_bdensity_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil bulk density(20_50)cm\isda_soil_bulk_density.nc'
ds_bdensity_20_50 = xr.open_dataset(nc_file_bdensity_20_50)

# Assuming the dataset contains 'mean_20_50', 'y', and 'x' variables
bulk_density_20_50 = ds_bdensity_20_50['mean_20_50']  # You might need to adjust variable names based on your NetCDF file structure
bulk_density_2d_20_50 =bulk_density_20_50.isel(time=0)
lat = ds_bdensity_20_50['y']
lon = ds_bdensity_20_50['x']

#1, 1.5

#Reclassifying according to the suitability levels
highly_suitable_threshold_bdensity_20_50 = 1.2
moderately_suitable_threshold_bdensity_20_50 = 1.4
marginally_suitable_threshold_bdensity_20_50 = 1.5


highly_suitable_bdensity_20_50 = np.where(bulk_density_2d_20_50 <= highly_suitable_threshold_bdensity_20_50, 1, 0)
moderately_suitable_bdensity_20_50 = np.where((bulk_density_2d_20_50 > highly_suitable_threshold_bdensity_20_50) & (bulk_density_2d_20_50 <= moderately_suitable_threshold_bdensity_20_50), 0.67, 0)
marginally_suitable_bdensity_20_50 = np.where((bulk_density_2d_20_50 > moderately_suitable_threshold_bdensity_20_50) & (bulk_density_2d_20_50 <= marginally_suitable_threshold_bdensity_20_50), 0.33, 0)
unsuitable_bdensity_20_50 = np.where(bulk_density_2d_20_50 > marginally_suitable_threshold_bdensity_20_50, 0, 0)

combined_classes_bdensity_20_50 = highly_suitable_bdensity_20_50 + moderately_suitable_bdensity_20_50 + marginally_suitable_bdensity_20_50 + unsuitable_bdensity_20_50

#SOIL NITROGEN 20-50

# Load the NetCDF file
nc_file_nitro_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil nitrogen total(20_50)cm\nitrogen_total.nc'
ds_nitro_20_50 = xr.open_dataset(nc_file_nitro_20_50)

# Assuming the dataset contains 'Nitrogen, total, predicted mean at 20-50 cm depth', 'y', and 'x' variables
nitrogen_total_20_50 = ds_nitro_20_50['Nitrogen, total, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_nitro_20_50['y']
lon = ds_nitro_20_50['x']

nitrogen_total_final_20_50 = nitrogen_total_20_50 / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_nitrogen_1_20_50, highly_suitable_threshold_nitrogen_2_20_50 = 1.5, 2
moderately_suitable_threshold_nitrogen_20_50 = 1
marginally_suitable_threshold_nitrogen_20_50 = 0.5


highly_suitable_nitrogen_20_50 = np.where((nitrogen_total_final_20_50 >= highly_suitable_threshold_nitrogen_1_20_50) & (nitrogen_total_final_20_50 <= highly_suitable_threshold_nitrogen_2_20_50), 1, 0)
moderately_suitable_nitrogen_20_50 = np.where((nitrogen_total_final_20_50 >= moderately_suitable_threshold_nitrogen_20_50) & (nitrogen_total_final_20_50 < highly_suitable_threshold_nitrogen_1_20_50), 0.67, 0)
marginally_suitable_nitrogen_20_50 = np.where((nitrogen_total_final_20_50 < moderately_suitable_threshold_nitrogen_20_50) & (nitrogen_total_final_20_50 >= marginally_suitable_threshold_nitrogen_20_50), 0.33, 0)
unsuitable_nitrogen_20_50 = np.where(nitrogen_total_final_20_50 < marginally_suitable_threshold_nitrogen_20_50, 0, 0)

combined_classes_nitrogen_20_50 = highly_suitable_nitrogen_20_50 + moderately_suitable_nitrogen_20_50 + marginally_suitable_nitrogen_20_50 + unsuitable_nitrogen_20_50

#SOIL CALCIUM 20-50

# Load the NetCDF file
nc_file_ca_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil calcium extractable(20_50)cm\calcium_extractable.nc'
ds_ca_20_50 = xr.open_dataset(nc_file_ca_20_50)

# Assuming the dataset contains 'Calcium, extractable, predicted mean at 20-50 cm depth', 'y', and 'x' variables
calcium_extractable_20_50 = ds_ca_20_50['Calcium, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_ca_20_50['y']
lon = ds_ca_20_50['x']

#163.02189636, 2207.34765625
#converting to g/kg
calcium_extractable_final_20_50 = calcium_extractable_20_50 / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_calcium_1_20_50, highly_suitable_threshold_calcium_2_20_50 = 0.6, 1.2
moderately_suitable_threshold_calcium_1_20_50, moderately_suitable_threshold_calcium_2_20_50 = 0.3, 1.5
marginally_suitable_threshold_calcium_1_20_50, marginally_suitable_threshold_calcium_2_20_50 = 0.15, 2


highly_suitable_calcium_20_50 = np.where((calcium_extractable_final_20_50 >= highly_suitable_threshold_calcium_1_20_50) & (calcium_extractable_final_20_50 <= highly_suitable_threshold_calcium_2_20_50), 1, 0)
moderately_suitable_calcium_20_50 = np.where((calcium_extractable_final_20_50 >= moderately_suitable_threshold_calcium_1_20_50) & (calcium_extractable_final_20_50 < highly_suitable_threshold_calcium_1_20_50) & (calcium_extractable_final_20_50 > highly_suitable_threshold_calcium_2_20_50) & (calcium_extractable_final_20_50 <= moderately_suitable_threshold_calcium_2_20_50), 0.67, 0)
marginally_suitable_calcium_20_50 = np.where((calcium_extractable_final_20_50 < moderately_suitable_threshold_calcium_1_20_50) & (calcium_extractable_final_20_50 >= marginally_suitable_threshold_calcium_1_20_50) & (calcium_extractable_final_20_50 > moderately_suitable_threshold_calcium_2_20_50) & (calcium_extractable_final_20_50 <= marginally_suitable_threshold_calcium_2_20_50), 0.33, 0)
unsuitable_calcium_20_50 = np.where((calcium_extractable_final_20_50 < marginally_suitable_threshold_calcium_1_20_50) & (calcium_extractable_final_20_50 > marginally_suitable_threshold_calcium_2_20_50), 0, 0)

combined_classes_calcium_20_50 = highly_suitable_calcium_20_50 + moderately_suitable_calcium_20_50 + marginally_suitable_calcium_20_50 + unsuitable_calcium_20_50

#SOIL CARBON ORGANIC 20-50

# Load the NetCDF file
nc_file_co_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil carbon organic(20_50)cm\carbon_organic.nc'
ds_co_20_50 = xr.open_dataset(nc_file_co_20_50)

# Assuming the dataset contains 'Carbon, organic, predicted mean at 20-50 cm depth', 'y', and 'x' variables
carbon_organic_20_50 = ds_co_20_50['Carbon, organic, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_co_20_50['y']
lon = ds_co_20_50['x']

#max_value = carbon_organic.max()
#print(max_value)

#8.02501392, 23.53253174

#Converting from g/kg to a percentage

carbon_organic_percentage_20_50 = (carbon_organic_20_50 / 1000) * 100

#Reclassifying according to the suitability levels
highly_suitable_threshold_co_20_50 = 3
moderately_suitable_threshold_co_20_50 = 2
marginally_suitable_threshold_co_20_50 = 1


highly_suitable_co_20_50 = np.where(carbon_organic_percentage_20_50 > highly_suitable_threshold_co_20_50, 1, 0)
moderately_suitable_co_20_50 = np.where((carbon_organic_percentage_20_50 <= highly_suitable_threshold_co_20_50) & (carbon_organic_percentage_20_50 >= moderately_suitable_threshold_co_20_50), 0.67, 0)
marginally_suitable_co_20_50 = np.where((carbon_organic_percentage_20_50 >= marginally_suitable_threshold_co_20_50) & (carbon_organic_percentage_20_50 < moderately_suitable_threshold_co_20_50), 0.33, 0)
unsuitable_co_20_50 = np.where(carbon_organic_percentage_20_50 < marginally_suitable_threshold_co_20_50, 0, 0)

combined_classes_co_20_50 = highly_suitable_co_20_50 + moderately_suitable_co_20_50 + marginally_suitable_co_20_50 + unsuitable_co_20_50

#SOIL CARBON TOTAL 20-50

# Load the NetCDF file
nc_file_ct_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil carbon total(20_50)cm\isda_soil_carbon_total.nc'
ds_ct_20_50 = xr.open_dataset(nc_file_ct_20_50)

# Assuming the dataset contains 'mean_20_50', 'y', and 'x' variables
carbon_total_20_50 = ds_ct_20_50['mean_20_50']  # You might need to adjust variable names based on your NetCDF file structure
carbon_total_2d_20_50 = carbon_total_20_50.isel(time=0)
lat = ds_ct_20_50['y']
lon = ds_ct_20_50['x']

#converting from g/kg to percentage 
carbon_total_percentage_20_50 = carbon_total_2d_20_50 / 10

#Reclassifying according to the suitability levels
highly_suitable_threshold_ct_20_50 = 2
moderately_suitable_threshold_ct_20_50 = 1.2
marginally_suitable_threshold_ct_20_50 = 0.6


highly_suitable_ct_20_50 = np.where(carbon_total_percentage_20_50 > highly_suitable_threshold_ct_20_50, 1, 0)
moderately_suitable_ct_20_50 = np.where((carbon_total_percentage_20_50 <= highly_suitable_threshold_ct_20_50) & (carbon_total_percentage_20_50 >= moderately_suitable_threshold_ct_20_50), 0.67, 0)
marginally_suitable_ct_20_50 = np.where((carbon_total_percentage_20_50 >= marginally_suitable_threshold_ct_20_50) & (carbon_total_percentage_20_50 < moderately_suitable_threshold_ct_20_50), 0.33, 0)
unsuitable_ct_20_50 = np.where(carbon_total_percentage_20_50 < marginally_suitable_threshold_ct_20_50, 0, 0)

combined_classes_ct_20_50 = highly_suitable_ct_20_50 + moderately_suitable_ct_20_50 + marginally_suitable_ct_20_50 + unsuitable_ct_20_50


#SOIL CATION EXCHANGE CAPACITY 20-50

nc_file_cec_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil cation exchange capacity(20_50)cm\cation_exchange_capacity.nc'
ds_cec_20_50 = xr.open_dataset(nc_file_cec_20_50)

# Assuming the dataset contains 'Effective Cation Exchange Capacity, predicted mean at 20-50 cm depth', 'y', and 'x' variables
cec_20_50 = ds_cec_20_50['Effective Cation Exchange Capacity, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_cec_20_50['y']
lon = ds_cec_20_50['x']

#Reclassifying according to the suitability levels
highly_suitable_threshold_cec_20_50 = 20
moderately_suitable_threshold_cec_20_50 = 10
marginally_suitable_threshold_cec_20_50 = 5


highly_suitable_cec_20_50 = np.where(cec_20_50 > highly_suitable_threshold_cec_20_50, 1, 0)
moderately_suitable_cec_20_50 = np.where((cec_20_50 <= highly_suitable_threshold_cec_20_50) & (cec_20_50 >= moderately_suitable_threshold_cec_20_50), 0.67, 0)
marginally_suitable_cec_20_50 = np.where((cec_20_50 >= marginally_suitable_threshold_cec_20_50) & (cec_20_50 < moderately_suitable_threshold_cec_20_50), 0.33, 0)
unsuitable_cec_20_50 = np.where(cec_20_50 < marginally_suitable_threshold_cec_20_50, 0, 0)

combined_classes_cec_20_50 = highly_suitable_cec_20_50 + moderately_suitable_cec_20_50 + marginally_suitable_cec_20_50 + unsuitable_cec_20_50

#SOIL CLAY CONTENT 20-50

# Load the NetCDF file
nc_file_clay_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil clay content(20_50)cm\isda_soil_clay_content.nc'
ds_clay_20_50 = xr.open_dataset(nc_file_clay_20_50)

# Assuming the dataset contains 'mean_20_50', 'y', and 'x' variables
clay_content_20_50 = ds_clay_20_50['mean_20_50']  # You might need to adjust variable names based on your NetCDF file structure
clay_content_2d_20_50 =clay_content_20_50.isel(time=0)
lat = ds_clay_20_50['y']
lon = ds_clay_20_50['x']

#Reclassifying according to the suitability levels
highly_suitable_threshold_clay_1_20_50, highly_suitable_threshold_clay_2_20_50 = 20, 35
moderately_suitable_threshold_clay_20_50 = 45
marginally_suitable_threshold_clay_1_20_50, marginally_suitable_threshold_clay_2_20_50 = 15, 60


highly_suitable_clay_20_50 = np.where((clay_content_2d_20_50 >= highly_suitable_threshold_clay_1_20_50) & (clay_content_2d_20_50 <= highly_suitable_threshold_clay_2_20_50), 1, 0)
moderately_suitable_clay_20_50 = np.where((clay_content_2d_20_50 > highly_suitable_threshold_clay_2_20_50) & (clay_content_2d_20_50 <= moderately_suitable_threshold_clay_20_50), 0.67, 0)
marginally_suitable_clay_20_50 = np.where((clay_content_2d_20_50 >= marginally_suitable_threshold_clay_1_20_50) & (clay_content_2d_20_50 < highly_suitable_threshold_clay_1_20_50) & (clay_content_2d_20_50 > moderately_suitable_threshold_clay_20_50) & (clay_content_2d_20_50 <= marginally_suitable_threshold_clay_2_20_50), 0.33, 0)
unsuitable_clay_20_50 = np.where((clay_content_2d_20_50 < marginally_suitable_threshold_clay_1_20_50) & (clay_content_2d_20_50 > marginally_suitable_threshold_clay_2_20_50), 0, 0)

combined_classes_clay_20_50 = highly_suitable_clay_20_50 + moderately_suitable_clay_20_50 + marginally_suitable_clay_20_50 + unsuitable_clay_20_50

#SOIL IRON EXTRACTABLE 20-50

# Load the NetCDF file
nc_file_fe_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil iron extractable(20_50)cm\iron_extractable.nc'
ds_fe_20_50 = xr.open_dataset(nc_file_fe_20_50)

# Assuming the dataset contains 'Iron, extractable, predicted mean at 20-50 cm depth', 'y', and 'x' variables
iron_extractable_20_50 = ds_fe_20_50['Iron, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_fe_20_50['y']
lon = ds_fe_20_50['x']

#converting from mg/kg to g/kg
iron_extractable_final_20_50 = iron_extractable_20_50 / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_iron_1_20_50, highly_suitable_threshold_iron_2_20_50 = 0.1, 0.3
moderately_suitable_threshold_iron_1_20_50, moderately_suitable_threshold_iron_2_20_50 = 0.05, 0.5
marginally_suitable_threshold_iron_1_20_50, marginally_suitable_threshold_iron_2_20_50 = 0.02, 0.8


highly_suitable_iron_20_50 = np.where((iron_extractable_final_20_50 >= highly_suitable_threshold_iron_1_20_50) & (iron_extractable_final_20_50 <= highly_suitable_threshold_iron_2_20_50), 1, 0)
moderately_suitable_iron_20_50 = np.where((iron_extractable_final_20_50 >= moderately_suitable_threshold_iron_1_20_50) & (iron_extractable_final_20_50 < highly_suitable_threshold_iron_1_20_50) & (iron_extractable_final_20_50 > highly_suitable_threshold_iron_2_20_50) & (iron_extractable_final_20_50 <= moderately_suitable_threshold_iron_2_20_50), 0.67, 0)
marginally_suitable_iron_20_50 = np.where((iron_extractable_final_20_50 < moderately_suitable_threshold_iron_1_20_50) & (iron_extractable_final_20_50 >= marginally_suitable_threshold_iron_1_20_50) & (iron_extractable_final_20_50 > moderately_suitable_threshold_iron_2_20_50) & (iron_extractable_final_20_50 <= marginally_suitable_threshold_iron_2_20_50), 0.33, 0)
unsuitable_iron_20_50 = np.where((iron_extractable_final_20_50 < marginally_suitable_threshold_iron_1_20_50) & (iron_extractable_final_20_50 > marginally_suitable_threshold_iron_2_20_50), 0, 0)

combined_classes_iron_20_50 = highly_suitable_iron_20_50 + moderately_suitable_iron_20_50 + marginally_suitable_iron_20_50 + unsuitable_iron_20_50

# SOIL MAGNESIUM 20-50

# Load the NetCDF file
nc_file_mg_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil magnesium extractable(20_50)cm\magnesium_extractable.nc'
ds_mg_20_50 = xr.open_dataset(nc_file_mg_20_50)

# Assuming the dataset contains 'Magnesium, extractable, predicted mean at 20-50 cm depth', 'y', and 'x' variables
magnesium_extractable_20_50 = ds_mg_20_50['Magnesium, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_mg_20_50['y']
lon = ds_mg_20_50['x']

#Converting from mg/kg to g/kg
magnesium_extractable_final_20_50 = magnesium_extractable_20_50 / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_mg_1_20_50, highly_suitable_threshold_mg_2_20_50 = 0.2, 0.4
moderately_suitable_threshold_mg_1_20_50 = 0.1
marginally_suitable_threshold_mg_1_20_50 = 0.05


highly_suitable_mg_20_50 = np.where((magnesium_extractable_final_20_50 >= highly_suitable_threshold_mg_1_20_50) & (magnesium_extractable_final_20_50 <= highly_suitable_threshold_mg_2_20_50), 1, 0)
moderately_suitable_mg_20_50 = np.where((magnesium_extractable_final_20_50 >= moderately_suitable_threshold_mg_1_20_50) & (magnesium_extractable_final_20_50 < highly_suitable_threshold_mg_1_20_50), 0.67, 0)
marginally_suitable_mg_20_50 = np.where((magnesium_extractable_final_20_50 < moderately_suitable_threshold_mg_1_20_50) & (magnesium_extractable_final_20_50 >= marginally_suitable_threshold_mg_1_20_50), 0.33, 0)
unsuitable_mg_20_50 = np.where((magnesium_extractable_final_20_50 < marginally_suitable_threshold_mg_1_20_50), 0, 0)

combined_classes_mg_20_50 = highly_suitable_mg_20_50 + moderately_suitable_mg_20_50 + marginally_suitable_mg_20_50 + unsuitable_mg_20_50

#SOIL PHOSPHOROUS 20-50

# Load the NetCDF file
nc_file_phosphorus_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil phosphorous extractable(20_50)cm\phosphorous_extractable.nc'
ds_phosphorus_20_50 = xr.open_dataset(nc_file_phosphorus_20_50)

# Assuming the dataset contains 'Phosphorus, extractable, predicted mean at 20-50 cm depth', 'y', and 'x' variables
phosphorus_extractable_20_50 = ds_phosphorus_20_50['Phosphorus, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_phosphorus_20_50['y']
lon = ds_phosphorus_20_50['x']

#Reclassifying according to the suitability levels
highly_suitable_threshold_phosphorus_1_20_50, highly_suitable_threshold_phosphorus_2_20_50 = 20, 50
moderately_suitable_threshold_phosphorus_20_50 = 10
marginally_suitable_threshold_phosphorus_20_50 = 5


highly_suitable_phosphorus_20_50 = np.where((phosphorus_extractable_20_50 >= highly_suitable_threshold_phosphorus_1_20_50) & (phosphorus_extractable_20_50 <= highly_suitable_threshold_phosphorus_2_20_50), 1, 0)
moderately_suitable_phosphorus_20_50 = np.where((phosphorus_extractable_20_50 >= moderately_suitable_threshold_phosphorus_20_50) & (phosphorus_extractable_20_50 < highly_suitable_threshold_phosphorus_1_20_50), 0.67, 0)
marginally_suitable_phosphorus_20_50 = np.where((phosphorus_extractable_20_50 < moderately_suitable_threshold_phosphorus_20_50) & (phosphorus_extractable_20_50 >= marginally_suitable_threshold_phosphorus_20_50), 0.33, 0)
unsuitable_phosphorus_20_50 = np.where(phosphorus_extractable_20_50 < marginally_suitable_threshold_phosphorus_20_50, 0, 0)

combined_classes_phosphorus_20_50 = highly_suitable_phosphorus_20_50 + moderately_suitable_phosphorus_20_50 + marginally_suitable_phosphorus_20_50 + unsuitable_phosphorus_20_50

#SOIL POTASSIUM 20-50

# Load the NetCDF file
nc_file_k_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil potassium extractable(20_50)cm\potassium_extractable.nc'
ds_k_20_50 = xr.open_dataset(nc_file_k_20_50)

# Assuming the dataset contains 'Potassium, extractable, predicted mean at 20-50 cm depth', 'y', and 'x' variables
potassium_extractable_20_50 = ds_k_20_50['Potassium, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_k_20_50['y']
lon = ds_k_20_50['x']

#Converting from mg/kg to g/kg
potassium_extractable_final_20_50 = potassium_extractable_20_50 / 1000

#Reclassifying according to the suitability levels
highly_suitable_threshold_k_1_20_50, highly_suitable_threshold_k_2_20_50 = 0.15, 0.25
moderately_suitable_threshold_k_20_50 = 0.10
marginally_suitable_threshold_k_20_50 = 0.05


highly_suitable_k_20_50 = np.where((potassium_extractable_final_20_50 >= highly_suitable_threshold_k_1_20_50) & (potassium_extractable_final_20_50 <= highly_suitable_threshold_k_2_20_50), 1, 0)
moderately_suitable_k_20_50 = np.where((potassium_extractable_final_20_50 >= moderately_suitable_threshold_k_20_50) & (potassium_extractable_final_20_50 < highly_suitable_threshold_k_1_20_50), 0.67, 0)
marginally_suitable_k_20_50 = np.where((potassium_extractable_final_20_50 < moderately_suitable_threshold_k_20_50) & (potassium_extractable_final_20_50 >= marginally_suitable_threshold_k_20_50), 0.33, 0)
unsuitable_k_20_50 = np.where(potassium_extractable_final_20_50 < marginally_suitable_threshold_k_20_50, 0, 0)

combined_classes_k_20_50 = highly_suitable_k_20_50 + moderately_suitable_k_20_50 + marginally_suitable_k_20_50 + unsuitable_k_20_50

#SOIL SAND CONTENT 20-50

# Load the NetCDF file
nc_file_sand_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil sand content(20_50)cm\isda_soil_sand_content.nc'
ds_sand_20_50 = xr.open_dataset(nc_file_sand_20_50)

# Assuming the dataset contains 'mean_20_50', 'y', and 'x' variables
sand_content_20_50 = ds_sand_20_50['mean_20_50']  # You might need to adjust variable names based on your NetCDF file structure
sand_content_2d_20_50 = sand_content_20_50.isel(time=0)
lat = ds_sand_20_50['y']
lon = ds_sand_20_50['x']


#Reclassifying according to the suitability levels
highly_suitable_threshold_sand_1_20_50, highly_suitable_threshold_sand_2_20_50 = 20, 40
moderately_suitable_threshold_sand_20_50  = 60
marginally_suitable_threshold_sand_1_20_50, marginally_suitable_threshold_sand_2_20_50 = 15, 70


highly_suitable_sand_20_50 = np.where((sand_content_2d_20_50 >= highly_suitable_threshold_sand_1_20_50) & (sand_content_2d_20_50 <= highly_suitable_threshold_sand_2_20_50), 1, 0)
moderately_suitable_sand_20_50  = np.where((sand_content_2d_20_50 > highly_suitable_threshold_sand_2_20_50) & (sand_content_2d_20_50 <= moderately_suitable_threshold_sand_20_50), 0.67, 0)
marginally_suitable_sand_20_50  = np.where((sand_content_2d_20_50 >= marginally_suitable_threshold_sand_1_20_50) & (sand_content_2d_20_50 < highly_suitable_threshold_sand_1_20_50) & (sand_content_2d_20_50 > moderately_suitable_threshold_sand_20_50) & (sand_content_2d_20_50 <= marginally_suitable_threshold_sand_2_20_50), 0.33, 0)
unsuitable_sand_20_50  = np.where((sand_content_2d_20_50 < marginally_suitable_threshold_sand_1_20_50) & (sand_content_2d_20_50 > marginally_suitable_threshold_sand_2_20_50), 0, 0)

combined_classes_sand_20_50 = highly_suitable_sand_20_50 + moderately_suitable_sand_20_50 + marginally_suitable_sand_20_50 + unsuitable_sand_20_50

#SOIL SILT CONTENT 20-50

# Load the NetCDF file
nc_file_silt_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil silt content(20_50)cm\isda_soil_silt_content.nc'
ds_silt_20_50 = xr.open_dataset(nc_file_silt_20_50)

# Assuming the dataset contains 'mean_20_50', 'y', and 'x' variables
silt_content_20_50 = ds_silt_20_50['mean_20_50']  # You might need to adjust variable names based on your NetCDF file structure
silt_content_2d_20_50 = silt_content_20_50.isel(time=0)
lat = ds_silt_20_50['y']
lon = ds_silt_20_50['x']

#Reclassifying according to the suitability levels
highly_suitable_threshold_silt_1_20_50, highly_suitable_threshold_silt_2_20_50 = 20, 40
moderately_suitable_threshold_silt_20_50 = 15
marginally_suitable_threshold_silt_20_50 = 10


highly_suitable_silt_20_50 = np.where((silt_content_2d_20_50 >= highly_suitable_threshold_silt_1_20_50) & (silt_content_2d_20_50 <= highly_suitable_threshold_silt_2_20_50), 1, 0)
moderately_suitable_silt_20_50 = np.where((silt_content_2d_20_50 >= moderately_suitable_threshold_silt_20_50) & (silt_content_2d_20_50 < highly_suitable_threshold_silt_1_20_50), 0.67, 0)
marginally_suitable_silt_20_50 = np.where((silt_content_2d_20_50 < moderately_suitable_threshold_silt_20_50) & (silt_content_2d_20_50 >= marginally_suitable_threshold_silt_20_50), 0.33, 0)
unsuitable_silt_20_50 = np.where(silt_content_2d_20_50 < marginally_suitable_threshold_silt_20_50, 0, 0)

combined_classes_silt_20_50 = highly_suitable_silt_20_50 + moderately_suitable_silt_20_50 + marginally_suitable_silt_20_50 + unsuitable_silt_20_50

#SOIL STONE CONTENT 20-50

# Load the NetCDF file
nc_file_stone_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil stone content(20_50)cm\stone_content.nc'
ds_stone_20_50 = xr.open_dataset(nc_file_stone_20_50)

# Assuming the dataset contains 'Stone content, predicted mean at 20-50 cm depth', 'y', and 'x' variables
stone_content_20_50 = ds_stone_20_50['Stone content, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_stone_20_50['y']
lon = ds_stone_20_50['x']

#Reclassifying according to the suitability levels
highly_suitable_threshold_stone_20_50 = 10
moderately_suitable_threshold_stone_20_50 = 20
marginally_suitable_threshold_stone_20_50 = 35


highly_suitable_stone_20_50 = np.where(stone_content_20_50 < highly_suitable_threshold_stone_20_50, 1, 0)
moderately_suitable_stone_20_50 = np.where((stone_content_20_50 >= highly_suitable_threshold_stone_20_50) & (stone_content_20_50 <= moderately_suitable_threshold_stone_20_50), 0.67, 0)
marginally_suitable_stone_20_50 = np.where((stone_content_20_50 <= marginally_suitable_threshold_stone_20_50) & (stone_content_20_50 > moderately_suitable_threshold_stone_20_50), 0.33, 0)
unsuitable_stone_20_50 = np.where(stone_content_20_50 > marginally_suitable_threshold_stone_20_50, 0, 0)

combined_classes_stone_20_50 = highly_suitable_stone_20_50 + moderately_suitable_stone_20_50 + marginally_suitable_stone_20_50 + unsuitable_stone_20_50

#SOIL SULPHUR EXTRACTABLE 20-50

# Load the NetCDF file
nc_file_sulphur_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil sulphur extractable(20_50)cm\sulphur_extractable.nc'
ds_sulphur_20_50 = xr.open_dataset(nc_file_sulphur_20_50)

# Assuming the dataset contains 'Sulphur, extractable, predicted mean at 20-50 cm depth', 'y', and 'x' variables
sulphur_extractable_20_50 = ds_sulphur_20_50['Sulphur, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_sulphur_20_50['y']
lon = ds_sulphur_20_50['x']

#Reclassifying according to the suitability levels
highly_suitable_threshold_sulphur_1_20_50, highly_suitable_threshold_sulphur_2_20_50 = 10, 30
moderately_suitable_threshold_sulphur_20_50 = 5
marginally_suitable_threshold_sulphur_20_50 = 3


highly_suitable_sulphur_20_50 = np.where((sulphur_extractable_20_50 >= highly_suitable_threshold_sulphur_1_20_50) & (sulphur_extractable_20_50 <= highly_suitable_threshold_sulphur_2_20_50), 1, 0)
moderately_suitable_sulphur_20_50 = np.where((sulphur_extractable_20_50 >= moderately_suitable_threshold_sulphur_20_50) & (sulphur_extractable_20_50 < highly_suitable_threshold_sulphur_1_20_50), 0.67, 0)
marginally_suitable_sulphur_20_50 = np.where((sulphur_extractable_20_50 < moderately_suitable_threshold_sulphur_20_50) & (sulphur_extractable_20_50 >= marginally_suitable_threshold_sulphur_20_50), 0.33, 0)
unsuitable_sulphur_20_50 = np.where(sulphur_extractable_20_50 < marginally_suitable_threshold_sulphur_20_50, 0, 0)

combined_classes_sulphur_20_50 = highly_suitable_sulphur_20_50 + moderately_suitable_sulphur_20_50 + marginally_suitable_sulphur_20_50 + unsuitable_sulphur_20_50

#SOIL ZINC EXTRACTABLE

# Load the NetCDF file
nc_file_zinc_20_50 = r'E:\project execution\geopandas\New thematic maps\soil_maps(20_50)cm\soil zinc extractable(20_50)cm\zinc_extractable.nc'
ds_zinc_20_50 = xr.open_dataset(nc_file_zinc_20_50)

# Assuming the dataset contains 'Zinc, extractable, predicted mean at 20-50 cm depth', 'y', and 'x' variables
zinc_extractable_20_50 = ds_zinc_20_50['Zinc, extractable, predicted mean at 20-50 cm depth']  # You might need to adjust variable names based on your NetCDF file structure
lat = ds_zinc_20_50['y']
lon = ds_zinc_20_50['x']

#Reclassifying according to the suitability levels
highly_suitable_threshold_zinc_1_20_50, highly_suitable_threshold_zinc_2_20_50 = 1.5, 10
moderately_suitable_threshold_zinc_20_50 = 1
marginally_suitable_threshold_zinc_20_50 = 0.5


highly_suitable_zinc_20_50 = np.where((zinc_extractable_20_50 >= highly_suitable_threshold_zinc_1_20_50) & (zinc_extractable_20_50 <= highly_suitable_threshold_zinc_2_20_50), 1, 0)
moderately_suitable_zinc_20_50 = np.where((zinc_extractable_20_50 >= moderately_suitable_threshold_zinc_20_50) & (zinc_extractable_20_50 < highly_suitable_threshold_zinc_1_20_50), 0.67, 0)
marginally_suitable_zinc_20_50 = np.where((zinc_extractable_20_50 < moderately_suitable_threshold_zinc_20_50) & (zinc_extractable_20_50 >= marginally_suitable_threshold_zinc_20_50), 0.33, 0)
unsuitable_zinc_20_50 = np.where(zinc_extractable_20_50 < marginally_suitable_threshold_zinc_20_50, 0, 0)

combined_classes_zinc_20_50 = highly_suitable_zinc_20_50 + moderately_suitable_zinc_20_50 + marginally_suitable_zinc_20_50 + unsuitable_zinc_20_50


#WEIGHTED OVERLAY OF SOIL MAPS(0_20)


soil_map_20_50 = (weight_ph * combined_classes_ph_20_50) + (weight_al * combined_classes_al_20_50) + (weight_bdepth * combined_classes_bd_20_50) + (weight_bdensity * combined_classes_bdensity_20_50) + (weight_nitrogen * combined_classes_nitrogen_20_50) + (weight_calcium * combined_classes_calcium_20_50) + (weight_carbon_organic * combined_classes_co_20_50) + (weight_ct * combined_classes_ct_20_50) + (weight_cec * combined_classes_cec_20_50) + (weight_clay * combined_classes_clay_20_50) + (weight_iron * combined_classes_iron_20_50) + (weight_mg * combined_classes_mg_20_50) + (weight_phosphorus * combined_classes_phosphorus_20_50) + (weight_potassium * combined_classes_k_20_50) + (weight_sand * combined_classes_sand_20_50) + (weight_silt * combined_classes_silt_20_50) + (weight_stone * combined_classes_stone_20_50) + (weight_sulphur * combined_classes_sulphur_20_50) + (weight_zinc * combined_classes_zinc_20_50)

#Combining 0_20 and 20_50
weight_soil_0_20 = 0.83333333
weight_soil_20_50 = 0.16666667

final_soil_map = (weight_soil_0_20 * soil_map_0_20) + (weight_soil_20_50  * soil_map_20_50)


# CLIMATE MAPs


#RAINFALL MAP

# Load your rainfall data from a NetCDF file
rainfall_data_path = r'E:\project execution\geopandas\New thematic maps\rainfall_map\2003_2023_rf_data.nc'
rainfall_data = xr.open_dataset(rainfall_data_path)

# Assuming your rainfall data has lat and lon dimensions, and a rainfall variable
# Adjust 'rainfall' below to match the variable name in your dataset
rainfall_values_2d = rainfall_data['rainfall']  # Replace 'rainfall' with your variable name
#lat = rainfall_data['y']
#lon = rainfall_data['x']

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
nc_file_mxtemp = r'E:\project execution\geopandas\New thematic maps\maximum temperature map\final_max_temperature.nc'
ds_mxtemp = xr.open_dataset(nc_file_mxtemp)

# Assuming the dataset contains 'air_temperature_at_2_metres', 'lat', and 'lon' variables
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
nc_file_mintemp = r'E:\project execution\geopandas\New thematic maps\minimum temperature map\final_min_temperature.nc'
ds_mintemp = xr.open_dataset(nc_file_mintemp)

# Assuming the dataset contains 'air_temperature_at_2_metres', 'lat', and 'lon' variables
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
nc_file_windspeed = r'E:\project execution\geopandas\New thematic maps\windspeed map\wind_speed_new.nc'
ds_windspeed = xr.open_dataset(nc_file_windspeed)

# Assuming the dataset contains '__xarray_dataarray_variable__', 'lat', and 'lon' variables
wind_speed = ds_windspeed['__xarray_dataarray_variable__']# You might need to adjust variable names based on your NetCDF file structure
wind_speed_2d = wind_speed.isel(time=0)
#lat = ds_windspeed['lat']
#lon = ds_windspeed['lon']

#Reclassifying according to the suitability levels
highly_suitable_threshold_windspeed_1, highly_suitable_threshold_windspeed_2 = 0, 2
moderately_suitable_threshold_windspeed = 4
marginally_suitable_threshold_windspeed = 6


highly_suitable_windspeed = np.where((wind_speed_2d >= highly_suitable_threshold_windspeed_1) & (wind_speed_2d <= highly_suitable_threshold_windspeed_2), 1, 0)
moderately_suitable_windspeed = np.where((wind_speed_2d <= moderately_suitable_threshold_windspeed) & (wind_speed_2d > highly_suitable_threshold_windspeed_2), 0.67, 0)
marginally_suitable_windspeed = np.where((wind_speed_2d > moderately_suitable_threshold_windspeed) & (wind_speed_2d <= marginally_suitable_threshold_windspeed), 0.33, 0)
unsuitable_windspeed = np.where(wind_speed_2d > marginally_suitable_threshold_windspeed, 0, 0)

combined_classes_windspeed = highly_suitable_windspeed + moderately_suitable_windspeed + marginally_suitable_windspeed + unsuitable_windspeed


#tentative weight till I use AHP to determine actual weights
weight_rf = 0.21676587
weight_mxtemp = 0.37103175
weight_mintemp = 0.37103175
weight_windspeed = 0.04117063


climate_map_1 = (weight_rf * combined_classes_rf) 
climate_map_2 = (weight_mxtemp * combined_classes_mxtemp) + (weight_mintemp * combined_classes_mintemp) + (weight_windspeed * combined_classes_windspeed)

#changing the shape of an array
new_array = np.zeros_like(climate_map_1)

new_array[:climate_map_2.shape[0], :climate_map_2.shape[1]] = climate_map_2

climate_data = climate_map_1 + new_array


# SLOPE MAP

# Load the NetCDF file
nc_file_slope = r'E:\project execution\geopandas\New thematic maps\slope map\slope_map.nc'
ds_slope = xr.open_dataset(nc_file_slope)


slope = ds_slope['slope']  # You might need to adjust variable names based on your NetCDF file structure
#lat_new = ds_slope['y']
#lon_new = ds_slope['x']

#reclassifying according to suitability classes
highly_suitable_threshold_slope = 4
moderately_suitable_threshold_slope = 8
marginally_suitable_threshold_slope = 20


highly_suitable_slope = np.where(slope < highly_suitable_threshold_slope, 1, 0)
moderately_suitable_slope = np.where((slope >= highly_suitable_threshold_slope) & (slope <= moderately_suitable_threshold_slope), 0.67, 0)
marginally_suitable_slope = np.where((slope > moderately_suitable_threshold_slope) & (slope <= marginally_suitable_threshold_slope), 0.33, 0)
unsuitable_slope = np.where(slope > marginally_suitable_threshold_slope, 0, 0)

combined_classes_slope = highly_suitable_slope + moderately_suitable_slope + marginally_suitable_slope + unsuitable_slope








#LAND USE LAND COVER 

# Load the NetCDF file
nc_file_lulc = r'E:\project execution\geopandas\New thematic maps\LULC\lulc_reclass.nc'
ds_lulc = xr.open_dataset(nc_file_lulc)

lulc_suitability = ds_lulc['lulc']  # Replace 'rainfall' with your variable name
lat_new = ds_lulc['lat']
lon_new = ds_lulc['lon']

threshold = 0.5

highly_suitable = np.where(lulc_suitability < threshold, 1, 0)
unsuitable = np.where(lulc_suitability > threshold, 0, 0)

combined_classes_lulc = highly_suitable + unsuitable



#changing shape of an array
new_final_soil_map = np.zeros_like(combined_classes_lulc)

new_final_soil_map[:final_soil_map.shape[0], :final_soil_map.shape[1]] = final_soil_map


#changing shape of an array
final_climate_data = np.zeros_like(combined_classes_lulc)

final_climate_data[:climate_data.shape[0], :climate_data.shape[1]] = climate_data 

#changing shape of an array
final_slope_map = np.zeros_like(combined_classes_lulc)

final_slope_map[:combined_classes_slope.shape[0], :combined_classes_slope.shape[1]] = combined_classes_slope

weight_lulc = 0.13003261
weight_climate = 0.47091777   
weight_soil = 0.31333322   
weight_slope = 0.0857164   

final_suitability_map = (weight_soil * new_final_soil_map) + (weight_climate * final_climate_data) + (weight_slope * final_slope_map) + (weight_lulc * combined_classes_lulc)

#Load study area shapefile
study_area = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp') 

#in mg/kg
boundaries = [0, 0.25, 0.5, 0.6, 1]  # Adjust these boundaries as needed
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
contour = ax.contourf(lon_new, lat_new, final_suitability_map, cmap=cmap, levels=boundaries, norm=norm)
study_area.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1.0)

# Label each polygon with the respective label from the metadata
for idx, row in study_area.iterrows():
    # Assuming the labels are in a column named 'label_column_name'
    label = row['UGA_adm4_9']
    # Get the centroid of the polygon to place the label
    centroid = row.geometry.centroid
    ax.text(centroid.x, centroid.y, label, ha='center', va='center', fontsize=5, bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

# Add a custom legend to the plot
legend = ax.legend(handles=patches, loc='upper left', title="Legend")


# Add colorbar, THE KEYWORD ARGUMENT FOR TICKS HAS BEEN CHANGED TO 'TICKS' FROM 'BOUNDARIES'
#cbar = plt.colorbar(contour, ax=ax, label='Suitability Map Values', norm=norm, boundaries=boundaries, ticks=ticks)
#cbar.set_ticklabels(labels)  # Apply the custom labels to the colorbar
ax.set_title('MAP OF BUIKWE DISTRICT SHOWING DIFFERENT CLASSES \n OF SUITABILITY FOR THE GROWTH OF COCOA')
ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.4)
plt.xlim(32.7, 33.4)
plt.ylim(0, 0.6)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()


 # Replace 'your_dem_file.tif' with the path to your DEM file
#dem_path = r'G:\soil aluminium map(0_20)cm\alu.tif'
#NB: DEM has to be in a projected CS for the slope calculations

#with rasterio.open(dem_path) as dem:
#    # Read the DEM data, resampling as needed
#    elevation = dem.read(1, resampling=Resampling.bilinear)

    # Get metadata for later use
#    meta = dem.meta

# Write the slope data to a new file
#slope_path = r'G:\soil aluminium map(0_20)cm\final_suitability.tif'
#with rasterio.open(slope_path, 'w', **meta) as dst:
#    dst.write(final_suitability_map.astype(rasterio.float32), 1)