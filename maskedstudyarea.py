import geopandas as gpd
import rasterio
from rasterio.plot import show
from rasterio.mask import mask
from rasterio.transform import from_origin
from shapely.geometry import box
from shapely.geometry import mapping
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

#Path to the tiff image
tiff_path = r'C:\Users\USER\Desktop\project execution\geopandas\masking the classified sentinel image\Classified_Sentinel2_Image.tif'

#Load the study area shapefile
study_area = gpd.read_file(r'C:\Users\USER\Desktop\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_admin.shp') 

#read the tiff file using rasterio
raster = rasterio.open(tiff_path)

#Create a Geodata frame for visualization
geometry = [box(*raster.bounds)]
gdf = gpd.GeoDataFrame(geometry=geometry, crs=raster.crs)

#Customizing the colors
custom_colors = ['#3cff13', '#5ad1d4', '#ffcf31', '#0c0c0c']
custom_cmap = ListedColormap(custom_colors)

#converting geometry of study_area to geojson format which is acceptable by the mask() function
study_area_geojson = mapping(study_area['geometry'][0])

#masking the tiff image to the study area 
masked_data, _ = mask(dataset=raster, shapes=[study_area_geojson], crop=True)

 # Assuming it's a single-band image
plt.imshow(masked_data[0], cmap=custom_cmap) 
plt.show()










