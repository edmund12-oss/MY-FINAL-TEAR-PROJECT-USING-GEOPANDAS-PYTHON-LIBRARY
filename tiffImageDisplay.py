import geopandas as gpd
import rasterio
from rasterio.plot import show
from shapely.geometry import box
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from rasterio.transform import from_origin
import matplotlib.patches as mpatches


#Path to the tiff image
tiff_path = r'E:\project execution\Datasets\classified sentinel image\mskdsentlimg.tif'


#read the tiff file using rasterio
raster = rasterio.open(tiff_path)

#Create a Geodata frame for visualization
geometry = [box(*raster.bounds)]
gdf = gpd.GeoDataFrame(geometry=geometry, crs=raster.crs)


#Load study area shapefile
#Load study area shapefile
cocoa_farm = gpd.read_file(r'E:\project execution\Datasets\Buikwe shapefile\Buikwe official admin boundaries shapefile\Buikwe_final.shp') 


#Display the tiff image using rasterio and matplotlib
fig, ax=plt.subplots(figsize=(8,8))
show(raster, ax=ax, cmap='viridis')
gdf.plot(ax=ax, facecolor='none', edgecolor='red')
cocoa_farm.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=0.6)


ax.set_title('A CLASSIFIED IMAGE OF BUIKWE DISTRICT')

#plt.xlim(32.7, 33.4)
#plt.ylim(0, 0.6)
plt.show()
