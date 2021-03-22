import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ncep_v1/datos/NCEP1_wind_daily_1948-2015.nc') #
A = dat['uwnd'] #.sel(time=slice('2009-01-01','2015-12-31')) 
AM = A.mean(dim='time')
ASTD = A.std(dim='time')

x, y = np.meshgrid(dat.lon, dat.lat)
lims = np.nanmax(np.abs(AM))

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([dat.lon[0]-360, dat.lon[-1], dat.lat[-1], dat.lat[0]],
              crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='grey',
         facecolor='lightgrey')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')

grid_style = dict(color='white', linestyle='--', linewidth=.3)
gl = ax.gridlines(dravw_labels=True, **grid_style)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

pc = ax.pcolormesh(x, y, AM, transform=ccrs.PlateCarree(), cmap='coolwarm', #
                   vmin=-lims, vmax=lims)
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('m/s') #
ax.set_title('Viento zonal medio NCEP vI 2009-2015') #
plt.savefig('/home/bock/Documents/tesis/vientos/ncep_v1/figuras/uwnd_medio2009-2015_ncepv1.png', bbox_inches='tight') #
plt.close()


fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([dat.lon[0]-360, dat.lon[-1], dat.lat[-1], dat.lat[0]],
              crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='grey',
         facecolor='lightgrey')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')

grid_style = dict(color='white', linestyle='--', linewidth=.3)
gl = ax.gridlines(draw_labels=True, **grid_style)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

pc = ax.pcolormesh(x, y, ASTD, transform=ccrs.PlateCarree(), cmap='afmhot') #
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('m/s') #
ax.set_title('STD del viento zonal NCEP vI 2009-2015') #
plt.savefig('/home/bock/Documents/tesis/vientos/ncep_v1/figuras/uwnd_std2009-2015_cfsr.png', bbox_inches='tight') #
plt.close()
