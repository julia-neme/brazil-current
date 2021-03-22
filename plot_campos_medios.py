import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

inp = input('Datos con path: ')
dat = xarray.open_dataset(inp).squeeze()
A = dat['curl'] #
AM = A.mean(dim='time')
ASTD = A.std(dim='time')

x, y = np.meshgrid(dat.lon, dat.lat)
lims = np.nanmax(np.abs(AM))*1e7

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([dat.lon[0], dat.lon[-1], dat.lat[-1], dat.lat[0]],
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

pc = ax.pcolormesh(x, y, AM*1e7, transform=ccrs.PlateCarree(), cmap='jet', #
                   vmin=-3, vmax=3)
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('10$^{-7}$ Pa/m') #
ax.set_title('Rotor medio ERA-Interim 2009-2015') #
plt.savefig('figuras/curl_medio2009-2015_eraint.png', bbox_inches='tight') #
plt.close()


fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([dat.lon[0], dat.lon[-1], dat.lat[-1], dat.lat[0]],
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

pc = ax.pcolormesh(x, y, ASTD*1e7, transform=ccrs.PlateCarree(), cmap='afmhot') #
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('10$^{-7}$ Pa/m') #
ax.set_title('STD del rotor ERA-Interim 2009-2015') #
plt.savefig('figuras/curl_std2009-2015_eraint.png', bbox_inches='tight') #
plt.close()
