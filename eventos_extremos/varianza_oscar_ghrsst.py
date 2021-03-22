import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_vel_atlsur_2009_2015.nc')
gh = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghsst_atlsur_2009_2015.nc')

x_o, y_o = np.meshgrid(osc.longitude, osc.latitude)
x_g, y_g = np.meshgrid(gh.lon, gh.lat)

sst = gh['analysed_sst'].sel(time = osc.time, method='nearest')

anom_osc = np.nanvar(osc['v'][:,0,:,:].values, axis=0)
anom_sst = np.nanvar(sst.values, axis=0)

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-65, -35, -45, -25],
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

pc = ax.pcolormesh(x_o, y_o, anom_osc, transform=ccrs.PlateCarree(), cmap='afmhot',
                   vmin=np.nanmin(anom_osc), vmax=.07)
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('m/s')
ax.set_title('Varianza de v (OSCAR)')
plt.savefig('/home/bock/Documents/tesis/oscar_varianza.png', bbox_inches='tight') #
plt.close()

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-65, -35, -45, -25],
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

pc = ax.pcolormesh(x_g, y_g, anom_sst, transform=ccrs.PlateCarree(), cmap='afmhot',
                   vmax=30)
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('$^{\circ}$ C')
ax.set_title('Varianza de SST (ghsst)')
plt.savefig('/home/bock/Documents/tesis/ghsst_varianza.png', bbox_inches='tight') #
plt.close()
