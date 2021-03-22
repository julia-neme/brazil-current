import numpy as np
import xarray
from scipy import stats
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

inp = input('Datos con path: ')
dat = xarray.open_dataset(inp).squeeze()

t = np.arange(0, len(dat.time), 1)
LT = np.empty([len(dat.lat), len(dat.lon)])
for i in range(0, len(dat.lat)):
    for j in range(0, len(dat.lon)):
        LF = stats.linregress(t, dat.curl[:,i,j]*1e7)
        if LF[3] <= 0.1:
            LT[i,j] = LF[0]*365
        else:
            LT[i,j] = np.nan

x, y = np.meshgrid(dat.lon, dat.lat)
lims = np.nanmax(np.abs(LT))

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

pc = ax.pcolormesh(x, y, LT, transform=ccrs.PlateCarree(), cmap='jet',
                   vmin=-.3, vmax=.3)
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('10$^{-7}$ Pa/m por año') #
ax.set_title('Tendencia del rotor CCMP 2009-2015') #
plt.savefig('figuras/trend_2009-2015_ccmp.png', bbox_inches='tight') #
plt.close()
