import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

gh = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghsst_atlsur_2009_2015.nc')

x = gh['analysed_sst'].sel(lat=slice(-20,-56), lon=slice(-80,-40))
msn = np.empty(np.shape(x[0,:,:]))
for i in range(0, len(x.lat)):
    for j in range(0, len(x.lon)):
        msn[i,j] = np.isnan(x[:,i,j]).sum()

x, y = np.meshgrid(x.lon, x.lat)
fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-80, -40, -56, -20], crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')

grid_style = dict(color='white', linestyle='--', linewidth=.3)
gl = ax.gridlines(draw_labels=True, **grid_style)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

clr = ax.pcolormesh(x, y, msn, transform=ccrs.PlateCarree())
cbar = fig.colorbar(clr, ax=ax, shrink=.7)

plt.show()
plt.close()
