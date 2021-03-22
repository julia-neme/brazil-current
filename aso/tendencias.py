import numpy as np
import xarray
from scipy import stats
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
gh = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrsst_atlsur_2009_2015.nc')

U = osc['u_anom'].sel(lat=slice(-20,-56), lon=slice(-80,-40))
V = osc['v_anom'].sel(lat=slice(-20,-56), lon=slice(-80,-40))
A = np.sqrt(U.values[:,:,:]**2 + V.values[:,:,:]**2)
mask = ~np.isnan(A)
t = np.arange(0, len(U.time), 1)
LT = np.empty([len(U.lat), len(U.lon)])
for i in range(0, len(U.lat)):
    for j in range(0, len(U.lon)):
        C = A[:, i, j]; C = C[mask[:,i,j]]
        if np.size(C) == 0:
            LT[i,j] = np.nan
        else:
            LF = stats.linregress(t[mask[:,i,j]], C)
            if LF[3] <= 0.05:
                LT[i,j] = LF[0]*72*10
            else:
                LT[i,j] = np.nan

x, y = np.meshgrid(U.lon, U.lat)

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

clr = ax.pcolormesh(x, y, LT, vmin=-0.4, vmax=0.4,
                  transform=ccrs.PlateCarree(), cmap='coolwarm')
cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('m/s por década')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/oscar_anom_trend.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

S = gh['analysed_sst'].sel(lat=slice(-20,-56), lon=slice(-80,-40))
mask = ~np.isnan(S)
t = np.arange(0, len(S.time), 1)
LT = np.empty([len(S.lat), len(S.lon)])
for i in range(0, len(S.lat)):
    for j in range(0, len(S.lon)):
        C = S[:, i, j]; C = C[mask[:,i,j]]
        if np.size(C) == 0:
            LT[i,j] = np.nan
        else:
            LF = stats.linregress(t[mask[:,i,j]], C)
            if LF[3] <= 0.05:
                LT[i,j] = LF[0]*365*10
            else:
                LT[i,j] = np.nan

x, y = np.meshgrid(S.lon, S.lat)

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

clr = ax.pcolormesh(x, y, LT, vmin=-2, vmax=2,
                  transform=ccrs.PlateCarree(), cmap='coolwarm')
cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('$^\circ$ C por década')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/ghrsst_trend.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()
