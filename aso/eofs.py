import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from mpl_toolkits.basemap import cm
from eofs.standard import Eof

def plotterchico(lat, lon, X, cmap, clevs, title, units):
    x, y = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(9,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-64, -22, -38, -22],
                  crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black',
             facecolor='white')
    ax.add_feature(tierra)
    ax.coastlines(resolution='50m')

    grid_style = dict(color='white', linestyle='--', linewidth=.3)
    gl = ax.gridlines(draw_labels=True, **grid_style)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    pc = ax.contourf(x, y, X, clevs, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
    b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3, linestyles='-',
                   transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax, shrink=0.9, orientation='horizontal')
    cbar.ax.set_ylabel(units)
    ax.set_title(title)
    return fig, ax
bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')
bati = bat.sel(y=slice(-56,-20), x=slice(-80,-40))
x_b, y_b = np.meshgrid(bati['x'], bati['y'])

dat = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_2009_2015.nc')
clim_nc = dat['curl'].groupby('time.month').mean('time').sel(lat=slice(-20, -40), lon=slice(-64,-22))

coslat = np.cos(np.deg2rad(clim_nc.lat.values)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(clim_nc.values, weights=wgts)
var = solver.varianceFraction()
plt.figure(1)
plt.bar(np.arange(0, len(var), 1), var*100)
plt.show()
plt.close()
n = input('Cuantos PC extraer: '); n = int(n)
eof = solver.eofs(neofs=n, eofscaling=2)
pc = solver.pcs(npcs=n, pcscaling=1)
vf = var[:n]

fig, ax = plotterchico(clim_nc.lat, clim_nc.lon, eof[0]*1e7, cm.GMT_no_green, np.arange(-0.5,0.51,.1), '', '10$^{-7}$ Pa/m')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/eof1_ncep.png', bbox_inches='tight'); plt.show()
fig, ax = plotterchico(clim_nc.lat, clim_nc.lon, eof[1]*1e7, cm.GMT_no_green, np.arange(-0.5,0.51,.1), '', '10$^{-7}$ Pa/m')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/eof2_ncep.png', bbox_inches='tight'); plt.show()

dat1 = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_2009_2015.nc')
clim_cf = dat1['curl'].groupby('time.month').mean('time').sel(lat=slice(-20, -40), lon=slice(-64,-22))

coslat1 = np.cos(np.deg2rad(clim_cf.lat.values)).clip(0., 1.)
wgts1 = np.sqrt(coslat1)[..., np.newaxis]
solver1 = Eof(clim_cf.values, weights=wgts1)
var1 = solver1.varianceFraction()
plt.figure(1)
plt.bar(np.arange(0, len(var1), 1), var1*100)
plt.show()
plt.close()
n = input('Cuantos PC extraer: '); n = int(n)
eof1 = solver1.eofs(neofs=n, eofscaling=2)
pc1 = solver1.pcs(npcs=n, pcscaling=1)
vf1 = var1[:n]

fig, ax = plotterchico(clim_cf.lat, clim_cf.lon, eof1[0]*1e7, cm.GMT_no_green, np.arange(-0.5,0.51,.1), '', '10$^{-7}$ Pa/m')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/eof1_cym_cfsr.png', bbox_inches='tight'); plt.show()
fig, ax = plotterchico(clim_cf.lat, clim_cf.lon, -eof1[1]*1e7, cm.GMT_no_green, np.arange(-0.5,0.51,.1), '', '10$^{-7}$ Pa/m')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/eof2_cym_cfsr.png', bbox_inches='tight'); plt.show()

f,ax = plt.subplots(1,figsize=(10,8))
ax.plot(np.arange(0,12,1), pc[:,0], color='lime', label='NCEP')
ax.plot(np.arange(0,12,1), pc1[:,0], color='orange', label='CFSR')
ax.grid(which="both", color='grey', linestyle=':')
ax.legend()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/pc1_cym.png', bbox_inches='tight'); plt.show()

f,ax = plt.subplots(1,figsize=(10,8))
ax.plot(np.arange(0,12,1), pc[:,1], color='lime', label='NCEP')
ax.plot(np.arange(0,12,1), -pc1[:,1], color='orange', label='CFSR')
ax.grid(which="both", color='grey', linestyle=':')
ax.legend()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/pc2_cym.png', bbox_inches='tight'); plt.show()
