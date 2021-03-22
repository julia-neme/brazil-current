import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from mpl_toolkits.basemap import cm

def plotter(lat, lon, X, u, v, cmap, clevs, title, units, txt):
    x, y = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-80, 30, -60, 0],
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
    gl.xlabel_style = {'size': 15}
    gl.ylabel_style = {'size': 15}
    pc = ax.contourf(x, y, X, clevs, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
    #qv = ax.quiver(x[::4, ::4], y[::4, ::4], u[::4, ::4], v[::4, ::4], transform=ccrs.PlateCarree())
    qv = ax.quiver(x, y, u, v, transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax, shrink=0.7)
    cbar.ax.set_ylabel(units)
    ax.set_title(title)
    ax.text(-0.1, 1, txt, transform=ax.transAxes, size=15)
    return fig, ax

ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')

crl = ncep['curl'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().mean(dim='time')
tx = ncep['taux'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().mean(dim='time')
ty = ncep['tauy'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().mean(dim='time')

crl1 = ncep['curl'].sel(time=slice('2009-01-01', '2015-12-31')).squeeze().mean(dim='time')
tx1 = ncep['taux'].sel(time=slice('2009-01-01', '2015-12-31')).squeeze().mean(dim='time')
ty1 = ncep['tauy'].sel(time=slice('2009-01-01', '2015-12-31')).squeeze().mean(dim='time')

fig, ax = plotter(ncep['lat'].values, ncep['lon'].values, crl.values*1e7, tx.values, ty.values, 'coolwarm', np.arange(-2.5, 2.6, .1), '', '10$^{-7}$Pa/m', '(a)')
plt.savefig('/home/bock/Documents/kk.png', bbox_inches='tight'); plt.show()

crl = cfsr['curl'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().mean(dim='time')
tx = cfsr['taux'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().mean(dim='time')
ty = cfsr['tauy'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().mean(dim='time')

crl1 = cfsr['curl'].sel(time=slice('2009-01-01', '2015-12-31')).squeeze().mean(dim='time')
tx1 = cfsr['taux'].sel(time=slice('2009-01-01', '2015-12-31')).squeeze().mean(dim='time')
ty1 = cfsr['tauy'].sel(time=slice('2009-01-01', '2015-12-31')).squeeze().mean(dim='time')

fig, ax = plotter(cfsr['lat'].values, cfsr['lon'].values, crl1.values*1e7, tx1.values, ty1.values, 'coolwarm', np.arange(-2.5, 2.6, .1), '', '10$^{-7}$Pa/m', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/curl_cfsr_medio_c.png', bbox_inches='tight'); plt.show()
