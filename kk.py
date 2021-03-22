import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from mpl_toolkits.basemap import cm

def plotter(lat, lon, X, cmap, clevs, title, units, txt):
    x, y = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-70, 30, -50, -15],
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
    #qv = ax.quiver(x, y, u, v, transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax, shrink=0.7)
    cbar.ax.set_ylabel(units)
    ax.set_title(title)
    #ax.text(-0.1, 1, txt, transform=ax.transAxes, size=15)
    return fig, ax

ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')

ncep_j = ncep['curl'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().groupby('time.month').mean(dim='time').sel(month=6)
ncep_e = ncep['curl'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().groupby('time.month').mean(dim='time').sel(month=1)
cfsr_j = cfsr['curl'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().groupby('time.month').mean(dim='time').sel(month=6)
cfsr_e = cfsr['curl'].sel(time=slice('1999-10-01', '2015-12-31')).squeeze().groupby('time.month').mean(dim='time').sel(month=1)

plotter(cfsr_e.lat, cfsr_e.lon, cfsr_j.values*1e7, 'coolwarm', np.arange(-2.5, 2.6, .1), ' ', ' ', ' ')
plt.savefig('/home/bock/cfsr_j.png', bbox_inches ='tight')
