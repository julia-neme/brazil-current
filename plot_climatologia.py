import numpy as np
import xarray
import os
from functions import climatologia_xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ncep_v1/datos/NCEP1_wind_daily_1948-2015.nc') #

A = climatologia_xarray(dat[''])
x, y = np.meshgrid(dat.lon, dat.lat)
tit = ['Enero', 'Febrero', 'Marzo', 'Abril', 'Mayo', 'Junio', 'Julio',
       'Agosto', 'Septiembre', 'Octubre', 'Noviembre', 'Diciembre']
sve = ['ncep_curl_01.png', 'ncep_curl_02.png', 'ncep_curl_03.png',
       'ncep_curl_04.png', 'ncep_curl_05.png', 'ncep_curl_06.png',
       'ncep_curl_07.png', 'ncep_curl_08.png', 'ncep_curl_09.png',
       'ncep_curl_10.png', 'ncep_curl_11.png', 'ncep_curl_12.png']

os.chdir('/home/bock/Documents/tesis/vientos/ncep_v2')
for i in range(0, len(A.month)):
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

    pc = ax.pcolormesh(x, y, A[i,:,:], transform=ccrs.PlateCarree(),
                       cmap='coolwarm', #
                       vmin=-lims, vmax=lims)
    cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
    cbar.ax.set_ylabel('') #
    ax.set_title('Viento zonal medio NCEP vI 2009-2015 - '+tit[i]) #
    plt.savefig(sve[i], bbox_inches='tight') #
    plt.close()
del clim, lat, lon
