import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
import shapely.geometry as sgeom
import cmocean
from mpl_toolkits.basemap import cm

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
osa = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')
bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')

def plotter(x, lat, lon, u, v, lat1, lon1, cmap, clevs, bati, lat_i, lat_f, lon_i, lon_f, cbar_label, title):

    x_o, y_o = np.meshgrid(lon1, lat1)
    x1, y1 = np.meshgrid(lon, lat)
    x_b, y_b = np.meshgrid(bati['x'], bati['y'])

    fig = plt.figure(figsize=(10,12))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([lon_i, lon_f, lat_f, lat_i], crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
    ax.add_feature(tierra)
    ax.coastlines(resolution='50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0),
                      draw_labels=True, color='white', linestyle='--', linewidth=.3)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    box = sgeom.box(minx=-51.67, maxx=-50.33, miny=-34.33, maxy=-34.33)
    ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none', linewidth=2.5,
                      edgecolor='black')

    clr = ax.contourf(x1, y1, x, clevs,
                      transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
    b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=2,
                   linestyles='-', transform=ccrs.PlateCarree())
    qvr = ax.quiver(x_o, y_o, u, v,
                    units='xy', scale=0.4/111139, transform=ccrs.PlateCarree())

    cbar = fig.colorbar(clr, ax=ax, shrink=.7)
    cbar.ax.set_ylabel(cbar_label)
    ax.quiverkey(qvr, .2, 0.8, .4, '40 cm/s', labelpos='E')

    plt.title(title)

    return fig, ax
bati = bat.sel(y=slice(-56,-20), x=slice(-80,-40))

d = ['2009-10-28', '2011-10-21', '2012-10-19', '2013-06-03', '2013-12-12', '2014-05-16',
     '2014-09-09', '2010-03-09', '2014-11-16', '2012-09-25', '2013-07-13', '2014-02-04',
     '2014-10-12']

for i in range(0, len(d)):
    U = osc['u'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44)).squeeze().sel(time=d[i], method='nearest')
    V = osc['v'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44)).squeeze().sel(time=d[i], method='nearest')
    SSTA = gha['sst_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).sel(time=d[i], method='nearest')
    AVM = osa['v_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).sel(time=d[i], method='nearest')
    AUM = osa['u_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).sel(time=d[i], method='nearest')

    fig, ax = plotter(AVM*100, AVM.lat, AVM.lon, U.values, V.values, U.latitude, V.longitude, 'coolwarm', np.arange(-40, 42, 2),
                    bati, -26, -40, -58, -44, 'cm/s', str(U['time.year'].values)+'-'+str(U['time.month'].values)+'-'+str(U['time.day'].values))
    plt.savefig('/home/bock/Documents/tesis/resultados/figs/comparacion_obs_mp/'+ d[i] + '_Va.png', bbox_inches='tight')

    fig, ax = plotter(SSTA, SSTA.lat, SSTA.lon, U.values, V.values, U.latitude, V.longitude, cm.GMT_no_green, np.arange(-2.3, 2.4, .1),
                    bati, -26, -40, -58, -44, '$^{\circ}$ C', d[i])
    plt.savefig('/home/bock/Documents/tesis/resultados/figs/comparacion_obs_mp/'+ d[i] + '_SSTA.png', bbox_inches='tight')

di = ['2014-12-01', '2009-09-22', '2012-05-19', '2012-12-11', '2013-10-31']
df = ['2014-12-05', '2009-10-17', '2012-05-29', '2012-12-12', '2013-11-10']

for i in range(0, len(di)):
    U = osc['u'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44)).squeeze().sel(time=slice(di[i], df[i])).mean(dim='time')
    V = osc['v'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44)).squeeze().sel(time=slice(di[i], df[i])).mean(dim='time')
    SSTA = gha['sst_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).sel(time=slice(di[i], df[i])).mean(dim='time')
    AVM = osa['v_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).sel(time=slice(di[i], df[i])).mean(dim='time')
    AUM = osa['u_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).sel(time=slice(di[i], df[i])).mean(dim='time')

    fig, ax = plotter(AVM*100, AVM.lat, AVM.lon, U.values, V.values, U.latitude, V.longitude, 'coolwarm', np.arange(-40, 42, 2),
                    bati, -26, -40, -58, -44, 'cm/s', di[i] + df[i])
    plt.savefig('/home/bock/Documents/tesis/resultados/figs/comparacion_obs_mp/'+ di[i] + '_Va.png', bbox_inches='tight')

    fig, ax = plotter(SSTA, SSTA.lat, SSTA.lon, U.values, V.values, U.latitude, V.longitude, cm.GMT_no_green, np.arange(-2.3, 2.4, .1),
                    bati, -26, -40, -58, -44, '$^{\circ}$ C', di[i] + df[i])
    plt.savefig('/home/bock/Documents/tesis/resultados/figs/comparacion_obs_mp/'+ di[i] + '_SSTA.png', bbox_inches='tight')
