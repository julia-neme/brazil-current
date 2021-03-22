import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
import shapely.geometry as sgeom
import imageio
import os

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')

# Cargo los datos #############################################################

serie = osc['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50)).mean(dim='longitude')[:,0]
serie1 = osc['v'].sel(latitude=-30, method='nearest').sel(longitude=slice(-49, -45)).sel(longitude=slice(-48, -47)).mean(dim='longitude')[:,0]
u = osc['u'].sel(latitude=slice(-25,-45), longitude=slice(-65,-35))
v = osc['v'].sel(latitude=slice(-25,-45), longitude=slice(-65,-35))
sp = np.sqrt(u**2 + v**2)
#ssta = gha['sst_anom'].sel(time = osc.time, method='nearest').sel(lon=slice(-65,-35), lat=slice(-25,-45))

x1, y1 = np.meshgrid(u.longitude, u.latitude)
#x2, y2 = np.meshgrid(ssta.lon, ssta.lat)
clevs_sp = np.arange(0, 0.55, 0.05)
clevs_vv = np.arange(-.3, .35, .05)
#clevs_ss = np.arange(-5, 5.5, 0.5)
kk = np.ones(len(osc.time))

def plotter(t):

    gs = gridspec.GridSpec(14, 8)
    fig = plt.figure(figsize=(16, 16))
    plt.clf()

    yr = str(osc['time.year'][t].values)
    mn = str(osc['time.month'][t].values).zfill(2)
    dy = str(osc['time.day'][t].values).zfill(2)
    plt.suptitle(yr + " " + mn + " " + dy)

    ax1 = plt.subplot(gs[0:9, 0:4], projection=ccrs.PlateCarree())
    ax1.set_title('Intensidad (colores)')
    ax1.set_extent([-58, -44, -37, -26], crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
    ax1.add_feature(tierra)
    ax1.coastlines(resolution='50m')
    gl = ax1.gridlines(color='white', linestyle='--', linewidth=.5, draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    box = sgeom.box(minx=-51.66, maxx=-50, miny=-34.33, maxy=-34.33)
    ax1.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='black', linewidth=2)
    box = sgeom.box(minx=-48, maxx=-47, miny=-30, maxy=-30)
    ax1.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='black', linewidth=2)

    cl = ax1.contourf(x1, y1, sp[t,0,:,:].values, clevs_sp, cmap='rainbow', transform=ccrs.PlateCarree(), extend='both')
    qvr = ax1.quiver(x1[::2,::2], y1[::2,::2], u[t,0,::2,::2].values, v[t,0,::2,::2].values,
                     transform=ccrs.PlateCarree())

    cbar = fig.colorbar(cl, ax=ax1, orientation='horizontal', shrink=.7, pad=.08, label='m/s')
    ax1.quiverkey(qvr, .2, 0.8, 1, '1 m/s', labelpos='E')

    ax2 = plt.subplot(gs[0:9, 4:8], projection=ccrs.PlateCarree())
    ax2.set_title('Vel. meridional (colores)')
    ax2.set_extent([-58, -44, -37, -26], crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
    ax2.add_feature(tierra)
    ax2.coastlines(resolution='50m')
    gl2 = ax2.gridlines(color='white', linestyle='--', linewidth=.5, draw_labels=True)
    gl2.xlabels_top = gl2.ylabels_left = gl2.ylabels_right = False
    gl2.xformatter = LONGITUDE_FORMATTER
    gl2.yformatter = LATITUDE_FORMATTER
    box = sgeom.box(minx=-51.66, maxx=-50, miny=-34.33, maxy=-34.33)
    ax2.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='black', linewidth=2)
    box = sgeom.box(minx=-48, maxx=-47, miny=-30, maxy=-30)
    ax2.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='black', linewidth=2)

    cl2 = ax2.contourf(x1, y1, v[t,0,:,:].values, clevs_vv, cmap='RdBu_r', transform=ccrs.PlateCarree(), extend='both')
    qvr2 = ax2.quiver(x1[::2,::2], y1[::2,::2], u[t,0,::2,::2].values, v[t,0,::2,::2].values,
                      transform=ccrs.PlateCarree())

    cbar2 = fig.colorbar(cl2, ax=ax2, orientation='horizontal', shrink=.7, pad=.08, label='ms/s')
    ax2.quiverkey(qvr2, .2, 0.8, 1, '1 m/s', labelpos='E')

    ax3 = plt.subplot(gs[9:11, 0:8])
    ax3.plot(osc['time'], serie1.values, 'k', linewidth=1)
    ax3.plot(osc['time'], serie1.mean().values*kk, 'grey', linewidth=1.5)
    ax3.plot(osc['time'], serie1.mean().values + serie1.std().values*kk, 'grey', linewidth=1.5, linestyle='--')
    ax3.plot(osc['time'], serie1.mean().values - serie1.std().values*kk, 'grey', linewidth=1.5, linestyle='--')
    ax3.plot(osc['time'][t], serie1[t].values, color='k', marker='o')
    plt.ylabel('Vel. meridional')
    ax3.set_title('30S')

    ax4 = plt.subplot(gs[12:14, 0:8])
    ax4.plot(osc['time'], serie.values, 'k', linewidth=1)
    ax4.plot(osc['time'], serie.mean().values*kk, 'grey', linewidth=1.5)
    ax4.plot(osc['time'], serie.mean().values + serie.std().values*kk, 'grey', linewidth=1.5, linestyle='--')
    ax4.plot(osc['time'], serie.mean().values - serie.std().values*kk, 'grey', linewidth=1.5, linestyle='--')
    ax4.plot(osc['time'][t], serie[t].values, color='k', marker='o')
    plt.ylabel('Vel. meridional')
    ax4.set_title('34.5S')

    plt.tight_layout()
    fig.subplots_adjust(hspace=.2, wspace=.2)
    plt.savefig('/home/bock/Documents/tesis/resultados/fig_for_an/ee_'+yr+mn+dy+'.png')

for i in range(0, len(osc['time'])):
    plotter(i)
    plt.close()

imagenes = []
os.chdir('/home/bock/Documents/tesis/resultados/fig_for_an/')
files = sorted(os.listdir())
for f in files:
    imagenes.append(imageio.imread(f))

imageio.mimsave('/home/bock/Documents/tesis/resultados/animaciones/ver_ee.gif', imagenes, fps=6)
