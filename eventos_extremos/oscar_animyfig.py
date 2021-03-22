import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.animation as anim
import imageio
import os

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_vel_atlsur_2009_2015.nc')

x, y = np.meshgrid(osc.longitude, osc.latitude)

def grafico_campos(i):

    fig = plt.figure(figsize=(9,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-65, -35, -45, -25], crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
    ax.add_feature(tierra)
    ax.coastlines(resolution='50m')

    grid_style = dict(color='white', linestyle='--', linewidth=.3)
    gl = ax.gridlines(draw_labels=True, **grid_style)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    clr = ax.pcolormesh(x, y, osc.v[i,0,:,:].values, vmin=-0.7, vmax=0.7,
                      transform=ccrs.PlateCarree(), cmap='RdBu_r')
    qvr = ax.quiver(x[::2,::2], y[::2,::2], osc.u[i,0,::2,::2].values, osc.v[i,0,::2,::2].values,
                    units='xy', scale=2/3/111139, transform=ccrs.PlateCarree())
    cbar = fig.colorbar(clr, ax=ax, shrink=.7)
    cbar.ax.set_ylabel('m/s')
    ax.quiverkey(qvr, 1.1, 1.1, 1, '1.5 m/s', labelpos='E')
    yr = str(osc['time.year'][i].values)
    mn = str(osc['time.month'][i].values).zfill(2)
    dy = str(osc['time.day'][i].values).zfill(2)
    ax.set_title(yr + " " + mn + " " + dy)

    salida = '/home/bock/Documents/tesis/resultados/figs/oscar_' + yr + mn + dy + '.png'
    plt.savefig(salida, dpi=250)

for i in range(0, len(osc.time)):
    grafico_campos(i)
    plt.close("all")

imagenes = []
os.chdir('/home/bock/Documents/tesis/resultados/figs/')
files = sorted(os.listdir())
for f in files:
    imagenes.append(imageio.imread(f))

imageio.mimsave('/home/bock/Documents/tesis/resultados/oscar_2009_2015.gif', imagenes, fps=3)
