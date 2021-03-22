import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.animation as anim
import cmocean
import imageio

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
gh = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')
u = osc['u_anom'].sel(lat=slice(-25,-45), lon=slice(-65,-35))
v = osc['v_anom'].sel(lat=slice(-25,-45), lon=slice(-65,-35))

x_o, y_o = np.meshgrid(u.lon, u.lat)
x_g, y_g = np.meshgrid(gh.lon.sel(lon=slice(-65,-35)), gh.lat.sel(lat=slice(-25,-45)))

sst = gh['sst_anom'].sel(time = osc.time, method='nearest').sel(lon=slice(-65,-35), lat=slice(-25,-45))

clevs = np.arange(-5, 6, 1)
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

    clr = ax.contourf(x_g, y_g, sst[i,:,:].values, clevs,
                      transform=ccrs.PlateCarree(), cmap='RdBu_r', extend='both')
    qvr = ax.quiver(x_o[::2,::2], y_o[::2,::2], u[i,::2,::2].values, v[i,::2,::2].values,
                    units='xy', scale=2/3/111139, transform=ccrs.PlateCarree())
    cbar = fig.colorbar(clr, ax=ax, shrink=.7)
    cbar.ax.set_ylabel('$^\circ$ C')
    ax.quiverkey(qvr, 1.1, 1.1, 1, '1.5 m/s', labelpos='E')
    yr = str(osc['time.year'][i].values)
    mn = str(osc['time.month'][i].values).zfill(2)
    dy = str(osc['time.day'][i].values).zfill(2)
    ax.set_title(yr + " " + mn + " " + dy)

    salida = '/home/bock/Documents/tesis/resultados/figs/osc/oscar_ghsst_anom_' + yr + mn + dy + '.png'
    plt.savefig(salida, dpi=250)

for i in range(0, len(osc.time)):
    grafico_campos(i)
    plt.close("all")

imagenes = []
os.chdir('/home/bock/Documents/tesis/resultados/figs/osc/')
files = sorted(os.listdir())
for f in files:
    imagenes.append(imageio.imread(f))

imageio.mimsave('/home/bock/Documents/tesis/resultados/animaciones/oscar_ghsst_anom_2009_2015.gif', imagenes, fps=3)
