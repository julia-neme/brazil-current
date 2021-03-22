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
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrsst_atlsur_2009_2015.nc')
bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')

def plotter(x, lat, lon, u, v, lat1, lon1, cmap, clevs, bati, lat_i, lat_f, lon_i, lon_f, cbar_label, title, kk):

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
    gl.xlabel_style = {'size': 15, 'rotation': 45}
    gl.ylabel_style = {'size': 15}
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
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_ylabel(cbar_label, size=18)
    ax.quiverkey(qvr, .1, 1.1, .4, '40 cm/s', labelpos='E', fontproperties={'size': 18})
    ax.text(0.08, .93, kk, transform=ax.transAxes, size=22)
    plt.title(title)

    return fig, ax

U = osc['u'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44)).squeeze()
V = osc['v'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44)).squeeze()
SSTA = gha['analysed_sst'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).sel(time=osc['time'], method='nearest').sel(lat=U.latitude, method='nearest').sel(lon=U.longitude, method='nearest')
AVM = osa['v_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44))
AUM = osa['u_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44))
bati = bat.sel(y=slice(-56,-20), x=slice(-80,-40))

###############################################################################
# CASO ANTICICLON
###############################################################################

AA = np.array([307, 308, 309])
aa = ['(a)', '(b)', '(c)', '(g)', '(h)', '(i)']
bb = ['(d)', '(e)', '(f)', '(j)', '(k)', '(l)']

for i in range(3, 6):
    idxED = AA+i
    n = str(i)
    yr = str(osc['time.year'][idxED[1]].values)
    mn = str(osc['time.month'][idxED[1]].values)
    dy = str(osc['time.day'][idxED[1]].values)

    SSTAm = SSTA.isel(time=idxED).mean(dim='time')
    AVMm = AVM.isel(time=idxED).mean(dim='time')
    AUMm = AUM.isel(time=idxED).mean(dim='time')

    #fig, ax = plotter(AVMm*100, AVM.lat, AVM.lon, AUMm.values, AVMm.values, AVM.lat, AVM.lon,
    #                'coolwarm', np.arange(-30, 31, 1), bati, -30.6, -39, -54, -44, 'cm/s', yr+mn+dy, aa[i-3])
    #plt.savefig('/home/bock/Documents/tesis/resultados/figs/casos_parti/AVM_' + '+'+n+ '.png', bbox_inches='tight');# plt.show()

    fig, ax = plotter(SSTAm-273, AVM.lat, AVM.lon, AUMm.values, AVMm.values, AVM.lat, AVM.lon,
                    cm.GMT_no_green, np.arange(10, 26.5, 0.5), bati, -30.6, -39, -54, -44, '$^{\circ}$ C', yr+mn+dy, bb[i-3])
    plt.savefig('/home/bock/Documents/tesis/resultados/figs/casos_parti/SST_' + '+'+n+ '.png', bbox_inches='tight'); #plt.show()
    plt.close('all')

###############################################################################
# CASO CICLON
###############################################################################

AA = np.array([354, 355, 356])
aa = ['(a)', '(b)', '(c)', '(g)', '(h)', '(i)']
bb = ['(d)', '(e)', '(f)', '(j)', '(k)', '(l)']

for i in range(6, 13):
    idxED = AA+i
    n = str(i)
    yr = str(osc['time.year'][idxED[1]].values)
    mn = str(osc['time.month'][idxED[1]].values)
    dy = str(osc['time.day'][idxED[1]].values)

    SSTAm = SSTA.isel(time=idxED).mean(dim='time')
    AVMm = AVM.isel(time=idxED).mean(dim='time')
    AUMm = AUM.isel(time=idxED).mean(dim='time')

    #fig, ax = plotter(AVMm*100, AVM.lat, AVM.lon, AUMm.values, AVMm.values, AVM.lat, AVM.lon,
    #                'coolwarm', np.arange(-30, 31, 1), bati, -30.6, -39, -54, -44, 'cm/s', yr+mn+dy, aa[i-6])
    #plt.savefig('/home/bock/Documents/tesis/resultados/figs/casos_parti/AVM_' + '+'+n+ '.png', bbox_inches='tight');# plt.show()

    fig, ax = plotter(SSTAm-273, AVM.lat, AVM.lon, AUMm.values, AVMm.values, AVM.lat, AVM.lon,
                    cm.GMT_no_green, np.arange(10, 26.5, 0.5), bati, -30.6, -39, -54, -44, '$^{\circ}$ C', yr+mn+dy, bb[i-6])
    plt.savefig('/home/bock/Documents/tesis/resultados/figs/casos_parti/SST_' + '+'+n+ '.png', bbox_inches='tight'); #plt.show()
    plt.close('all')
