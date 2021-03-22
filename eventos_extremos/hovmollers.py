import xarray
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import cm

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')

v = osc['v_anom'].sel(lat=-34.5, method='nearest').sel(lon=slice(-52, -30)).sel(lon=slice(-51.66, -45))
ssta = gha['sst_anom'].sel(lon=slice(-51.67,-45), lat=-34.5).sel(time=osc['time'], method='nearest')

idxED = np.array([18, 24, 70, 75, 93, 101, 172, 177, 185, 190, 196, 201, 326, 331, 366, 375, 417, 423, 461, 467])
idxEF = np.array([63, 67, 104, 112, 272, 278, 318, 325, 359, 363, 378, 391, 426, 433, 439, 444, 470, 475, 489, 497])

x = np.arange(0, len(osc['time']), 1)
xx, yy = np.meshgrid(v.lon, x)
xxb, yyb = np.meshgrid(ssta.lon, x)

n = 106; v_torcido = []; v_torcido_lat = []; v_torcido_lon = []
for i in range(84, 84+14):
    v_torcido.append(osc['v_anom'].isel(lon=i, lat=n).values)
    v_torcido_lon.append(osc['lon'].isel(lon=i).values)
    v_torcido_lat.append(osc['lat'].isel(lat=n).values)
    n = n-1
v_torcido = np.array(v_torcido)
v_torcido_lon = np.array(v_torcido_lon)
v_torcido_lat = np.array(v_torcido_lat)

def plotter(x, lat, lon, u, v, cmap, clevs, bati, lat_i, lat_f, lon_i, lon_f, cbar_label):
    import cartopy.crs as ccrs
    from cartopy.feature import NaturalEarthFeature
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.gridspec as gridspec
    import shapely.geometry as sgeom
    import numpy as np
    import xarray

    x_o, y_o = np.meshgrid(u.longitude, v.latitude)
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
    l1 = sgeom.LineString([(-51.67,-34.33), (-45, -34.33)])
    l2 = sgeom.LineString([(-52,-35), (-47.66, -31)])
    ax.add_geometries([l1,l2], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='orange', linewidth=2)

    clr = ax.contourf(x1, y1, np.sqrt(u**2+v**2)*100, clevs,
                      transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
    cbloc = ax.contour(x1, y1, np.sqrt(u**2+v**2), levels = 0.1, transform=ccrs.PlateCarree(),
                       linestyles='--', colors='k', linewidths=1.2)
    b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3,
                    transform=ccrs.PlateCarree())
    cbar = fig.colorbar(clr, ax=ax, shrink=.7)
    cbar.ax.set_ylabel('m/s')

    return fig, ax

os = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
U = os['u'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44))[:,0,:,:].mean(dim='time')
V = os['v'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44))[:,0,:,:].mean(dim='time')
bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')
bati = bat.sel(y=slice(-56,-20), x=slice(-80,-40))

fig, ax = plotter(V.values, V.latitude, V.longitude, U, V , 'rainbow', np.arange(0, 41, 1), bati, -26,-40,-58,-44, 'cm/s')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/hovmoller_ubicacion.png', bbox_inches='tight')


xx, yy = np.meshgrid(x, v_torcido_lat)
xx1, yy1 = np.meshgrid(v.lon, x)

f, ax = plt.subplots(1, figsize=(6,10))
cl = ax.contourf(xx1, yy1, v*100, np.arange(-30,35,5), cmap='coolwarm', extend='both')
ze = ax.contour(xx1, yy1, v*100, levels=[0], colors='white')

for i in range(0, len(idxED)):
    ax.plot((xx1[0,0], xx1[0,3]), (idxED[i], idxED[i]), color='g')

for i in range(0, len(idxEF)):
    ax.plot((xx1[0,0], xx1[0,3]), (idxEF[i], idxEF[i]), color='k')
cbar = f.colorbar(cl, ax=ax, shrink=.7)
cbar.ax.set_ylabel('cm/s')
ax.set_xlabel('Longitud')
#ax.set_yticks([0,100,200,300,400,500])
#ax.set_yticklabels(['2009-01-01', '2010-05-22', '2011-10-11', '2013-03-02', '2014-07-22', '2015-12-11'], rotation=45)
ax.text(-0.3, .95, '(a)', transform=ax.transAxes, size=15)
#plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/hovmoller_loncorto.png', bbox_inches='tight')
plt.show()
plt.close()

f1, ax1 = plt.subplots(1, figsize=(15,5))
cl1 = ax1.contourf(xx, yy, v_torcido*100, np.arange(-30,35,5), cmap='coolwarm', extend='both')
ze1 = ax1.contour(xx, yy, v_torcido*100, levels=[0], colors='white')

for i in range(0, len(idxED)):
    ax1.plot((idxED[i], idxED[i]), (yy[2,0], yy[3,0]), color='g')

for i in range(0, len(idxEF)):
    ax1.plot((idxEF[i], idxEF[i]), (yy[2,0], yy[3,0]), color='k')
cbar1 = f1.colorbar(cl1, ax=ax1, shrink=.7)
cbar1.ax.set_ylabel('cm/s')
ax1.set_ylabel('Latitud')
ax1.set_xticks([0,100,200,300,400,500])
ax1.set_xticklabels(['2009-01-01', '2010-05-22', '2011-10-11', '2013-03-02', '2014-07-22', '2015-12-11'])
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/hovmoller_torcido.png', bbox_inches='tight')
plt.show()
