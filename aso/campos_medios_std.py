import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from mpl_toolkits.basemap import cm

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
gh = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrsst_atlsur_2009_2015.nc')
bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')

# Calculo de campos medios ####################################################

u_mean = osc['u'].sel(latitude=slice(-20,-56), longitude=slice(-80,-40)).mean(dim='time')
v_mean = osc['v'].sel(latitude=slice(-20,-56), longitude=slice(-80,-40)).mean(dim='time')
sst_mean = gh['analysed_sst'].sel(lat=slice(-20,-56), lon=slice(-80,-40)).mean(dim='time')
speed = np.sqrt(u_mean**2 + v_mean**2)
bati = bat.sel(y=slice(-56,-20), x=slice(-80,-40))

x_o, y_o = np.meshgrid(u_mean.longitude, u_mean.latitude)
x_g, y_g = np.meshgrid(sst_mean.lon, sst_mean.lat)
x_b, y_b = np.meshgrid(bati['x'], bati['y'])

# Figura ASO ##################################################################

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-80, -40, -56, -20], crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')
gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0),
                  draw_labels=True, color='white', linestyle='--', linewidth=.3)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
box = sgeom.box(minx=-58, maxx=-44, miny=-40, maxy=-26)
ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                  edgecolor='black')

clr = ax.contourf(x_b, y_b, speed[0,:,:], np.arange(0, .575, .025),
                  transform=ccrs.PlateCarree(), cmap='rainbow', extend='both')
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3, linestyles='-',
               transform=ccrs.PlateCarree())
qvr = ax.quiver(x_o[::2,::2], y_o[::2,::2], u_mean[0,::2,::2].values, v_mean[0,::2,::2].values,
                units='xy', scale=0.2/111139, transform=ccrs.PlateCarree())

cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('m/s')
ax.quiverkey(qvr, 1.1, 0.97, .2, '0.2 m/s', labelpos='E')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/oscar_medio.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-80, -40, -56, -20], crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')
gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0),
                  draw_labels=True, color='white', linestyle='--', linewidth=.3)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
box = sgeom.box(minx=-58, maxx=-44, miny=-40, maxy=-26)
ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                  edgecolor='black')

clr = ax.contourf(x_g, y_g, sst_mean[:,:]-273, np.arange(3, 27, 1),
                  transform=ccrs.PlateCarree(), cmap=cm.GMT_no_green, extend='both')
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3, linestyles='-',
             transform=ccrs.PlateCarree())
qvr = ax.quiver(x_o[::2,::2], y_o[::2,::2], u_mean[0,::2,::2].values, v_mean[0,::2,::2].values,
                units='xy', scale=0.2/111139, transform=ccrs.PlateCarree())

cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('$^{\circ}$ C')
ax.quiverkey(qvr, 1.1, 0.97, .2, '0.2 m/s', labelpos='E')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/oscar_ghrsst_medio.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

# Busco la CB en distintas latitudes. Criterio: supere la velocidad 0.1 m/s S #

def find_CB(lat, lon_i):
    x = speed.sel(latitude=lat, method='nearest').sel(longitude=slice(lon_i, speed.longitude[-1].values))
    idx = x>0.1; idx = idx[0].values
    i = 0
    while idx[i] != True:
        i = i+1
    else:
        n = nn = i
        while idx[nn] == True:
            nn = nn+1
    wre = np.arange(n, nn+1, 1)
    CB = x.isel(longitude=wre)
    return(CB)

lt_cb = [-36, -35, -34.5, -33, -31.5, -30, -28.5]
lg_cb = [-53.3, -52.7, -52.46, -51.05, -49.4, -48.7, -48.1]

CBs = []
for i in range(0, len(lt_cb)):
    CBs.append(find_CB(lt_cb[i], lg_cb[i]))

cb1 = xarray.DataArray(CBs[0][0,:], coords=(CBs[0].longitude,),
                       dims=['lon1'], name='CB_loc1').astype('float32')
cb1 = cb1.assign_coords(lat1= CBs[0].latitude); cb1 = cb1.expand_dims('lat1')
cb2 = xarray.DataArray(CBs[1][0,:], coords=(CBs[1].longitude,),
                       dims=['lon2'], name='CB_loc2').astype('float32')
cb2 = cb2.assign_coords(lat2= CBs[1].latitude); cb2 = cb2.expand_dims('lat2')
cb3 = xarray.DataArray(CBs[2][0,:], coords=(CBs[2].longitude,),
                       dims=['lon3'], name='CB_loc3').astype('float32')
cb3 = cb3.assign_coords(lat3= CBs[2].latitude); cb3 = cb3.expand_dims('lat3')
cb4 = xarray.DataArray(CBs[3][0,:], coords=(CBs[3].longitude,),
                       dims=['lon4'], name='CB_loc4').astype('float32')
cb4 = cb4.assign_coords(lat4= CBs[3].latitude); cb4 = cb4.expand_dims('lat4')
cb5 = xarray.DataArray(CBs[4][0,:], coords=(CBs[4].longitude,),
                       dims=['lon5'], name='CB_loc5').astype('float32')
cb5 = cb5.assign_coords(lat5= CBs[4].latitude); cb5 = cb5.expand_dims('lat5')
cb6 = xarray.DataArray(CBs[5][0,:], coords=(CBs[5].longitude,),
                       dims=['lon6'], name='CB_loc6').astype('float32')
cb6 = cb6.assign_coords(lat6= CBs[5].latitude); cb6 = cb6.expand_dims('lat6')
cb7 = xarray.DataArray(CBs[6][0,:], coords=(CBs[6].longitude,),
                       dims=['lon7'], name='CB_loc7').astype('float32')
cb7 = cb7.assign_coords(lat7= CBs[6].latitude); cb7 = cb7.expand_dims('lat7')

res = cb1.to_dataset(name='Loc_7')
res['Loc_6'] = cb2
res['Loc_5'] = cb3
res['Loc_4'] = cb4
res['Loc_3'] = cb5
res['Loc_2'] = cb6
res['Loc_1'] = cb7
res.to_netcdf('/home/bock/Documents/tesis/datos/CB_osc_medio.nc', mode='w')

# Figura regional #############################################################

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-58, -44, -40, -26], crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')
gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0),
                  draw_labels=True, color='white', linestyle='--', linewidth=.3)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

box = sgeom.box(minx=-51.67, maxx=-50.33, miny=-34.33, maxy=-34.33)
ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                  edgecolor='darkorange', linestyle='--', linewidth=2.7)

#clr = ax.contourf(x_o, y_o, speed[0,:,:], np.arange(0, 0.401, 0.01),
#                  transform=ccrs.PlateCarree(), cmap='rainbow', extend='both')
clr = ax.contourf(x_b, y_b, bati['z'].values, np.arange(-7000,200,200), transform=ccrs.PlateCarree(), cmap='bone', extend='both')
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='r', linewidths=1.6,
             linestyles='-', transform=ccrs.PlateCarree())
#cbloc = ax.contour(x_o, y_o, speed[0,:,:], levels = 0.1, transform=ccrs.PlateCarree(),
#                   linestyles='--', colors='r', linewidths=1.2)
#qvr = ax.quiver(x_o[::1,::1], y_o[::1,::1], u_mean[0,::1,::1].values, v_mean[0,::1,::1].values,
#                units='xy', scale=0.4/111139, transform=ccrs.PlateCarree())

cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('m')
#ax.quiverkey(qvr, .2, 0.8, .4, '0.4 m/s', labelpos='E')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/batimetria.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-58, -44, -40, -26], crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')
gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0),
                  draw_labels=True, color='white', linestyle='--', linewidth=.3)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

for i in range(2, 3):
    a = CBs[i].latitude; aa = CBs[i].longitude[0]; aaa = CBs[i].longitude[-2]
    box = sgeom.box(minx=aa, maxx=aaa, miny=a, maxy=a)
    ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='magenta', linestyle='--', linewidth=2.7)

clr = ax.contourf(x_g, y_g, sst_mean[:,:]-273, np.arange(10, 26, 0.5),
                  transform=ccrs.PlateCarree(), cmap=cm.GMT_no_green, extend='both')
cbloc = ax.contour(x_o, y_o, speed[0,:,:], levels = 0.1, transform=ccrs.PlateCarree(),
                   linestyles='--', colors='k', linewidths=1.2)
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='lime', linewidths=1.4,
             linestyles='-', transform=ccrs.PlateCarree())
qvr = ax.quiver(x_o[::1,::1], y_o[::1,::1], u_mean[0,::1,::1].values, v_mean[0,::1,::1].values,
                units='xy', scale=0.4/111139, transform=ccrs.PlateCarree())

cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('$^{\circ}$ C')
ax.quiverkey(qvr, .2, 0.8, .4, '0.4 m/s', labelpos='E')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/oscar_ghrsst_medio_regional.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

# Calculo de campos de STD ####################################################

u_std = osc['u'].sel(latitude=slice(-20,-56), longitude=slice(-80,-40)).std(dim='time')
v_std = osc['v'].sel(latitude=slice(-20,-56), longitude=slice(-80,-40)).std(dim='time')
sst_std = gh['analysed_sst'].sel(lat=slice(-20,-56), lon=slice(-80,-40)).std(dim='time')
speed_std = np.sqrt(u_std**2 + v_std**2)

# Figura ASO ##################################################################

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-80, -40, -56, -20], crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')
gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0),
                  draw_labels=True, color='white', linestyle='--', linewidth=.3)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
box = sgeom.box(minx=-58, maxx=-44, miny=-40, maxy=-26)
ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                  edgecolor='black')

clr = ax.contourf(x_o, y_o, speed_std[0,:,:], np.arange(0, .52, .02),
                  transform=ccrs.PlateCarree(), cmap='cubehelix_r', extend='both')
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3, linestyles='-',
               transform=ccrs.PlateCarree())
cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('m/s')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/oscar_std.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-80, -40, -56, -20], crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')
gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0),
                  draw_labels=True, color='white', linestyle='--', linewidth=.3)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
box = sgeom.box(minx=-58, maxx=-44, miny=-40, maxy=-26)
ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                  edgecolor='black')

clr = ax.contourf(x_g, y_g, sst_std[:,:], np.arange(0, 6.05, .05),
                  transform=ccrs.PlateCarree(), cmap='cubehelix_r', extend='both')
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3, linestyles='-',
               transform=ccrs.PlateCarree())
cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('$^{\circ}$ C')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/ghrsst_std.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()
