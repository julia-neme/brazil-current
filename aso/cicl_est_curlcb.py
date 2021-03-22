import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import cm

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')
u = osc['u'].isel(time=slice(501, -1)).sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()
v = osc['v'].isel(time=slice(501, -1)).sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()

# Construyo la serie y la climatologia de la velocidad proyectada

d = np.rad2deg(np.arctan2(v.mean(dim='time').item(), u.mean(dim='time').item()))
n = np.sqrt(u.mean(dim='time').item()**2 + v.mean(dim='time').item()**2)
u_B = u.mean(dim='time').item()/n
v_B = v.mean(dim='time').item()/n
P = np.empty(len(u.time))
for i in range(0, len(u.time)):
    P[i] = np.dot([u_B, v_B], [u[i].item(), v[i].item()])

u_clim = u.groupby('time.month').mean('time')
v_clim = v.groupby('time.month').mean('time')

d = np.rad2deg(np.arctan2(v_clim.mean().item(), u_clim.mean().item()))
n = np.sqrt(u_clim.mean().item()**2 + v_clim.mean().item()**2)
u_clim_B = u_clim.mean().item()/n
v_clim_B = v_clim.mean().item()/n
P_clim = np.empty(len(u_clim))
for i in range(0, len(u_clim)):
    P_clim[i] = np.dot([u_clim_B, v_clim_B], [u_clim[i].item(), v_clim[i].item()])

# EOFs de la climatologia del rotor del viento

ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
curl_ncep = ncep['curl'].sel(time=slice('1999-01-01', '2015-12-31')).sel(lat=slice(-15, -50)).sel(lon=slice(-70, 30))
curl_ncep_clim = curl_ncep.groupby('time.month').mean('time')

cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')
curl_cfsr = cfsr['curl'].sel(time=slice('1999-01-01', '2015-12-31')).sel(lat=slice(-15, -50)).sel(lon=slice(-70, 30))
curl_cfsr_clim = curl_cfsr.groupby('time.month').mean('time')

from eofs.standard import Eof

# ncep

coslat = np.cos(np.deg2rad(curl_ncep_clim.lat.values)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(curl_ncep_clim.values, weights=wgts)
var = solver.varianceFraction()

plt.bar(np.arange(0, len(var), 1), var*100); plt.show()

n = 1
eof_ncep = solver.eofs(neofs=n, eofscaling=2)
pc_ncep = solver.pcs(npcs=n, pcscaling=1)
vf_ncep = var[:n]

# cfsr

coslat = np.cos(np.deg2rad(curl_cfsr_clim.lat.values)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(curl_cfsr_clim.values, weights=wgts)
var = solver.varianceFraction()

plt.bar(np.arange(0, len(var), 1), var*100); plt.show()

n = 1
eof_cfsr = solver.eofs(neofs=n, eofscaling=2)
pc_cfsr = solver.pcs(npcs=n, pcscaling=1)
vf_cfsr = var[:n]

# Comparo los componentes principales y el ciclo de la CB

# Figuras

def plotter(lat, lon, C, serie):
    x, y = np.meshgrid(lon, lat)

    gs = gridspec.GridSpec(10, 5)
    fig = plt.figure(figsize=(10, 7))
    plt.clf()

    ax1 = plt.subplot(gs[1:8, :], projection=ccrs.PlateCarree())
    ax1.set_extent([-70, 30, -50, -15], crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
    ax1.add_feature(tierra)
    ax1.coastlines(resolution='50m')
    gl = ax1.gridlines(color='white', linestyle='--', linewidth=.5, draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    cl = ax1.contourf(x, y, C*1e7, np.arange(-1,1.05,.05), cmap='coolwarm', transform=ccrs.PlateCarree(), extend='both')
    cbar = fig.colorbar(cl, ax=ax1, orientation='vertical', shrink=.7, pad=.08, label='[10$^{-7}$ Pa/m]')
    ax1.text(-0.15, 1, '(a)', transform=ax1.transAxes, size=20)

    ax2 = plt.subplot(gs[8:10, :-1])
    ax2.plot(-P_clim*100, linewidth = 1.5, color='k')
    ax2.tick_params(axis='y', labelcolor='k')
    ax3 = ax2.twinx()
    ax3.plot(serie, linewidth = 1.5, color='r')
    ax3.tick_params(axis='y', labelcolor='r')
    ax2.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11])
    ax2.set_xticklabels(['Ene', 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.', 'Dic.'])
    ax3.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
    ax2.set_ylabel('Vel.proy. [cm/s]', color='k')
    ax2.grid(linestyle = ':', color='grey')
    ax2.text(-0.15, 1, '(b)', transform=ax2.transAxes, size=20)

    plt.tight_layout()
    fig.subplots_adjust(hspace=.2, wspace=.2)

plotter(curl_ncep_clim.lat, curl_ncep_clim.lon, eof_ncep1[0,:,:], pc_ncep1)
plt.savefig('/home/bock/Documents/tesis/resultados_1/ncep_cb.png', bbox_inches='tight'); plt.show()
plotter(curl_cfsr_clim.lat, curl_cfsr_clim.lon, eof_cfsr1[0,:,:], pc_cfsr1)
plt.savefig('/home/bock/Documents/tesis/resultados_1/cfsr_cb.png', bbox_inches='tight'); plt.show()

def plotterchico(lat, lon, X, cmap, clevs, units, txt):
    x, y = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(9,7))
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

    pc = ax.contourf(x, y, X, clevs, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax, shrink=0.9, orientation='horizontal')
    cbar.ax.set_title(units)
    ax.text(-0.1, 1, txt, transform=ax.transAxes, size=15)
    return fig, ax

plotterchico(curl_ncep.lat, curl_ncep.lon, eof_ncep1[0,:,:]*1e7, 'coolwarm', np.arange(-1, 1.1, 0.1), '[10$^{-7}$ Pa/m]', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/eof1_ncep_l.png', bbox_inches='tight'); plt.show()
plotterchico(curl_cfsr.lat, curl_cfsr.lon, eof_cfsr1[0,:,:]*1e7, 'coolwarm', np.arange(-1, 1.1, 0.1), '[10$^{-7}$ Pa/m]', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/eof1_cfsr_l.png', bbox_inches='tight'); plt.show()

mth = np.arange(0, 12, 1)
f, ax = plt.subplots(figsize=(15,5))
#ax.plot(mth, pc_ncep[:,0], color='k', linewidth=1.5, linestyle='--', label='2009-2015')
#ax.plot(mth, pc_cfsr[:,0], color='r', linewidth=1.5, linestyle='--', label='2009-2015')
ax.plot(mth, pc_ncep[:,0], color='k', linewidth=1.5, label='1999-2015')
ax.plot(mth, pc_cfsr[:,0], color='r', linewidth=1.5, label='1999-2015')
ax1 = ax.twinx()
ax1.set_ylabel('[cm/s]', size=18)
ax1.plot(mth, -P_clim*100, color='b', linewidth=2)
#ax1.plot(mth, -P_clim1*100, color='b', linewidth=2, linestyle='--')
ax.grid(color='grey', linestyle=':')
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11])
ax.set_xticklabels(['Ene.', 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.', 'Dic.'])
f.subplots_adjust(hspace=0.05)
ax.tick_params(labelsize=18)
ax1.tick_params(labelsize=18)
plt.savefig('/home/bock/Documents/k.png', bbox_inches='tight')
plt.show()

# Correlaciones

R_ncep = np.empty(6); R_ncep[0] = np.corrcoef(-P_clim1, pc_ncep[:,0])[0,1]
R_cfsr = np.empty(6); R_cfsr[0] = np.corrcoef(-P_clim1, pc_cfsr[:,0])[0,1]
for i in range(1, 6):
    kk = np.empty(12)
    kk[0:i] = pc_ncep[-i:, 0]
    kk[i:] = pc_ncep[0:-i, 0]
    R_ncep[i] = np.corrcoef(-P_clim1, kk)[0,1]
    kk = np.empty(12)
    kk[0:i] = pc_cfsr[-i:, 0]
    kk[i:] = pc_cfsr[0:-i, 0]
    R_cfsr[i] = np.corrcoef(-P_clim1, kk)[0,1]

mth = np.arange(0, 12, 1)

f, ax = plt.subplots(figsize=(10,5))
ax.plot(mth, pc_ncep[:,0], color='k', linestyle='--', label='1979-2015')
ax.plot(mth, pc_cfsr[:,0], color='r', linestyle='--', label='1979-2015')
ax.plot(mth, pc_ncep1[:,0], color='k', label='2009-2015')
ax.plot(mth, pc_cfsr1[:,0], color='r', label='2009-2015')
ax.grid(color='grey', linestyle=':')
ax.set_xticks([1,2,3,4,5,6,7,8,9,10])
ax.set_xticklabels([ 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.'], size=15)
f.subplots_adjust(hspace=0.05)
plt.savefig('/home/bock/Documents/tesis/resultados_1/pcs_ambosper.png', bbox_inches='tight')
plt.show()
