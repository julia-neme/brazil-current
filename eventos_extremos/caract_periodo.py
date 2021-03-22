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
osa = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')
bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')

# Campo medio OSCAR ###########################################################

#U = osa['u'].squeeze().isel(time=slice(501,-1)).mean(dim='time')
#V = osa['v'].squeeze().isel(time=slice(501,-1)).mean(dim='time')
U = osc['u'].squeeze().mean(dim='time')
V = osc['v'].squeeze().mean(dim='time')
S = np.sqrt(U**2+V**2)
bati = bat.sel(y=slice(-40,-20), x=slice(-58,-44))
x_o, y_o = np.meshgrid(U.longitude, U.latitude)
x_b, y_b = np.meshgrid(bati['x'], bati['y'])

fig = plt.figure(figsize=(10,12))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([-55, -48, -36, -32], crs=ccrs.PlateCarree())
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

clr = ax.contourf(x_o, y_o, S*100, np.arange(0,42, 2),
                  transform=ccrs.PlateCarree(), cmap='rainbow', extend='both')
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3,
              linestyles='-', transform=ccrs.PlateCarree())
cbloc = ax.contour(x_o, y_o, S, levels = 0.1, transform=ccrs.PlateCarree(),
                   linestyles='--', colors='w', linewidths=1.5)
qvr = ax.quiver(x_o[::1,::1], y_o[::1,::1], U.values, V.values,
                units='xy', scale=0.4/111139, transform=ccrs.PlateCarree())
box = sgeom.box(minx=-51.67, maxx=-50.33, miny=-34.33, maxy=-34.33)
ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none', linewidth=2.5,
                  edgecolor='black')
ax.scatter(Vc.longitude, [Vc.latitude, Vc.latitude, Vc.latitude, Vc.latitude], marker='x', color='k', transform=ccrs.PlateCarree())
cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('cm/s', size=15)
cbar.ax.tick_params(labelsize=15)
ax.quiverkey(qvr, .18, 0.75, .4, '40 cm/s', labelpos='E', fontproperties={'size': 15})
#ax.text(-0.15, 1, '(a)', transform=ax.transAxes, size=15)

#plt.savefig('/home/bock/Documents/kk.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

# STD de OSCAR ################################################################

Sstd = np.sqrt(osa['u'].squeeze()**2 + osa['v'].squeeze()**2).std(dim='time')

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
gl.xlabel_style = {'size': 15, 'rotation': 45}
gl.ylabel_style = {'size': 15}

clr = ax.contourf(x_o, y_o, Sstd*100, np.arange(0,32, 2),
                  transform=ccrs.PlateCarree(), cmap='cubehelix_r', extend='both')
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3,
              linestyles='-', transform=ccrs.PlateCarree())

cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('cm/s', size=15)
cbar.ax.tick_params(labelsize=15)
#ax.text(-0.15, 1, '(a)', transform=ax.transAxes, size=15)

plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/oscarstd_1999-2015.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

# Ciclo anual OSCAR ###########################################################

Vc = osc['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()
Vl = osa['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze().isel(time=slice(500, -1))

Vc_clim = Vc.groupby('time.month').mean('time')
Vl_clim = Vl.groupby('time.month').mean('time')
Vc_sem = Vc.groupby('time.month').std('time')/np.sqrt(len(Vc.time)-1)
Vl_sem = Vl.groupby('time.month').std('time')/np.sqrt(len(Vl.time)-1)
mth = np.arange(0, 12, 1)

fig, ax = plt.subplots(1, figsize=(10,3))
ax.plot(mth, Vl_clim*100, color='k', linewidth=1.5, label='1999-2015')
ax.plot(mth, Vc_clim*100, color='b', linewidth=1.5, label='2009-2015')
ax.fill_between(mth, (Vl_clim-Vl_sem)*100, (Vl_clim+Vl_sem)*100, color='k', alpha=0.3)
ax.fill_between(mth, (Vc_clim-Vc_sem)*100, (Vc_clim+Vc_sem)*100, color='k', alpha=0.3)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11])
ax.set_xticklabels(['Ene', 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.', 'Dic.'])
ax.set_ylabel('V [cm/s]', size=15)
ax.tick_params(labelsize=15)
plt.grid(color='grey', linestyle=':')
#plt.legend()
#ax.text(-0.15, 1, '(b)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/fig_3_4a.png', dpi=250, bbox_inches='tight')
plt.show()

# Espectro ####################################################################
from scipy import signal
from scipy.stats import chi2
from matplotlib import ticker as mticker
def welch_spectrum(x, fs, ws, alfa):

    """
    Calcula el espectro de Welch. ws es el largo de la ventana dado en cantidad
    de datos y alfa es el nivel de significancia para los intervalos de confia-
    nza.
    """

    Fs = fs
    win = ws*Fs
    ovlp = win/2
    nfft = win

    f, Pxx = signal.welch(x, Fs, nperseg=win, noverlap=ovlp)

    # Calcula los niveles de confianza

    ddof = np.round((8/3)*len(x)/win)
    c = chi2.ppf([1-alfa/2, alfa/2], ddof)
    c = ddof/c
    CI_dw = Pxx*c[0]
    CI_up = Pxx*c[1]

    plt.figure()
    plt.plot(f, f*Pxx, color='k')
    plt.fill_between(f, f*CI_dw, f*CI_up, color='k', alpha=.5)
    plt.xscale('log')
    plt.xlabel('log(f)')
    plt.ylabel('f*PSD')
    plt.show()
    plt.close()

    PSD = [f, Pxx, CI_dw, CI_up]
    return PSD

#Hago la señal teorica de ruido blanco para ver la significancia

noise_c = np.random.normal(0, np.nanstd(Vc.values)*2., len(Vc))
noise_l = np.random.normal(0, np.nanstd(Vl.values)*2., len(Vl))

PSD_vc = welch_spectrum(Vc.values, 1, 292, 0.05)
PSD_vl = welch_spectrum(Vl.values, 1, 292, 0.05)
PSD_vc_n = welch_spectrum(noise_c, 1, 292, 0.05)
PSD_vl_n = welch_spectrum(noise_l, 1, 292, 0.05)
SL_vc = np.ones(len(PSD_vc_n[0]))*(np.sqrt(PSD_vc_n[1]/2.)*2.*noise_c.std()).mean()
SL_vl = np.ones(len(PSD_vl_n[0]))*(np.sqrt(PSD_vl_n[1]/2.)*2.*noise_l.std()).mean()

fig, ax = plt.subplots(1, figsize=(13,6))
ax.plot(5/PSD_vl[0], PSD_vl[0]*PSD_vl[1]*100, color='k', linewidth=1.5, label='1992-2015')
ax.fill_between(5/PSD_vl[0], PSD_vl[0]*PSD_vl[2]*100, PSD_vl[0]*PSD_vl[3]*100, color='k', linestyle='--', linewidth=1.5, alpha=.3)
ax.plot(5/PSD_vc[0], PSD_vc[0]*PSD_vc[1]*100, color='b', linewidth=1.5, label='2009-2015')
ax.fill_between(5/PSD_vc[0], PSD_vc[0]*PSD_vc[2]*100, PSD_vc[0]*PSD_vc[3]*100, color='b', linestyle='--', linewidth=1.5, alpha=.3)
#plt.plot(5/PSD_vl_n[0], PSD_vl_n[0]*SL_vl*100, color='k', linewidth=1.5, linestyle='--')
#plt.plot(5/PSD_vc_n[0], PSD_vl_n[0]*SL_vc*100, color='b', linewidth=1.5, linestyle='--')
ax.set_ylabel('f*PSD [cm$^{2}$/s$^{2}$]', size=15)
ax.set_xlabel('Período [días]', size=15)
plt.xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
ax.tick_params(which='both', labelsize=15)
ax.grid(which="both", color='grey', linestyle=':')
#plt.legend()
#ax.text(-0.1, 1, '(b)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/fig_3_5a.png', dpi=250, bbox_inches='tight')
plt.show()

# Campo medio y std GHRSST ####################################################

ghr = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrsst_atlsur_2009_2015.nc')
SST = ghr['analysed_sst'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).mean(dim='time')
x_g, y_g = np.meshgrid(SST.lon, SST.lat)

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
gl.xlabel_style = {'size': 15, 'rotation': 45}
gl.ylabel_style = {'size': 15}

clr = ax.contourf(x_g, y_g, SST-273, np.arange(12, 26.5, .5),
                  transform=ccrs.PlateCarree(), cmap=cm.GMT_no_green, extend='both')
b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3,
              linestyles='-', transform=ccrs.PlateCarree())

cbar = fig.colorbar(clr, ax=ax, shrink=.7)
cbar.ax.set_ylabel('$^{\circ}$ C', size=15)
cbar.ax.tick_params(labelsize=15)
#ax.text(-0.15, 1, '(b)', transform=ax.transAxes, size=15)

plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/sst_mean_2009-2015.png', dpi=250, bbox_inches='tight')
plt.show()
plt.close()

sst = ghr['analysed_sst'].sel(lat=-34.5, method='nearest').sel(lon=slice(-51.66, -50.33)).mean(dim='lon')
sstclim = sst.groupby('time.month').mean('time')
sstsem = sst.std('time')/np.sqrt(len(ghr.time)-1)

fig, ax = plt.subplots(1, figsize=(10,3))
ax.plot(mth, sstclim-273, color='r', linewidth=1.5)
ax.fill_between(mth, sstclim-273-sstsem, sstclim-273+sstsem, color='r', alpha=0.3)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11])
ax.set_xticklabels(['Ene', 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.', 'Dic.'])
ax.tick_params(labelsize=15)
ax.set_ylabel('$^{\circ}$ C', size=15)
plt.grid(color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/sst_clim.png', dpi=250, bbox_inches='tight')
plt.show()

# Espectro de ssta
ssta = sst.groupby('time.month') - sstclim-273

noise_sst = np.random.normal(0, np.nanstd(ssta.values)*2., len(ssta))
PSD_sst_n = welch_spectrum(noise_sst, 1, 730, 0.05)
SL_sst = np.ones(len(PSD_sst_n[0]))*(np.sqrt(PSD_sst_n[1]/2.)*2.*noise_sst.std()).mean()
PSD_sst = welch_spectrum(ssta, 1, 730, 0.05)

fig, ax = plt.subplots(1, figsize=(12,6))
ax.plot(1/PSD_sst[0], PSD_sst[0]*PSD_sst[1], color='r', linewidth=1.5)
ax.fill_between(1/PSD_sst[0], PSD_sst[0]*PSD_sst[2], PSD_sst[0]*PSD_sst[3], color='r', linestyle='--', linewidth=1.5, alpha=.3)
plt.plot(1/PSD_sst_n[0], PSD_sst_n[0]*SL_sst, color='r', linewidth=1.5, linestyle='--')
ax.set_ylabel('f*PSD [$^{\circ}$ C$^{2}$]', size=15)
ax.set_xlabel('Período [días]', size=15)
plt.xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
ax.tick_params(which='both', labelsize=15)
ax.grid(which="both", color='grey', linestyle=':')
plt.legend()
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/sstpsd.png', dpi=250, bbox_inches='tight')
plt.show()
