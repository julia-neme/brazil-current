import numpy as np
import xarray
from scipy import stats
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def plotter(lat, lon, X, cmap, vmi, vma, units, txt):
    x, y = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-80, 30, -60, 0],
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
    pc = ax.pcolormesh(x, y, X, vmin=vmi, vmax=vma, cmap=cmap, transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax, shrink=0.7)
    cbar.ax.set_ylabel(units)
    ax.text(-0.1, 1, txt, transform=ax.transAxes, size=15)
    return fig, ax

ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')

u = ncep['uwnd'].sel(lat=slice(0, -60)).squeeze()
t = np.arange(0, len(u.time), 1)
uncep_trend = np.empty(np.shape(u[0,:,:]))
for i in range(0, len(u.lat)):
    for j in range(0, len(u.lon)):
        kk = stats.linregress(t, u[:,i,j])
        if kk.pvalue <= 0.05:
            uncep_trend[i,j] = kk.slope*365*10
        else:
            uncep_trend[i,j] = np.nan
fig, ax = plotter(u.lat, u.lon, uncep_trend, 'coolwarm', -.8, .8, 'm/s por década', '(a)'); plt.show()
#plt.savefig('/home/bock/Documents/tesis/resultados_1/uncep_trend.png', bbox_inches='tight'); plt.close()

u = cfsr['uwnd'].sel(lat=slice(0, -60)).squeeze()
t = np.arange(0, len(u.time), 1)
uncep_trend1 = np.empty(np.shape(u[0,:,:]))
for i in range(0, len(u.lat)):
    for j in range(0, len(u.lon)):
        kk = stats.linregress(t, u[:,i,j])
        if kk.pvalue <= 0.05:
            uncep_trend1[i,j] = kk.slope*365*10
        else:
            uncep_trend1[i,j] = np.nan
fig, ax = plotter(u.lat, u.lon, uncep_trend1, 'coolwarm', -.8, .8, 'm/s por década', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/ucfsr_trend.png', bbox_inches='tight'); plt.close()

u = ncep['vwnd'].sel(lat=slice(0, -60)).squeeze()
t = np.arange(0, len(u.time), 1)
uncep_trend = np.empty(np.shape(u[0,:,:]))
for i in range(0, len(u.lat)):
    for j in range(0, len(u.lon)):
        kk = stats.linregress(t, u[:,i,j])
        if kk.pvalue <= 0.05:
            uncep_trend[i,j] = kk.slope*365*10
        else:
            uncep_trend[i,j] = np.nan
fig, ax = plotter(u.lat, u.lon, uncep_trend, 'coolwarm', -.8, .8, 'm/s por década', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/vncep_trend.png', bbox_inches='tight'); plt.close()

u = cfsr['vwnd'].sel(lat=slice(0, -60)).squeeze()
t = np.arange(0, len(u.time), 1)
uncep_trend2 = np.empty(np.shape(u[0,:,:]))
for i in range(0, len(u.lat)):
    for j in range(0, len(u.lon)):
        kk = stats.linregress(t, u[:,i,j])
        if kk.pvalue <= 0.05:
            uncep_trend2[i,j] = kk.slope*365*10
        else:
            uncep_trend2[i,j] = np.nan
fig, ax = plotter(u.lat, u.lon, uncep_trend2, 'coolwarm', -.8, .8, 'm/s por década', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/vcfsr_trend.png', bbox_inches='tight'); plt.close()

u = ncep['curl'].sel(lat=slice(0, -60)).squeeze()
t = np.arange(0, len(u.time), 1)
mask = ~np.isnan(u)
uncep_trend = np.empty(np.shape(u[0,:,:]))
for i in range(0, len(u.lat)):
    for j in range(0, len(u.lon)):
        C = u[:, i, j]; C = C[mask[:,i,j]]
        kk = stats.linregress(t[mask[:,i,j]], C)
        if kk.pvalue <= 0.05:
            uncep_trend[i,j] = kk.slope*365*10
        else:
            uncep_trend[i,j] = np.nan
fig, ax = plotter(u.lat, u.lon, uncep_trend*1e7, 'coolwarm', -.8, .8, '10$^{-7}$Pa/m por década', '(e)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/crlncep_trend.png', bbox_inches='tight'); plt.close()

u = cfsr['curl'].sel(lat=slice(0, -60)).squeeze()
t = np.arange(0, len(u.time), 1)
mask = ~np.isnan(u)
uncep_trend3 = np.empty(np.shape(u[0,:,:]))
for i in range(0, len(u.lat)):
    for j in range(0, len(u.lon)):
        C = u[:, i, j]; C = C[mask[:,i,j]]
        kk = stats.linregress(t[mask[:,i,j]], C)
        if kk.pvalue <= 0.05:
            uncep_trend3[i,j] = kk.slope*365*10
        else:
            uncep_trend3[i,j] = np.nan
fig, ax = plotter(u.lat, u.lon, uncep_trend3*1e7, 'coolwarm', -.8, .8, '10$^{-7}$Pa/m por década', '(f)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/crlcfsr_trend.png', bbox_inches='tight'); plt.close()

# tendencias del promedio zonal

t = np.arange(0,len(ncep.time), 1)
U = ncep['uwnd'].squeeze().mean(dim='lon')
Ut_ncep = np.empty(len(ncep.lat))
for i in range(0, len(ncep.lat)):
    kk = stats.linregress(t, U[:,i])
    if kk.pvalue<0.05:
        Ut_ncep[i] = kk.slope*3650
    else:
        Ut_ncep[i] = np.nan
U = cfsr['uwnd'].squeeze().mean(dim='lon')
Ut_cfsr = np.empty(len(cfsr.lat))
for i in range(0, len(cfsr.lat)):
    kk = stats.linregress(t, U[:,i])
    if kk.pvalue<0.05:
        Ut_cfsr[i] = kk.slope*3650
    else:
        Ut_cfsr[i] = np.nan

fig, ax = plt.subplots(1, figsize=(7, 10))
ax.plot(Ut_ncep, ncep['lat'], color='k', marker='*')
ax.plot(Ut_cfsr, cfsr['lat'],  color='r', marker='*')
ax.grid(color='grey', linestyle=':')
ax.set_yticks([-70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0])
ax.set_yticklabels(['70$^{\circ}$S', '65$^{\circ}$S', '60$^{\circ}$S', '55$^{\circ}$S', '50$^{\circ}$S', '45$^{\circ}$S', '40$^{\circ}$S', '35$^{\circ}$S', '30$^{\circ}$S', '25$^{\circ}$S', '20$^{\circ}$S', '15$^{\circ}$S', '10$^{\circ}$S', '5$^{\circ}$S', '0$^{\circ}$S'])
ax.tick_params(labelsize=15)
ax.set_xlabel('Tendencias de U [m/s por década]', size=15)
ax.set_ylabel('Latitud', size=15)
ax.text(-0.1, 1, '(a)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/umediocuenca_trend_b.png', dpi=250, bbox_inches='tight')
plt.show()

U = ncep['vwnd'].squeeze().mean(dim='lon')
Ut_ncep = np.empty(len(ncep.lat))
for i in range(0, len(ncep.lat)):
    kk = stats.linregress(t, U[:,i])
    if kk.pvalue<0.05:
        Ut_ncep[i] = kk.slope*3650
    else:
        Ut_ncep[i] = np.nan
U = cfsr['vwnd'].squeeze().mean(dim='lon')
Ut_cfsr = np.empty(len(cfsr.lat))
for i in range(0, len(cfsr.lat)):
    kk = stats.linregress(t, U[:,i])
    if kk.pvalue<0.05:
        Ut_cfsr[i] = kk.slope*3650
    else:
        Ut_cfsr[i] = np.nan

fig, ax = plt.subplots(1, figsize=(7, 10))
ax.plot(Ut_ncep, ncep['lat'], color='k', marker='*')
ax.plot(Ut_cfsr, cfsr['lat'],  color='r', marker='*')
ax.grid(color='grey', linestyle=':')
ax.set_yticks([-70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0])
ax.set_yticklabels(['70$^{\circ}$S', '65$^{\circ}$S', '60$^{\circ}$S', '55$^{\circ}$S', '50$^{\circ}$S', '45$^{\circ}$S', '40$^{\circ}$S', '35$^{\circ}$S', '30$^{\circ}$S', '25$^{\circ}$S', '20$^{\circ}$S', '15$^{\circ}$S', '10$^{\circ}$S', '5$^{\circ}$S', '0$^{\circ}$S'])
ax.tick_params(labelsize=15)
ax.set_xlabel('Tendencias de V [m/s por década]', size=15)
ax.text(-0.1, 1, '(b)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/vmediocuenca_trend_b.png', dpi=250, bbox_inches='tight')
plt.show()

U = ncep['curl'].squeeze().mean(dim='lon')
Ut_ncep = np.empty(len(ncep.lat))
for i in range(0, len(ncep.lat)):
    kk = stats.linregress(t, U[:,i])
    if kk.pvalue<0.05:
        Ut_ncep[i] = kk.slope*3650
    else:
        Ut_ncep[i] = np.nan
U = cfsr['curl'].squeeze().mean(dim='lon')
Ut_cfsr = np.empty(len(cfsr.lat))
for i in range(0, len(cfsr.lat)):
    kk = stats.linregress(t, U[:,i])
    if kk.pvalue<0.05:
        Ut_cfsr[i] = kk.slope*3650
    else:
        Ut_cfsr[i] = np.nan

fig, ax = plt.subplots(1, figsize=(7, 10))
ax.plot(Ut_ncep*1e7, ncep['lat'], color='k', marker='*')
ax.plot(Ut_cfsr*1e7, cfsr['lat'],  color='r', marker='*')
ax.grid(color='grey', linestyle=':')
ax.set_yticks([-70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0])
ax.set_yticklabels(['70$^{\circ}$S', '65$^{\circ}$S', '60$^{\circ}$S', '55$^{\circ}$S', '50$^{\circ}$S', '45$^{\circ}$S', '40$^{\circ}$S', '35$^{\circ}$S', '30$^{\circ}$S', '25$^{\circ}$S', '20$^{\circ}$S', '15$^{\circ}$S', '10$^{\circ}$S', '5$^{\circ}$S', '0$^{\circ}$S'])
ax.tick_params(labelsize=15)
ax.set_xlabel('Tendencias del rotor [10$^{-7}$ Pa/m por década]', size=15)
ax.text(-0.1, 1, '(c)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/cmediocuenca_trend_b.png', dpi=250, bbox_inches='tight')
plt.show()
