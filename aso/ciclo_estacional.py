import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import cm
from mpl_toolkits.axes_grid1 import AxesGrid

def get_eddof(x):
    import numpy as np
    from scipy import integrate
    x = x[~np.isnan(x)]
    if np.size(x) == 0:
        edof = np.nan
    else:
        N = np.nanmax(len(x))
        N1 = N-1
        x = x - np.nanmean(x)

        # Calcula la funcion de autocorrelacion para lags desde
        # -N1 a N1. Por lo tanto, la correlacion sin lag (la
        # varianza de la serie) estara en c[N1] = c[N-1].

        c = np.correlate(x, x, 'full')

        # Normaliza la funcion de autocorrelacion segun N-1-k donde
        # k es el numero de lags, positivo.

        lags = np.abs(np.arange(-N1+1, N1, 1))
        cn = c[1:-1]/(N-1-lags)
        Var = cn[N1-1]

        # Busca el primer zero-crossing

        n = 0
        while (cn[N1+n] > 0) and (n < N1):
            n = n+1

        # Calcula el tiempo integral y los EDoF

        T = integrate.simps(cn[N1-1-n:N1+n])/Var

        edof = N/T
        if (np.isnan(edof) == False) and (np.isinf(edof) == False):
            edof = int(edof)
        else:
            edof = np.nan
        return edof

dat1 = xarray.open_dataset('/home/bock/Documents/tesis/datos/reynolds_sst_1982_2015.nc')
sst = dat1['sst'].sel(lat=slice(-60,0), lon=slice(-80,30)).groupby('time.month').mean('time')
sst2 = dat1['sst'].sel(time=slice('2009-01-01', '2015-12-31')).sel(lat=slice(-60,0), lon=slice(-80,30)).groupby('time.month').mean('time')

dat2 = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')
cfsr = dat2['curl'].sel(time=slice('1982-01-01', '2015-12-31'), lat=slice(0,-60), lon=slice(-80,30)).groupby('time.month').mean('time')
cfsr2 = dat2['curl'].sel(time=slice('2009-01-01', '2015-12-31'), lat=slice(0,-60), lon=slice(-80,30)).groupby('time.month').mean('time')

dat3 = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
ncep = dat3['curl'].sel(time=slice('1982-01-01', '2015-12-31'), lat=slice(0,-60), lon=slice(-80,30)).groupby('time.month').mean('time')
ncep2 = dat3['curl'].sel(time=slice('2009-01-01', '2015-12-31'), lat=slice(0,-60), lon=slice(-80,30)).groupby('time.month').mean('time')

x1, y1 = np.meshgrid(sst.lon, sst.lat)
x2, y2 = np.meshgrid(cfsr.lon, cfsr.lat)
x3, y3 = np.meshgrid(ncep.lon, ncep.lat)

def plotter(V, x, y, cmap, clevs, units, save):
    fig, axes = plt.subplots(4, 3, sharex=True, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(15,10))
    n = 0
    txt = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']
    for j in range(0,4):
        for i in range(0,3):
            axes[j,i].set_extent([-80, 30, -60, 0], crs=ccrs.PlateCarree())
            tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
            axes[j,i].add_feature(tierra)
            axes[j,i].coastlines(resolution='50m')

            cl = axes[j,i].contourf(x, y, V[n,:,:], clevs, cmap=cmap, transform=ccrs.PlateCarree(), extend='both')

            axes[j,i].text(0.1, 0.86, txt[n], transform=axes[j,i].transAxes, size=15)
            n = n+1

            grid_style = dict(color='white', linestyle='--', linewidth=0)
            gl = axes[j,0].gridlines(draw_labels=True, **grid_style)
            gl.ylabels_right = gl.xlabels_top = gl.xlabels_bottom = False
            gl.yformatter = LATITUDE_FORMATTER
            gl1 = axes[-1,i].gridlines(draw_labels=True, **grid_style)
            gl1.xlabels_top = gl1.ylabels_right = gl1.ylabels_left = False
            gl1.xformatter = LONGITUDE_FORMATTER
    fig.subplots_adjust(wspace=.02, hspace=.02)
    cbar = fig.colorbar(cl, ax=axes.ravel().tolist(), shrink=0.6)
    cbar.ax.set_ylabel(units)
    plt.savefig(save, bbox_inches='tight')

plotter(sst, x1, y1, cm.GMT_no_green, np.arange(0,30,1), '$^{\circ}$ C', '/home/bock/Documents/tesis/documentos/template/figuras2/sst_clim.png')
plotter(cfsr*1e7, x2, y2, 'coolwarm', np.arange(-3,3.2,.2), '10$^{-7}$ Pa/m', '/home/bock/Documents/tesis/documentos/template/figuras2/cfsr_clim.png')
plotter(ncep*1e7, x3, y3, 'coolwarm', np.arange(-3,3.2,.2), '10$^{-7}$ Pa/m', '/home/bock/Documents/tesis/documentos/template/figuras2/ncep_clim.png')

sst_c1 = sst.sel(lat=slice(-40, -20))
sst_c2 = sst2.sel(lat=slice(-40, -20))
for i in range(0, 12):
    for j in range(0, len(sst_c2.lat)):
        for k in range(0, len(sst_c2.lon)):
            if sst_c2[i,j,k] == -np.inf:
                sst_c2[i,j,k] = np.nan
            if sst_c1[i,j,k] == -np.inf:
                sst_c1[i,j,k] = np.nan

sst_s1 = sst_c1.mean(dim='lon').mean(dim='lat')
sst_std1 = sst_c1.std(dim='lon').std(dim='lat')/get_eddof(sst_s1)
sst_s2 = sst_c2.mean(dim='lon').mean(dim='lat')
sst_std2 = sst_c2.std(dim='lon').std(dim='lat')/get_eddof(sst_s2)

mth = np.arange(0, 12, 1)
fig, ax = plt.subplots(1, figsize=(10,3))
ax.plot(mth, sst_s1, color='r', linewidth=1.5, label='1999-2015')
ax.plot(mth, sst_s2, color='k', linewidth=1.5, label='2009-2015')
ax.fill_between(mth, sst_s1-sst_std1, sst_s1+sst_std1, color='r', alpha=0.3)
ax.fill_between(mth, sst_s2-sst_std2, sst_s2+sst_std2, color='k', alpha=0.3)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11])
ax.set_xticklabels(['Ene', 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.', 'Dic.'])
ax.set_ylabel('SST [$^{\circ}$]')
plt.grid(color='grey', linestyle=':')
plt.legend()
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras2/sst_clim_serie.png', dpi=250, bbox_inches='tight')
plt.show()

cfsr_s1 = cfsr.sel(lat=-34.5, method='nearest').mean(dim='lon')*1e7
cfsr_std1 = cfsr.std(dim='lon').std(dim='lat')/get_eddof(cfsr_s1)*1e7
cfsr_s2 = cfsr2.sel(lat=-34.5, method='nearest').mean(dim='lon')*1e7
cfsr_std2 = cfsr2.std(dim='lon').std(dim='lat')/get_eddof(cfsr_s2)*1e7

fig, ax = plt.subplots(1, figsize=(10,3))
ax.plot(mth, cfsr_s1, color='lime', linewidth=1.5, label='1999-2015')
ax.plot(mth, cfsr_s2, color='k', linewidth=1.5, label='2009-2015')
ax.fill_between(mth, cfsr_s1-cfsr_std1, cfsr_s1+cfsr_std1, color='lime', alpha=0.3)
ax.fill_between(mth, cfsr_s2-cfsr_std2, cfsr_s2+cfsr_std2, color='k', alpha=0.3)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11])
ax.set_xticklabels(['Ene', 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.', 'Dic.'])
ax.set_ylabel(r'$\nabla$ x $\tau$ [10$^{-7} Pa/m]$')
plt.grid(color='grey', linestyle=':')
plt.legend()
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras2/cfsr_clim_serie.png', dpi=250, bbox_inches='tight')
plt.show()

ncep_s1 = ncep.sel(lat=-34.5, method='nearest').mean(dim='lon')*1e7
ncep_std1 = ncep.std(dim='lon').std(dim='lat')/get_eddof(ncep_s1)*1e7
ncep_s2 = ncep2.sel(lat=-34.5, method='nearest').mean(dim='lon')*1e7
ncep_std2 = ncep2.std(dim='lon').std(dim='lat')/get_eddof(ncep_s2)*1e7

fig, ax = plt.subplots(1, figsize=(10,3))
ax.plot(mth, ncep_s1, color='orange', linewidth=1.5, label='1999-2015')
ax.plot(mth, ncep_s2, color='k', linewidth=1.5, label='2009-2015')
ax.fill_between(mth, ncep_s1-ncep_std1, ncep_s1+ncep_std1, color='orange', alpha=0.3)
ax.fill_between(mth, ncep_s2-ncep_std2, ncep_s2+ncep_std2, color='k', alpha=0.3)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11])
ax.set_xticklabels(['Ene', 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.', 'Dic.'])
ax.set_ylabel(r'$\nabla$ x $\tau$ [10$^{-7} Pa/m]$')
plt.grid(color='grey', linestyle=':')
plt.legend()
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras2/ncep_clim_serie.png', dpi=250, bbox_inches='tight')
plt.show()
