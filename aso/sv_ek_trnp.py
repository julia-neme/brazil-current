import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from mpl_toolkits.basemap import cm

def plotter(lat, lon, X, cmap, clevs, title, units):
    x, y = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(9,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-80, 30, -70, 0],
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
    box = sgeom.box(minx=-80, maxx=-40, miny=-56, maxy=-20)
    ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='black')
    pc = ax.contourf(x, y, X, clevs, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
    cbar.ax.set_ylabel(units)
    ax.set_title(title)
    return fig, ax

dat = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_2009_2015.nc')
C_nc = dat['curl'].mean(dim='time').squeeze()
Tx_nc = dat['taux'].mean(dim='time').squeeze()
Ty_nc = dat['tauy'].mean(dim='time').squeeze()
x, y = np.meshgrid(C_nc.lon, C_nc.lat)
fig, ax = plotter(C_nc.lat, C_nc.lon, C_nc.values*1e7, 'coolwarm', np.arange(-3, 3, 0.1), '', '10$^{-7}$Pa/m')
ax.quiver(x, y, Tx_nc.values, Ty_nc.values, transform=ccrs.PlateCarree())
plt.savefig('/home/bock/Documents/tesis/resultados/figs/crlncep_medio.png', bbox_inches='tight')

dat1 = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_2009_2015.nc')
C_cf = dat1['curl'].mean(dim='time').squeeze()
Tx_cf = dat1['taux'].mean(dim='time').squeeze()
Ty_cf = dat1['tauy'].mean(dim='time').squeeze()
x1, y1 = np.meshgrid(C_cf.lon, C_cf.lat)
fig, ax = plotter(C_cf.lat, C_cf.lon, C_cf.values*1e7, 'coolwarm', np.arange(-3, 3, 0.1), '', '10$^{-7}$Pa/m')
ax.quiver(x1[::3,::3], y1[::3,::3], Tx_cf[::3,::3].values, Ty_cf[::3,::3].values, transform=ccrs.PlateCarree())
plt.savefig('/home/bock/Documents/tesis/resultados/figs/crlcfsr_medio.png', bbox_inches='tight')

# Sverdrup ####################################################################

bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')
bati = bat.sel(y=slice(-56,-20), x=slice(-80,-40))
x_b, y_b = np.meshgrid(bati['x'], bati['y'])

def plotterchico(lat, lon, X, cmap, clevs, title, units):
    x, y = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(9,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-80, -40, -56, -20],
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
    b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.3, linestyles='-',
                   transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
    cbar.ax.set_ylabel(units)
    ax.set_title(title)
    return fig, ax

fnc = 2 * 7.29e-5 * np.sin(np.deg2rad(C_nc.lat))
betanc = 2 * 7.29e-5 / 6.3781e6 * np.cos(np.deg2rad(dat.lat))
svmync = C_nc / (betanc * 1027)

fig, ax = plotterchico(C_nc.lat, C_nc.lon, svmync.values, 'coolwarm', np.arange(-10, 10.1, 0.5), '', '1/s')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/svmyncep_medio.png', bbox_inches='tight')

fcf = 2 * 7.29e-5 * np.sin(np.deg2rad(C_cf.lat))
betacf = 2 * 7.29e-5 / 6.3781e6 * np.cos(np.deg2rad(dat1.lat))
svmycf = C_cf / (betacf * 1027)

fig, ax = plotterchico(C_cf.lat, C_cf.lon, svmycf.values, 'coolwarm', np.arange(-10, 10.1, 0.5), '', '1/s')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/svmycfsr_medio.png', bbox_inches='tight')

def get_coasts(lat, lon):
    import numpy as np
    from mpl_toolkits.basemap import Basemap

    lon1 = np.empty(np.shape(lon))
    for i in range(0, len(lon)):
        if lon[i] > 30:
            lon1[i] = -360+lon[i]
        else:
            lon1[i] = lon[i]

    llclon = np.min(lon1); llclat = np.min(lat)
    urclon = np.max(lon1); urclat = np.max(lat)
    m = Basemap(projection='merc', area_thresh=10000, llcrnrlat=llclat,
                urcrnrlat=urclat, llcrnrlon=llclon, urcrnrlon=urclon,
                resolution='i')

    land = np.empty([len(lat), len(lon1)])
    for i in range(0, len(lat)):
        for j in range(0, len(lon1)):
            x, y = m(lon1[j],lat[i])
            land[i,j] = m.is_land(x,y)

    EC = np.empty(len(lat))
    WC = np.empty(len(lat))
    ss = int(len(lon1)/2)
    kk = np.diff(land, axis=1)
    for i in range(0, len(lat)):
        if any(kk[i,:] == -1):
            WC[i] = int(np.where(kk[i, 1:ss] == -1)[0][0]) + 2
        else:
            WC[i] = 0

    for i in range(0, len(lat)):
        if any(kk[i,ss:] == 1):
            EC[i] = int(np.where(kk[i, ss:] == 1)[0][0]) + 1 + ss
        else:
            EC[i] = len(lon1)
    return EC, WC, land
def zonal_integration(x, lat, lon):
    import numpy as np
    from scipy import integrate

    EC, WC, land = get_coasts(lat, lon); del land
    if np.ndim(x) == 2:
        x_int = np.empty(len(lat))
        for i in range(1, len(lat)):
            n = len(lon[int(WC[i]):int(EC[i])])
            h = np.abs(lon[1]-lon[0])*60*1.852*1000*np.cos(lat[i]*np.pi/180)
            xx = np.empty(n); xx[0] = 0
            for j in range(0, n-1):
                xx[j+1] = xx[j] + h
            x_int[i] = integrate.simps(x[i, int(WC[i]):int(EC[i])], xx)
    else:
        x_int = np.empty([len(x[:,0,0]), len(lat)])
        for i in range(1, len(lat)):
            n = len(lon[int(WC[i]):int(EC[i])])
            h = np.abs(lon[1]-lon[0])*60*1.852*1000*np.cos(lat[i]*np.pi/180)
            xx = np.empty(n); xx[0] = 0
            for j in range(0, n-1):
                xx[j+1] = xx[j] + h
            x_int[:, i] = integrate.simps(x[:, i, int(WC[i]):int(EC[i])], xx)
    return x_int

svmync_int = zonal_integration(svmync, svmync.lat, svmync.lon)
svmycf_int = zonal_integration(svmycf, svmycf.lat, svmycf.lon)

f,ax = plt.subplots(1,figsize=(10,10))
ax.plot(svmync_int[0:22]/1e6, svmync.lat[0:22], color='lime', label='NCEP')
ax.plot(svmycf_int[0:86]/1e6, svmycf.lat[0:86], color='orange', label='CFSR')
ax.plot(svmycf_int[69]/1e6, svmycf.lat[69], color='k', marker='*')
ax.plot(svmync_int[18]/1e6, svmync.lat[18], color='k', marker='*')
ax.grid(which="both", color='grey', linestyle=':')
ax.set_xlabel('My [Sv/m]')
ax.legend()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/sv_int_345.png', bbox_inches='tight')

Svnc = dat['curl'].sel(time=osc['time'], method='nearest').squeeze()/(betanc*1027)
Svcf = dat1['curl'].sel(time=osc['time'], method='nearest').squeeze()/(betacf*1027)
Svnc_int = zonal_integration(Svnc[:,17:20,:].values, Svnc.lat[17:20], Svnc.lon)
Svcf_int = zonal_integration(Svcf[:,68:71,:].values, Svcf.lat[68:71], Svcf.lon)

from scipy import signal
s1, s2 = signal.butter(6, 1/5, 'low')
Svnc_intf = signal.filtfilt(s1, s2, Svnc_int[:,1])
Svcf_intf = signal.filtfilt(s1, s2, Svcf_int[:,1])

f,ax = plt.subplots(1,figsize=(10,10))
ax.plot(osc['time'], Svnc_int[:,1]/1e6, color='k', linewidth=.8, label='NCEP')
ax.plot(osc['time'], Svcf_int[:,1]/1e6, color='dimgray', linewidth=.8, label='CFSR')
ax.plot(osc['time'], Svnc_intf/1e6, color='lime', linewidth=1.3, label='NCEP')
ax.plot(osc['time'], Svcf_intf/1e6, color='red', linewidth=1.3, label='CFSR')
ax.grid(which="both", color='grey', linestyle=':')
ax.set_ylabel('My [Sv/m]')
ax.xaxis.set_minor_locator(mdates.MonthLocator())
ax.grid(which="both", color='grey', linestyle=':')
ax.legend()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/sv_int_345_serie.png', bbox_inches='tight')

Svnc = dat['curl'].squeeze()/(betanc*1027)
Svcf = dat1['curl'].squeeze()/(betacf*1027)
Svnc_int = zonal_integration(Svnc[:,:19,:].values, Svnc.lat[:19], Svnc.lon)
Svcf_int = zonal_integration(Svcf[:,:71,:].values, Svcf.lat[:71], Svcf.lon)

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

eddof_nc = np.empty(np.shape(Svnc_int[0,:]))
for i in range(0, len(Svnc_int[0,:])):
    eddof_nc[i] = get_eddof(Svnc_int[:,i])
eddof_cf = np.empty(np.shape(Svcf_int[0,:]))
for i in range(0, len(Svcf_int[0,:])):
    eddof_cf[i] = get_eddof(Svcf_int[:,i])

Svnc_sem = np.std(Svnc_int, axis=0)/np.sqrt(eddof_nc-1)
Svcf_sem = np.std(Svcf_int, axis=0)/np.sqrt(eddof_cf-1)

f,ax = plt.subplots(1,figsize=(10,10))
ax.plot(np.mean(Svnc_int, axis=0)/1e6, Svnc.lat[:19], color='lime', label='NCEP')
ax.plot(np.mean(Svcf_int, axis=0)/1e6, Svcf.lat[:71], color='orange', label='CFSR')
ax.fill_betweenx(Svnc.lat[:19], (np.mean(Svnc_int, axis=0)-Svnc_sem)/1e6, (np.mean(Svnc_int, axis=0)+Svnc_sem)/1e6, color='lime', alpha=.3)
ax.fill_betweenx(Svcf.lat[:71], (np.mean(Svcf_int, axis=0)-Svcf_sem)/1e6, (np.mean(Svcf_int, axis=0)+Svcf_sem)/1e6, color='orange', alpha=.3)
plt.axhline(-34.5, color='k', linestyle=':')
ax.grid(which="both", color='grey', linestyle=':')
ax.set_xlabel('My [Sv]')
ax.set_ylabel('Latitud')
ax.legend()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/sv_int.png', bbox_inches='tight')
plt.show()

# Ekman #######################################################################

ekmync = -Tx_nc / (fnc * 1027)
fig, ax = plotterchico(C_nc.lat, C_nc.lon, ekmync.values, 'coolwarm', np.arange(-1.5, 1.6, 0.1), '', '1/s')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/ekmyncep_medio.png', bbox_inches='tight')

ekmycf = -Tx_cf / (fcf * 1027)
fig, ax = plotterchico(C_cf.lat, C_cf.lon, ekmycf.values, 'coolwarm', np.arange(-1.5, 1.6, 0.1), '', '1/s')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/ekmycfsr_medio.png', bbox_inches='tight')

# En la transecta de SAMOC para comparar con OSCAR ############################

from scipy import signal
s1, s2 = signal.butter(6, 1/5, 'low')

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
dx = (1.34 * 60.0 * 1.8520 * np.cos(-34.33*(np.pi/180)) * 1000)
v = osc['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -45)).sel(longitude=slice(-51.66, -45)).mean(dim='longitude')[:,0]

A = dat['curl'].sel(time=osc['time'], method='nearest').sel(lat=slice(-33,-35), lon=slice(-52,-50))
v_svnc = A / (1027 * 2 * 7.29e-5 / 6.3781e6 * np.cos(np.deg2rad(-34.33)))
v_svnc = v_svnc[:,0,0]*dx

A = dat1['curl'].sel(time=osc['time'], method='nearest').sel(lat=slice(-33,-35), lon=slice(-52,-50)).mean(dim='lat').mean(dim='lon')
v_svcf = A * dx / (1027 * 2 * 7.29e-5 / 6.3781e6 * np.cos(np.deg2rad(-34.33)))

A = dat['taux'].sel(time=osc['time'], method='nearest').sel(lat=slice(-33,-35), lon=slice(-52,-50))
v_eknc = -A / (1027 * 2 * 7.29e-5  * np.sin(np.deg2rad(-34.33)))
v_eknc = v_eknc[:,0,0]*dx

A = dat1['taux'].sel(time=osc['time'], method='nearest').sel(lat=slice(-33,-35), lon=slice(-52,-50)).mean(dim='lat').mean(dim='lon')
v_ekcf = -A * dx / (1027 * 2 * 7.29e-5 * np.sin(np.deg2rad(-34.33)))

f,ax = plt.subplots(1,figsize=(15,10))
ax.plot(osc['time'], signal.filtfilt(s1,s2,v_svnc)/1e6, label='Sv. NCEP')
ax.plot(osc['time'], signal.filtfilt(s1,s2,v_svcf)/1e6, label='Sv. CFSR')
ax.plot(osc['time'], signal.filtfilt(s1,s2,v_eknc)/1e6, label='Ek. NCEP')
ax.plot(osc['time'], signal.filtfilt(s1,s2,v_ekcf)/1e6, label='Ek. CFSR')
ax.xaxis.set_minor_locator(mdates.MonthLocator())
ax.grid(which="both", color='grey', linestyle=':')
ax.set_ylabel('My [Sv/m]')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/my_samoc.png', bbox_inches='tight')

from scipy.stats import chi2

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

PSD_nc = welch_spectrum(signal.filtfilt(s1,s2,v_svnc), 1, 146, 0.05)
PSD_cf = welch_spectrum(signal.filtfilt(s1,s2,v_svcf), 1, 146, 0.05)

f, ax = plt.subplots(2, sharex=True, figsize=(15,10))
ax[0].plot(5/PSD_nc[0], PSD_nc[0]*PSD_nc[1], color='k', linewidth=1.5)
ax[0].fill_between(5/PSD_nc[0], PSD_nc[0]*PSD_nc[2], PSD_nc[0]*PSD_nc[3], color='k', alpha=.5)
ax[1].plot(5/PSD_cf[0], PSD_cf[0]*PSD_cf[1], color='k', linewidth=1.5)
ax[1].fill_between(5/PSD_cf[0], PSD_cf[0]*PSD_cf[2], PSD_cf[0]*PSD_cf[3], color='k', alpha=.5)
ax[0].set_ylabel('f*PSD [m/s]')
ax[1].set_ylabel('f*PSD [$^{\circ}$ C]')
ax[1].set_xlabel('Período [días]')
f.subplots_adjust(hspace=0.05)
plt.xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
plt.gca().invert_xaxis()
ax[0].grid(which="major", color='grey', linestyle=':')
ax[1].grid(which="major", color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/my_samoc_welch.png', bbox_inches='tight')
plt.show()
