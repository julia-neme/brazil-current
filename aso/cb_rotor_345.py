import xarray
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')
u = osc['u'].isel(time=slice(501,-1)).sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()
v = osc['v'].isel(time=slice(501,-1)).sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()

# Construyo la serie y la climatologia de la velocidad proyectada

d = np.rad2deg(np.arctan2(v.mean(dim='time').item(), u.mean(dim='time').item()))
n = np.sqrt(u.mean(dim='time').item()**2 + v.mean(dim='time').item()**2)
u_B = u.mean(dim='time').item()/n
v_B = v.mean(dim='time').item()/n
P = np.empty(len(u.time))
for i in range(0, len(u.time)):
    P[i] = np.dot([u_B, v_B], [u[i].item(), v[i].item()])

ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
curl_ncep = ncep['curl'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(-30, -40))
cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')
curl_cfsr = cfsr['curl'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(-30, -40))

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
    return EC, WC

ncep_ec, ncep_wc = get_coasts(curl_ncep['lat'].values, curl_ncep['lon'].values)
cfsr_ec, cfsr_wc = get_coasts(curl_cfsr['lat'].values, curl_cfsr['lon'].values)

ncep_34 = curl_ncep.isel(lat=2).isel(lon=slice(int(ncep_wc[2]), int(ncep_ec[2]))).mean(dim='lon').sel(time=u['time'], method='nearest')
cfsr_34 = curl_cfsr.isel(lat=9).isel(lon=slice(int(cfsr_wc[9]), int(cfsr_ec[9]))).mean(dim='lon').sel(time=u['time'], method='nearest')

# Espectros

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

PSD_ncep = welch_spectrum(ncep_34.values*1e7, 1, 730, 0.05)
PSD_cfsr = welch_spectrum(cfsr_34.values*1e7, 1, 730, 0.05)
PSD_P = welch_spectrum(P, 1, 146, 0.05)

fig, ax = plt.subplots(1, figsize=(10,6))
ax.plot(1/PSD_ncep[0], PSD_ncep[0]*PSD_ncep[1], color='k', linewidth=1.5)
ax.fill_between(1/PSD_ncep[0], PSD_ncep[0]*PSD_ncep[2], PSD_ncep[0]*PSD_ncep[3], color='k', linestyle='--', linewidth=1.5, alpha=.3)
ax.set_ylabel('f*PSD [(10$^{-7}$ Pa/m)^{2}]')
ax.set_xlabel('Período [días]')
plt.xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
#plt.gca().invert_xaxis()
ax.grid(which="both", color='grey', linestyle=':')
ax.text(-0.1, 1, '(b)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/ncep34_psd.png', dpi=250, bbox_inches='tight')
plt.show()

# Correlaciones con distintos filtros

from scipy import stats
from scipy import signal

kk1 = signal.detrend(ncep_34.values*1e7, type='linear')
kk2 = signal.detrend(cfsr_34.values*1e7, type='linear')
kk3 = signal.detrend(-P*100, type='linear')
reg_ncep = stats.linregress(kk1, kk3)
reg_cfsr = stats.linregress(kk2, kk3)

# figura
gs = gridspec.GridSpec(4, 8)
fig = plt.figure(figsize=(16, 9))

ax1 = plt.subplot(gs[0:2, 0:6])
ax1.plot(u['time'], kk3, color='k')
ax1.set_ylim([-30,30])
ax1.set_yticks([-30, -20, -10, 0, 10, 20, 30])
ax1.tick_params(axis='y', labelcolor='k')
ax1.set_ylabel('Vel. proy. [cm/s]', color='k')
ax11 = ax1.twinx()
ax11.plot(u['time'], kk1, color='r')
ax11.set_ylim([-3,3])
ax11.set_yticks([-3, -2, -1, 0, 1, 2, 3])
ax11.tick_params(axis='y', labelcolor='r')
ax11.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
ax1.grid(which="major", color='grey', linestyle=':')
ax2 = plt.subplot(gs[2:4, 0:6])
ax2.plot(u['time'], kk3, color='k')
ax2.set_ylim([-30,30])
ax2.set_yticks([-30, -20, -10, 0, 10, 20, 30])
ax2.tick_params(axis='y', labelcolor='k')
ax2.set_ylabel('Vel. proy. [cm/s]', color='k')
ax22 = ax2.twinx()
ax22.plot(u['time'], kk2, color='r')
ax22.set_ylim([-3,3])
ax22.set_yticks([-3, -2, -1, 0, 1, 2, 3])
ax22.tick_params(axis='y', labelcolor='r')
ax22.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
ax2.xaxis.set_minor_locator(mdates.MonthLocator())
ax2.grid(which="major", color='grey', linestyle=':')

ax3 = plt.subplot(gs[0:2, 6:8])
ax3.scatter(kk1, kk3, color='k', s=1)
ax3.plot(kk1, reg_ncep.slope*kk1+reg_ncep.intercept, color='b')
ax3.tick_params(axis='y', labelcolor='k')
ax3.tick_params(axis='x', labelcolor='r', labelrotation=60)
ax3.grid(which="major", color='grey', linestyle=':')
ax4 = plt.subplot(gs[2:4, 6:8])
ax4.scatter(kk2, kk3, color='k', s=1)
ax4.plot(kk2, reg_cfsr.slope*kk2+reg_cfsr.intercept, color='b')
ax4.tick_params(axis='y', labelcolor='k')
ax4.tick_params(axis='x', labelcolor='r', labelrotation=60)
ax4.grid(which="major", color='grey', linestyle=':')
plt.subplots_adjust(wspace=1.2)
ax1.text(0, 1.1, 'Corr. sin filtro (R_ncep = '+str(round(reg_ncep.rvalue, 2))+') (R_cfsr = '+str(round(reg_cfsr.rvalue, 2))+')',
         transform=ax1.transAxes)
plt.savefig('/home/bock/Documents/tesis/resultados_1/corr_sinfilt.png', dpi=250, bbox_inches='tight')
#plt.show()

fc = np.arange(6, 74, 2)
for i in range(0, len(fc)):
    s1, s2 = signal.butter(6, 1/fc[i], 'low')
    ncep_f = signal.filtfilt(s1, s2, kk1)
    cfsr_f = signal.filtfilt(s1, s2, kk2)
    cb_f = signal.filtfilt(s1, s2, kk3)

    reg_ncep = stats.linregress(ncep_f, cb_f)
    reg_cfsr = stats.linregress(cfsr_f, cb_f)
    gs = gridspec.GridSpec(4, 8)
    fig = plt.figure(figsize=(16, 9))

    ax1 = plt.subplot(gs[0:2, 0:6])
    ax1.plot(u['time'], cb_f, color='k')
    ax1.set_ylim([-30,30])
    ax1.set_yticks([-30, -20, -10, 0, 10, 20, 30])
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.set_ylabel('Vel. proy. [cm/s]', color='k')
    ax11 = ax1.twinx()
    ax11.plot(u['time'], ncep_f, color='r')
    ax11.set_ylim([-3,3])
    ax11.set_yticks([-3, -2, -1, 0, 1, 2, 3])
    ax11.tick_params(axis='y', labelcolor='r')
    ax11.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
    ax1.xaxis.set_minor_locator(mdates.MonthLocator())
    ax1.grid(which="major", color='grey', linestyle=':')
    ax2 = plt.subplot(gs[2:4, 0:6])
    ax2.plot(u['time'], cb_f, color='k')
    ax2.set_ylim([-30,30])
    ax2.set_yticks([-30, -20, -10, 0, 10, 20, 30])
    ax2.tick_params(axis='y', labelcolor='k')
    ax2.set_ylabel('Vel. proy. [cm/s]', color='k')
    ax22 = ax2.twinx()
    ax22.plot(u['time'], cfsr_f, color='r')
    ax22.set_ylim([-3,3])
    ax22.set_yticks([-3, -2, -1, 0, 1, 2, 3])
    ax22.tick_params(axis='y', labelcolor='r')
    ax22.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
    ax2.xaxis.set_minor_locator(mdates.MonthLocator())
    ax2.grid(which="major", color='grey', linestyle=':')

    ax3 = plt.subplot(gs[0:2, 6:8])
    ax3.scatter(ncep_f[::fc[i]], cb_f[::fc[i]], color='k', s=10)
    ax3.plot(ncep_f, reg_ncep.slope*ncep_f+reg_ncep.intercept, color='b')
    ax3.tick_params(axis='y', labelcolor='k')
    ax3.tick_params(axis='x', labelcolor='r', labelrotation=60)
    ax3.grid(which="major", color='grey', linestyle=':')
    ax4 = plt.subplot(gs[2:4, 6:8])
    ax4.scatter(cfsr_f[::fc[i]], cb_f[::fc[i]], color='k', s=10)
    ax4.plot(cfsr_f, reg_cfsr.slope*cfsr_f+reg_cfsr.intercept, color='b')
    ax4.tick_params(axis='y', labelcolor='k')
    ax4.tick_params(axis='x', labelcolor='r', labelrotation=60)
    ax4.grid(which="major", color='grey', linestyle=':')
    plt.subplots_adjust(wspace=1.2)
    ax1.text(0, 1.1, 'Corr. filtro'+str(5*fc[i])+'(R_ncep = '+str(round(reg_ncep.rvalue, 2))+') (R_cfsr = '+str(round(reg_cfsr.rvalue, 2))+')',
             transform=ax1.transAxes)
    plt.savefig('/home/bock/Documents/tesis/resultados_1/corr_filt'+str(5*fc[i])+'.png', dpi=250, bbox_inches='tight')
    plt.show()

# Mejor correlacion: la de f_130.

s1, s2 = signal.butter(6, 1/26, 'low')
ncep_f = signal.filtfilt(s1, s2, ncep_34.values*1e7)
cfsr_f = signal.filtfilt(s1, s2, cfsr_34.values*1e7)
cb_f = signal.filtfilt(s1, s2, -P*100)

PSD_ncepf = welch_spectrum(ncep_f*1e7, 1, 584, 0.05)
PSD_cfsrf = welch_spectrum(cfsr_f*1e7, 1, 584, 0.05)
PSD_Pf = welch_spectrum(cb_f, 1, 584, 0.05)

fig, ax = plt.subplots(1, figsize=(10,6))
ax.plot(5/PSD_Pf[0], PSD_Pf[0]*PSD_Pf[1], color='k', linewidth=1.5)
ax.fill_between(5/PSD_Pf[0], PSD_Pf[0]*PSD_Pf[2], PSD_Pf[0]*PSD_Pf[3], color='k', linestyle='--', linewidth=1.5, alpha=.3)
ax.set_ylabel('f*PSD [(10$^{-7}$ Pa/m)^{2}]')
ax.set_xlabel('Período [días]')
plt.xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
#plt.gca().invert_xaxis()
ax.grid(which="both", color='grey', linestyle=':')
ax.text(-0.1, 1, '(b)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/cb34f_psd.png', dpi=250, bbox_inches='tight')
plt.show()

# tendencias (ojo que cambie la def. de de las series filtradas y les saque la tendencia, x si lo quiero hacer de vuelta)
x = np.arange(0, len(ncep_34.values), 1)
ncep_t = stats.linregress(x, signal.filtfilt(s1, s2, ncep_34.values*1e7))
cfsr_t = stats.linregress(x, signal.filtfilt(s1, s2, cfsr_34.values*1e7))
cb_t = stats.linregress(x, signal.filtfilt(s1, s2, -P*100))

gs = gridspec.GridSpec(4, 1)
fig = plt.figure(figsize=(13, 9))
ax1 = plt.subplot(gs[0:2, :])
ax1.plot(u['time'], cb_f, color='k', linewidth=1)
ax1.plot(u['time'],cb_t.slope*x + cb_t.intercept, color='k', linewidth=2)
ax1.set_ylabel('Vel. proy. [cm/s]', color='k')
ax11 = ax1.twinx()
ax11.plot(u['time'], ncep_f, color='r', linewidth=1)
ax11.plot(u['time'], ncep_t.slope*x + ncep_t.intercept, color='r', linewidth=2)
ax11.set_ylim([-0.2,1.4])
ax11.set_yticks([0, 0.4, 0.8, 1.2])
ax11.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
ax11.tick_params(axis='y', labelcolor='r')
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
ax1.grid(which="major", color='grey', linestyle=':')
ax1.text(0, 1.1, 'Tendencia: '+str(round(ncep_t.slope*730, 2))+' 10$^{-7}$ Pa/m por década y '+str(round(cb_t.slope*730, 2))+' cm/s por década.',
         transform=ax1.transAxes)

ax2 = plt.subplot(gs[2:4, :])
ax2.plot(u['time'], cb_f, color='k', linewidth=1)
ax2.plot(u['time'],cb_t.slope*x + cb_t.intercept, color='k', linewidth=2)
ax2.set_ylabel('Vel. proy. [cm/s]', color='k')
ax22 = ax2.twinx()
ax22.plot(u['time'], cfsr_f, color='r', linewidth=1)
ax22.plot(u['time'], cfsr_t.slope*x + cfsr_t.intercept, color='r', linewidth=2)
ax22.set_ylim([-0.2,1.4])
ax22.set_yticks([0, 0.4, 0.8, 1.2])
ax22.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
ax22.tick_params(axis='y', labelcolor='r')
ax2.xaxis.set_minor_locator(mdates.MonthLocator())
ax2.grid(which="major", color='grey', linestyle=':')
ax2.text(0, 1.1, 'Tendencia: '+str(round(cfsr_t.slope*730, 2))+' 10$^{-7}$ Pa/m por década y '+str(round(cb_t.slope*730, 2))+' cm/s por década.',
         transform=ax2.transAxes)
plt.subplots_adjust(hspace=1.2)
plt.savefig('/home/bock/Documents/tesis/resultados_1/tendencias.png', dpi=250, bbox_inches='tight')
plt.show()

#tendencia para el periodo corto
x = np.arange(0, len(ncep_34.sel(time=slice('2009-01-01', '2015-12-31')).values), 1)
ncep_tc = stats.linregress(x, signal.filtfilt(s1, s2, ncep_34.sel(time=slice('2009-01-01', '2015-12-31')).values*1e7))
cfsr_tc = stats.linregress(x, signal.filtfilt(s1, s2, cfsr_34.sel(time=slice('2009-01-01', '2015-12-31')).values*1e7))
cb_tc = stats.linregress(x, signal.filtfilt(s1, s2, -P[-len(x):]*100))

# espectros
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

noise_ncep = np.random.normal(0, np.nanstd(ncep_f)*2., len(ncep_f))
noise_cfsr = np.random.normal(0, np.nanstd(cfsr_f)*2., len(cfsr_f))
noise_P = np.random.normal(0, np.nanstd(cb_f)*2., len(cb_f))

PSD_ncep = welch_spectrum(ncep_f, 1, 292, 0.05)
PSD_cfsr = welch_spectrum(cfsr_f, 1, 292, 0.05)
PSD_P = welch_spectrum(cb_f, 1, 292, 0.05)

PSD_ncep_n = welch_spectrum(ncep_f, 1, 292, 0.05)
PSD_cfsr_n = welch_spectrum(cfsr_f, 1, 292, 0.05)
PSD_P_n = welch_spectrum(cb_f, 1, 292, 0.05)
SL_ncep = np.ones(len(PSD_ncep_n[0]))*(np.sqrt(PSD_ncep_n[1]/2.)*2.*noise_ncep.std()).mean()
SL_cfsr = np.ones(len(PSD_cfsr_n[0]))*(np.sqrt(PSD_cfsr_n[1]/2.)*2.*noise_cfsr.std()).mean()
SL_P = np.ones(len(PSD_P_n[0]))*(np.sqrt(PSD_P_n[1]/2.)*2.*noise_P.std()).mean()

fig, ax = plt.subplots(1, figsize=(10,6))
ax.plot(5/PSD_P[0], PSD_P[0]*PSD_P[1], color='k', linewidth=1.5)
ax.fill_between(5/PSD_P[0], PSD_P[0]*PSD_P[2], PSD_P[0]*PSD_P[3],
                color='k', linestyle='--', linewidth=1.5, alpha=.3)
ax.plot(5/PSD_P_n[0], PSD_P_n[0]*SL_P, color='b', linewidth=1.5, linestyle='--')
ax.set_ylabel('f*PSD [(cm$^{2}$/s$^{2}$]')
ax.set_xlabel('Período [días]')
plt.xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
#plt.gca().invert_xaxis()
ax.grid(which="both", color='grey', linestyle=':')
#ax.text(-0.1, 1, '(b)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/cb_psd.png', dpi=250, bbox_inches='tight')
plt.show()

# Correlaciones laggeadas con la serie filtrada

Rlag_ncep = [np.corrcoef(signal.detrend(cb_f), signal.detrend(ncep_f))[0][1]]
N = int(len(u['time'])/3)
CF_ncep = [1.96/np.sqrt(N-3)]
for i in range(1, N):
    Rlag_ncep.append(np.corrcoef(signal.detrend(cb_f[0:-i]), signal.detrend(ncep_f[i:]))[0][1])
    CF_ncep.append(1.96/np.sqrt(N-i-3))

Rlag_cfsr = [np.corrcoef(signal.detrend(cb_f), signal.detrend(cfsr_f))[0][1]]
N = int(len(u['time'])/3)
CF_cfsr = [1.96/np.sqrt(N-3)]
for i in range(1, N):
    Rlag_cfsr.append(np.corrcoef(signal.detrend(cb_f[0:-i]), signal.detrend(cfsr_f[i:]))[0][1])
    CF_cfsr.append(1.96/np.sqrt(N-i-3))

Rlag_ncep = np.array(Rlag_ncep)
Rlag_cfsr = np.array(Rlag_cfsr)
CF_ncep = np.array(CF_ncep)
CF_cfsr = np.array(CF_cfsr)

f, ax = plt.subplots(2, sharex=True, figsize=(15,10))
ax[0].bar(np.arange(0, len(Rlag_ncep),1), Rlag_ncep, width=.3, color='k')
ax[1].bar(np.arange(0, len(Rlag_cfsr),1), Rlag_cfsr, width=.3, color='k')
ax[1].plot(CF_cfsr, color='k', linestyle='--')
ax[1].plot(-CF_cfsr, color='k', linestyle='--')
ax[0].plot(CF_ncep, color='k', linestyle='--')
ax[0].plot(-CF_ncep, color='k', linestyle='--')
f.subplots_adjust(hspace=0.05)
ax[0].grid(which="major", color='grey', linestyle=':')
ax[1].grid(which="major", color='grey', linestyle=':')
ax[1].set_xlabel('Lags [días]')
ax[0].set_ylim([-0.55, 0.55])
ax[1].set_ylim([-0.55, 0.55])
#ax[0].set_xlim([0, 73])
ax[0].set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
ax[0].set_xticklabels([0, 50, 100, 150, 200, 250, 300, 350])
ax[1].set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
ax[1].set_xticklabels([0, 50, 100, 150, 200, 250, 300, 350])
#plt.savefig('/home/bock/Documents/tesis/resultados_1/corr_lag.png', bbox_inches='tight')
plt.show()

gs = gridspec.GridSpec(4, 1)
fig = plt.figure(figsize=(13, 9))
ax1 = plt.subplot(gs[0:2, :])
ax1.plot(u['time'][:-3], cb_f[3:], color='k', linewidth=1)
ax1.set_ylabel('Vel. proy. [cm/s]', color='k')
ax11 = ax1.twinx()
ax11.plot(u['time'][:-3], ncep_f[:-3], color='r', linewidth=1)
ax11.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
ax11.tick_params(axis='y', labelcolor='r')
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
ax1.grid(which="major", color='grey', linestyle=':')

ax2 = plt.subplot(gs[2:4, :])
ax2.plot(u['time'], cb_f, color='k', linewidth=1)
ax2.set_ylabel('Vel. proy. [cm/s]', color='k')
ax22 = ax2.twinx()
ax22.plot(u['time'][:-3], cfsr_f[:-3], color='r', linewidth=1)
ax22.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='r')
ax22.tick_params(axis='y', labelcolor='r')
ax2.xaxis.set_minor_locator(mdates.MonthLocator())
ax2.grid(which="major", color='grey', linestyle=':')
plt.subplots_adjust(hspace=1.2)
plt.savefig('/home/bock/Documents/tesis/resultados_1/lag_15dias.png', dpi=250, bbox_inches='tight')
plt.show()
