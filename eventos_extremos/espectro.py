import xarray
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from scipy import signal
from scipy import stats
import matplotlib.ticker as mticker
from scipy.stats import chi2

osca = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')

v = osca['v_anom'].sel(lat=-34.5, method='nearest').sel(lon=slice(-52, -48)).sel(lon=slice(-51.66, -50.33)).mean(dim='lon').squeeze()
s1, s2 = signal.butter(6, 1/5, 'low')
vf = signal.filtfilt(s1, s2, v)

ssta = gha['sst_anom'].sel(lat=-34.5, method='nearest').sel(lon=slice(-52, -48)).sel(lon=slice(-51.66, -50.33)).mean(dim='lon').squeeze()
s1, s2 = signal.butter(6, 1/30, 'low')
sstaf = signal.filtfilt(s1, s2, ssta)

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

PSD_v = welch_spectrum(v.values, 1, 146, 0.05)
PSD_s = welch_spectrum(ssta.values, 1, 730, 0.05)

f, ax = plt.subplots(1, sharex=True, figsize=(10,5))
ax.plot(5/PSD_v[0], PSD_v[0]*PSD_v[1]*100, marker='*', color='b', linewidth=1.5)
ax.fill_between(5/PSD_v[0], PSD_v[0]*PSD_v[2]*100, PSD_v[0]*PSD_v[3]*100, color='k', alpha=.4)
ax.set_ylabel('f*PSD [(cm$^{2}$/s$^{2}$]')
ax.set_xlabel('Período [días]')
plt.xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
plt.gca().invert_xaxis()
ax.grid(which="both", color='grey', linestyle=':')
#plt.savefig('/home/bock/Documents/tesis/resultados/figs/psd_avm.png', bbox_inches='tight')
plt.show()

f, ax = plt.subplots(1, sharex=True, figsize=(10,5))
ax.plot(1/PSD_s[0], PSD_s[0]*PSD_s[1]*100, color='r', linewidth=1.5)
ax.fill_between(1/PSD_s[0], PSD_s[0]*PSD_s[2]*100, PSD_s[0]*PSD_s[3]*100, color='k', alpha=.4)
ax.set_ylabel('f*PSD [($^{\circ}$ C$^{2}$]')
ax.set_xlabel('Período [días]')
plt.xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
plt.gca().invert_xaxis()
ax.grid(which="both", color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/psd_ssta.png', bbox_inches='tight')
plt.show()

f, ax = plt.subplots(2, sharex=True, figsize=(12,7))
ax[0].plot(5/PSD_v[0], PSD_v[0]*PSD_v[1]*100, color='b', linewidth=1.5)
ax[0].fill_between(5/PSD_v[0], PSD_v[0]*PSD_v[2]*100, PSD_v[0]*PSD_v[3]*100, color='k', alpha=.4)
ax[1].plot(1/PSD_s[0], PSD_s[0]*PSD_s[1], color='r', linewidth=1.5)
ax[1].fill_between(1/PSD_s[0], PSD_s[0]*PSD_s[2], PSD_s[0]*PSD_s[3], color='k', alpha=.4)
ax[0].grid(which="both", color='grey', linestyle=':')
ax[1].grid(which="both", color='grey', linestyle=':')
ax[0].set_ylabel('f*DEP [(cm$^{2}$/s$^{2}$]')
ax[1].set_ylabel('f*DEP [($^{\circ}$ C$^{2}$]')
ax[1].set_xlabel('Período [días]')
plt.gca().set_xscale('log')
plt.gca().xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().xaxis.set_tick_params(rotation=60, which='both')
plt.gca().invert_xaxis()
f.subplots_adjust(hspace=0.05)
ax[0].text(-0.1, 0.9, '(a)', transform=ax[0].transAxes, size=15)
ax[1].text(-0.1, 0.9, '(b)', transform=ax[1].transAxes, size=15)

plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/psd_avm_ssta.png', bbox_inches='tight')
plt.show()

# Porcentaje por banda ########################################################

vvar = integrate.simps(PSD_v[0][1:]*PSD_v[1][1:]*100, np.log(PSD_v[0][1:]))
