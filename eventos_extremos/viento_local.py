import numpy as np
import xarray
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from scipy import stats
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import cm
from scipy import signal

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')
cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_curl_atlsur_2009_2015.nc')
ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_curl_atlsur_2009_2015.nc')

v = osc['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -45)).sel(longitude=slice(-51.66, -45)).mean(dim='longitude')[:,0]
ssta = gha['sst_anom'].sel(lon=slice(-51.67,-45), lat=-34.5).sel(time=osc['time'], method='nearest').mean(dim='lon')

crl_cfsr = cfsr['curl'].sel(time=osc['time'], method='nearest').sel(lat=slice(-33,-35), lon=slice(-52,-50))
crl_ncep = ncep['curl'].sel(time=osc['time'], method='nearest').sel(lat=slice(-33,-35), lon=slice(-52,-50))
crl_cfsr_m = crl_cfsr.mean(dim='lon').mean(dim='lat')
crl_ncep_m = crl_ncep.mean(dim='lon').mean(dim='lat')

s1, s2 = signal.butter(6, 1/5, 'low')
crl_ncep_f = signal.filtfilt(s1, s2, crl_ncep_m)
crl_cfsr_f = signal.filtfilt(s1, s2, crl_cfsr_m)

f, ax = plt.subplots(2, sharex=True, figsize=(15,10))
#ax[0].plot(osc['time'], crl_ncep[:,0,0]*1e7, color='gray', linewidth=1)
#for i in range(0, 5):
#    ax[1].plot(osc['time'], crl_cfsr[:,i,:]*1e7, color='gray', linewidth=1)
ax[1].plot(osc['time'], crl_cfsr_m*1e4, color='k', linewidth=1)
ax[0].plot(osc['time'], crl_ncep_m*1e4, color='k', linewidth=1)
#ax[1].plot(osc['time'], crl_cfsr_f*1e7, color='r', linewidth=1.5)
#ax[0].plot(osc['time'], crl_ncep_f*1e7, color='r', linewidth=1.5)
ax[0].set_ylabel('Rotor [10$^{-4}$ Pa/m]')
ax[1].set_ylabel('Rotor [10$^{-4}$ Pa/m]')
f.subplots_adjust(hspace=0.05)
#ax[0].set_ylim([-10,10])
#ax[1].set_ylim([-10,10])
ax[0].xaxis.set_minor_locator(mdates.MonthLocator())
ax[0].grid(which="major", color='grey', linestyle=':')
ax[1].xaxis.set_minor_locator(mdates.MonthLocator())
ax[1].grid(which="major", color='grey', linestyle=':')
#plt.savefig('/home/bock/Documents/tesis/resultados/figs/serie_rotor.png', bbox_inches='tight')
plt.show()

# Correlaciones ###############################################################

cfsr_msk = np.ma.array(crl_cfsr_f, mask=np.abs(crl_cfsr_f*1e7)>3)
ncep_msk = np.ma.array(crl_ncep_f, mask=np.abs(crl_ncep_f*1e7)>3)

LR_cfsr = stats.linregress(v.values, cfsr_msk*1e7)
LR_ncep = stats.linregress(v.values, ncep_msk*1e7)

gs = gridspec.GridSpec(4, 8)
fig = plt.figure(figsize=(15, 10))

ax1 = plt.subplot(gs[0:2, 0:6])
ax1.plot(osc['time'], v.values, color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.set_ylabel('V [m/s]', color='blue')
ax11 = ax1.twinx()
ax11.plot(osc['time'], crl_ncep_f*1e7, color='lime')
ax11.tick_params(axis='y', labelcolor='lime')
ax11.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='lime')
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
ax1.grid(which="major", color='grey', linestyle=':')

ax2 = plt.subplot(gs[2:4, 0:6], sharex=ax1)
ax2.plot(osc['time'], v.values, color='blue')
ax2.tick_params(axis='y', labelcolor='blue')
ax2.set_ylabel('V [m/s]', color='blue')
ax22 = ax2.twinx()
ax22.plot(osc['time'], crl_cfsr_f*1e7, color='lime')
ax22.tick_params(axis='y', labelcolor='lime')
ax22.set_ylabel('Rotor [10$^{-7}$ Pa/m]', color='lime')
ax2.xaxis.set_minor_locator(mdates.MonthLocator())
ax2.grid(which="major", color='grey', linestyle=':')

ax3 = plt.subplot(gs[0:2, 6:8])
ax3.scatter(v.values, crl_ncep_f*1e7, color='k', s=1)
ax3.plot(v.values, LR_ncep.slope*v.values+LR_ncep.intercept, color='red')
ax3.tick_params(axis='y', labelcolor='lime')
ax3.tick_params(axis='x', labelcolor='blue', labelrotation=60)
ax3.grid(which="major", color='grey', linestyle=':')

ax4 = plt.subplot(gs[2:4, 6:8], sharex=ax3)
ax4.scatter(v.values, crl_cfsr_f*1e7, color='k', s=1)
ax4.plot(v.values, LR_cfsr.slope*v.values+LR_cfsr.intercept, color='red')
ax4.tick_params(axis='y', labelcolor='lime')
ax4.tick_params(axis='x', labelcolor='blue', labelrotation=60)
ax4.grid(which="major", color='grey', linestyle=':')

plt.tight_layout()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/correlacion_rotor.png', bbox_inches='tight')
plt.show()

# Idem con u (por geostrofia) #################################################

cfsr_u = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_wind_atlsur_2009_2015.nc')
ncep_u = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_wind_atlsur_2009_2015.nc')

u_cfsr = cfsr_u['uwnd'].sel(time=osc['time'], method='nearest').sel(lat=slice(-33,-35), lon=slice(-52,-50))
u_ncep = ncep_u['uwnd'].sel(time=osc['time'], method='nearest').sel(lat=slice(-33,-35), lon=slice(-52,-50)).squeeze()
u_cfsr_m = u_cfsr.mean(dim='lon').mean(dim='lat')
u_ncep_m = u_ncep
u_ncep_f = signal.filtfilt(s1, s2, u_ncep_m)
u_cfsr_f = signal.filtfilt(s1, s2, u_cfsr_m)

f, ax = plt.subplots(2, sharex=True, figsize=(15,10))
for i in range(0, 5):
    ax[1].plot(osc['time'], u_cfsr[:,i,:], color='gray', linewidth=1)
ax[1].plot(osc['time'], u_cfsr_m, color='k', linewidth=1)
ax[0].plot(osc['time'], u_ncep_m, color='k', linewidth=1)
ax[1].plot(osc['time'], u_cfsr_f, color='r', linewidth=1.5)
ax[0].plot(osc['time'], u_ncep_f, color='r', linewidth=1.5)
ax[0].set_ylabel('U [m/s]')
ax[1].set_ylabel('U [m/s]')
f.subplots_adjust(hspace=0.05)
#ax[0].set_ylim([-10,10])
#ax[1].set_ylim([-10,10])
ax[0].xaxis.set_minor_locator(mdates.MonthLocator())
ax[0].grid(which="major", color='grey', linestyle=':')
ax[1].xaxis.set_minor_locator(mdates.MonthLocator())
ax[1].grid(which="major", color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/serie_u.png', bbox_inches='tight')
plt.show()

LR_cfsru = stats.linregress(v.values, u_cfsr_f)
LR_ncepu = stats.linregress(v.values, u_ncep_f)

gs = gridspec.GridSpec(4, 8)
fig = plt.figure(figsize=(15, 10))

ax1 = plt.subplot(gs[0:2, 0:6])
ax1.plot(osc['time'], v.values, color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.set_ylabel('V [m/s]', color='blue')
ax11 = ax1.twinx()
ax11.plot(osc['time'], u_ncep_f, color='lime')
ax11.tick_params(axis='y', labelcolor='lime')
ax11.set_ylabel('U [m/s]', color='lime')
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
ax1.grid(which="major", color='grey', linestyle=':')

ax2 = plt.subplot(gs[2:4, 0:6], sharex=ax1)
ax2.plot(osc['time'], v.values, color='blue')
ax2.tick_params(axis='y', labelcolor='blue')
ax2.set_ylabel('V [m/s]', color='blue')
ax22 = ax2.twinx()
ax22.plot(osc['time'], u_cfsr_f, color='lime')
ax22.tick_params(axis='y', labelcolor='lime')
ax22.set_ylabel('U [m/s]', color='lime')
ax2.xaxis.set_minor_locator(mdates.MonthLocator())
ax2.grid(which="major", color='grey', linestyle=':')

ax3 = plt.subplot(gs[0:2, 6:8])
ax3.scatter(v.values, u_ncep_f, color='k', s=1)
ax3.plot(v.values, LR_ncepu.slope*v.values+LR_ncepu.intercept, color='red')
ax3.tick_params(axis='y', labelcolor='lime')
ax3.tick_params(axis='x', labelcolor='blue', labelrotation=60)
ax3.grid(which="major", color='grey', linestyle=':')

ax4 = plt.subplot(gs[2:4, 6:8], sharex=ax3)
ax4.scatter(v.values, u_cfsr_f, color='k', s=1)
ax4.plot(v.values, LR_cfsru.slope*v.values+LR_cfsru.intercept, color='red')
ax4.tick_params(axis='y', labelcolor='lime')
ax4.tick_params(axis='x', labelcolor='blue', labelrotation=60)
ax4.grid(which="major", color='grey', linestyle=':')

plt.tight_layout()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/correlacion_u.png', bbox_inches='tight')
plt.show()

# Correlacion laggeada con u ##################################################

R_lags = [np.corrcoef(v.values, u_cfsr_f)[0][1]]
N = int(len(osc['time'])/3)
CF = [1.96/np.sqrt(N*3-3)]
for i in range(1, N):
    R_lags.append(np.corrcoef(v.values[0:-i], u_cfsr_f[i:])[0][1])
    CF.append(1.96/np.sqrt(N*3-i-3))
CFm = -np.ones(len(CF))*CF
R_lags1 = [np.corrcoef(v.values, u_ncep_f)[0][1]]
CF1 = [1.96/np.sqrt(N*3-3)]
for i in range(1, N):
    R_lags1.append(np.corrcoef(v.values[0:-i], u_ncep_f[i:])[0][1])
    CF1.append(1.96/np.sqrt(N*3-i-3))
CFm1 = -np.ones(len(CF1))*CF1

f, ax = plt.subplots(2, sharex=True, figsize=(15,10))
ax[0].bar(np.arange(0, len(R_lags),1), R_lags, width=.3, color='k')
ax[1].bar(np.arange(0, len(R_lags1),1), R_lags1, width=.3, color='k')
ax[1].plot(CF1, color='k', linestyle='--')
ax[1].plot(CFm1, color='k', linestyle='--')
ax[0].plot(CF, color='k', linestyle='--')
ax[0].plot(CFm, color='k', linestyle='--')
f.subplots_adjust(hspace=0.05)
ax[0].grid(which="major", color='grey', linestyle=':')
ax[1].grid(which="major", color='grey', linestyle=':')
ax[1].set_xlabel('Lags (1 lag son 5 días)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/autocorr_u_v.png', bbox_inches='tight')
plt.show()
