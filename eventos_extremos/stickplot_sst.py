import xarray
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from scipy import signal
from scipy import stats

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')

u = osc['u'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude')[:,0]
v = osc['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude')[:,0]

ssta = gha['sst_anom'].sel(lon=slice(-51.67,-50.33), lat=-34.5).mean(dim='lon')
ssta = ssta.sel(time=osc['time'], method='nearest')

# Chequeo si tiene una tendencia importante y me conviene restar la regresion #
# lineal o solo con sacar la media estamos ####################################

x = np.arange(0, len(osc['time']), 1)
LR_v = stats.linregress(x, v)

plt.figure(1, figsize=(12,6))
plt.plot(osc['time'], v.values, color='k')
plt.plot(osc['time'], np.mean(v.values)*np.ones(len(osc['time'])), color='darkgreen', linewidth=2.5, label='Media')
plt.plot(osc['time'], LR_v.slope*x+LR_v.intercept, color='darkred', linewidth=2.5, label='Lin. Regress')
plt.gca().xaxis.set_minor_locator(mdates.MonthLocator())
plt.grid(color='grey', linestyle=':')
plt.ylabel('m/s')
plt.title('Velocidad meridional en sección 34.33S')
plt.legend()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/vserie_tendylinreg.png', bbox_inches='tight')
plt.show()
plt.close()

# Tendencia es significativa, le saco la regresion a la serie #################

v_detrend = v.values - (LR_v.slope*x+LR_v.intercept)
v_std = v_detrend.std()

# Busco los eventos extremos ##################################################

idx_ED = np.where(v.values > (LR_v.slope*x+LR_v.intercept+v_std))
idx_EF = np.where(v.values > (LR_v.slope*x+LR_v.intercept-v_std))

y1 = LR_v.slope*x+LR_v.intercept+v_std
y2 = LR_v.slope*x+LR_v.intercept-v_std

plt.figure(2, figsize=(12,6))
plt.plot(osc['time'], v.values, color='k')
plt.plot(osc['time'], LR_v.slope*x+LR_v.intercept, color='darkred', linewidth=2.5)
plt.plot(osc['time'], y1, color='darkred', linewidth=2.5, linestyle='--')
plt.plot(osc['time'], y2, color='darkred', linewidth=2.5, linestyle='--')
plt.fill_between(osc['time'].values, v.values, y1, where=v.values >= y1, color='darkred', alpha=.5)
plt.fill_between(osc['time'].values, v.values, y2, where=v.values <= y2, color='darkred', alpha=.5)
plt.gca().xaxis.set_minor_locator(mdates.MonthLocator())
plt.grid(color='grey', linestyle=':')
plt.ylabel('m/s')
plt.title('Velocidad meridional en sección 34.33S')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/vserie_constd.png', bbox_inches='tight')
plt.show()
plt.close()

# Pruebo distintos stick plots ################################################

fig, axs = plt.subplots(figsize=(12,6))
ax2 = axs.twinx()
q = axs.quiver(osc['time'].values, [[0]*len(osc['time'])], u, v, scale=.25, scale_units='inches', units='y',
              width=0.0005, headwidth=0, headlength=0, headaxislength=0)
axs.quiverkey(q, 0.075, 0.925, 0.25, label='0.25 m/s', labelpos='N')
axs.axes.get_yaxis().set_visible(False)
ax2.plot(osc['time'], ssta, color='darkred', label='SSTA')
ax2.legend()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position('left')
ax2.set_ylabel('$^{\circ}$ C')
ax2.xaxis.set_minor_locator(mdates.MonthLocator())
ax2.grid(color='grey', linestyle=':', which='both')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/stickplt_sinfiltrar.png', bbox_inches='tight')
plt.show()
plt.close()

s1, s2 = signal.butter(6, 1/4, 'low')
u_f = signal.filtfilt(s1, s2, u.values)
v_f = signal.filtfilt(s1, s2, v.values)
ssta_f = signal.filtfilt(s1, s2, ssta.values)

fig, axs = plt.subplots(figsize=(12,6))
ax2 = axs.twinx()
q = axs.quiver(osc['time'].values, [[0]*len(osc['time'])], u_f, v_f, scale=.25, scale_units='inches', units='y',
              width=0.0005, headwidth=0, headlength=0, headaxislength=0)
axs.quiverkey(q, 0.075, 0.925, 0.25, label='0.25 m/s', labelpos='N')
axs.axes.get_yaxis().set_visible(False)
ax2.plot(osc['time'].values, ssta_f, color='darkred', label='SSTA')
ax2.legend()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position('left')
ax2.set_ylabel('$^{\circ}$ C')
ax2.grid(color='grey', linestyle=':')
ax2.xaxis.set_minor_locator(mdates.MonthLocator())
plt.savefig('/home/bock/Documents/tesis/resultados/figs/stickplt_lowpass_20d.png', bbox_inches='tight')
plt.show()
plt.close()

LR_ssta = stats.linregress(x, ssta_f)
y11 = LR_ssta.slope*x+LR_ssta.intercept + np.std(ssta_f)
y22 = LR_ssta.slope*x+LR_ssta.intercept - np.std(ssta_f)

fig, axs = plt.subplots(figsize=(12,6))
axs.plot(osc['time'], v.values, color='k', label='V')
axs.plot(osc['time'], LR_v.slope*x+LR_v.intercept, color='k')
axs.plot(osc['time'], y1, color='k', linestyle='--')
axs.plot(osc['time'], y2, color='k', linestyle='--')
axs.fill_between(osc['time'].values, v.values, y1, where=v.values >= y1, color='k', alpha=.5)
axs.fill_between(osc['time'].values, v.values, y2, where=v.values <= y2, color='k', alpha=.5)
axs.xaxis.set_minor_locator(mdates.MonthLocator())
axs.grid(color='grey', linestyle=':')
axs.set_ylabel('m/s')

ax2 = axs.twinx()
ax2.plot(osc['time'].values, ssta_f, color='darkred', label='SSTA')
ax2.plot(osc['time'], y11-np.std(ssta_f), color='darkred')
ax2.plot(osc['time'], y11, color='darkred', linestyle='--')
ax2.plot(osc['time'], y22, color='darkred', linestyle='--')
ax2.fill_between(osc['time'].values, ssta_f, y11, where=ssta_f >= y11, color='darkred', alpha=.5)
ax2.fill_between(osc['time'].values, ssta_f, y22, where=ssta_f <= y22, color='darkred', alpha=.5)
ax2.set_ylabel('$^{\circ}$ C')
ax2.xaxis.set_minor_locator(mdates.MonthLocator())
plt.legend()
plt.savefig('/home/bock/Documents/tesis/resultados/figs/ssta_lowpass_20d_v_ee.png', bbox_inches='tight')
plt.show()
plt.close()

fig, axs = plt.subplots(2,1, figsize=(12,9))

axs[0].plot(osc['time'], v.values, color='k')
axs[0].plot(osc['time'], LR_v.slope*x+LR_v.intercept, color='k', linewidth=2.5)
axs[0].plot(osc['time'], y1, color='k', linewidth=2.5, linestyle='--')
axs[0].plot(osc['time'], y2, color='k', linewidth=2.5, linestyle='--')
axs[0].fill_between(osc['time'].values, v.values, y1, where=v.values >= y1, color='k', alpha=.5)
axs[0].fill_between(osc['time'].values, v.values, y2, where=v.values <= y2, color='k', alpha=.5)
axs[0].xaxis.set_minor_locator(mdates.MonthLocator())
axs[0].grid(color='grey', linestyle=':')
axs[0].set_ylabel('m/s')
axs[0].set_title('Vel. meridional' + "\n" + 'Slope: ' + str(LR_v.slope*73) + '$^{\circ}$ C por año')

axs[1].plot(osc['time'], ssta_f, color='k')
axs[1].plot(osc['time'], y11-np.std(ssta_f), color='darkred', linewidth=2.5)
axs[1].plot(osc['time'], y11, color='darkred', linewidth=2.5, linestyle='--')
axs[1].plot(osc['time'], y22, color='darkred', linewidth=2.5, linestyle='--')
axs[1].fill_between(osc['time'].values, ssta_f, y11, where=ssta_f >= y11, color='darkred', alpha=.5)
axs[1].fill_between(osc['time'].values, ssta_f, y22, where=ssta_f <= y22, color='darkred', alpha=.5)
axs[1].xaxis.set_minor_locator(mdates.MonthLocator())
axs[1].grid(color='grey', linestyle=':')
axs[1].set_ylabel('$^{\circ}$ C')
axs[1].set_title('SSTA' + "\n" + 'Slope: ' + str(LR_ssta.slope*73) + 'm/s por año')

plt.savefig('/home/bock/Documents/tesis/resultados/figs/ssta_lowpass_20d_v_ee_subplt.png', bbox_inches='tight')
plt.show()
plt.close()

R_lags = [np.corrcoef(v, ssta_f)[0][1]]
N = int(len(osc['time'])/3)
CF = [1.96/np.sqrt(N*3-3)]
for i in range(1, N):
    R_lags.append(np.corrcoef(v[0:-i], ssta_f[i:])[0][1])
    CF.append(1.96/np.sqrt(N*3-i-3))
CFm = -np.ones(len(CF))*CF

plt.figure(3, figsize=(12,6))
plt.bar(np.arange(0, len(R_lags),1), R_lags, width=.3, color='k')
plt.plot(CF, color='k', linestyle='--')
plt.plot(CFm, color='k', linestyle='--')
plt.ylim([-.5,.5])
plt.xlabel('Lags (1unidad = 5 días)')
plt.grid('grey', linestyle=':')
plt.title('Correlación entre vel. meridional y SSTA')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/corr_ssta_v.png', bbox_inches='tight')
plt.show()
