import xarray
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from scipy import signal
from scipy import stats
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
u = osc['u_anom'].sel(lat=-34.5, method='nearest').sel(lon=slice(-52, -48)).sel(lon=slice(-51.66, -50.33)).mean(dim='lon').squeeze()
v = osc['v_anom'].sel(lat=-34.5, method='nearest').sel(lon=slice(-52, -48)).sel(lon=slice(-51.66, -50.33)).mean(dim='lon').squeeze()
osc1 = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
uu = osc1['u'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()
vv = osc1['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()

d = np.rad2deg(np.arctan2(vv.mean(dim='time').item(), uu.mean(dim='time').item()))
n = np.sqrt(uu.mean(dim='time').item()**2 + vv.mean(dim='time').item()**2)

u_B = uu.mean(dim='time').item()/n
v_B = vv.mean(dim='time').item()/n

P = np.empty(len(uu.time))
for i in range(0, len(uu.time)):
    P[i] = np.dot([u_B, v_B], [u[i].item(), v[i].item()])

s1, s2 = signal.butter(6, 1/5, 'low')
vf = signal.filtfilt(s1, s2, v.values)
Pf = signal.filtfilt(s1, s2, P)

plt.figure(figsize=(10,3))
plt.plot(osc['time'], Pf*100, color='k', linewidth=1.2)
plt.plot(osc['time'], np.mean(Pf)*np.ones(len(Pf))*100, color='k', linewidth=2)
plt.plot(osc['time'], (np.mean(Pf)+np.std(Pf))*np.ones(len(Pf))*100, color='k', linestyle='--', linewidth=1.7)
plt.plot(osc['time'], (np.mean(Pf)-np.std(Pf))*np.ones(len(Pf))*100, color='k', linestyle='--', linewidth=1.7)
plt.plot(osc['time'], vf*100, color='b', linewidth=1.2)
plt.plot(osc['time'], np.mean(vf)*np.ones(len(vf))*100, color='b', linewidth=2)
plt.plot(osc['time'], (np.mean(vf)+np.std(vf))*np.ones(len(vf))*100, color='b', linestyle='--', linewidth=1.7)
plt.plot(osc['time'], (np.mean(vf)-np.std(vf))*np.ones(len(vf))*100, color='b', linestyle='--', linewidth=1.7)
plt.ylabel('AVM [cm/s]', size=15)
plt.gca().xaxis.set_minor_locator(mdates.MonthLocator())
plt.grid(which="major", color='grey', linestyle=':')
plt.tick_params(labelsize=15)
#plt.savefig('/home/bock/Documents/tesis/resultados/figs/eestd.png', bbox_inches='tight')
plt.show()

LR = stats.linregress(vf*100, Pf*100)

gs = gridspec.GridSpec(1, 8)
fig = plt.figure(figsize=(13, 4))
ax1 = plt.subplot(gs[0, 0:6])
ax1.plot(osc['time'], vf*100, color='blue', label='Vel. mer.')
ax1.set_ylabel('V [cm/s]', size=15)
ax1.plot(osc['time'], -Pf*100, color='black', label='Vel. proy.')
#ax1.plot(osc['time'], (np.mean(vf)+np.std(vf))*np.ones(len(vf))*100, color='b', linestyle='--', linewidth=1.7)
#ax1.plot(osc['time'], (np.mean(vf)-np.std(vf))*np.ones(len(vf))*100, color='b', linestyle='--', linewidth=1.7)
#plt.plot(osc['time'], (np.mean(-Pf)+np.std(-Pf))*np.ones(len(-Pf))*100, color='lime', linestyle='--', linewidth=1.7)
#plt.plot(osc['time'], (np.mean(-Pf)-np.std(-Pf))*np.ones(len(-Pf))*100, color='lime', linestyle='--', linewidth=1.7)
ax1.xaxis.set_minor_locator(mdates.MonthLocator())
ax1.grid(which="major", color='grey', linestyle=':')
ax3 = plt.subplot(gs[0, 6:8])
ax3.scatter(vf*100, -Pf*100, color='k', s=1)
ax3.plot(vf*100, -LR.slope*100*vf-LR.intercept*100, color='red')
ax3.grid(which="major", color='grey', linestyle=':')
ax1.tick_params(labelsize=15)
ax3.tick_params(labelsize=15)
plt.tight_layout()
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/fig_3_0_b.png', bbox_inches='tight')
plt.show()

# Como afectan los eventos a la CB ############################################

osa = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')
u = osa['u'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze().isel(time=slice(501, -1))
v = osa['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze().isel(time=slice(501, -1))
vv = v.groupby('time.month') - v.groupby('time.month').mean('time')

d = np.rad2deg(np.arctan2(v.mean(dim='time').item(), u.mean(dim='time').item()))
n = np.sqrt(u.mean(dim='time').item()**2 + v.mean(dim='time').item()**2)

u_B = u.mean(dim='time').item()/n
v_B = v.mean(dim='time').item()/n

P = np.empty(len(u.time))
for i in range(0, len(u.time)):
    P[i] = np.dot([u_B, v_B], [u[i].item(), v[i].item()])

s1, s2 = signal.butter(6, 1/5, 'low')
vf = signal.filtfilt(s1, s2, vv.values)
Pf = signal.filtfilt(s1, s2, P)

idxED = np.array([1, 5, 18, 26, 52, 60, 230, 238, 464, 469, 558, 569, 683, 689,
                   759, 766, 992, 996, 1031, 1040, 1082, 1088, 1127, 1132])

idxEF = np.array([127, 132, 215, 220, 240, 247, 455, 459, 526, 533, 595, 601, 609,
                  613, 728, 733, 983, 990, 1043, 1057, 1091, 1098, 1134, 1141])

cc_ed1 = np.empty([12])
cc_ef1 = np.empty([12])
for i in range(0, 12):
    cc_ed1[i] = np.mean(Pf[idxED[2*i]:idxED[2*i+1]])#/np.mean(Pf)
    cc_ef1[i] = np.mean(Pf[idxEF[2*i]:idxEF[2*i+1]])#/np.mean(Pf)
