import numpy as np
import xarray
from matplotlib import pyplot as plt
from scipy import signal

def series_temp(dat, latitud, longitud):
    A = dat['curl'].sel(lat=-30, lon=-45, method='nearest')
    s1, s2 = signal.butter(6, 1/30, 'low')
    B = signal.filtfilt(s1, s2, A)
    t = dat.time
    return A, B, t

inp = ['/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/cfsr/CFSR_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ccmp/CCMP_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_2009-2015.nc']
lbl = ['NCEP vI', 'NCEP vII', 'CFSR', 'CCMP', 'ERA-Interim']

fig1 = plt.figure(1, figsize=(10,7))
fig2 = plt.figure(2, figsize=(10,7))
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
for i in range(0, 5):
    dat = xarray.open_dataset(inp[i]).squeeze()
    A, B, t = series_temp(dat, -30, -45)
    ax1.plot(t, A*1e7, label=lbl[i], linewidth=.5)
    ax2.plot(t, B*1e7, label=lbl[i])
    del A, B, t

ax1.legend(loc='upper right')
ax1.grid(linestyle=':')
ax1.set_ylabel('10$^{-7}$ Pa/m')
ax1.set_ylim(-25, 25)
ax1.set_title('Rotor diario - 30S 45W')
fig1.savefig('figuras/curl_3045.png', bbox_inches='tight')

ax2.legend(loc='upper right')
ax2.grid(linestyle=':')
ax2.set_ylabel('10$^{-7}$ Pa/m')
ax2.set_title('Rotor filtrado - 30S 45W')
fig2.savefig('figuras/curl_filt_3045.png', bbox_inches='tight')
plt.show()
