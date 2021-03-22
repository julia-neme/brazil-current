import numpy as np
import xarray
from matplotlib import pyplot as plt
from scipy import signal

inp = ['/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_stress_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_stress_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ccmp/CCMP_stress_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/cfsr/CFSR_stress_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/era_int/ERA_stress_daily_2009-2015.nc']
lbl = ['NCEP vI', 'NCEP vII', 'CCMP', 'CFSR', 'ERA-Interim']

def ekman_trn(dat):
    tx = dat['taux'].sel(lat=-30, lon=-45, method='nearest')
    ty = dat['tauy'].sel(lat=-30, lon=-45, method='nearest')
    f = 2*7.29e-5*np.sin(np.deg2rad(tx.lat))
    Mx = ty/f
    My = -tx/f
    return dat.time, Mx, My

fig1 = plt.figure(1, figsize=(10,7))
fig2 = plt.figure(2, figsize=(10,7))
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
fig3 = plt.figure(3, figsize=(10,7))
fig4 = plt.figure(4, figsize=(10,7))
ax3 = fig3.add_subplot(111)
ax4 = fig4.add_subplot(111)
for i in range(0, 5):
    dat = xarray.open_dataset(inp[i]).squeeze()
    t, Mx, My = ekman_trn(dat)
    ax1.plot(t, Mx/1e6, label=lbl[i])
    ax2.plot(t, My/1e6, label=lbl[i])

    s1, s2 = signal.butter(6, 1/30, 'low')
    x = signal.filtfilt(s1, s2, Mx)
    xx = signal.filtfilt(s1, s2, My)
    ax3.plot(t, x/1e6, label=lbl[i])
    ax4.plot(t, xx/1e6, label=lbl[i])

    del t, Mx, My, x, xx
ax1.legend(loc='upper right')
ax1.grid(linestyle=':')
ax1.set_ylabel('Sv')
ax1.set_title('Mx Ekman - 30S 45W')
fig1.savefig('figuras/mxek_3045.png', bbox_inches='tight')

ax2.legend(loc='upper right')
ax2.grid(linestyle=':')
ax2.set_ylabel('Sv')
ax2.set_title('My Ekman - 30S 45W')
fig2.savefig('figuras/myek_3045.png', bbox_inches='tight')

ax3.legend(loc='upper right')
ax3.grid(linestyle=':')
ax3.set_ylabel('Sv')
ax3.set_title('Mx Ekman filtrado - 30S 45W')
fig3.savefig('figuras/mxek_filt_3045.png', bbox_inches='tight')

ax4.legend(loc='upper right')
ax4.grid(linestyle=':')
ax4.set_ylabel('Sv')
ax4.set_title('My Ekman filtrado - 30S 45W')
fig4.savefig('figuras/myek_filt_3045.png', bbox_inches='tight')
plt.show()
