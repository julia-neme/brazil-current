import numpy as np
import xarray
from functions import get_coasts
from matplotlib import pyplot as plt
from scipy import signal
from scipy import stats

def latmax_clim(dat):
    A = dat['curl'].resample('1MS', dim='time', how='mean')
    B = np.empty([len(A.lat), len(A.time)])
    EC, WC, land = get_coasts(A.lat, A.lon); del land

    for i in range(0, len(A.lat)):
        B[i, :] = A[:, i, int(WC[i]):int(EC[i])].mean(dim='lon')

    C = B.argmax(axis=0)
    LM = np.empty(np.shape(C))
    for i in range(0, len(C)):
        LM[i] = A.lat[C[i]].item()
    return LM, A.time

inp = ['/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_1979-2015.nc',
       '/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_1979-2015.nc',
       '/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_1979-2015.nc']

lbl = ['NCEP vI', 'NCEP vII', 'ERA-Interim']
cls = ['k', 'r', 'g']

fig1 = plt.figure(1, figsize=(10,7))
fig2 = plt.figure(2, figsize=(10,7))
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
for i in range(0, 3):
    dat = xarray.open_dataset(inp[i]).squeeze()
    LM, t = latmax_clim(dat)
    ax1.plot(t, LM, label=lbl[i], linewidth=.5, color=cls[i])

    s1, s2 = signal.butter(6, 1/12, 'low')
    B = signal.filtfilt(s1, s2, LM)

    k = t.sel(time='1993-01-01')
    n = np.where(t==k)[0][0]
    tt = np.arange(0, len(t), 1)
    LF1 = stats.linregress(tt[:n], B[:n])
    LF2 = stats.linregress(tt[n:], B[n:])

    ax2.plot(t, B, label=lbl[i], color=cls[i])
    ax2.plot(t[:n], tt[:n]*LF1[0] + LF1[1], color=cls[i], linestyle='--')
    ax2.plot(t[n:], tt[n:]*LF2[0] + LF2[1], color=cls[i], linestyle='--')
    print(LF1[0]*120)
    print(LF2[0]*120)
    print(LF1[-1]*120)
    print(LF2[-1]*120)
    del LM, B, t, LF1, LF2

ax1.legend(loc='upper right')
ax1.grid(linestyle=':')
ax1.set_title('Latitud máximo del rotor promediado (val. mensuales)')
fig1.savefig('figuras/latmax.png', bbox_inches='tight')

ax2.legend(loc='upper right')
ax2.grid(linestyle=':')
ax2.set_title('Latitud máximo del rotor promediado (filtrado)')
fig2.savefig('figuras/latmax_filt.png', bbox_inches='tight')
plt.show()
