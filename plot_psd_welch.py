import numpy as np
import xarray
from matplotlib import pyplot as plt
from matplotlib import mlab
from scipy import signal

def welch_puntos(dat, latitud, longitud):
    A = dat['curl'].sel(lat=latitud, lon=longitud, method='nearest')
    ps, f = mlab.psd(signal.detrend(A), NFFT=256, Fs=1)
    return ps, f

inp = ['/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ccmp/CCMP_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/cfsr/CFSR_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_2009-2015.nc']

la = [-30, -35, -40]
lo = [-45, -15, -10]

PSD = []; Fr = []
for t in range(0, len(inp)):
    for i in range(0, len(la)):
        for j in range(0, len(lo)):
            dat = xarray.open_dataset(inp[t])
            A, B = welch_puntos(dat, la[i], lo[j])
            PSD.append(A); Fr.append(B)
            del A, B
tit = ['PSD 30S 45W', 'PSD 30S 15W', 'PSD 30S 10E',
       'PSD 35S 45W', 'PSD 35S 15W', 'PSD 35S 10E',
       'PSD 40S 45W', 'PSD 40S 15W', 'PSD 40S 10E',]
sve = ['figuras/psd_3045b.png', 'figuras/psd_3015b.png', 'figuras/psd_3010b.png',
       'figuras/psd_3545b.png', 'figuras/psd_3515b.png', 'figuras/psd_3510b.png',
       'figuras/psd_4045b.png', 'figuras/psd_4015b.png', 'figuras/psd_4010b.png',]
for i in range(0, 9):
    plt.figure(i, figsize=(10,5))
    plt.plot(Fr[5*i], PSD[5*i], label='NCEP vI', linewidth=.8)
    plt.plot(Fr[5*i+1], PSD[5*i+1], label='NCEP vII', linewidth=.8)
    plt.plot(Fr[5*i+2], PSD[5*i+2], label='CCMP', linewidth=.8)
    plt.plot(Fr[5*i+3], PSD[5*i+3], label='CFSR', linewidth=.8)
    plt.plot(Fr[5*i+4], PSD[5*i+4], label='ERA-Interim', linewidth=.8)
    plt.legend()
    plt.grid(linestyle=':')
    plt.xlabel('f')
    plt.ylabel('PSD')
    plt.title(tit[i])
    plt.savefig(sve[i])
    plt.show()
    plt.close()
