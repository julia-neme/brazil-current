import numpy as np
import xarray
from matplotlib import pyplot as plt
from scipy import signal
from DoFs_CIs import get_eddof
from DoFs_CIs import get_chi2conf

inp = ['/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ccmp/CCMP_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/cfsr/CFSR_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_2009-2015.nc']
lbl = ['NCEP vI', 'NCEP vII', 'CCMP', 'CFSR', 'ERA-Interim']

plt.figure(1)
for i in range(0, len(inp)):
    dat = xarray.open_dataset(inp[i])
    A = dat['curl'].sel(lat=-30, lon=-45, method='nearest')
    B = signal.detrend(A)
    s1, s2 = signal.butter(6, 1/30, 'low')
    x = signal.filtfilt(s1, s2, B)

    f, Pxx = signal.welch(x, fs=1, nperseg=600, noverlap=300)

    ddof = np.round((8/3)*len(x)/600)
    CI = get_chi2conf(Pxx, 0.05, ddof)

    plt.plot(f, f*Pxx, label=lbl[i])
    plt.fill_between(f, f*CI[0], f*CI[1], alpha=.5)
    del dat, A, B, s1, s2, x, f, Pxx, ddof, CI
plt.xscale( "log" )
plt.grid(linestyle=':')
plt.legend()
plt.xlabel('log(f)')
plt.ylabel('f*PSD')
plt.savefig('figuras/psd_welch_ci.png', bbox_inches='tight')
plt.show()
