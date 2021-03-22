import numpy as np
import xarray
from functions import get_coasts
from matplotlib import pyplot as plt

def latmax_clim(dat):
    A = dat['curl'].resample('1MS', dim='time', how='mean')
    B = np.empty([len(dat.lat), len(A.time)])
    EC, WC = get_coasts(dat.lat, dat.lon)
    for i in range(0, len(dat.lat)):
        B[i, :] = A[:, i, int(WC[i]):int(EC[i])].mean(dim='lon')
    C = B.argmax(axis=0)
    LM = np.empty(np.shape(C))
    for i in range(0, len(C)):
        LM[i] = A.lat[C[i]].item()
    return LM, A.time

inp = ['/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ccmp/CCMP_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/cfsr/CFSR_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/viesntos/era_int/ERA_curl_daily_2009-2015.nc']

lm = np.empty([84, 5])
for i in range(0, 5):
    dat = xarray.open_dataset(inp[i])
    lm[:,i], t = latmax_clim(dat)


plt.figure(figsize=(10,7))
plt.plot(t, lm[:,0], label='NCEP vI - 40.95$^{\circ}$S')
plt.plot(t, lm[:,1], label='NCEP vII - 40.95$^{\circ}$S')
plt.plot(t, lm[:,2], label='CCMP - 39.875$^{\circ}$S')
plt.plot(t, lm[:,3], label='CFSR - 40$^{\circ}$S')
plt.plot(t, lm[:,4], label='ERA-Interim - 39.75$^{\circ}$S')
plt.legend()
plt.grid(linestyle=':')
plt.xlabel('Mes')
plt.ylabel('Latitud')
plt.title('Latitud del máximo del rotor del viento promediado en la cuenca')
plt.savefig('figuras/curl_latmax_promediado.png', bbox_inches='tight')
plt.close()
