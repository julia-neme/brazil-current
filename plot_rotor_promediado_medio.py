import numpy as np
import xarray
from functions import get_coasts
from matplotlib import pyplot as plt

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_2009-2015.nc')
NC1 = dat['curl'].mean(dim='time')
NC1_lat = dat.lat;
EC,WC = get_coasts(NC1_lat, dat.lon)
NC1_za = np.empty(len(NC1_lat))
for i in range(0, len(NC1_lat)):
    NC1_za[i] = np.mean(NC1[i, int(WC[i]):int(EC[i])])
del EC, WC, NC1, dat

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_2009-2015.nc')
NC2 = dat['curl'].mean(dim='time')
NC2_lat = dat.lat
EC,WC = get_coasts(NC2_lat, dat.lon)
NC2_za = np.empty(len(NC2_lat))
for i in range(0, len(NC2_lat)):
    NC2_za[i] = np.mean(NC2[i, int(WC[i]):int(EC[i])])
del EC, WC, NC2, dat

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ccmp/CCMP_curl_daily_2009-2015.nc')
NC2 = dat['curl'].mean(dim='time')
CCM_lat = dat.lat
EC,WC = get_coasts(CCM_lat, dat.lon)
CCM_za = np.empty(len(CCM_lat))
for i in range(0, len(CCM_lat)):
    CCM_za[i] = np.mean(NC2[i, int(WC[i]):int(EC[i])])
del EC, WC, NC2, dat

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/cfsr/CFSR_curl_daily_2009-2015.nc')
NC2 = dat['curl'].mean(dim='time')
CFS_lat = dat.lat
EC,WC = get_coasts(CFS_lat, dat.lon)
CFS_za = np.empty(len(CFS_lat))
for i in range(0, len(CFS_lat)):
    CFS_za[i] = np.mean(NC2[i, int(WC[i]):int(EC[i])])
del EC, WC, NC2, dat

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_2009-2015.nc')
NC2 = dat['curl'].mean(dim='time')
ERA_lat = dat.lat
EC,WC = get_coasts(ERA_lat, dat.lon)
ERA_za = np.empty(len(ERA_lat))
for i in range(0, len(ERA_lat)):
    ERA_za[i] = np.mean(NC2[i, int(WC[i]):int(EC[i])])
del EC, WC, NC2, dat

plt.figure(figsize=(7,10))
plt.plot(NC1_za*1e7, NC1_lat, label='NCEP vI')
plt.plot(NC2_za*1e7, NC2_lat, label='NCEP vII')
plt.plot(CCM_za*1e7, CCM_lat, label='CCMP')
plt.plot(CFS_za*1e7, CFS_lat, label='CFSR')
plt.plot(ERA_za*1e7, ERA_lat, label='ERA-Interim')
plt.legend()
plt.grid(linestyle=':')
plt.xlabel('10$^{-7}$ Pa/m')
plt.title('Rotor medio promediado en la cuenca')
plt.savefig('figuras/curl_promediado_cuenca.png')
plt.close()
