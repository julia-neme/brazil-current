import numpy as np
import xarray
from functions import climatologia_xarray
from functions import zonal_integration
from matplotlib import pyplot as plt
from functions import get_coasts
from scipy import integrate

def mysv_clim(dat):
    A = climatologia_xarray(dat['curl'])
    B = A.sel(lat=-40, method='nearest')

    EC, WC = get_coasts(A.lat, A.lon)
    C = np.empty(len(A[:,0,0]))
    kk = np.where(A.lat==B.lat)[0][0]
    n = len(B.lon[int(WC[kk]):int(EC[kk])])
    h = np.abs(B.lon[1]-B.lon[0])*60*1.852*1000*np.cos(B.lat*np.pi/180)
    xx = np.empty(n); xx[0] = 0
    for j in range(0, n-1):
        xx[j+1] = xx[j] + h
    C = integrate.simps(B[:, int(WC[kk]):int(EC[kk])], xx)
    D = B.lat.item()
    CC = C/(1027*2*7.29e-5*np.cos(np.deg2rad(D))/6371000)
    return CC, D

inp = ['/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/ccmp/CCMP_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/cfsr/CFSR_curl_daily_2009-2015.nc',
       '/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_2009-2015.nc']

mysv = np.empty([12, 5]); lat = np.empty(5)
for i in range(0, len(inp)):
    d = xarray.open_dataset(inp[i]).squeeze()
    CC, D = mysv_clim(d)
    mysv[:,i] = CC
    lat[i] = D
    del CC, D

t = np.arange(1,13,1)
plt.figure(figsize=(10,7))
plt.plot(t, mysv[:,0]/1e6, label='NCEP vI - 40.95$^{\circ}$S')
plt.plot(t, mysv[:,1]/1e6, label='NCEP vII - 40.95$^{\circ}$S')
plt.plot(t, mysv[:,2]/1e6, label='CCMP - 39.875$^{\circ}$S')
plt.plot(t, mysv[:,3]/1e6, label='CFSR - 40$^{\circ}$S')
plt.plot(t, mysv[:,4]/1e6, label='ERA-Interim - 39.75$^{\circ}$S')
plt.legend()
plt.grid(linestyle=':')
plt.xlabel('Mes')
plt.ylabel('Sv')
plt.title('Climatologia del transporte meridional de Sv. integrado')
plt.savefig('figuras/curl_integrado_40s.png', bbox_inches='tight')
plt.close()
