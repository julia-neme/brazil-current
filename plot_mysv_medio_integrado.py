import numpy as np
import xarray
from functions import get_coasts
from functions import zonal_integration
from matplotlib import pyplot as plt

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_2009-2015.nc').squeeze()
NC1 = dat['curl'].mean(dim='time')
NC1_lat = dat.lat; NC1_lon = dat.lon; del dat
dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_2009-2015.nc').squeeze()
NC2 = dat['curl'].mean(dim='time')
NC2_lat = dat.lat; NC2_lon = dat.lon; del dat
dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ccmp/CCMP_curl_daily_2009-2015.nc').squeeze()
CCM = dat['curl'].mean(dim='time')
CCM_lat = dat.lat; CCM_lon = dat.lon; del dat
dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/cfsr/CFSR_curl_daily_2009-2015.nc').squeeze()
CFS = dat['curl'].mean(dim='time')
CFS_lat = dat.lat; CFS_lon = dat.lon; del dat
dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_2009-2015.nc').squeeze()
ERA = dat['curl'].mean(dim='time')
ERA_lat = dat.lat; ERA_lon = dat.lon; del dat

NC1_int = zonal_integration(NC1, NC1_lat, NC1_lon)/(1027*2*7.29e-5*np.cos(np.deg2rad(NC1_lat))/6371000)
NC2_int = zonal_integration(NC2, NC2_lat, NC2_lon)/(1027*2*7.29e-5*np.cos(np.deg2rad(NC2_lat))/6371000)
CCM_int = zonal_integration(CCM, CCM_lat, CCM_lon)/(1027*2*7.29e-5*np.cos(np.deg2rad(CCM_lat))/6371000)
CFS_int = zonal_integration(CFS, CFS_lat, CFS_lon)/(1027*2*7.29e-5*np.cos(np.deg2rad(CFS_lat))/6371000)
ERA_int = zonal_integration(ERA, ERA_lat, ERA_lon)/(1027*2*7.29e-5*np.cos(np.deg2rad(ERA_lat))/6371000)

plt.figure(figsize=(7,10))
plt.plot(NC1_int[1:]/1e6, NC1_lat[1:], label='NCEP vI')
plt.plot(NC2_int[1:]/1e6, NC2_lat[1:], label='NCEP vII')
plt.plot(CCM_int[1:]/1e6, CCM_lat[1:], label='CCMP')
plt.plot(CFS_int[1:]/1e6, CFS_lat[1:], label='CFSR')
plt.plot(ERA_int[1:]/1e6, ERA_lat[1:], label='ERA-Interim')
plt.ticklabel_format(useOffset=False)
plt.legend()
plt.ylim(-37, -32)
plt.xlim(18, 43)
plt.grid(linestyle=':')
plt.xlabel('Sv')
plt.title('Transporte meridional medio de Sverdrup integrado en la cuenca')
plt.savefig('figuras/mysv_medio_integrado_zoom.png')
plt.close()
