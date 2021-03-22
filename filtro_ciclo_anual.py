import xarray
import numpy as np

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')

clim_v = osc['v'].groupby('time.month').mean('time')
anom_v = osc['v'].groupby('time.month') - clim_v

clim_u = osc['u'].groupby('time.month').mean('time')
anom_u = osc['u'].groupby('time.month') - clim_u

u_anom = xarray.DataArray(anom_u[:,0,:,:], coords=(osc.time, osc.latitude, osc.longitude),
       dims=['time', 'lat', 'lon'], name='u_anom').astype('float32')
v_anom = xarray.DataArray(anom_v[:,0,:,:], coords=(osc.time, osc.latitude, osc.longitude),
       dims=['time', 'lat', 'lon'], name='v_anom').astype('float32')
u_clim = xarray.DataArray(clim_u[:,0,:,:], coords=(clim_u.month, osc.latitude, osc.longitude),
       dims=['month', 'lat', 'lon'], name='u_clim').astype('float32')
v_clim = xarray.DataArray(clim_v[:,0,:,:], coords=(clim_v.month, osc.latitude, osc.longitude),
       dims=['month', 'lat', 'lon'], name='v_clim').astype('float32')

oup = input('Archivo del stress (output) c/path: ')
res = u_anom.to_dataset(name='u_anom')
res['v_anom'] = v_anom
res['u_clim'] = u_clim
res['v_clim'] = v_clim
res.to_netcdf(oup, mode='w')

gh = xarray.open_dataset('/home/bock/Documents/tesis/datos/3.nc')
sst = gh['analysed_sst']

clim_sst = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghsst_clim_atlsur_2009_2015.nc')
anom_sst = sst.groupby('time.month') - clim_sst['sst_clim']

sst_anom = xarray.DataArray(anom_sst, coords=(sst.time, sst.lat, sst.lon),
       dims=['time', 'lat', 'lon'], name='sst_anom').astype('float32')
#sst_clim = xarray.DataArray(clim_sst, coords=(clim_sst.month, sst.lat, sst.lon),
#       dims=['month', 'lat', 'lon'], name='sst_clim').astype('float32')

oup = input('Archivo del stress (output) c/path: ')
res = sst_anom.to_dataset(name='sst_anom')
#res['sst_clim'] = sst_clim
res.to_netcdf(oup, mode='w')
