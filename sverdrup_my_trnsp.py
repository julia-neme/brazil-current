"""
Calcula el transporte meridional de Sverdrup

@author: Julia Neme
"""


import numpy as np
import xarray

inp = input('File to process with path: ')
dat = xarray.open_dataset(inp).squeeze(drop=True)
dat = dat.astype('float64')

My = np.empty(np.shape(dat.curl))
beta = 2*7.29e-5*np.cos(np.deg2rad(dat.lat))/6371000
for i in range(0, len(dat.lat)):
    My[:,i,:] = dat.curl[:,i,:]/(beta[i]*1027)

my = xarray.DataArray(My, coords=(dat.time, dat.lat, dat.lon),
       dims=['time', 'lat', 'lon'], name='My').astype('float32')
oup = input('Archivo del curl (output) c/path: ')
res = my.to_dataset(name='My')
res.to_netcdf(oup, mode='w')
