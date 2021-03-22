import numpy as np
import xarray

inp = input('Archivo de vientos (input) c/path: ')
dat = xarray.open_dataset(inp, autoclose=True).squeeze()
dat = dat.astype('float64')

tx = np.empty(np.shape(dat.uwnd))
ty = np.empty(np.shape(dat.uwnd))

for i in range(0, len(dat.latitude)):
    for j in range(0, len(dat.longitude)):
        WS = np.sqrt(dat.uwnd[:,i,j]**2 + dat.vwnd[:,i,j]**2)
        Cd = (2.70/WS + 0.142 + WS/13.09)/1000
        tx[:,i,j] = 1.22 * Cd * dat.uwnd[:,i,j] * WS
        ty[:,i,j] = 1.22 * Cd * dat.vwnd[:,i,j] * WS
del Cd, WS

taux = xarray.DataArray(tx, coords=(dat.time, dat.latitude, dat.longitude),
       dims=['time', 'lat', 'lon'], name='taux').astype('float32')
tauy = xarray.DataArray(ty, coords=(dat.time, dat.latitude, dat.longitude),
       dims=['time', 'lat', 'lon'], name='tauy').astype('float32')

oup = input('Archivo del stress (output) c/path: ')
res = taux.to_dataset(name='taux')
res['tauy'] = tauy
res.to_netcdf(oup, mode='w')
