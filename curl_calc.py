import numpy as np
import xarray

inp = input('Archivo de vientos (input) c/path: ')
dat = xarray.open_dataset(inp).squeeze(drop=True)
dat = dat.astype('float64')

x_lon, y_lat = np.meshgrid(dat.lon, dat.lat)
x = (x_lon-x_lon[0,0]) * 60.0 * 1.8520 * np.cos(y_lat*(np.pi/180)) * 1000
y = (y_lat-y_lat[0,0]) * 60.0 * 1.8520 * 1000

rotor = np.empty(np.shape(dat.taux))
for i in range(1, len(dat.lat)-1):
    for j in range(1, len(dat.lon)-1):
        A1 = (dat.tauy[:,i,j+1]-dat.tauy[:,i,j]) / (x[i,j+1]-x[i,j])
        A2 = (dat.tauy[:,i,j]-dat.tauy[:,i,j-1]) / (x[i,j]-x[i,j-1])
        B1 = (dat.taux[:,i+1,j]-dat.taux[:,i,j]) / (y[i+1,j]-y[i,j])
        B2 = (dat.taux[:,i,j]-dat.taux[:,i-1,j]) / (y[i,j]-y[i-1,j])
        rotor[:, i, j] = 0.5*(A1+A2) - 0.5*(B1+B2)

curl = xarray.DataArray(rotor, coords=(dat.time, dat.lat, dat.lon),
       dims=['time', 'lat', 'lon'], name='curl').astype('float32')

oup = input('Archivo del curl (output) c/path: ')
res = curl.to_dataset(name='curl')
res.to_netcdf(oup, mode='w')
