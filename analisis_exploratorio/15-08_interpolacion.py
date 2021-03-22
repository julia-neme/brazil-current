import numpy as np
import netCDF4 as nc
from scipy import interpolate

data = nc.Dataset('/home/bock/Documents/tesis/vientos/ncep_v1/datos/NCEP_curl_daily_1948-2015.nc', 'r')
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
kkk = data.variables['curl'][:,:,:]
tim = nc.num2date(data.variables['time'][:],
       units=data.variables['time'].units,
       calendar=data.variables['time'].calendar)

h=0
while tim[h].year != 2009:
    h = h+1
crl = kkk[h:,:,:]; del kkk
ff = data.variables['time'][h:]

data1 = nc.Dataset('/home/bock/Documents/tesis/vientos/cfsr/datos/CFSR_curl_daily_1979-2015.nc', 'r')
lat1 = data.variables['lat'][:]
lon1 = data.variables['lon'][:]
data1.close()

for i in range(0, len(lon)):
    if lon[i] > 30:
        lon[i] = -360+lon[i]
x, y = np.meshgrid(lon, lat)
x1, y1 = np.meshgrid(lon1, lat1)

zz = np.empty([len(crl[:,0,0]), len(x1[:,0]), len(x1[0,:])])
ww = np.transpose([x.flatten(), y.flatten()])
for i in range(0, len(crl[:,0,0])):
    kk = crl[i,:,:].flatten()
    zz[i,:,:] = interpolate.griddata(ww, kk, (x1,y1), method='cubic')
    del kk

out = nc.Dataset('/home/bock/Documents/tesis/vientos/ncep_v1/datos/NCEP_curl_interp_2009-2015.nc', 'w')

out.createDimension('lon', len(lon1))
out.createDimension('lat', len(lat1))
out.createDimension('time', len(crl[:,0,0]))

longitud = out.createVariable('lon', 'f4', 'lon')
latitud = out.createVariable('lat', 'f4', 'lat')
cc = out.createVariable('curl', 'f4', ('time', 'lat', 'lon'))
tt = out.createVariable('time', 'f4', 'time')

longitud[:] = lon1
latitud[:] = lat1
cc[:,:,:] = zz
tt[:] = ff

tt.units = data.variables['time'].units
tt.calendar = data.variables['time'].calendar

data.close()
out.close()
