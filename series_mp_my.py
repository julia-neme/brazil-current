import xarray
import numpy as np
from matplotlib import pyplot as plt

data = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')

curl = data['curl'].groupby('time.year').mean(dim = 'time')

def get_coasts(lat, lon):
    import numpy as np
    from mpl_toolkits.basemap import Basemap

    lon1 = np.empty(np.shape(lon))
    for i in range(0, len(lon)):
        if lon[i] > 30:
            lon1[i] = -360+lon[i]
        else:
            lon1[i] = lon[i]

    llclon = np.min(lon1); llclat = np.min(lat)
    urclon = np.max(lon1); urclat = np.max(lat)
    m = Basemap(projection='merc', area_thresh=10000, llcrnrlat=llclat,
                urcrnrlat=urclat, llcrnrlon=llclon, urcrnrlon=urclon,
                resolution='i')

    land = np.empty([len(lat), len(lon1)])
    for i in range(0, len(lat)):
        for j in range(0, len(lon1)):
            x, y = m(lon1[j],lat[i])
            land[i,j] = m.is_land(x,y)

    EC = np.empty(len(lat))
    WC = np.empty(len(lat))
    ss = int(len(lon1)/2)
    kk = np.diff(land, axis=1)
    for i in range(0, len(lat)):
        if any(kk[i,:] == -1):
            WC[i] = int(np.where(kk[i, 1:ss] == -1)[0][0]) + 2
        else:
            WC[i] = 0

    for i in range(0, len(lat)):
        if any(kk[i,ss:] == 1):
            EC[i] = int(np.where(kk[i, ss:] == 1)[0][0]) + 1 + ss
        else:
            EC[i] = len(lon1)
    return EC, WC
def zonal_integration(x, lat, lon, EC, WC):
    import numpy as np
    from scipy import integrate

    x_int = np.empty([len(lat)])
    for i in range(1, len(lat)):
        n = len(lon[int(WC[i]):int(EC[i])])
        h = np.abs(lon[1]-lon[0])*60*1.852*1000*np.cos(lat[i]*np.pi/180)
        xx = np.empty(n); xx[0] = 0
        for j in range(0, n-1):
            xx[j+1] = xx[j] + h
        x_int[i] = integrate.simps(x[i, int(WC[i]):int(EC[i])], xx)

    return x_int

ec, wc = get_coasts(curl['lat'], curl['lon'])
ec[0] = ec[1]; wc[0] = wc[1]
my = np.empty([len(curl['year']), len(curl['lat'][0:70])])
beta = 2 * 7.29e-5 / 6.3781e6 * np.cos(np.deg2rad(curl['lat'][0:70]))
for i in range(0, len(curl['year'])):
    my[i, :] = zonal_integration(curl[i, 0:70, :].values, curl['lat'][0:70], curl['lon'], ec, wc)
    my[i, :] = my[i, :] / (beta*1e6*1027)

sv = xarray.DataArray(my[:, 1:70], coords=(curl['year'], curl['lat'][1:70]), dims=['year', 'lat'], name='my_sv').astype('float32')
res = sv.to_dataset(name='my_sv')
res.to_netcdf('/home/bock/Documents/cfsr.nc', mode='w')

for i in range(0, len(curl['year'])):
    plt.plot(my[i, 1:70], curl['lat'][1:70])
plt.xlabel('My [Sv]')
plt.ylabel('Latitud')
plt.savefig('/home/bock/Documents/cfsr.png', bbox_inches = 'tight')
