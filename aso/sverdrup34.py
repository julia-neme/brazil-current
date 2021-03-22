import numpy as np
import xarray
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')
curl_cfsr = cfsr['curl'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(0, -40))

curl_ncep = ncep['curl'].sel(lat=slice(0, -40)).mean(dim='time').squeeze()
curl_cfsr = cfsr['curl'].sel(lat=slice(0, -40)).mean(dim='time').squeeze()
curl_ncep1 = ncep['curl'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(0, -40)).mean(dim='time').squeeze()
curl_cfsr1 = cfsr['curl'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(0, -40)).mean(dim='time').squeeze()
curl_ncep2 = ncep['curl'].sel(time=slice('2009-01-01', '2015-12-31')).sel(lat=slice(0, -40)).mean(dim='time').squeeze()
curl_cfsr2 = cfsr['curl'].sel(time=slice('2009-01-01', '2015-12-31')).sel(lat=slice(0, -40)).mean(dim='time').squeeze()

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
def zonal_integration(x, lat, lon):
    import numpy as np
    from scipy import integrate

    EC, WC = get_coasts(lat, lon)

    x_int = np.empty([len(lat)])
    for i in range(1, len(lat)):
        n = len(lon[int(WC[i]):int(EC[i])])
        h = np.abs(lon[1]-lon[0])*60*1.852*1000*np.cos(lat[i]*np.pi/180)
        xx = np.empty(n); xx[0] = 0
        for j in range(0, n-1):
            xx[j+1] = xx[j] + h
        x_int[i] = integrate.simps(x[i, int(WC[i]):int(EC[i])], xx)

    return x_int

ncep_int = zonal_integration(curl_ncep.values, curl_ncep.lat, curl_ncep.lon)
cfsr_int = zonal_integration(curl_cfsr.values, curl_cfsr.lat, curl_cfsr.lon)
ncep_int1 = zonal_integration(curl_ncep1.values, curl_ncep.lat, curl_ncep.lon)
cfsr_int1 = zonal_integration(curl_cfsr1.values, curl_cfsr.lat, curl_cfsr.lon)
ncep_int2 = zonal_integration(curl_ncep2.values, curl_ncep.lat, curl_ncep.lon)
cfsr_int2 = zonal_integration(curl_cfsr2.values, curl_cfsr.lat, curl_cfsr.lon)

betancep = 2 * 7.29e-5 / 6.3781e6 * np.cos(np.deg2rad(curl_ncep.lat))
betacfsr = 2 * 7.29e-5 / 6.3781e6 * np.cos(np.deg2rad(curl_cfsr.lat))

sv_ncep = ncep_int/(1e6*betancep*1027)
sv_ncep1 = ncep_int1/(1e6*betancep*1027)
sv_ncep2 = ncep_int2/(1e6*betancep*1027)
sv_cfsr = cfsr_int/(1e6*betacfsr*1027)
sv_cfsr1 = cfsr_int1/(1e6*betacfsr*1027)
sv_cfsr2 = cfsr_int2/(1e6*betacfsr*1027)

f,ax = plt.subplots(1,figsize=(10,10))
ax.plot(sv_ncep, curl_ncep.lat, color='k', linestyle=':')
ax.plot(sv_cfsr, curl_cfsr.lat, color='r', linestyle=':')
ax.plot(sv_ncep1, curl_ncep.lat, color='k', linestyle='--')
ax.plot(sv_cfsr1, curl_cfsr.lat, color='r', linestyle='--')
ax.plot(sv_ncep2, curl_ncep.lat, color='k')
ax.plot(sv_cfsr2, curl_cfsr.lat, color='r')
ax.grid(which="both", color='grey', linestyle=':')
ax.set_yticks([-40, -35, -30, -25, -20, -15, -10, -5, 0])
ax.set_yticklabels(['40$^{\circ}$S', '35$^{\circ}$S', '30$^{\circ}$S', '25$^{\circ}$S',
                    '20$^{\circ}$S', '15$^{\circ}$S', '10$^{\circ}$S', '5$^{\circ}$S', '0$^{\circ}$S'])
ax.set_xlabel('My [Sv]')
ax.set_ylabel('Latitud')
plt.savefig('/home/bock/Documents/tesis/resultados_1/sv_int.png', bbox_inches='tight')
plt.show()

sv34_ncep = ncep_int[:,-3]
sv34_cfsr = cfsr_int[:,-12]

from scipy import signal

s1, s2 = signal.butter(6, 1/130, 'low')
sv34_ncep_f = signal.filtfilt(s1, s2, sv34_ncep)
sv34_cfsr_f = signal.filtfilt(s1, s2, sv34_cfsr)

f,ax = plt.subplots(1,figsize=(10,5))
ax.plot(curl_ncep['time'], sv34_ncep_f/1e6/betancep[2].values/1027, color='k')
ax.plot(curl_cfsr['time'], sv34_cfsr_f/1e6/betacfsr[9].values/1027, color='r')
ax.grid(which="both", color='grey', linestyle=':')
ax.set_ylabel('My [Sv/m]')
ax.grid(which="major", color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/resultados_1/sv_int_345_serie.png', bbox_inches='tight')
plt.show()
