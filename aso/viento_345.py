import xarray
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')
u = osc['u'].isel(time=slice(501,-1)).sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()
v = osc['v'].isel(time=slice(501,-1)).sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()

# Construyo las series

d = np.rad2deg(np.arctan2(v.mean(dim='time').item(), u.mean(dim='time').item()))
n = np.sqrt(u.mean(dim='time').item()**2 + v.mean(dim='time').item()**2)
u_B = u.mean(dim='time').item()/n
v_B = v.mean(dim='time').item()/n
P = np.empty(len(u.time))
for i in range(0, len(u.time)):
    P[i] = np.dot([u_B, v_B], [u[i].item(), v[i].item()])

ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
curl_ncep = ncep['curl'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(-30, -40)).squeeze()
u_ncep = ncep['uwnd'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(-30, -40)).squeeze()
v_ncep = ncep['vwnd'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(-30, -40)).squeeze()
cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')
curl_cfsr = cfsr['curl'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(-30, -40)).squeeze()
u_cfsr = cfsr['uwnd'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(-30, -40)).squeeze()
v_cfsr = cfsr['vwnd'].sel(time=slice('1999-10-06', '2015-12-31')).sel(lat=slice(-30, -40)).squeeze()

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

ncep_ec, ncep_wc = get_coasts(curl_ncep['lat'].values, curl_ncep['lon'].values)
cfsr_ec, cfsr_wc = get_coasts(curl_cfsr['lat'].values, curl_cfsr['lon'].values)

ncep_c34 = curl_ncep.isel(lat=2).isel(lon=slice(int(ncep_wc[2]), int(ncep_ec[2]))).mean(dim='lon').sel(time=u['time'], method='nearest')
cfsr_c34 = curl_cfsr.isel(lat=9).isel(lon=slice(int(cfsr_wc[9]), int(cfsr_ec[9]))).mean(dim='lon').sel(time=u['time'], method='nearest')
ncep_u34 = u_ncep.isel(lat=2).isel(lon=slice(int(ncep_wc[2]), int(ncep_ec[2]))).mean(dim='lon').sel(time=u['time'], method='nearest')
cfsr_u34 = u_cfsr.isel(lat=9).isel(lon=slice(int(cfsr_wc[9]), int(cfsr_ec[9]))).mean(dim='lon').sel(time=u['time'], method='nearest')
ncep_v34 = v_ncep.isel(lat=2).isel(lon=slice(int(ncep_wc[2]), int(ncep_ec[2]))).mean(dim='lon').sel(time=u['time'], method='nearest')
cfsr_v34 = v_cfsr.isel(lat=9).isel(lon=slice(int(cfsr_wc[9]), int(cfsr_ec[9]))).mean(dim='lon').sel(time=u['time'], method='nearest')

from scipy import stats
from scipy import signal

s1, s2 = signal.butter(6, 1/26, 'low')
ncep_u34_f = signal.filtfilt(s1, s2, ncep_u34)
ncep_v34_f = signal.filtfilt(s1, s2, ncep_v34)
ncep_c34_f = signal.filtfilt(s1, s2, ncep_c34)
cfsr_u34_f = signal.filtfilt(s1, s2, cfsr_u34)
cfsr_v34_f = signal.filtfilt(s1, s2, cfsr_v34)
cfsr_c34_f = signal.filtfilt(s1, s2, cfsr_c34)
P_f = signal.filtfilt(s1, s2, P)

# Tendencias

t = np.arange(0, len(P), 1)
lg1 = stats.linregress(t, ncep_u34_f)
lg2 = stats.linregress(t, cfsr_u34_f)
fig, ax = plt.subplots(1, figsize=(12, 5))
ax.plot(u.time, ncep_u34_f, color='k', linewidth=1)
ax.plot(u.time, lg1.intercept + lg1.slope*t, color='k', linewidth=1.5)
ax.plot(u.time, cfsr_u34_f, color='r', linewidth=1)
ax.plot(u.time, lg2.intercept + lg2.slope*t, color='r', linewidth=1.5)
ax.grid(which="both", color='grey', linestyle=':')
ax.set_ylabel('U [m/s]')
ax.text(-0.1, 1, '(a)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/u34_trend.png', dpi=250, bbox_inches='tight')
plt.show()

lg1 = stats.linregress(t, ncep_v34_f)
lg2 = stats.linregress(t, cfsr_v34_f)
fig, ax = plt.subplots(1, figsize=(12, 5))
ax.plot(u.time, ncep_v34_f, color='k', linewidth=1)
ax.plot(u.time, lg1.intercept + lg1.slope*t, color='k', linewidth=1.5)
ax.plot(u.time, cfsr_v34_f, color='r', linewidth=1)
ax.plot(u.time, lg2.intercept + lg2.slope*t, color='r', linewidth=1.5)
ax.grid(which="both", color='grey', linestyle=':')
ax.set_ylabel('V [m/s]')
ax.text(-0.1, 1, '(b)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/v34_trend.png', dpi=250, bbox_inches='tight')
plt.show()

lg1 = stats.linregress(t, ncep_c34_f*1e7)
lg2 = stats.linregress(t, cfsr_c34_f*1e7)
fig, ax = plt.subplots(1, figsize=(12, 5))
ax.plot(u.time, ncep_c34_f*1e7, color='k', linewidth=1)
ax.plot(u.time, lg1.intercept + lg1.slope*t, color='k', linewidth=1.5)
ax.plot(u.time, cfsr_c34_f*1e7, color='r', linewidth=1)
ax.plot(u.time, lg2.intercept + lg2.slope*t, color='r', linewidth=1.5)
ax.grid(which="both", color='grey', linestyle=':')
ax.set_ylabel('Rotor [10$^{-7}$ Pa/m]')
ax.text(-0.1, 1, '(c)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/c34_trend.png', dpi=250, bbox_inches='tight')
plt.show()

lg1 = stats.linregress(t, -P_f*100)
fig, ax = plt.subplots(1, figsize=(12, 5))
ax.plot(u.time, -P_f*100, color='b', linewidth=1)
ax.plot(u.time, lg1.intercept + lg1.slope*t, color='b', linewidth=1.5)
ax.grid(which="both", color='grey', linestyle=':')
ax.set_ylabel('$V_{P}$ [cm/s]')
ax.text(-0.1, 1, '(d)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/cb34_trend.png', dpi=250, bbox_inches='tight')
plt.show()
