import numpy as np
import netCDF4 as nc
import os
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy.ma as ma

os.chdir('/home/bock/Documents/tesis/vientos/validacion')
data = nc.Dataset('scow.nc', 'r')
lon = data.variables['longitude'][:]
lat = data.variables['latitude'][:]
mps = list(data.variables.keys())
scow = np.empty([12, len(lat), len(lon)])
for i in range(0, len(mps)-2):
    scow[i,:,:] = data.variables[mps[i+2]][:,:]
data.close()

kk = scow
scow[scow==-9999]=np.nan
scow = ma.masked_where(np.isnan(scow),scow)
lon_ticks = [280, 300, 320, 340, 0, 20]
lat_ticks = [-60, -40, -20, 0]
sve = ['scow_curl_01.png', 'scow_curl_02.png', 'scow_curl_03.png',
       'scow_curl_04.png', 'scow_curl_05.png', 'scow_curl_06.png',
       'scow_curl_07.png', 'scow_curl_08.png', 'scow_curl_09.png',
       'scow_curl_10.png', 'scow_curl_11.png', 'scow_curl_12.png']
for i in range(0, 12):
    fig = plt.figure()
    m = Basemap(projection='merc', llcrnrlat=-70.0, urcrnrlat=0.0,
                llcrnrlon=-80.0, urcrnrlon=30.0, resolution='i')
    m.fillcontinents()
    m.drawparallels(lat_ticks, labels=[1,0,0,0], linewidth=0)
    m.drawmeridians(lon_ticks, labels=[0,0,0,1], linewidth=0)
    x, y = m(*np.meshgrid(lon, lat))
    g = m.pcolormesh(x, y, scow[i,:,:], cmap='jet', vmin=-3, vmax=3)
    cb = fig.colorbar(g, fraction=0.046, pad=0.04, orientation='vertical')
    cb.set_label('10$^{-7}$ Pa m$^{-1}$')
    plt.title('scow')
    plt.tight_layout()
    plt.savefig(sve[i])
    plt.close()
