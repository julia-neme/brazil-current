import numpy as np
import xarray
from functions import get_coasts
from functions import climatologia_xarray
from matplotlib import pyplot as plt

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_2009-2015.nc').squeeze()
A = climatologia_xarray(dat['curl'])
B = np.empty([len(dat.lat), 12])
EC, WC = get_coasts(dat.lat, dat.lon)
for i in range(0, len(dat.lat)):
    B[i, :] = np.mean(A[:, i, int(WC[i]):int(EC[i])], axis=1)

t = np.arange(1,13,1)
x, y = np.meshgrid(t, dat.lat)

fig = plt.figure(figsize=(7,9))
ax = fig.add_subplot(111)
pc = plt.pcolormesh(x, y, B*1e7, cmap='jet', vmin=-1.5, vmax=1.5)
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('10$^{-7}$ Pa/m')
plt.xlabel('Mes')
plt.ylabel('Latitud')
plt.title('Rotor promediado en la cuenca ERA')
plt.savefig('figuras/curl_promediado_climatologia_eraint.png', bbox_inches='tight')
plt.close()
