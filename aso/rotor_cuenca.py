import numpy as np
import xarray
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import cm
from mpl_toolkits.basemap import cm

dat = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_2009_2015.nc')
clim_nc = dat['curl'].groupby('time.month').mean('time').mean(dim='lon').squeeze()

dat1 = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_2009_2015.nc')
clim_cf = dat1['curl'].groupby('time.month').mean('time').mean(dim='lon').squeeze()

x1, y1 = np.meshgrid(np.arange(0, 12, 1), clim_nc.lat.values)
x2, y2 = np.meshgrid(np.arange(0, 12, 1), clim_cf.lat.values)

f, ax = plt.subplots(1,2, sharey=True, figsize=(10,8))
cl = ax[0].contourf(x1, y1, np.transpose(clim_nc.values*1e7), np.arange(-1,1.1,.1), cmap=cm.GMT_no_green, extend='both')
ze = ax[1].contourf(x2, y2, np.transpose(clim_cf.values*1e7), np.arange(-1,1.1,.1), cmap=cm.GMT_no_green, extend='both')
cbar=f.colorbar(cl,  ax=ax[1], shrink=.7)
cbar.ax.set_ylabel('10$^{-7}$ Pa/m')
f.subplots_adjust(hspace=0.05)
ax[0].set_ylabel('Latitud')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/clim_rotor_medio.png', bbox_inches='tight')
plt.show()
