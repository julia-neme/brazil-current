import numpy as np
import xarray
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from mpl_toolkits.basemap import cm

def plotter(lat, lon, X, cmap, clevs, units, txt):
    x, y = np.meshgrid(lon, lat)
    fig = plt.figure(figsize=(9,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-80, 30, -60, 0],
                  crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black',
             facecolor='white')
    ax.add_feature(tierra)
    ax.coastlines(resolution='50m')

    grid_style = dict(color='white', linestyle='--', linewidth=.3)
    gl = ax.gridlines(draw_labels=True, **grid_style)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    pc = ax.contourf(x, y, X, clevs, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax, shrink=0.7)
    cbar.ax.set_ylabel(units)
    ax.text(-0.1, 1, txt, transform=ax.transAxes, size=15)
    return fig, ax

ncep = xarray.open_dataset('/home/bock/Documents/tesis/datos/ncep2_atlsur_1979_2015.nc')
cfsr = xarray.open_dataset('/home/bock/Documents/tesis/datos/cfsr_atlsur_1979_2015.nc')

# U y V e intensidad

u_ncep = ncep['uwnd'].sel(time=slice('1999-10-06', '2015-12-31')).mean(dim='time').squeeze()
u_cfsr = cfsr['uwnd'].sel(time=slice('1999-10-06', '2015-12-31')).mean(dim='time').squeeze()
fig, ax = plotter(u_ncep.lat, u_ncep.lon, u_ncep.values, 'coolwarm', np.arange(-9.5, 10, .5), 'm/s', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/uncep_medio.png', bbox_inches='tight'); plt.show()
fig, ax = plotter(u_cfsr.lat, u_cfsr.lon, u_cfsr.values, 'coolwarm', np.arange(-9.5, 10, .5), 'm/s', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/ucfsr_medio.png', bbox_inches='tight'); plt.show()

v_ncep = ncep['vwnd'].sel(time=slice('1999-10-06', '2015-12-31')).mean(dim='time').squeeze()
v_cfsr = cfsr['vwnd'].sel(time=slice('1999-10-06', '2015-12-31')).mean(dim='time').squeeze()
fig, ax = plotter(v_ncep.lat, v_ncep.lon, v_ncep.values, 'coolwarm', np.arange(-4.5, 10, .5), 'm/s', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/vncep_medio.png', bbox_inches='tight'); plt.show()
fig, ax = plotter(v_cfsr.lat, v_cfsr.lon, v_cfsr.values, 'coolwarm', np.arange(-4.5, 10, .5), 'm/s', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/vcfsr_medio.png', bbox_inches='tight'); plt.show()

fig, ax = plotter(v_ncep.lat, v_ncep.lon, np.sqrt(u_ncep**2+v_ncep**2), 'viridis_r', np.arange(0, 10, .5), 'm/s', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/intncep_medio.png', bbox_inches='tight'); plt.show()
fig, ax = plotter(v_cfsr.lat, v_cfsr.lon, np.sqrt(u_cfsr**2+v_cfsr**2), 'viridis_r', np.arange(0, 10, .5), 'm/s', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/intcfsr_medio.png', bbox_inches='tight'); plt.show()

u_ncep1 = ncep['uwnd'].sel(time=slice('2009-01-01', '2015-12-31')).mean(dim='time').squeeze()
u_cfsr1 = cfsr['uwnd'].sel(time=slice('2009-01-01', '2015-12-31')).mean(dim='time').squeeze()
fig, ax = plotter(u_ncep.lat, u_ncep.lon, u_ncep.values, 'coolwarm', np.arange(-9.5, 10, .5), 'm/s', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/uncep_medio_c.png', bbox_inches='tight'); plt.show()
fig, ax = plotter(u_cfsr.lat, u_cfsr.lon, u_cfsr.values, 'coolwarm', np.arange(-9.5, 10, .5), 'm/s', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/ucfsr_medio_c.png', bbox_inches='tight'); plt.show()

v_ncep1 = ncep['vwnd'].sel(time=slice('2009-01-01', '2015-12-31')).mean(dim='time').squeeze()
v_cfsr1 = cfsr['vwnd'].sel(time=slice('2009-01-01', '2015-12-31')).mean(dim='time').squeeze()
fig, ax = plotter(v_ncep.lat, v_ncep.lon, v_ncep.values, 'coolwarm', np.arange(-4.5, 10, .5), 'm/s', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/vncep_medio_c.png', bbox_inches='tight'); plt.show()
fig, ax = plotter(v_cfsr.lat, v_cfsr.lon, v_cfsr.values, 'coolwarm', np.arange(-4.5, 10, .5), 'm/s', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/vcfsr_medio_c.png', bbox_inches='tight'); plt.show()

fig, ax = plotter(v_ncep.lat, v_ncep.lon, np.sqrt(u_ncep**2+v_ncep**2), 'viridis_r', np.arange(0, 10, .5), 'm/s', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/intncep_medio_c.png', bbox_inches='tight'); plt.show()
fig, ax = plotter(v_cfsr.lat, v_cfsr.lon, np.sqrt(u_cfsr**2+v_cfsr**2), 'viridis_r', np.arange(0, 10, .5), 'm/s', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/intcfsr_medio_c.png', bbox_inches='tight'); plt.show()

# Diferencias

fig, ax = plotter(u_ncep.lat, u_ncep.lon, u_ncep1-u_ncep, 'coolwarm', np.arange(-1, 1.1, .1), 'm/s', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/uncep_diff.png', bbox_inches='tight'); plt.show()
fig, ax = plotter(u_cfsr.lat, u_cfsr.lon, u_cfsr1-u_cfsr, 'coolwarm', np.arange(-1, 1.1, .1), 'm/s', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados_1/ucfsr_diff.png', bbox_inches='tight'); plt.show()

# Viento zonal promediado en la cuenca

u_ncep2 = ncep['uwnd'].mean(dim='time').squeeze()
u_cfsr2 = cfsr['uwnd'].mean(dim='time').squeeze()
v_ncep2 = ncep['vwnd'].mean(dim='time').squeeze()
v_cfsr2 = cfsr['vwnd'].mean(dim='time').squeeze()

f,ax = plt.subplots(1,figsize=(7,10))
ax.plot(u_ncep.mean(dim='lon'), ncep.lat, color='k', linestyle='--')
ax.plot(u_ncep1.mean(dim='lon'), ncep.lat, color='k')
ax.plot(u_ncep2.mean(dim='lon'), ncep.lat, color='k', linestyle=':')
ax.plot(u_cfsr.mean(dim='lon'), cfsr.lat, color='r', linestyle='--')
ax.plot(u_cfsr1.mean(dim='lon'), cfsr.lat, color='r')
ax.plot(u_cfsr2.mean(dim='lon'), cfsr.lat, color='r', linestyle=':')
ax.grid(which="both", color='grey', linestyle=':')
ax.set_yticks([-70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0])
ax.set_yticklabels(['70$^{\circ}$S', '65$^{\circ}$S', '60$^{\circ}$S', '55$^{\circ}$S', '50$^{\circ}$S', '45$^{\circ}$S',
                    '40$^{\circ}$S', '35$^{\circ}$S', '30$^{\circ}$S', '25$^{\circ}$S',
                    '20$^{\circ}$S', '15$^{\circ}$S', '10$^{\circ}$S', '5$^{\circ}$S', '0$^{\circ}$S'])
ax.set_xlabel('U [m/s]')
ax.set_ylabel('Latitud')
ax.text(-0.1, 1, '(a)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/u_medio_cuenca.png', bbox_inches='tight')
plt.show()

f,ax = plt.subplots(1,figsize=(7,10))
ax.plot(v_ncep.mean(dim='lon'), ncep.lat, color='k', linestyle='--')
ax.plot(v_ncep1.mean(dim='lon'), ncep.lat, color='k')
ax.plot(v_ncep2.mean(dim='lon'), ncep.lat, color='k', linestyle=':')
ax.plot(v_cfsr.mean(dim='lon'), cfsr.lat, color='r', linestyle='--')
ax.plot(v_cfsr1.mean(dim='lon'), cfsr.lat, color='r')
ax.plot(v_cfsr2.mean(dim='lon'), cfsr.lat, color='r', linestyle=':')
ax.grid(which="both", color='grey', linestyle=':')
ax.set_yticks([-70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0])
ax.set_yticklabels(['70$^{\circ}$S', '65$^{\circ}$S', '60$^{\circ}$S', '55$^{\circ}$S', '50$^{\circ}$S', '45$^{\circ}$S',
                    '40$^{\circ}$S', '35$^{\circ}$S', '30$^{\circ}$S', '25$^{\circ}$S',
                    '20$^{\circ}$S', '15$^{\circ}$S', '10$^{\circ}$S', '5$^{\circ}$S', '0$^{\circ}$S'])
ax.set_xlabel('V [m/s]')
ax.set_ylabel('Latitud')
ax.text(-0.1, 1, '(b)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/v_medio_cuenca.png', bbox_inches='tight')
plt.show()

# Rotor

u_ncep = ncep['curl'].sel(time=slice('1999-10-06', '2015-12-31')).mean(dim='time').squeeze()
u_cfsr = cfsr['curl'].sel(time=slice('1999-10-06', '2015-12-31')).mean(dim='time').squeeze()
u_ncep1 = ncep['curl'].sel(time=slice('2009-01-01', '2015-12-31')).mean(dim='time').squeeze()
u_cfsr1 = cfsr['curl'].sel(time=slice('2009-01-01', '2015-12-31')).mean(dim='time').squeeze()
u_ncep2 = ncep['curl'].mean(dim='time').squeeze()
u_cfsr2 = cfsr['curl'].mean(dim='time').squeeze()

f,ax = plt.subplots(1,figsize=(7,10))
ax.plot(u_ncep.mean(dim='lon')*1e7, ncep.lat, color='k', linestyle='--')
ax.plot(u_ncep1.mean(dim='lon')*1e7, ncep.lat, color='k')
ax.plot(u_ncep2.mean(dim='lon')*1e7, ncep.lat, color='k', linestyle=':')
ax.plot(u_cfsr.mean(dim='lon')*1e7, cfsr.lat, color='r', linestyle='--')
ax.plot(u_cfsr1.mean(dim='lon')*1e7, cfsr.lat, color='r')
ax.plot(u_cfsr2.mean(dim='lon')*1e7, cfsr.lat, color='r', linestyle=':')
ax.grid(which="both", color='grey', linestyle=':')
ax.set_yticks([-70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0])
ax.set_yticklabels(['70$^{\circ}$S', '65$^{\circ}$S', '60$^{\circ}$S', '55$^{\circ}$S', '50$^{\circ}$S', '45$^{\circ}$S',
                    '40$^{\circ}$S', '35$^{\circ}$S', '30$^{\circ}$S', '25$^{\circ}$S',
                    '20$^{\circ}$S', '15$^{\circ}$S', '10$^{\circ}$S', '5$^{\circ}$S', '0$^{\circ}$S'])
ax.set_xlabel('Rotor [10$^{-7}$ Pa/m]')
ax.set_ylabel('Latitud')
ax.text(-0.1, 1, '(c)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/resultados_1/crl_medio_cuenca.png', bbox_inches='tight')
plt.show()
