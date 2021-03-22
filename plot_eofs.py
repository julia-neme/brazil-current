import numpy as np
import xarray
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from functions import climatologia_xarray
from functions import get_coasts
from eofs.standard import Eof
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def eofs_as(dat):
    A = climatologia_xarray(dat['curl']).values
    global land
    EC, WC, land = get_coasts(dat.lat, dat.lon)

    msk = np.empty(np.shape(A))
    for i in range(0, len(A[:,0,0])):
        msk[i,:,:] = land
        B = np.ma.array(A, mask=msk)
    from get_eddof import get_eddof
    edof = np.empty([len(dat.lat), len(dat.lon)])
    for i in range(0, len(dat.lat)):
        for j in range(0, len(dat.lon)):
            if msk[0,i,j] == False:
                edof[i,j] = get_eddof(B[:,i,j])
            else:
                edof[i,j] = np.nan

    dof = int(np.nanmean(edof))
    coslat = np.cos(np.deg2rad(dat.lat.values)).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = Eof(B, center=True, weights=wgts, ddof=dof)

    eof = solver.eofs(neofs=10, eofscaling=2)
    pc = solver.pcs(npcs=10, pcscaling=1)
    varfrac = solver.varianceFraction()
    eigvals = solver.eigenvalues()

    x, y = np.meshgrid(dat.lon, dat.lat)

    return eof, pc, varfrac, x, y, edof

dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ncep_v1/NCEP1_curl_daily_2009-2015.nc').squeeze()
NC1_eof, NC1_pc, NC1_vf, NC1_x, NC1_y = eofs_as(dat); del dat
dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ncep_v2/NCEP2_curl_daily_2009-2015.nc').squeeze()
NC2_eof, NC2_pc, NC2_vf, NC2_x, NC2_y = eofs_as(dat); del dat
dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/ccmp/CCMP_curl_daily_2009-2015.nc').squeeze()
CCM_eof, CCM_pc, CCM_vf, CCM_x, CCM_y = eofs_as(dat); del dat
dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/cfsr/CFSR_curl_daily_2009-2015.nc').squeeze()
CFS_eof, CFS_pc, CFS_vf, CFS_x, CFS_y = eofs_as(dat); del dat
dat = xarray.open_dataset('/home/bock/Documents/tesis/vientos/era_int/ERA_curl_daily_2009-2015.nc').squeeze()
ERA_eof, ERA_pc, ERA_vf, ERA_x, ERA_y = eofs_as(dat); del dat


fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection=ccrs.Mercator())

ax.set_extent([NC1_x[0,0], NC1_x[0,-1], NC1_y[0,0], NC1_y[-1,0]],
              crs=ccrs.PlateCarree())
tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='grey',
         facecolor='lightgrey')
ax.add_feature(tierra)
ax.coastlines(resolution='50m')

grid_style = dict(color='white', linestyle='--', linewidth=.3)
gl = ax.gridlines(draw_labels=True, **grid_style)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

pc = ax.pcolormesh(NC1_x, NC1_y, NC1_eof[1,:,:]*1e7, transform=ccrs.PlateCarree(),
                   cmap='jet', vmin=-.75, vmax=.75)
cbar = fig.colorbar(pc, ax=ax, shrink=0.9)
cbar.ax.set_ylabel('10$^{-7}$ Pa/m') #
ax.set_title('NCEP vI - EOF 2 - 17.88%') #
plt.savefig('figuras/eof2_ncepv1.png', bbox_inches='tight') #
plt.close()

t = np.arange(1,13,1)
plt.figure(figsize=(10,7))
plt.plot(t, NC1_pc[:,1], label='NCEP vI')
plt.plot(t, NC2_pc[:,1], label='NCEP vII')
plt.plot(t, CCM_pc[:,1], label='CCMP')
plt.plot(t, CFS_pc[:,1], label='CFSR')
plt.plot(t, ERA_pc[:,1], label='ERA-Interim')
plt.legend()
plt.grid(linestyle=':')
plt.xlabel('Mes')
plt.title('Componentes principales del segundo EOF')
plt.savefig('figuras/pc2_as.png')
plt.close()
