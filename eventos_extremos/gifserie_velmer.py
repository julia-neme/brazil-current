import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
import shapely.geometry as sgeom
import matplotlib.dates as mdates
from scipy import stats
from scipy import signal
import imageio
import os

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')

osca = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
U = osca['u_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44))
V = osca['v_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44))
serie = osca['v_anom'].sel(lat=-34.5, method='nearest').sel(lon=slice(-52, -48)).sel(lon=slice(-51.66, -50.33)).mean(dim='lon')
x1, y1 = np.meshgrid(V.lon, V.lat)
bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')
bati = bat.sel(y=slice(-56,-20), x=slice(-80,-40))
x_b, y_b = np.meshgrid(bati['x'], bati['y'])

s1, s2 = signal.butter(6, 1/5, 'low')
vf = signal.filtfilt(s1, s2, serie)
yM = np.mean(vf)*np.ones(len(vf))
ySM = np.mean(vf)*np.ones(len(vf)) + np.std(vf)
ySm = np.mean(vf)*np.ones(len(vf)) - np.std(vf)

idx_ED = np.where(vf>np.mean(vf)+np.std(vf))[0]
kk = np.diff(idx_ED)
ED = []; n=0; nn=0
for i in range(0, len(kk)):
    if kk[i]==1:
        nn = nn+1
    else:
        d = (nn-n+1)*5
        f_i = osc['time'][idx_ED[n]].values
        v_mean = np.mean(vf[idx_ED[n]:idx_ED[nn]+1])
        v_max = np.max(vf[idx_ED[n]:idx_ED[nn]+1])
        n = nn+1
        nn = nn+1
        ED.append([f_i, d, v_mean, v_max])
del kk, f_i, d, v_mean, v_max
ED.append([osc['time'][idx_ED[n]].values,
           (idx_ED[-1]-idx_ED[n]+1)*5,
           np.mean(vf[idx_ED[n]:idx_ED[-1]+1]),
           np.max(vf[idx_ED[n]:idx_ED[-1]+1])])
ED = np.array(ED)
k = np.where(ED[:,1]>=20)
EDp = ED[k]; del ED

idx_EF = np.where(vf<np.mean(vf)-np.std(vf))[0]
kk = np.diff(idx_EF)
EF = []; n=0; nn=0
for i in range(0, len(kk)):
    if kk[i]==1:
        nn = nn+1
    else:
        d = (nn-n+1)*5
        f_i = osc['time'][idx_EF[n]].values
        v_mean = np.mean(vf[idx_EF[n]:idx_EF[nn]+1])
        v_max = np.min(vf[idx_EF[n]:idx_EF[nn]+1])
        n = nn+1
        nn = nn+1
        EF.append([f_i, d, v_mean, v_max])
del kk, f_i, d, v_mean, v_max
EF.append([osc['time'][idx_EF[n]].values,
           (idx_EF[-1]-idx_EF[n]+1)*5,
           np.mean(vf[idx_EF[n]:idx_EF[-1]+1]),
           np.min(vf[idx_EF[n]:idx_EF[-1]+1])])
EF = np.array(EF)
k = np.where(EF[:,1]>=20)
EFp = EF[k]; del EF
del idx_ED, idx_EF

idxED = []
for i in range(0, len(EDp[:,0])):
    t = np.where(osc['time']==EDp[i,0])[0][0]
    idxED.append(t); idxED.append(int(EDp[i,1]/5+t))
idxEF = []
for i in range(0, len(EFp[:,0])):
    t = np.where(osc['time']==EFp[i,0])[0][0]
    idxEF.append(t); idxEF.append(int(EFp[i,1]/5+t))

UU = np.empty(np.shape(U.values))
VV = np.empty(np.shape(V.values))
for i in range(0, len(U.lat)):
    for j in range(0, len(V.lon)):
        UU[:, i, j] = signal.filtfilt(s1, s2, U[:,i,j].values)
        VV[:, i, j] = signal.filtfilt(s1, s2, V[:,i,j].values)

def plotter(t):

    gs = gridspec.GridSpec(10, 1)
    fig = plt.figure(figsize=(12, 12))
    plt.clf()

    yr = str(osc['time.year'][t].values)
    mn = str(osc['time.month'][t].values).zfill(2)
    dy = str(osc['time.day'][t].values).zfill(2)
    plt.suptitle(yr + " " + mn + " " + dy)

    ax1 = plt.subplot(gs[1:7, 0], projection=ccrs.PlateCarree())
    ax1.set_title('AVM (colores)')
    ax1.set_extent([-58, -44, -40, -26], crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
    ax1.add_feature(tierra)
    ax1.coastlines(resolution='50m')
    gl = ax1.gridlines(color='white', linestyle='--', linewidth=.5, draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    box = sgeom.box(minx=-51.66, maxx=-50.33, miny=-34.33, maxy=-34.33)
    ax1.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='lime', linewidth=2.7)

    cl = ax1.contourf(x1, y1, VV[t,:,:]*100, np.arange(-50,55,5), cmap='RdBu_r', transform=ccrs.PlateCarree(), extend='both')
    qvr = ax1.quiver(x1, y1, UU[t,::1,::1], VV[t,::1,::1], units='inches', scale=.8,
                     transform=ccrs.PlateCarree())
    b = ax1.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1.4,
                 linestyles='-', transform=ccrs.PlateCarree())
    cbar = fig.colorbar(cl, ax=ax1, orientation='vertical', shrink=.7, pad=.08, label='m/s')
    ax1.quiverkey(qvr, .2, 0.8, .2, '20 cm/s', labelpos='E')

    ax3 = plt.subplot(gs[8:10, 0])
    ax3.plot(serie['time'], vf*100, 'k', linewidth=1)
    ax3.plot(serie['time'], yM*100, 'k', linewidth=2.5)
    ax3.plot(serie['time'], ySM*100, 'k', linewidth=2.5, linestyle='--')
    ax3.plot(serie['time'], ySm*100, 'k', linewidth=2.5, linestyle='--')
    for i in range(0, int(len(idxED)/2)):
        ax3.fill_between(serie['time'][idxED[2*i]:idxED[2*i+1]].values, vf[idxED[2*i]:idxED[2*i+1]]*100,
                           (np.mean(vf)+np.std(vf))*100, color='blue', alpha=.7)
    for i in range(0, int(len(idxEF)/2)):
        ax3.fill_between(serie['time'][idxEF[2*i]:idxEF[2*i+1]].values, vf[idxEF[2*i]:idxEF[2*i+1]]*100,
                           (np.mean(vf)-np.std(vf))*100, color='red', alpha=.7)
    ax3.plot(serie['time'][t].values, vf[t]*100, marker='o', color='k')
    ax3.xaxis.set_minor_locator(mdates.MonthLocator())
    ax3.xaxis.grid(which="minor", color='grey', linestyle=':')
    plt.ylabel('AVM [cm/s]')

    plt.tight_layout()
    fig.subplots_adjust(hspace=.2, wspace=.2)
    plt.savefig('/home/bock/Documents/tesis/resultados/fig_for_an/ee_'+yr+mn+dy+'.png', bbox_inches='tight')

for i in range(0, len(osca['time'])):
    plotter(i)
    plt.close()

imagenes = []
os.chdir('/home/bock/Documents/tesis/resultados/fig_for_an/')
files = sorted(os.listdir())
for f in files:
    imagenes.append(imageio.imread(f))

imageio.mimsave('/home/bock/Documents/tesis/resultados/animaciones/new.gif', imagenes, fps=6)
