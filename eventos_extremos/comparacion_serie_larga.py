import xarray
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from scipy import signal
from scipy import stats
import matplotlib.ticker as mticker
from scipy.stats import chi2

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')
v = osc['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()
v = v.isel(time=slice(501,-1))
v = v.groupby('time.month') - v.groupby('time.month').mean('time')

s1, s2 = signal.butter(6, 1/5, 'low')
vf = signal.filtfilt(s1, s2, v.values)

idx_ED = np.where(vf>np.mean(vf)+np.std(vf))[0]
kk = np.diff(idx_ED)
ED = []; n=0; nn=0
for i in range(0, len(kk)):
    if kk[i]==1:
        nn = nn+1
    else:
        d = (nn-n+1)*5
        f_i = v['time'][idx_ED[n]].values
        v_mean = np.mean(vf[idx_ED[n]:idx_ED[nn]+1])
        v_max = np.max(vf[idx_ED[n]:idx_ED[nn]+1])
        n = nn+1
        nn = nn+1
        ED.append([f_i, d, v_mean, v_max])
del kk, f_i, d, v_mean, v_max
ED.append([v['time'][idx_ED[n]].values,
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
        f_i = v['time'][idx_EF[n]].values
        v_mean = np.mean(vf[idx_EF[n]:idx_EF[nn]+1])
        v_max = np.min(vf[idx_EF[n]:idx_EF[nn]+1])
        n = nn+1
        nn = nn+1
        EF.append([f_i, d, v_mean, v_max])
del kk, f_i, d, v_mean, v_max
EF.append([v['time'][idx_EF[n]].values,
           (idx_EF[-1]-idx_EF[n]+1)*5,
           np.mean(vf[idx_EF[n]:idx_EF[-1]+1]),
           np.min(vf[idx_EF[n]:idx_EF[-1]+1])])
EF = np.array(EF)
k = np.where(EF[:,1]>=20)
EFp = EF[k]; del EF
del idx_ED, idx_EF

idxED = []
for i in range(0, len(EDp[:,0])):
    t = np.where(v['time']==EDp[i,0])[0][0]
    idxED.append(t); idxED.append(int(EDp[i,1]/5+t))
idxEF = []
for i in range(0, len(EFp[:,0])):
    t = np.where(v['time']==EFp[i,0])[0][0]
    idxEF.append(t); idxEF.append(int(EFp[i,1]/5+t))

f, ax = plt.subplots(1, sharex=True, figsize=(10,3))
ax.set_xlim([v['time'][0].values, v['time'][-1].values])
ax.plot(v['time'], vf*100, color='k', linewidth=1)
ax.plot(v['time'], np.mean(vf)*np.ones(len(vf))*100, color='k', linewidth=2)
ax.plot(v['time'], (np.mean(vf)+np.std(vf))*np.ones(len(vf))*100, color='r', linestyle='--', linewidth=2)
ax.plot(v['time'], (np.mean(vf)-np.std(vf))*np.ones(len(vf))*100, color='b', linestyle='--', linewidth=2)
for i in range(0, int(len(idxED)/2)):
    ax.fill_between(v['time'][idxED[2*i]:idxED[2*i+1]].values, vf[idxED[2*i]:idxED[2*i+1]]*100,
                       np.mean(vf)*100+np.std(vf)*100, color='r', alpha=.5)
for i in range(0, int(len(idxEF)/2)):
    ax.fill_between(v['time'][idxEF[2*i]:idxEF[2*i+1]].values, vf[idxEF[2*i]:idxEF[2*i+1]]*100,
                       np.mean(vf)*100-np.std(vf)*100, color='b', alpha=.5)
ax.set_ylabel('$V_A$ [cm/s]')
ax.xaxis.set_minor_locator(mdates.YearLocator())
ax.grid(which="major", color='grey', linestyle=':')
#plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/identificacion_eventos_larga.png', bbox_inches='tight')
plt.show()

# Hovmoller ###################################################################

va = osc['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -30)).squeeze()
va = va.groupby('time.month') - va.groupby('time.month').mean('time')
x = np.arange(0, len(va['time']), 1)
xx, yy = np.meshgrid(va.longitude, x)

f, ax = plt.subplots(1, figsize=(6,10))
cl = ax.contourf(xx, yy, va*100, np.arange(-30,35,5), cmap='coolwarm', extend='both')
ze = ax.contour(xx, yy, va*100, levels=[0], colors='white')
for i in range(0, int(len(idxED)/2)):
    ax.plot((xx[0,1], xx[0,6]), (idxED[2*i+1]+500, idxED[2*i+1]+500), color='g')
    ax.plot((xx[0,1], xx[0,6]), (idxED[2*i]+500, idxED[2*i]+500), color='g')
for i in range(0, int(len(idxEF)/2)):
    ax.plot((xx[0,1], xx[0,6]), (idxEF[2*i+1]+500, idxEF[2*i+1]+500), color='k')
    ax.plot((xx[0,1], xx[0,6]), (idxEF[2*i]+500, idxEF[2*i]+500), color='k')

cbar = f.colorbar(cl, ax=ax, shrink=.7)
cbar.ax.set_ylabel('cm/s')
ax.set_xlabel('Longitud')
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/hovmoller_serielarga.png', bbox_inches='tight')
plt.show()
plt.close()

# Ahora que elegi los eventos hago el grafiquito y saco sus datas #############

idxED = np.array(idxED); idxEF = np.array(idxEF)

EDPE = np.array([0, 1, 2, 3, 4, 5, 8, 9, 22, 23, 28, 29, 32, 33, 34, 35, 38, 39, 40, 41, 42, 43, 44, 45])
kk = np.array([0, 1, 2, 4, 11, 14, 16, 17, 19, 20, 21, 22])

EFPE = np.array([4, 5, 8, 9, 10, 11, 14, 15, 20, 21, 26, 27, 28, 29, 30, 31, 36, 37, 40, 41, 42, 43, 44, 45])
kk2 = np.array([2, 4, 5, 7, 10, 13, 14, 15, 18, 20, 21, 22])

f, ax = plt.subplots(1, sharex=True, figsize=(10,3))
ax.set_xlim([v['time'][0].values, v['time'][-1].values])
ax.plot(v['time'], vf*100, color='k', linewidth=1)
ax.plot(v['time'], np.mean(vf)*np.ones(len(vf))*100, color='k', linewidth=2)
ax.plot(v['time'], (np.mean(vf)+np.std(vf))*np.ones(len(vf))*100, color='r', linestyle='--', linewidth=2)
ax.plot(v['time'], (np.mean(vf)-np.std(vf))*np.ones(len(vf))*100, color='b', linestyle='--', linewidth=2)
for i in range(0, int(len(idxED[EDPE])/2)):
    ax.fill_between(v['time'][idxED[EDPE[2*i]]:idxED[EDPE[2*i+1]]].values, vf[idxED[EDPE[2*i]]:idxED[EDPE[2*i+1]]]*100,
                       np.mean(vf)*100+np.std(vf)*100, color='r', alpha=.5)
for i in range(0, int(len(idxEF[EFPE])/2)):
    ax.fill_between(v['time'][idxEF[EFPE[2*i]]:idxEF[EFPE[2*i+1]]].values, vf[idxEF[EFPE[2*i]]:idxEF[EFPE[2*i+1]]]*100,
                       np.mean(vf)*100-np.std(vf)*100, color='b', alpha=.5)
ax.set_ylabel('$V_A$ [cm/s]')
ax.xaxis.set_minor_locator(mdates.YearLocator())
ax.grid(which="major", color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/identificacion_eventos_larga.png', bbox_inches='tight')
plt.show()

# COMPOSITES ##################################################################

import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
import shapely.geometry as sgeom
import cmocean
from mpl_toolkits.basemap import cm

osa = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anomclim_1992_2015.nc')
bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')

def plotter(x, lat, lon, u, v, lat1, lon1, cmap, clevs, bati, lat_i, lat_f, lon_i, lon_f, cbar_label, title, kk):

    x_o, y_o = np.meshgrid(lon1, lat1)
    x1, y1 = np.meshgrid(lon, lat)
    x_b, y_b = np.meshgrid(bati['x'], bati['y'])

    fig = plt.figure(figsize=(10,12))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([lon_i, lon_f, lat_f, lat_i], crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
    ax.add_feature(tierra)
    ax.coastlines(resolution='50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0),
                      draw_labels=True, color='white', linestyle='--', linewidth=.3)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'rotation': 45}
    gl.ylabel_style = {'size': 15}
    box = sgeom.box(minx=-51.67, maxx=-50.33, miny=-34.33, maxy=-34.33)
    ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none', linewidth=2.5,
                      edgecolor='black')

    clr = ax.contourf(x1, y1, x, clevs,
                      transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
    b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=2,
                   linestyles='-', transform=ccrs.PlateCarree())
    qvr = ax.quiver(x_o, y_o, u, v,
                    units='xy', scale=0.2/111139, transform=ccrs.PlateCarree())

    cbar = fig.colorbar(clr, ax=ax, shrink=.7)
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_ylabel(cbar_label, size=18)
    ax.quiverkey(qvr, .2, 0.8, .2, '20 cm/s', labelpos='E', fontproperties={'size': 18})
    ax.text(0.08, .93, kk, transform=ax.transAxes, size=22)
    plt.title(title)

    return fig, ax

AVM = osa['v_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).isel(time=slice(500, -1))
AUM = osa['u_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).isel(time=slice(500, -1))
bati = bat.sel(y=slice(-40,-26), x=slice(-58,-44))

###############################################################################
# EVENTOS DEBILES
###############################################################################

# DURANTE #####################################################################

idxEDd = np.array([1, 2, 3, 4, 5, 18, 19, 20, 21, 22, 23, 24, 25, 26, 52, 53, 54,
                   55, 56, 57, 58, 59, 60, 230, 231, 232, 233, 234, 235, 236, 237,
                   238, 464, 465, 466, 467, 468, 469, 558, 559, 560, 561, 562, 563,
                   564, 565, 566, 567, 568, 569, 683, 684, 685, 686, 687, 688, 689,
                   759, 760, 761, 762, 763, 764, 765, 766, 992, 993, 994, 995, 996,
                   1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1082,
                   1083, 1084, 1085, 1086, 1087, 1088, 1127, 1128, 1129, 1130, 1131,
                   1132])

AVMm = AVM.isel(time=idxEDd).mean(dim='time')
AUMm = AUM.isel(time=idxEDd).mean(dim='time')
N = len(idxEDd)

msk = (AVMm/(AVM.isel(time=idxEDd).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(BB*100, AVM.lat, AVM.lon, AA, BB, AVM.lat, AUM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composites_serie_larga/compositeL_ED_AVM.png', bbox_inches='tight')

# ANTES #####################################################################

idxA = np.array([16, 17, 15, 50, 51, 49, 228, 229, 227, 462, 463, 461, 556, 557, 556,
                 681, 682, 680, 756, 757, 755, 990, 991, 989, 1029, 1030, 1028, 1080,
                 1081, 1079, 1125, 1126, 1124])

AVMm = AVM.isel(time=idxA).mean(dim='time')
AUMm = AUM.isel(time=idxA).mean(dim='time')
N = len(idxA)

msk = (AVMm/(AVM.isel(time=idxA).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(BB*100, AVM.lat, AVM.lon, AA, BB, AVM.lat, AUM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composites_serie_larga/compositeL_ED_AVM_-15.png', bbox_inches='tight'); plt.show()

# DESPUES #####################################################################

idxD = np.array([6, 7, 8, 28, 29, 27, 62, 63, 61, 240, 241, 239, 471, 472, 470,
                 571, 572, 570, 691, 692, 690, 767, 768, 766, 998, 999, 997,
                 1042, 1043, 1041, 1090, 1091, 1089, 1134, 1135, 1133])

AVMm = AVM.isel(time=idxD).mean(dim='time')
AUMm = AUM.isel(time=idxD).mean(dim='time')
N = len(idxD)

msk = (AVMm/(AVM.isel(time=idxD).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(BB*100, AVM.lat, AVM.lon, AA, BB, AVM.lat, AUM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composites_serie_larga/compositeL_ED_AVM_+10.png', bbox_inches='tight')

###############################################################################
# EVENTOS FUERTES
###############################################################################

# DURANTE #####################################################################

idxEFd = np.array([128, 129, 130, 131, 132, 127, 216, 217, 218, 219, 220, 215,
                   241, 242, 243, 244, 245, 246, 247, 240, 456, 457, 458, 459,
                   455, 527, 528, 529, 530, 531, 532, 533, 526, 596, 597, 598,
                   599, 600, 601, 595, 610, 611, 612, 613, 609, 729, 730, 731,
                   732, 733, 728, 984, 985, 986, 987, 988, 989, 990, 983, 1044,
                   1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054,
                   1055, 1056, 1057, 1043, 1092, 1093, 1094, 1095, 1096, 1097,
                   1098, 1091, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1134])

AVMm = AVM.isel(time=idxEFd).mean(dim='time')
AUMm = AUM.isel(time=idxEFd).mean(dim='time')
N = len(idxEFd)

msk = (AVMm/(AVM.isel(time=idxEFd).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(BB*100, AVM.lat, AVM.lon, AA, BB, AVM.lat, AUM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composites_serie_larga/compositeL_EF_AVM.png', bbox_inches='tight')

# ANTES #######################################################################

idxA = np.array([125, 126, 124, 213, 214, 213, 238, 239, 237, 453, 454, 452, 524,
                 525, 523, 593, 594, 592, 607, 608, 606, 726, 727, 725, 981, 982,
                 980, 1041, 1042, 1040, 1089, 1090, 1088, 1132, 1133, 1131])

AVMm = AVM.isel(time=idxA).mean(dim='time')
AUMm = AUM.isel(time=idxA).mean(dim='time')
N = len(idxA)

msk = (AVMm/(AVM.isel(time=idxA).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(BB*100, AVM.lat, AVM.lon, AA, BB, AVM.lat, AUM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composites_serie_larga/compositeL_EF_AVM_-15.png', bbox_inches='tight'); plt.show()

# DESPUES #####################################################################

idxD = np.array([134, 135, 133, 222, 223, 221, 249, 250, 248, 461, 462, 460, 535,
                 536, 534, 603, 604, 602, 615, 616, 614, 735, 736, 734, 992, 993,
                 991, 1059, 1060, 1058, 1100, 1101, 1099, 1143, 1144, 1142])

AVMm = AVM.isel(time=idxD).mean(dim='time')
AUMm = AUM.isel(time=idxD).mean(dim='time')
N = len(idxD)

msk = (AVMm/(AVM.isel(time=idxD).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(BB*100, AVM.lat, AVM.lon, AA, BB, AVM.lat, AUM.lon,
              'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composites_serie_larga/compositeL_EF_AVM_+10.png', bbox_inches='tight')

# Distribucion temporal de eventos ############################################

n = np.array([0, 1, 2, 3])
d = np.array([3, 2, 4, 3])*100/12
ff = np.array([5, 2, 2, 3])*100/12

f, ax = plt.subplots(1, sharex=True, figsize=(6, 5))
ax.bar(n, d, color='darkgrey', width=0.6, edgecolor='b')
ax.set_ylabel('% del total')
ax.set_xticks([0,1,2,3])
ax.set_xticklabels(['EFM', 'AMJ', 'JAS', 'OND'])
ax.text(-0.15, 1, '(a)', transform=ax.transAxes, size=15)
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/estacionalidad_ed.png', bbox_inches='tight')
plt.show()

y = np.arange(1999, 2016, 1)
k = np.arange(0, 24, 1)
y1 = np.array([1, 18, 52, 230, 464, 558, 683, 759, 992, 1031, 1082, 1127])
d1 = np.array([20, 40, 40, 40, 25, 55, 30, 35, 20, 45, 30, 25])
y2 = np.array([127, 215, 240, 455, 526, 595, 609, 728, 983, 1043, 1091, 1134])
d2 = np.array([25, 25, 35, 20, 35, 30, 20, 25, 35, 70, 35, 35])
v1 = np.array([9.71, 11.21, 12.05, 20.62, 10.22, 19.09, 13.12, 13.90, 11.08, 16.11, 15.40, 12.94])
v2 = np.array([-14.91, -12.46, -10.60, -10.50, -13.86, -12.26, -12.53, -10.61, -15.48, -13.2, -15.88, -14.95])

f, ax = plt.subplots(1, sharex=True, figsize=(6, 5))
ax.plot(k, osa.time[np.sort(np.concatenate((y1, y2)), axis=0)], color='darkgrey', marker='o', linestyle='None')
ax.set_ylabel('% del total')
ax.text(-0.15, 1, '(a)', transform=ax.transAxes, size=15)
#plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/estacionalidad_ed.png', bbox_inches='tight')
plt.show()

# Veo en que porcentaje intensificaron o disminuyeron la CB en la dir media ###
