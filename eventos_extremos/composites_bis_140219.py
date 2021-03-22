import xarray
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
import shapely.geometry as sgeom
import cmocean
from mpl_toolkits.basemap import cm

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_2009_2015.nc')
osa = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')
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
    gl.xlabel_style = {'size': 20, 'rotation': 45}
    gl.ylabel_style = {'size': 20}
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
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_ylabel(cbar_label, size=20)
    ax.quiverkey(qvr, .2, 0.8, .2, '20 cm/s', labelpos='E', fontproperties={'size': 20})
    ax.text(0.08, .93, kk, transform=ax.transAxes, size=22)
    plt.title(title)

    return fig, ax

U = osc['u'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44)).squeeze()
V = osc['v'].sel(latitude=slice(-26,-40), longitude=slice(-58,-44)).squeeze()
SSTA = gha['sst_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44)).sel(time=osc['time'], method='nearest').sel(lat=U.latitude, method='nearest').sel(lon=U.longitude, method='nearest')
AVM = osa['v_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44))
AUM = osa['u_anom'].sel(lat=slice(-26,-40), lon=slice(-58,-44))
bati = bat.sel(y=slice(-56,-20), x=slice(-80,-40))

###############################################################################
# EVENTOS DEBILES
###############################################################################

# DURANTE #####################################################################

idxED = np.array([18, 19, 20, 21, 22, 23, 24, 93, 94, 95,
                  96, 97, 98, 99, 100, 101, 172, 173, 174, 175, 176, 177, 326, 327,
                  328, 329, 330, 331, 366, 367, 368, 369, 370, 371, 372, 373, 375,
                  417, 418, 419, 420, 421, 422, 423, 461, 462, 463, 464, 465, 466, 467])

SSTAm = SSTA.isel(time=idxED).mean(dim='time')
AVMm = AVM.isel(time=idxED).mean(dim='time')
AUMm = AUM.isel(time=idxED).mean(dim='time')

N = len(idxED)

msk = (AVMm/(AVM.isel(time=idxED).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(h)')
plt.savefig('/home/bock/Documents/kk.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxED).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(k)')
plt.savefig('/home/bock/Documents/kk1.png', bbox_inches='tight'); plt.show()

# ANTES #######################################################################

idxA = np.array([15, 16, 17, 90, 91, 92, 169, 170, 171, 323, 324, 325, 363, 364,
                 365, 414, 415, 416, 458, 459, 460])

SSTAm = SSTA.isel(time=idxA).mean(dim='time')
AVMm = AVM.isel(time=idxA).mean(dim='time')
AUMm = AUM.isel(time=idxA).mean(dim='time')

N = len(idxA)

msk = (AVMm/(AVM.isel(time=idxA).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(g)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_AVM_antes.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6669)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6669)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6669)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(j)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_SSTA_antes.png', bbox_inches='tight'); plt.show()

# DESPUES #####################################################################

idxD = np.array([25, 26, 27, 102, 103, 104, 178, 179, 180, 332, 333, 334, 376,
                 377, 378, 424, 425, 426, 468, 469, 470])

SSTAm = SSTA.isel(time=idxD).mean(dim='time')
AVMm = AVM.isel(time=idxD).mean(dim='time')
AUMm = AUM.isel(time=idxD).mean(dim='time')

N = len(idxD)

msk = (AVMm/(AVM.isel(time=idxD).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon,  BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(i)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_AVM_despues.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxD).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6933)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(l)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_SSTA_despues.png', bbox_inches='tight'); plt.show()

# ANTES LAG 1 #################################################################

idxA1 = idxA-1

SSTAm = SSTA.isel(time=idxA1).mean(dim='time')
AVMm = AVM.isel(time=idxA1).mean(dim='time')
AUMm = AUM.isel(time=idxA1).mean(dim='time')

N = len(idxA1)

msk = (AVMm/(AVM.isel(time=idxA1).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_AVM_antes_1.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA1).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6933)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(f)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_SSTA_antes_1.png', bbox_inches='tight'); plt.show()

# ANTES LAG 2 #################################################################

idxA2 = idxA-2

SSTAm = SSTA.isel(time=idxA2).mean(dim='time')
AVMm = AVM.isel(time=idxA2).mean(dim='time')
AUMm = AUM.isel(time=idxA2).mean(dim='time')

N = len(idxA2)

msk = (AVMm/(AVM.isel(time=idxA2).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_AVM_antes_2.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA2).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6933)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(e)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_SSTA_antes_2.png', bbox_inches='tight'); plt.show()

# ANTES LAG 3 #################################################################

idxA3 = idxA-3

SSTAm = SSTA.isel(time=idxA3).mean(dim='time')
AVMm = AVM.isel(time=idxA3).mean(dim='time')
AUMm = AUM.isel(time=idxA3).mean(dim='time')

N = len(idxA3)

msk = (AVMm/(AVM.isel(time=idxA3).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_AVM_antes_3.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA3).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6933)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_SSTA_antes_3.png', bbox_inches='tight'); plt.show()

# ANTES LAG 4 #################################################################

idxA4 = idxA-4

Um = U.isel(time=idxA4).mean(dim='time')
Vm = V.isel(time=idxA4).mean(dim='time')
SSTAm = SSTA.isel(time=idxA4).mean(dim='time')
AVMm = AVM.isel(time=idxA4).mean(dim='time')
AUMm = AUM.isel(time=idxA4).mean(dim='time')

N = len(idxA4)

msk = (AVMm/(AVM.isel(time=idxA4).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.6933)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_AVM_antes_4.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA4).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.6933)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_SSTA_antes_4.png', bbox_inches='tight'); plt.show()

# ANTES LAG 5 #################################################################

idxA5 = idxA-5

Um = U.isel(time=idxA5).mean(dim='time')
Vm = V.isel(time=idxA5).mean(dim='time')
SSTAm = SSTA.isel(time=idxA5).mean(dim='time')
AVMm = AVM.isel(time=idxA5).mean(dim='time')
AUMm = AUM.isel(time=idxA5).mean(dim='time')

N = len(idxA5)

msk = (AVMm/(AVM.isel(time=idxA5).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.6933)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_AVM_antes_5.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA5).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.6933)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_SSTA_antes_5.png', bbox_inches='tight'); plt.show()

# ANTES LAG 6 #################################################################

idxA6 = idxA-6

Um = U.isel(time=idxA6).mean(dim='time')
Vm = V.isel(time=idxA6).mean(dim='time')
SSTAm = SSTA.isel(time=idxA6).mean(dim='time')
AVMm = AVM.isel(time=idxA6).mean(dim='time')
AUMm = AUM.isel(time=idxA6).mean(dim='time')

N = len(idxA6)

msk = (AVMm/(AVM.isel(time=idxA6).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.6933)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_AVM_antes_6.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA6).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.6933)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.6933)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6933)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_ED_SSTA_antes_6.png', bbox_inches='tight'); plt.show()


###############################################################################
# EVENTOS FUERTES
###############################################################################

# Durante #####################################################################

idxEF = np.array([63, 64, 65, 66, 67, 318, 319, 320, 321, 322, 323, 324, 325,
                  378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389,
                  390, 391, 426, 427, 428, 429, 430, 431, 432, 433, 470, 471,
                  472, 473, 474, 475])

SSTAm = SSTA.isel(time=idxEF).mean(dim='time')
AVMm = AVM.isel(time=idxEF).mean(dim='time')
AUMm = AUM.isel(time=idxEF).mean(dim='time')

N = len(idxEF)

msk = (AVMm/(AVM.isel(time=idxEF).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.6839)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.6839)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(h)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxEF).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.6839)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.6839)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.6839)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', 'k)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA.png', bbox_inches='tight'); plt.show()

# ANTES LAG 0 #################################################################

idxA = np.array([60, 61, 62, 315, 316, 317, 375, 376, 377, 423, 424, 425, 467, 468, 469])

SSTAm = SSTA.isel(time=idxA).mean(dim='time')
AVMm = AVM.isel(time=idxA).mean(dim='time')
AUMm = AUM.isel(time=idxA).mean(dim='time')

N = len(idxA)

msk = (AVMm/(AVM.isel(time=idxA).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(g)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM_antes.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(j)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA_antes.png', bbox_inches='tight'); plt.show()

# DESPUES #####################################################################

idxD = np.array([68, 69, 70, 326, 327, 328, 392, 393, 394, 434, 435, 436, 476, 477, 478])

SSTAm = SSTA.isel(time=idxD).mean(dim='time')
AVMm = AVM.isel(time=idxD).mean(dim='time')
AUMm = AUM.isel(time=idxD).mean(dim='time')

N = len(idxD)

msk = (AVMm/(AVM.isel(time=idxD).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(i)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM_despues.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxD).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(l)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA_despues.png', bbox_inches='tight'); plt.show()

# ANTES LAG 1 #################################################################

idxA1 = idxA-1

SSTAm = SSTA.isel(time=idxA1).mean(dim='time')
AVMm = AVM.isel(time=idxA1).mean(dim='time')
AUMm = AUM.isel(time=idxA1).mean(dim='time')

N = len(idxA1)

msk = (AVMm/(AVM.isel(time=idxA1).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(c)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM_antes_1.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA1).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(f)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA_antes_1.png', bbox_inches='tight'); plt.show()

# ANTES LAG 2 #################################################################

idxA2 = idxA-2

SSTAm = SSTA.isel(time=idxA2).mean(dim='time')
AVMm = AVM.isel(time=idxA2).mean(dim='time')
AUMm = AUM.isel(time=idxA2).mean(dim='time')

N = len(idxA2)

msk = (AVMm/(AVM.isel(time=idxA2).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM_antes_2.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA2).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(e)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA_antes_2.png', bbox_inches='tight'); plt.show()

# ANTES LAG 3 #################################################################

idxA3 = idxA-3

SSTAm = SSTA.isel(time=idxA3).mean(dim='time')
AVMm = AVM.isel(time=idxA3).mean(dim='time')
AUMm = AUM.isel(time=idxA3).mean(dim='time')

N = len(idxA3)

msk = (AVMm/(AVM.isel(time=idxA3).std(dim='time'))*np.sqrt(N-1))
BB = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, BB, CC, AVM.lat, AVM.lon,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM_antes_3.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA3).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(AUMm, mask=np.abs(msk)<1.7531)
BB = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB,AVM.lat, AVM.lon,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(d)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA_antes_3.png', bbox_inches='tight'); plt.show()

# ANTES LAG 4 #################################################################

idxA4 = idxA-4

Um = U.isel(time=idxA4).mean(dim='time')
Vm = V.isel(time=idxA4).mean(dim='time')
SSTAm = SSTA.isel(time=idxA4).mean(dim='time')
AVMm = AVM.isel(time=idxA4).mean(dim='time')
AUMm = AUM.isel(time=idxA4).mean(dim='time')

N = len(idxA4)

msk = (AVMm/(AVM.isel(time=idxA4).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.7531)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM_antes_4.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA4).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.7531)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA_antes_4.png', bbox_inches='tight'); plt.show()

# ANTES LAG 5 #################################################################

idxA5 = idxA-5

Um = U.isel(time=idxA5).mean(dim='time')
Vm = V.isel(time=idxA5).mean(dim='time')
SSTAm = SSTA.isel(time=idxA5).mean(dim='time')
AVMm = AVM.isel(time=idxA5).mean(dim='time')
AUMm = AUM.isel(time=idxA5).mean(dim='time')

N = len(idxA5)

msk = (AVMm/(AVM.isel(time=idxA5).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.7531)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM_antes_5.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA5).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.7531)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA_antes_5.png', bbox_inches='tight'); plt.show()

# ANTES LAG 6 #################################################################

idxA6 = idxA-6

Um = U.isel(time=idxA6).mean(dim='time')
Vm = V.isel(time=idxA6).mean(dim='time')
SSTAm = SSTA.isel(time=idxA6).mean(dim='time')
AVMm = AVM.isel(time=idxA6).mean(dim='time')
AUMm = AUM.isel(time=idxA6).mean(dim='time')

N = len(idxA6)

msk = (AVMm/(AVM.isel(time=idxA6).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.7531)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(AVMm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC*100, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 'coolwarm', np.arange(-20, 21, 1), bati, -26, -40, -58, -44, 'cm/s', '', '(a)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_AVM_antes_6.png', bbox_inches='tight'); plt.show()

msk = (SSTAm/(SSTA.isel(time=idxA6).std(dim='time'))*np.sqrt(N-1))
AA = np.ma.array(Um, mask=np.abs(msk)<1.7531)
BB = np.ma.array(Vm, mask=np.abs(msk)<1.7531)
CC = np.ma.array(SSTAm, mask=np.abs(msk)<1.7531)

fig, ax = plotter(CC, AVM.lat, AVM.lon, AA, BB, Um.latitude, Um.longitude,
                 cm.GMT_no_green, np.arange(-1, 1.01, .01), bati, -26, -40, -58, -44, '$^{\circ}$ C', '', '(b)')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/composite_EF_SSTA_antes_6.png', bbox_inches='tight'); plt.show()
