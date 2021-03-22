import xarray
import numpy as np
import datetime as dt
import pandas as pd
import pytz
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import imageio
import shapely.geometry as sgeom
import matplotlib.patches as mpatches
import os

raw = xarray.open_dataset('/home/bock/Documents/tesis/datos/tracks_AVISO_DT2014_daily_web.nc', decode_cf=False)
del raw.j1.attrs['units']
ds = xarray.decode_cf(raw)

def julian2date(julian):
    """
    Calendar date from julian date.
    Works only for years past 1582!

    Parameters
    ----------
    julian : numpy.ndarray or double
        Julian day.

    Returns
    -------
    year : numpy.ndarray or int32
        Year.
    month : numpy.ndarray or int32
        Month.
    day : numpy.ndarray or int32
        Day.
    hour : numpy.ndarray or int32
        Hour.
    minute : numpy.ndarray or int32
        Minute.
    second : numpy.ndarray or int32
        Second.
    """
    min_julian = 2299160
    max_julian = 1827933925

    julian = np.atleast_1d(np.array(julian, dtype=float))

    if np.min(julian) < min_julian or np.max(julian) > max_julian:
        raise ValueError("Value of Julian date is out of allowed range.")

    jn = (np.round(julian + 0.0000001)).astype(np.int32)

    jalpha = (((jn - 1867216) - 0.25) / 36524.25).astype(np.int32)
    ja = jn + 1 + jalpha - (np.int32(0.25 * jalpha))
    jb = ja + 1524
    jc = (6680.0 + ((jb - 2439870.0) - 122.1) / 365.25).astype(np.int32)
    jd = (365.0 * jc + (0.25 * jc)).astype(np.int32)
    je = ((jb - jd) / 30.6001).astype(np.int32)

    day = jb - jd - np.int64(30.6001 * je)
    month = je - 1
    month = (month - 1) % 12 + 1
    year = jc - 4715
    year = year - (month > 2)

    fraction = (julian + 0.5 - jn).astype(np.float64)
    eps = (np.float64(1e-12) * np.abs(jn)).astype(np.float64)
    eps.clip(min=np.float64(1e-12), max=None)
    hour = (fraction * 24. + eps).astype(np.int64)
    hour.clip(min=0, max=23)
    fraction -= hour / 24.
    minute = (fraction * 1440. + eps).astype(np.int64)
    minute = minute.clip(min=0, max=59)
    second = (fraction - minute / 1440.) * 86400.
    second = second.clip(min=0, max=None)
    microsecond = ((second - np.int32(second)) * 1e6).astype(np.int32)
    microsecond = microsecond.clip(min=0, max=999999)
    second = second.astype(np.int32)

    return year, month, day, hour, minute, second, microsecond
def julian2datetime(julian, tz=None):
    """
    converts julian date to python datetime
    default is not time zone aware

    Parameters
    ----------
    julian : float
        julian date
    """
    year, month, day, hour, minute, second, microsecond = julian2date(julian)
    if type(julian) == np.array or type(julian) == np.memmap or \
            type(julian) == np.ndarray or type(julian) == np.flatiter:
        return np.array([dt.datetime(y, m, d, h, mi, s, ms, tz)
                         for y, m, d, h, mi, s, ms in
                         zip(np.atleast_1d(year),
                             np.atleast_1d(month),
                             np.atleast_1d(day),
                             np.atleast_1d(hour),
                             np.atleast_1d(minute),
                             np.atleast_1d(second),
                             np.atleast_1d(microsecond))])

    return dt.datetime(year, month, day, hour, minute, second, microsecond, tz)
f = julian2datetime(ds.j1.values)

ind_lat = np.where(np.abs(ds.lat.values+34.5) < 3)
ind_lon = np.where(np.abs(ds.lon.values-315) < 15)
ind_zon = np.intersect1d(ind_lat, ind_lon)

yr = [2009, 2009, 2010, 2011, 2011, 2011, 2013, 2014, 2014, 2015]
mt = [4, 12, 4, 5, 7, 9, 7, 1, 10, 5]
ind_ed = []
for j in range(0, len(yr)):
    ind_f = []
    for i in range(0, len(ind_zon)):
        if f[ind_zon[i]].year == yr[j] and f[ind_zon[i]].month == mt[j] and ds['cyc'][ind_zon[i]] == -1:
            ind_f.append(ind_zon[i])
    if len(ind_f) != 0:
        ind_ed.append(np.unique(ds.track[ind_f]))
    else:
        ind_ed.append([np.nan])

yr = [2009, 2010, 2012, 2013, 2013, 2014, 2014, 2015, 2015, 2015]
mt = [11, 6, 10, 6, 12, 4, 12, 2, 7, 10]
ind_ef = []
for j in range(0, len(yr)):
    ind_f = []
    for i in range(0, len(ind_zon)):
        if f[ind_zon[i]].year == yr[j] and f[ind_zon[i]].month == mt[j] and ds['cyc'][ind_zon[i]] == 1:
            ind_f.append(ind_zon[i])
    if len(ind_f) != 0:
        ind_ef.append(np.unique(ds.track[ind_f]))
    else:
        ind_ef.append([np.nan])

bat = xarray.open_dataset('/home/bock/Documents/tesis/batimetria/ETOPO1_Bed_g_gmt4.grd')
bati = bat.sel(y=slice(-40,-25), x=slice(-60,-35))
x, y = np.meshgrid(bati.x, bati.y)

def circle_degrees(radius, lat_i, lon_i):
    lon = np.empty([13])
    lat = np.empty([13])
    i = 0; n = 0
    while i < 360:
        lat[n] = radius*np.cos(np.deg2rad(i))
        lon[n] = radius*np.sin(np.deg2rad(i))
        n += 1
        i += 30
    lat[-1] = lat[0]; lon[-1] = lon[0]
    lat = lat_i + lat/(60*1.852)
    lon = lon_i + lon/(60*1.852*np.cos(np.deg2rad(lat_i)))
    return lat, lon

def grafico_campos(nro_ed, mnth, dy, output):

    fig = plt.figure(figsize=(9,7))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())

    ax.set_extent([-60, -35, -40, -25], crs=ccrs.PlateCarree())
    tierra = NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='white')
    ax.add_feature(tierra)
    ax.coastlines(resolution='50m')

    grid_style = dict(color='white', linestyle='--', linewidth=.3)
    gl = ax.gridlines(draw_labels=True, **grid_style)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15}
    gl.ylabel_style = {'size': 15}
    box = sgeom.box(minx=-51.67, maxx=-50.33, miny=-34.33, maxy=-34.33)
    ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none', edgecolor='yellow', linewidth=3)

    ax.contourf(x, y, bati['z'].values, cmap='binary_r', extend='both', transform=ccrs.PlateCarree())
    ax.contour(x, y, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1, transform=ccrs.PlateCarree())

    kk = np.where(ds['track'] == nro_ed)[0]
    ax.plot(ds['lon'][kk]-360, ds['lat'][kk], color='lime', linewidth=3, transform=ccrs.PlateCarree())
    for i in range(0, len(kk)):
        if f[kk[i]].month == mnth and f[kk[i]].day == dy:

            ax.plot(ds['lon'][kk[i]]-360, ds['lat'][kk[i]], color='red', linewidth=3, marker="^", transform=ccrs.PlateCarree())
            llat, llon = circle_degrees(ds['L'][kk[i]].item(), ds['lat'][kk[i]].item(), ds['lon'][kk[i]].item()-360)
            ax.plot(llon, llat, color='red', linewidth=3, transform=ccrs.PlateCarree())
            radius = ds['L'][kk[i]]
            print(radius)
    salida = '/home/bock/Documents/' + output + '.png'
    plt.savefig(salida, bbox_inches='tight', dpi=250)

grafico_campos(ind_ed[4][5], 7, 27, 'ed5')
plt.close('all')

for i in range(0, len(ind_ef[8])):
    grafico_campos(ind_ef[8][i], 2, 5, 'ef9')
plt.show()

"""
D1: 2; 117.42km
D2: 9; 97.16km
D3: -
D4: 4; 72.94km
D5: 5; 57.08km
D6: 1; 101.28km
D7: 12; 87.83km
D8: 8; 106.45km
D9: 4; 81.78km
F1: 9; 79.33km
F2: -
F3: 6; 109.42km
F4: 0; 119.39km
F5: 6; 97.72km
F6: 8; 135.46km
F7: 1; 77.90km
F8: 4; 64.33km
F9: -
F10: -
"""
