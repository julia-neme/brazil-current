"""
Función que grafica las velocidades de OSCAR en quiver, con algún campo de en
contornos detrás. Incluye la batimetria y la transecta de SAMOC en la CB.

def plotter(x, lat, lon, u, v, cmap, clevs, bati, lat_i, lat_f, lon_i, lon_f,
            cbar_label)

    x: campo a realizar en contornos.
    lat, lon: latitud y longitud del campo x.
    u, v = componentes de la velocidad de OSCAR, en formato xarray. Si tiene la
           dimensión extra (depth), sacarsela.
    cmap = colormap para x.
    clevs = niveles para x dados en un numpy array.
    bati = el archivo de batimetria dado en formato xarray. Esta función solo
           va a graficar los niveles -200m y -1000m.
    lat_i = latitud norte del mapa.
    lat_f = latitud sur del mapa.
    lon_i = longitud oeste del mapa.
    lon_f = longitud este del mapa.
    cbar_label = unidades de x.

@author: JB
"""


def plotter(x, lat, lon, u, v, cmap, clevs, bati, lat_i, lat_f, lon_i, lon_f,
            cbar_label):
    import cartopy.crs as ccrs
    from cartopy.feature import NaturalEarthFeature
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.gridspec as gridspec
    import shapely.geometry as sgeom
    import numpy as np
    import xarray

    x_o, y_o = np.meshgrid(u.longitude, v.latitude)
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
    box = sgeom.box(minx=-51.67, maxx=-50.33, miny=-34.33, maxy=-34.33)
    ax.add_geometries([box], ccrs.PlateCarree(), facecolor='none',
                      edgecolor='black')

    clr = ax.contourf(x1, y1, x, clevs,
                      transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
    b = ax.contour(x_b, y_b, bati['z'].values, levels=[-1000,-200], colors='k', linewidths=1,
                   linestyles, transform=ccrs.PlateCarree())
    qvr = ax.quiver(x_o, y_o, u.values, v.values,
                    units='xy', scale=0.2/111139, transform=ccrs.PlateCarree())

    cbar = fig.colorbar(clr, ax=ax, shrink=.7)
    cbar.ax.set_ylabel('m/s')
    ax.quiverkey(qvr, .2, 0.8, .2, '0.2 m/s', labelpos='E')

    return fig, ax
