# get_coasts : los puntos de grilla más cercanos a las costas este y oeste ####

def get_coasts(lat, lon):
    import numpy as np
    from mpl_toolkits.basemap import Basemap

    lon1 = np.empty(np.shape(lon))
    for i in range(0, len(lon)):
        if lon[i] > 30:
            lon1[i] = -360+lon[i]
        else:
            lon1[i] = lon[i]

    llclon = np.min(lon1); llclat = np.min(lat)
    urclon = np.max(lon1); urclat = np.max(lat)
    m = Basemap(projection='merc', area_thresh=10000, llcrnrlat=llclat,
                urcrnrlat=urclat, llcrnrlon=llclon, urcrnrlon=urclon,
                resolution='i')

    land = np.empty([len(lat), len(lon1)])
    for i in range(0, len(lat)):
        for j in range(0, len(lon1)):
            x, y = m(lon1[j],lat[i])
            land[i,j] = m.is_land(x,y)

    EC = np.empty(len(lat))
    WC = np.empty(len(lat))
    ss = int(len(lon1)/2)
    kk = np.diff(land, axis=1)
    for i in range(0, len(lat)):
        if any(kk[i,:] == -1):
            WC[i] = int(np.where(kk[i, 1:ss] == -1)[0][0]) + 2
        else:
            WC[i] = 0

    for i in range(0, len(lat)):
        if any(kk[i,ss:] == 1):
            EC[i] = int(np.where(kk[i, ss:] == 1)[0][0]) + 1 + ss
        else:
            EC[i] = len(lon1)
    return EC, WC, land

# zonal_integration : integra zonalmente entre costa este y oeste con una regla
# de simpson.

def zonal_integration(x, lat, lon):
    import numpy as np
    from functions import get_coasts
    from scipy import integrate

    EC, WC = get_coasts(lat, lon)

    x_int = np.empty([len(x[:,0,0]), len(lat)])
    for i in range(1, len(lat)):
        n = len(lon[int(WC[i]):int(EC[i])])
        h = np.abs(lon[1]-lon[0])*60*1.852*1000*np.cos(lat[i]*np.pi/180)
        xx = np.empty(n); xx[0] = 0
        for j in range(0, n-1):
            xx[j+1] = xx[j] + h
        x_int[:, i] = integrate.simps(x[:, i, int(WC[i]):int(EC[i])], xx)

    return x_int

# climatologia : hace la climatologia del campo. El tiempo tiene que estar en #
# un xarray.

def climatologia_xarray(x):
    import xarray

    CL = x.groupby('time.month').mean('time')

    return CL

# climatologia : hace la climatologia del campo. El tiempo tiene que estar en #
# numpyarray.

def climatologia_ncdf4(x, time):
    import numpy as np

    m = np.arange(1,13,1)
    mths = np.empty(len(time))
    for i in range(0, len(time)):
        mths[i] = time[i].month

    CL = np.empty([12, len(lat), len(lon)])
    for i in range(1, 13):
        a = np.where(mths == i)
        CL[i-1,:,:] = np.mean(x[a,:,:], axis=1)[0,:,:]
        del a

    return CL
