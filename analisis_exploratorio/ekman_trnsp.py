"""
Calcula el tranporte zonal y meridional de Ekman a partir de la tension del
vientos.

@author: Julia Neme
"""

def ek_trnsp(data):
    import xarray
    import numpy as np

    tx = dat['taux']; ty = dat['tauy']
    f = 2*7.29e-5*np.sin(np.deg2rad(dat.lat))
    Mx = ty/f
    My = -tx/f

    ET = [Mx, My]
    return ET
