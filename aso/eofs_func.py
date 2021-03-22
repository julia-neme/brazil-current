"""
Calcula los eofs como campos espaciales y las componentes principales como
series temporales. Al correr se selecciona la cantidad de eofs a extraer.

Los datos deben ser procesados previamente (sacarles tendencias y ciclos)

@author: Julia Neme
"""

def get_eofs(x):
    import numpy as np
    import xarray
    from eofs.xarray import Eof
    from matplotlib import pyplot as plt

    coslat = np.cos(np.deg2rad(x.lat.values)).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]

    from eddof import get_eddof
    DF = np.empty(np.shape(x[0,:,:]))
    for i in range(0, len(x.lat)):
        for j in range(0, len(x.lon)):
            DF[i, j] = get_eddof(x[:, i, j].values)
    edof = np.nanmin(DF)
    print(edof)
    solver = Eof(x, weights=wgts, ddof=edof)

    var = solver.varianceFraction()
    plt.figure(1)
    plt.bar(np.arange(0, len(var), 1), var*100)
    plt.show()
    plt.close()

    n = input('Cuantos PC extraer: '); n = int(n)
    eof = solver.eofs(neofs=n, eofscaling=2)
    pc = solver.pcs(npcs=n, pcscaling=1)
    vf = var[:n]

    EOFs = [eof, pc, vf]
    return EOFs
