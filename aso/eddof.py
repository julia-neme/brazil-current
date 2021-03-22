"""
Calcula los grados de libertad efectivos (EDoF) segun la expresion Ndt/T, donde
N es la cantidad de datos, dt es el intervalo muestral y T es la integral time
scale.

T se calcula integrando el area bajo la curva de la funcion de autocorrelacion
normalizada hasta el primer zero-crossing y dividiendo por la varianza de la
serie.

Ref.: Emery & Thomson
@author: Julia Neme
"""

def get_eddof(x):
    import numpy as np
    from scipy import integrate
    x = x[~np.isnan(x)]
    if np.size(x) == 0:
        edof = np.nan
    else:
        N = np.nanmax(len(x))
        N1 = N-1
        x = x - np.nanmean(x)

        # Calcula la funcion de autocorrelacion para lags desde
        # -N1 a N1. Por lo tanto, la correlacion sin lag (la
        # varianza de la serie) estara en c[N1] = c[N-1].

        c = np.correlate(x, x, 'full')

        # Normaliza la funcion de autocorrelacion segun N-1-k donde
        # k es el numero de lags, positivo.

        lags = np.abs(np.arange(-N1+1, N1, 1))
        cn = c[1:-1]/(N-1-lags)
        Var = cn[N1-1]

        # Busca el primer zero-crossing

        n = 0
        while (cn[N1+n] > 0) and (n < N1):
            n = n+1

        # Calcula el tiempo integral y los EDoF

        T = integrate.simps(cn[N1-1-n:N1+n])/Var

        edof = N/T
        if (np.isnan(edof) == False) and (np.isinf(edof) == False):
            edof = int(edof)
        else:
            edof = np.nan
        return edof
