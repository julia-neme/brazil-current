"""
Calcula el espectro de Welch y sus intervalos de confianza segun la dist.
chi2. Se usa un overlap de 50% y se supone que no se hace ningun frequency-
band averaging.

La frecuencia de muestreo es Fs = 1.

@author: Julia Neme
"""

def welch_spectrum(x, npseg, alfa):
    from scipy import signal
    from scipy.stats import chi2
    from matplotlib import pyplot as plt
    import numpy as np
    # Calcula el espectro

    nfft = npseg
    ovlp = npseg/2
    f, Pxx = signal.welch(x, fs=1, nperseg=npseg, noverlap=ovlp)

    # Calcula los niveles de confianza

    ddof = np.round((8/3)*len(x)/npseg)
    c = chi2.ppf([1-alfa/2, alfa/2], ddof)
    c = ddof/c
    CI_dw = Pxx*c[0]
    CI_up = Pxx*c[1]

    plt.figure()
    plt.plot(f, f*Pxx, color='k')
    plt.fill_between(f, f*CI_dw, f*CI_up, color='k', alpha=.5)
    plt.xscale('log')
    plt.xlabel('log(f)')
    plt.ylabel('f*PSD')
    plt.show()
    plt.close()

    PSD = [f, Pxx, CI_dw, CI_up]
    return PSD
