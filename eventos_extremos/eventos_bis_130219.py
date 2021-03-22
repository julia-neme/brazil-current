import xarray
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from scipy import signal
from scipy import stats
import matplotlib.ticker as mticker
from scipy.stats import chi2

osc = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_clim_atlsur_2009_2015.nc')
osca = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_anom_atlsur_2009_2015.nc')
gha = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrssta_atlsur_2009_2015.nc')
gh = xarray.open_dataset('/home/bock/Documents/tesis/datos/ghrsst_clim_atlsur_2009_2015.nc')
oscl = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')

clim_vl = oscl['vm'].groupby('time.month').mean('time')

# Ciclo anual de V ############################################################

clim_v = osc['v_clim']
clim_sst = gh['sst_clim']

A = clim_v.sel(lat=-34.5, method='nearest').sel(lon=slice(-52, -48)).sel(lon=slice(-51.66, -50.33)).mean(dim='lon')
B = clim_sst.sel(lon=slice(-51.67,-50.33), lat=-34.5).mean(dim='lon')
C = clim_vl.sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude')
mth = np.arange(0, 12, 1)

f, ax = plt.subplots(2, sharex=True, figsize=(10,5))
ax[0].set_xlim([0, 11])

ax[0].plot(mth, A*100, color='b')
ax[1].plot(mth, B-273, color='r')
ax[0].grid(which="major", color='grey', linestyle=':')
ax[1].grid(which="major", color='grey', linestyle=':')
ax[0].set_ylabel('V [cm/s]')
ax[1].set_ylabel('SST [$^{\circ}$ C]')
ax[0].set_xticks([1,2,3,4,5,6,7,8,9,10])
ax[0].set_xticklabels([ 'Feb.', 'Mar.', 'Abr.', 'May.', 'Jun.', 'Jul.', 'Ago.', 'Sep.', 'Oct.', 'Nov.'], size=15)
f.subplots_adjust(hspace=0.05)
ax[0].text(-0.1, 0.9, '(a)', transform=ax[0].transAxes, size=15)
ax[1].text(-0.1, 0.9, '(b)', transform=ax[1].transAxes, size=15)

plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/ciclo_climatologico_vsst.png', bbox_inches='tight')
plt.show()

# Calculo eventos sin el ciclo anual ##########################################

v = osca['v_anom'].sel(lat=-34.5, method='nearest').sel(lon=slice(-52, -48)).sel(lon=slice(-51.66, -50.33)).mean(dim='lon').squeeze()
s1, s2 = signal.butter(6, 1/5, 'low')
vf = signal.filtfilt(s1, s2, v)

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

eddof = get_eddof(vf)
sem = np.nanstd(vf)/np.sqrt(eddof-1)

plt.figure(figsize=(10,3))
plt.xlim([osca['time'][0].values, osca['time'][-1].values])
plt.plot(osca['time'], np.zeros(len(v)), color='k', linestyle=':', linewidth=.9)
plt.plot(osca['time'], vf*100, color='b', linewidth=1.2)
plt.plot(osca['time'], np.mean(vf)*np.ones(len(vf))*100, color='k', linewidth=2)
plt.plot(osca['time'], (np.mean(vf)+np.std(vf))*np.ones(len(vf))*100, color='k', linestyle='--', linewidth=1.7)
plt.plot(osca['time'], (np.mean(vf)-np.std(vf))*np.ones(len(vf))*100, color='k', linestyle='--', linewidth=1.7)
plt.ylabel('AVM [cm/s]')
plt.gca().xaxis.set_minor_locator(mdates.MonthLocator())
plt.grid(which="major", color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/resultados/figs/eestd.png', bbox_inches='tight')
plt.show()

idx_ED = np.where(vf>np.mean(vf)+np.std(vf))[0]
kk = np.diff(idx_ED)
ED = []; n=0; nn=0
for i in range(0, len(kk)):
    if kk[i]==1:
        nn = nn+1
    else:
        d = (nn-n+1)*5
        f_i = osca['time'][idx_ED[n]].values
        v_mean = np.mean(vf[idx_ED[n]:idx_ED[nn]+1])
        v_max = np.max(vf[idx_ED[n]:idx_ED[nn]+1])
        n = nn+1
        nn = nn+1
        ED.append([f_i, d, v_mean, v_max])
del kk, f_i, d, v_mean, v_max
ED.append([osca['time'][idx_ED[n]].values,
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
        f_i = osca['time'][idx_EF[n]].values
        v_mean = np.mean(vf[idx_EF[n]:idx_EF[nn]+1])
        v_max = np.min(vf[idx_EF[n]:idx_EF[nn]+1])
        n = nn+1
        nn = nn+1
        EF.append([f_i, d, v_mean, v_max])
del kk, f_i, d, v_mean, v_max
EF.append([osca['time'][idx_EF[n]].values,
           (idx_EF[-1]-idx_EF[n]+1)*5,
           np.mean(vf[idx_EF[n]:idx_EF[-1]+1]),
           np.min(vf[idx_EF[n]:idx_EF[-1]+1])])
EF = np.array(EF)
k = np.where(EF[:,1]>=20)
EFp = EF[k]; del EF
del idx_ED, idx_EF

idxED = []
for i in range(0, len(EDp[:,0])):
    t = np.where(osca['time']==EDp[i,0])[0][0]
    idxED.append(t); idxED.append(int(EDp[i,1]/5+t))
idxEF = []
for i in range(0, len(EFp[:,0])):
    t = np.where(osca['time']==EFp[i,0])[0][0]
    idxEF.append(t); idxEF.append(int(EFp[i,1]/5+t))

f, ax = plt.subplots(1, sharex=True, figsize=(10,3))
ax.set_xlim([osca['time'][0].values, osca['time'][-1].values])
ax.plot(osca['time'], vf*100, color='k', linewidth=1)
ax.plot(osca['time'], np.mean(vf)*np.ones(len(vf))*100, color='k', linewidth=2)
ax.plot(osca['time'], (np.mean(vf)+np.std(vf))*np.ones(len(vf))*100, color='r', linestyle='--', linewidth=2)
ax.plot(osca['time'], (np.mean(vf)-np.std(vf))*np.ones(len(vf))*100, color='b', linestyle='--', linewidth=2)
for i in range(0, int(len(idxED)/2)):
    ax.fill_between(osca['time'][idxED[2*i]:idxED[2*i+1]].values, vf[idxED[2*i]:idxED[2*i+1]]*100,
                       np.mean(vf)*100+np.std(vf)*100, color='r', alpha=.5)
for i in range(0, int(len(idxEF)/2)):
    ax.fill_between(osca['time'][idxEF[2*i]:idxEF[2*i+1]].values, vf[idxEF[2*i]:idxEF[2*i+1]]*100,
                       np.mean(vf)*100-np.std(vf)*100, color='b', alpha=.5)
ax.set_ylabel('$V_A$ [cm/s]')
ax.xaxis.set_minor_locator(mdates.MonthLocator())
ax.grid(which="major", color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/identificacion_eventos.png', bbox_inches='tight')
plt.show()

# SSTA ########################################################################

ssta = gha['sst_anom'].sel(lon=slice(-51.67,-50.33), lat=-34.5).mean(dim='lon')
s1, s2 = signal.butter(6, 1/30, 'low')
sstaf = signal.filtfilt(s1, s2, ssta)

ED = np.array([18, 24, 70, 75, 93, 101, 172, 177, 185, 190, 196, 201, 326, 331, 366,
               375, 417, 423, 461, 467])
EF = np.array([63, 67, 104, 112, 272, 278, 318, 325, 359, 363, 378, 391, 426, 433,
      439, 444, 470, 475, 489, 497])

f, ax = plt.subplots(1, sharex=True, figsize=(10,3))
ax.set_xlim([gha['time'][0].values, gha['time'][-1].values])
ax.plot(gha['time'], sstaf, color='k', linewidth=1)
ax.plot(gha['time'], np.mean(sstaf)*np.ones(len(sstaf)), color='k', linewidth=2)
ax.plot(gha['time'], (np.mean(sstaf)+np.std(sstaf))*np.ones(len(sstaf)), color='k', linestyle='--', linewidth=2)
ax.plot(gha['time'], (np.mean(sstaf)-np.std(sstaf))*np.ones(len(sstaf)), color='k', linestyle='--', linewidth=2)
for i in range(0, int(len(ED)/2)):
    kk = osca.isel(time=ED[2*i]+1)
    kk1 = osca.isel(time=ED[2*i+1]+1)
    kk2 = ssta.sel(time=slice(kk['time'].values, kk1['time'].values))
    k1 = np.where(ssta['time']==kk2.time[0])[0][0]
    k2 = np.where(ssta['time']==kk2.time[-1])[0][0]
    ax.plot(kk2['time'], sstaf[k1:k2+1], color='r', linewidth=2.5)
for i in range(0, int(len(EF)/2)):
    kk = osca.isel(time=EF[2*i]+1)
    kk1 = osca.isel(time=EF[2*i+1]+1)
    kk2 = ssta.sel(time=slice(kk['time'].values, kk1['time'].values))
    k1 = np.where(ssta['time']==kk2.time[0])[0][0]
    k2 = np.where(ssta['time']==kk2.time[-1])[0][0]
    ax.plot(kk2['time'], sstaf[k1:k2+1], color='b', linewidth=2.5)
ax.set_ylabel('SSTA [$^{\circ}$ C]')
ax.xaxis.set_minor_locator(mdates.MonthLocator())
ax.grid(which="major", color='grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/ssta_eventos.png', bbox_inches='tight')
plt.show()

# Correlacion V y SSTA ########################################################

ssta_o = ssta.sel(time=osca['time'], method='nearest')
ise = []
for i in range(0, len(ssta_o.time)):
    ise.append(np.where(ssta.time==ssta_o.time[i])[0][0])
ise = np.array(ise)

R_lags = [np.corrcoef(vf, sstaf[ise])[0][1]]
N = 60
CF = [1.96/np.sqrt(N*3-3)]
for i in range(1, N-1):
    R_lags.append(np.corrcoef(vf[0:-i], sstaf[ise[:-i]+i])[0][1])
    CF.append(1.96/np.sqrt(N*3-i-3))
CFm = -np.ones(len(CF))*CF

f, ax = plt.subplots(1, sharex=True, figsize=(10,3))
ax.bar(np.arange(0, len(R_lags),1), R_lags, width=.3, color='k')
ax.plot(CF, color='k', linestyle='--')
ax.plot(CFm, color='k', linestyle='--')
ax.set_xlabel('Lags [días]')
ax.set_ylabel('Coeficiente de correlación')
ax.grid('grey', linestyle=':')
plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/corr_ssta_v.png', bbox_inches='tight')
plt.show()

# Eventos serie larga oscar ###################################################

osc1 = xarray.open_dataset('/home/bock/Documents/tesis/datos/oscar_atlsur_1992_2015.nc')

serie = osc1['v'].sel(latitude=-34.5, method='nearest').sel(longitude=slice(-52, -48)).sel(longitude=slice(-51.66, -50.33)).mean(dim='longitude').squeeze()
serie_a = serie.groupby('time.month') - serie.groupby('time.month').mean('time')

s1, s2 = signal.butter(6, 1/5, 'low')
serie_f = signal.filtfilt(s1, s2, serie_a.values)

idx_ED = np.where(serie_f>np.mean(serie_f)+np.std(serie_f))[0]
kk = np.diff(idx_ED)
ED = []; n=0; nn=0
for i in range(0, len(kk)):
    if kk[i]==1:
        nn = nn+1
    else:
        d = (nn-n+1)*5
        f_i = osc1['time'][idx_ED[n]].values
        v_mean = np.mean(serie_f[idx_ED[n]:idx_ED[nn]+1])
        v_max = np.max(serie_f[idx_ED[n]:idx_ED[nn]+1])
        n = nn+1
        nn = nn+1
        ED.append([f_i, d, v_mean, v_max])
del kk, f_i, d, v_mean, v_max
ED.append([osc1['time'][idx_ED[n]].values,
           (idx_ED[-1]-idx_ED[n]+1)*5,
           np.mean(serie_f[idx_ED[n]:idx_ED[-1]+1]),
           np.max(serie_f[idx_ED[n]:idx_ED[-1]+1])])
ED = np.array(ED)
k = np.where(ED[:,1]>=20)
EDp = ED[k]; del ED

idx_EF = np.where(serie_f<np.mean(serie_f)-np.std(serie_f))[0]
kk = np.diff(idx_EF)
EF = []; n=0; nn=0
for i in range(0, len(kk)):
    if kk[i]==1:
        nn = nn+1
    else:
        d = (nn-n+1)*5
        f_i = osc1['time'][idx_EF[n]].values
        v_mean = np.mean(serie_f[idx_EF[n]:idx_EF[nn]+1])
        v_max = np.min(serie_f[idx_EF[n]:idx_EF[nn]+1])
        n = nn+1
        nn = nn+1
        EF.append([f_i, d, v_mean, v_max])
del kk, f_i, d, v_mean, v_max
EF.append([osc1['time'][idx_EF[n]].values,
           (idx_EF[-1]-idx_EF[n]+1)*5,
           np.mean(serie_f[idx_EF[n]:idx_EF[-1]+1]),
           np.min(serie_f[idx_EF[n]:idx_EF[-1]+1])])
EF = np.array(EF)
k = np.where(EF[:,1]>=20)
EFp = EF[k]; del EF
del idx_ED, idx_EF

idxED = []
for i in range(0, len(EDp[:,0])):
    t = np.where(osc1['time']==EDp[i,0])[0][0]
    idxED.append(t); idxED.append(int(EDp[i,1]/5+t))
idxEF = []
for i in range(0, len(EFp[:,0])):
    t = np.where(osc1['time']==EFp[i,0])[0][0]
    idxEF.append(t); idxEF.append(int(EFp[i,1]/5+t))


f, ax = plt.subplots(1, sharex=True, figsize=(10,3))
ax.set_xlim([osc1['time'][0].values, osc1['time'][-1].values])
ax.plot(osc1['time'], serie_f*100, color='k', linewidth=1)
ax.plot(osc1['time'], np.mean(serie_f)*np.ones(len(serie_f))*100, color='k', linewidth=2)
ax.plot(osc1['time'], (np.mean(serie_f)+np.std(serie_f))*np.ones(len(serie_f))*100, color='r', linestyle='--', linewidth=2)
ax.plot(osc1['time'], (np.mean(serie_f)-np.std(serie_f))*np.ones(len(serie_f))*100, color='b', linestyle='--', linewidth=2)
for i in range(0, int(len(idxED)/2)):
    ax.fill_between(osc1['time'][idxED[2*i]:idxED[2*i+1]].values, serie_f[idxED[2*i]:idxED[2*i+1]]*100,
                       np.mean(serie_f)*100+np.std(serie_f)*100, color='r', alpha=.5)
for i in range(0, int(len(idxEF)/2)):
    ax.fill_between(osc1['time'][idxEF[2*i]:idxEF[2*i+1]].values, serie_f[idxEF[2*i]:idxEF[2*i+1]]*100,
                       np.mean(serie_f)*100-np.std(serie_f)*100, color='b', alpha=.5)
ax.set_ylabel('AVM [cm/s]')
ax.xaxis.set_minor_locator(mdates.YearLocator())
ax.grid(which="major", color='grey', linestyle=':')
#plt.savefig('/home/bock/Documents/tesis/documentos/template/figuras/identificacion_eventos.png', bbox_inches='tight')
plt.show()
