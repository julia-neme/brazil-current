import numpy as np
import xarray
from matplotlib import pyplot as plt
from functions import get_coasts
from functions import climatologia_xarray

inp = input('File to process with path: ')
dat = xarray.open_dataset(inp).squeeze(drop=True)
dat = dat.astype('float64')

EC, WC = get_coasts(dat.lat,dat.lon)

# Valores diarios
A = np.empty([len(dat.time), len(dat.lat)])
for i in range(0,len(dat.lat)):
    A[:,i] = np.mean(dat.curl[:,i,int(WC[i]):int(EC[i])], axis=1)

LM = np.empty([len(dat.time), 2])
for i in range(0, len(dat.time)):
    xx = np.where(A==np.max(A[i,:]))[1][0]
    LM[i,0] = dat.lat[xx]
    LM[i,1] = A[i,xx]

# Valores mensuales

B = dat.curl.resample('1MS', dim='time', how='mean')
AA = np.empty([len(B[:,0,0]), len(dat.lat)])
for i in range(0,len(dat.lat)):
    AA[:,i] = np.mean(B[:,i,int(WC[i]):int(EC[i])], axis=1)

LMM = np.empty([len(AA[:,0]), 2])
for i in range(0, len(AA[:,0])):
    xx = np.where(AA==np.max(AA[i,:]))[1][0]
    LMM[i,0] = dat.lat[xx]
    LMM[i,1] = AA[i,xx]

# Climatologia

C = climatologia_xarray(dat.curl, dat.time)
AAA = np.empty([len(C[:,0,0]), len(dat.lat)])
for i in range(0,len(dat.lat)):
    AAA[:,i] = np.mean(C[:,i,int(WC[i]):int(EC[i])], axis=1)

LMC = np.empty([len(AAA[:,0]), 2])
for i in range(0, len(AAA[:,0])):
    xx = np.where(AAA==np.max(AAA[i,:]))[1][0]
    LMC[i,0] = dat.lat[xx]
    LMC[i,1] = AAA[i,xx]

# Ploteo

plt.figure(figsize=(10,7))
plt.plot(dat.time, LM[:,0], color='grey', label='Valores diarios')
plt.plot(dat.time, LMM[:,0], color='olive', label='Valores mensuales')
plt.grid(color='lighygrey', linestyle=':')
plt.show()

A = np.empty([len(dat.time), len(dat.lat)])
for i in range(0,len(dat.lat)):
    A[:,i] = np.mean(dat.curl[:,i,int(WC[i]):int(EC[i])], axis=1)

LM = np.empty([len(dat.time), 2])
for i in range(0, len(dat.time)):
    xx = np.where(A==np.max(A[i,:]))[1][0]
    LM[i,0] = dat.lat[xx]
    LM[i,1] = A[i,xx]

# Valores mensuales

B = dat.curl.resample('1MS', dim='time', how='mean')
AA = np.empty([len(B[:,0,0]), len(dat.lat)])
for i in range(0,len(dat.lat)):
    AA[:,i] = np.mean(B[:,i,int(WC[i]):int(EC[i])], axis=1)

LMM = np.empty([len(AA[:,0]), 2])
for i in range(0, len(AA[:,0])):
    xx = np.where(AA==np.max(AA[i,:]))[1][0]
    LMM[i,0] = dat.lat[xx]
    LMM[i,1] = AA[i,xx]

# Climatologia

C = climatologia_xarray(dat.curl)
AAA = np.empty([len(C[:,0,0]), len(dat.lat)])
for i in range(0,len(dat.lat)):
    AAA[:,i] = np.mean(C[:,i,int(WC[i]):int(EC[i])], axis=1)

LMC = np.empty([len(AAA[:,0]), 2])
for i in range(0, len(AAA[:,0])):
    xx = np.where(AAA==np.max(AAA[i,:]))[1][0]
    LMC[i,0] = dat.lat[xx]
    LMC[i,1] = AAA[i,xx]

# Ploteo

plt.figure(figsize=(16,5))
plt.plot(dat.time, LM[:,0], color='lightgrey', label='Valores diarios',
         linewidth=.5, linestyle='--')
plt.plot(B.time, LMM[:,0], color='olive', label='Valores mensuales',
         linewidth=2)
plt.margins(0,.2)
plt.grid(color='grey', linestyle=':')
plt.legend()
plt.title('Latitud del máximo del rotor promediado en la cuenca - NCEP I')
plt.show()
plt.close()

plt.figure(figsize=(10,5))
plt.plot(C.month, LMC[:,0], color='navy', linewidth=2)
plt.margins(0,.2)
plt.grid(color='grey', linestyle=':')
plt.xlabel('Mes')
plt.title('Climatología de la latitud del máximo del rotor promediado en la cuenca - NCEP I')
plt.show()
plt.close()
