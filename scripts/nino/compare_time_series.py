import xarray as xr
import sys
sys.path.append('../nino')
from extract_nino import read_index
import numpy as np
import matplotlib.pyplot as plt
import apecosm.ts as ts
import scipy.signal as sig

dnino, nino = read_index('../data/index/oni.data')
iok = np.nonzero(np.isnan(nino) == False)[0]
dnino = dnino[iok]
nino = nino[iok]
time = np.arange(len(dnino))
ynino = dnino // 100
mnino = dnino - 100 * ynino

labels = np.array(['%.4d-%.2d' %(y, m) for y,m in zip(ynino, mnino)])

data = xr.open_dataset("data/simulated_enso_index.nc")
years = data['time.year'].values
months = data['time.month'].values
timemod = np.arange(len(years)) + 8 * 12
nmod = len(timemod)
enso = data['enso'].values
ensof = np.zeros(enso.shape)
clim, enso = ts.get_monthly_clim(enso)
enso = sig.detrend(enso)

index = np.arange(3)
for i in range(1, nmod - 1):
    ensof[i] = enso[index].mean()
    index += 1

ensof[ensof == 0] = np.nan

istart = np.nonzero(dnino == 195802)[0][0]
iend = np.nonzero(dnino == 201811)[0][0] + 1

test = np.corrcoef(ensof[1:-1], nino[istart:iend])

istart = np.nonzero(dnino == 195801)[0][0]
iend = np.nonzero(dnino == 201812)[0][0] + 1
test = np.corrcoef(enso, nino[istart:iend])
print(test[0, 1])

xticks = np.arange(0, len(time), 3 * 12)
print(xticks)
print(dnino[8*12:])

plt.figure(figsize=(12, 8))
ax = plt.gca()
plt.fill_between(time, 0, nino, where=(nino>0), color='firebrick', interpolate=True)
plt.fill_between(time, 0, nino, where=(nino<0), color='steelblue', interpolate=True)
plt.plot(timemod, ensof, 'k', label='Sim.')
#plt.legend(loc=0)
ax.set_xticks(time[xticks])
ax.set_xticklabels(labels[xticks], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(time.min(), time.max())
ax.set_ylim(-4, 4)
plt.savefig('compare_timeseries.png', bbox_inches='tight')

