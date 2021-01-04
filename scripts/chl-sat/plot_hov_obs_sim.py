import sys
sys.path.append('../nino')
from extract_nino import read_index

daten, nino = read_index(filename='/home/barrier/Work/apecosm/ORCA1/DATA_APE_ORCA1/apecosm/index/oni.data')

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

plt.figure()
plt.subplots_adjust(left=0.09, top = 0.95, right=0.98, bottom=0.2, wspace=0.1)

# processing simulated data

data = xr.open_dataset('data/anom_chl_simulated.nc')
date = data['time_counter.year'] * 100 + data['time_counter.month']
date = date.values
time = np.arange(len(date))
lon = data['x'].values

year = date//100
month = date - year *100
iticks = np.nonzero(month == 1)[0]
iticks = iticks[::2]

datestr = ['%.4d-%.2d' %(y, m) for y, m in zip(year, month)]
datestr = np.array(datestr)

inino = np.nonzero((daten <= date.max()) & (daten >= date.min()))[0]
daten = daten[inino]
nino = nino[inino]

ccc = 0.2

ax1 = plt.subplot(131)
cs = plt.pcolormesh(lon, time, data['chl'])
cs.set_clim(-ccc, ccc)
plt.xlabel('Lon')
ax1.set_yticks(time[iticks])
ax1.set_yticklabels(datestr[iticks], rotation=45, va='top')
plt.grid(True)
plt.title('Sim.')

# processing obs data

lonmod = lon

data = xr.open_dataset('data/anom_chl_obs.nc')
date = data['time.year'] * 100 + data['time.month']
date = date.values
time = np.arange(len(date))
lon = data['lon'].values
chl_anom = data['chl_anom'].to_masked_array()

delta = lon[1] - lon[0]
nint = int(np.round(1. / delta))

lonint = []
chlint = []
off = 0
iii = np.arange(nint) + off
iii = iii.astype(np.int)

while(iii[-1] < len(lon)):
    lontemp = np.mean(lon[iii])
    chltemp = np.mean(chl_anom[:, iii], axis=1)
    lonint.append(lontemp)
    chlint.append(chltemp)
    iii += nint

lonint = np.array(lonint)
chlint = np.array(chlint).T

lon = lonint
chl = chlint

ax2 = plt.subplot(132, sharey=ax1)
cs = plt.pcolormesh(lon, time, chlint)
plt.xlabel('Lon')
cs.set_clim(-ccc, ccc)
plt.title('Sat.')

cax = plt.axes([0.17, 0.09, 0.4, 0.02])
cb = plt.colorbar(cs, cax, orientation='horizontal')
cb.set_label('Chl anom (mg/m3)')
plt.grid(True)

# processing nino index

ax3 = plt.subplot(133, sharey=ax1)
plt.plot(nino, time, lw=1)
plt.fill_betweenx(time, nino, 0, where=(nino>0), color='firebrick', interpolate=True)
plt.fill_betweenx(time, nino, 0, where=(nino<0), color='steelblue', interpolate=True)
ax3.set_ylim(time.min(), time.max())
ax3.set_xlim(-3, 3)
plt.title('ONI index')
plt.grid(True)

ax3.tick_params(axis='y', labelleft=False)
ax2.tick_params(axis='y', labelleft=False)

plt.savefig('hov_chl.pdf', bbox_inches='tight')

