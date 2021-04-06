import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import sys
sys.path.append('../../nino')
from extract_nino import read_index
import scipy.signal as sig
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import string
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.ticker as mticker
from cycler import cycler

data = xr.open_dataset('data/ORCA1_JRAC02_CORMSK_CYC1_FINAL_ConstantFields.nc')
data = data.isel(w=[14, 45, 80])
wstep = data['weight_step'].values
print(wstep)

dnino, nino = read_index(filename='data/index/oni.data')
ynino = dnino // 100
iok = np.nonzero((ynino <= 2018) & (ynino >= 1958))
nino = nino[iok]
dnino = dnino[iok]

data = xr.open_dataset("data/eof_oni_annual_density_20.nc")
eof = data['eofmap'].to_masked_array()
pc = data['eofpc'].values
var = data['eofvar'].values * 100
time = data['time'].values
ntime = len(time)

years = np.array([t.year for t in time])
month = np.array([t.month for t in time])
time = np.arange(ntime)
date = ['%.4d-%.2d' %(y, m) for y,m in zip(years, month)]
time = np.arange(ntime)

ncom, nbins, neof = var.shape
coms = ['Epi.']
sizes = ['3cm', '20cm', '90cm']

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = 120
lattext = 30
            
fig = plt.figure(figsize=(12, 8))
axgr = ImageGrid(fig, 111,  nrows_ncols=(3, 2), label_mode='L', aspect=False, share_all=True, axes_pad=[0.5, 0.5])
cbar_axes = axgr.cbar_axes
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

cpt = 1

letters = list(string.ascii_lowercase)

clim = [50, 5, 3]

for c in range(1):
    for s in range(nbins):
        for e in range(2):

            print("Comm = %d, Size = %d, EOF = %d" %(c, s, e + 1))

            corrcoef = np.corrcoef(pc[c, s, e, :], nino)[0, 1]

            if(corrcoef < 0):
                eof[c, s, e, :, :] *= -1
                pc[c, s, e, : ] *= -1

            eof[c, s, e]  *= wstep[s]
            corr = sig.correlate(pc[c, s, e, :], nino) / ntime
            lags = sig.correlation_lags(ntime, ntime)
            ilags = np.nonzero((lags >= 0) & (lags <= 3*12))[0]

            lags = lags[ilags]
            corr = corr[ilags]
            print(corr)

            test = np.abs(corr)
            iok = np.nonzero(test == test.max())[0][0]

            ax = axout[cpt - 1]
            ax.plot(lags, corr, color='black')
            xmin = 3 * 12
            ymin = 0.8
            print(lags[iok], corr[iok])
            ax.plot(lags[iok], corr[iok], marker='o', color='red', linestyle='none', zorder=1000, markersize=10)
            ax.set_title('Epi., %s, PC %d (%.2f' %(sizes[s], e + 1, var[c, s, e]) + '%)')
            ax.set_ylim(-ymin, ymin)
            ax.grid(True, linestyle='--', linewidth=0.5)
            ax.axvline(0, -1, 1, color='red')

            cpt += 1

plt.savefig('correlation_full.png', bbox_inches='tight')
plt.close(fig)

