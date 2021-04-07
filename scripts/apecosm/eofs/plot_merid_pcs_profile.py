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
if 'w' in data.dims:
    data = data.isel(w=[14, 45, 80])
else:
    data = data.isel(wpred=[14, 45, 80])

wstep = data['weight_step'].values
print(wstep)

dnino, nino = read_index(filename='data/index/oni.data')
ynino = dnino // 100
iok = np.nonzero((ynino <= 2018) & (ynino >= 1958))
nino = nino[iok]
dnino = dnino[iok]

lonmax = -150

data = xr.open_dataset("data/eof_merid_profile_%.f.nc" %(lonmax))
pc = data['eofpc'].values
var = data['eofvar'].values * 100
time = data['time'].values
ntime = len(time)

years = np.array([t.year for t in time])
month = np.array([t.month for t in time])
time = np.arange(ntime)
date = ['%.4d-%.2d' %(y, m) for y,m in zip(years, month)]
time = np.arange(ntime)

ndn, ncom, nbins, neof = var.shape
coms = ['Epi.']

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = 120
lattext = 30
            

letters = list(string.ascii_lowercase)

classes = ['0-3cm', '3cm-20cm', '20cm-90cm', '90cm-200cm']
comm = ['Epi.', 'Mig.', 'Meso.']
dn = ['Day', 'Night']

cmin = 999
cmax = 0

for d in range(2):

    fig = plt.figure(figsize=(12, 8))
    axgr = ImageGrid(fig, 111,  nrows_ncols=(4, 2), label_mode='L', aspect=False, share_all=True, axes_pad=[0.5, 0.5])
    cbar_axes = axgr.cbar_axes
    axout = list(enumerate(axgr))
    axout = [p[1] for p in axout]

    print("@@@@@@@@@@@@@@@@@@@ d = ", d)

    cpt = 1

    for c in range(1):
        for s in range(nbins):
            for e in range(2):

                print("Comm = %d, Size = %d, EOF = %d" %(c, s, e + 1))

                corrcoef = np.corrcoef(nino, pc[d, c, s, e, :])[0, 1]
                print(corrcoef)
                if(corrcoef < 0):
                    pc[d, c, s, e, : ] *= -1

                if(e == 0):
                    cmin = min(cmin, abs(corrcoef))
                    cmax = max(cmax, abs(corrcoef))

                ax = axout[cpt - 1]
                ax.plot(time, pc[d, c, s, e, :], color='black')
                ax.fill_between(time, 0, nino, where=nino>0, interpolate=True, color='firebrick')
                ax.fill_between(time, 0, nino, where=nino<0, interpolate=True, color='steelblue')
                ax.set_xlim(time.min(), time.max())
                stride =  5 * 12
                ax.set_xticks(time[::stride])
                ax.set_xticklabels(date[::stride], rotation=45, ha='right')
                ax.set_title('Epi., %s, PC %d (%.2f' %(classes[s], e + 1, var[d, c, s, e]) + '\%)')
                ax.set_ylim(-3, 3)
                ax.grid(True, linestyle='--', linewidth=0.5)

                cpt += 1

    plt.savefig('pc_merid_profile_%d_%.f.png' %(d, lonmax), bbox_inches='tight')
    plt.close(fig)

print(cmin, cmax)
