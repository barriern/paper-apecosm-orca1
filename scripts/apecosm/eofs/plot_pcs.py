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
data = data.isel(wpred=[14, 45, 80])
wstep = data['weight_step'].values
print(wstep)

dnino, nino = read_index(filename='data/index/oni.data')
ynino = dnino // 100
iok = np.nonzero((ynino <= 2018) & (ynino >= 1958))
nino = nino[iok]
dnino = dnino[iok]

mesh = xr.open_dataset("/Users/Nicolas/Work/sent/apecosm/ORCA1/mesh_mask_eORCA1_v2.2.nc")
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

data = xr.open_dataset("data/eof_oni_annual_density.nc")
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

            corrcoef = np.corrcoef(nino, pc[c, s, e, :])[0, 1]
            if(corrcoef < 0):
                eof[c, s, e, :, :] *= -1
                pc[c, s, e, : ] *= -1

            eof[c, s, e]  *= wstep[s]

            ax = axout[cpt - 1]
            ax.plot(time, pc[c, s, e, :], color='black')
            ax.fill_between(time, 0, nino, where=nino>0, interpolate=True, color='firebrick')
            ax.fill_between(time, 0, nino, where=nino<0, interpolate=True, color='steelblue')
            ax.set_xlim(time.min(), time.max())
            stride =  5 * 12
            ax.set_xticks(time[::stride])
            ax.set_xticklabels(date[::stride], rotation=45, ha='right')
            ax.set_title('Epi., %s, PC %d (%.2f' %(sizes[s], e + 1, var[c, s, e]) + '%)')
            ax.set_ylim(-3, 3)
            ax.grid(True, linestyle='--', linewidth=0.5)

            cpt += 1

plt.savefig('pc_full.png', bbox_inches='tight')
plt.close(fig)

