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

wpred = [14, 45, 80]

tpidata = xr.open_dataset("../../nino/filt_tpi.nc")
tpi = tpidata['tpi_filt'].to_masked_array()
#tpi = tpidata['tpi_raw'].to_masked_array()
tpi = (tpi - np.mean(tpi)) / np.std(tpi)

data = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
data = data.isel(wpred=wpred)
wstep = data['weight_step'].values
print(wstep)

dnino, nino = read_index(filename='../../data/index/oni.data')
ynino = dnino // 100
iok = np.nonzero((ynino <= 2018) & (ynino >= 1958))
nino = nino[iok]
dnino = dnino[iok]

for lll in range(20, 60, 10):

    data = xr.open_dataset("data/eof_full_density_%d.nc" %lll)
    data = data.isel(bins=wpred)
    eof = data['eofmap'].to_masked_array()
    pc = data['eofpc'].values
    var = data['eofvar'].values * 100
    time = data['time'].values
    ntime = len(time)

    years = np.array([t.year for t in time])
    month = np.array([t.month for t in time])
    time = np.arange(ntime)
    date = ['%.4d-%.2d' %(y, m) for y, m in zip(years, month)]
    time = np.arange(ntime)

    nbins, neof = var.shape
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

    for s in range(nbins):
        for e in range(2):

            print("Comm = %d, Size = %d, EOF = %d" %(0, s, e + 1))

            corrcoef = np.corrcoef(nino, pc[s, e, :])[0, 1]
            if(corrcoef < 0):
                eof[s, e, :, :] *= -1
                pc[s, e, : ] *= -1

            print('Lalala')
            print(eof.shape)
            print(eof[s, e].shape)
            print(wstep.shape)
            eof[s, e]  *= wstep[s]

            ax = axout[cpt - 1]
            ax.plot(time, pc[s, e, :], color='black')
            ax.plot(time, tpi, color='gold')
            ax.fill_between(time, 0, nino, where=nino>0, interpolate=True, color='firebrick')
            ax.fill_between(time, 0, nino, where=nino<0, interpolate=True, color='steelblue')
            ax.set_xlim(time.min(), time.max())
            stride =  5 * 12
            ax.set_xticks(time[::stride])
            ax.set_xticklabels(date[::stride], rotation=45, ha='right')
            ax.set_title('Epi., %s, PC %d (%.2f' %(sizes[s], e + 1, var[s, e]) + '\%)')
            ax.set_ylim(-3, 3)
            ax.grid(True, linestyle='--', linewidth=0.5)

            cpt += 1

    plt.savefig('fig7_%d.png' %lll, bbox_inches='tight')
    plt.close(fig)
