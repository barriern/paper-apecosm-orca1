# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import string
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
sys.path.append('../nino')
from extract_nino import read_index
import numpy as np

latmax = 20
lonmin = 150
lonmax = -120
# -

pci = xr.open_dataset('../data/filt_tpi.nc')
pci = pci['tpi_filt'].to_masked_array()
pci = pci / pci.std()
ipci = np.nonzero(~np.ma.getmaskarray(pci))
ipci

dnino, nino = read_index(filename='../data/index/oni.data')

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred' : 'l'})
const['l'] = const['length'].values * 100
const

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'
dirin = 'data/'

eof = xr.open_dataset('%s/full_eof_pacific_OOPE_latmax_%d_lonmin_%d_lonmax_%d.nc' %(dirin, latmax, lonmin, lonmax))
eof = eof.rename({'w': 'l'})
eof['l'] = const['l']

var = eof['eofvar']
var

pc = eof['eofpcs']
pc

dates = ['%.4d-%.2d' %(d.year, d.month) for d in pc['time'].values]
dates

time = np.arange(nino.shape[0])

# +
letters = list(string.ascii_lowercase)
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = time[-1] - 50
lattext = 2.7

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
plt.rcParams['font.size'] = 15

fig = plt.figure(figsize=(12, 10), facecolor='white')
axgr = AxesGrid(fig, 111,  aspect=False, nrows_ncols=(3, 2), axes_pad=(0.6, 0.4), cbar_pad=0.05, share_all=True)
cbar_axes = axgr.cbar_axes

cpt = 0
for l in [3, 20, 90]:
    for e in range(2):
        
        pctemp = pc.sel(l=l, method='nearest').isel(eof=e)
        vartemp = var.sel(l=l, method='nearest').isel(eof=e)

        corr = np.corrcoef(pctemp, nino)[0, 1]
        if(corr < 0):
            pctemp *= -1
            corr *= -1
        print('L=%.fcm, EOF %d, corr_oni = %f' %(l, e + 1, corr))
        print('L=%.fcm, EOF %d, corr_pci = %f' %(l, e + 1, np.corrcoef(pci[ipci], pctemp[ipci])[0, 1]))
        pciroll = pctemp.rolling(time=7*12, center=True).mean()
            
        ax = axgr[cpt]
        ax.fill_between(time, 0, nino, where=(nino > 0), facecolor='firebrick', interpolate=True, label='ONI index')
        ax.fill_between(time, 0, nino, where=(nino < 0), facecolor='steelblue', interpolate=True)
        cs = ax.plot(time, pctemp, color='k', label='PC')
        ax.plot(time, pciroll, color='plum', linewidth=4, label='Rolling mean PC', linestyle='-')
        ax.plot(pci, color='Gold', linewidth=4, label='PCI index')
        title = 'L=%.fcm, EOF %d (%.f' %(l, e + 1, vartemp) + '\%' + ')'
        ax.set_title(title)
        ax.set_xlim(0, time.max())
        if(cpt == 0): 
            ax.legend(loc='upper left', ncol=2, fontsize=10)
        ax.grid()
        ax.text(lontext, lattext, letters[cpt] + ')', ha='center', va='center', bbox=dicttext)
        stride= 5 * 12
        ax.set_xticks(time[2*12::stride])
        ax.set_xticklabels(dates[2*12::stride], ha='right', rotation=45)
        
        cpt += 1
plt.savefig('pcs_latmax_%d_lonmin_%d_lonmax_%d.png' %(latmax, lonmin, lonmax), bbox_inches='tight')
# -


