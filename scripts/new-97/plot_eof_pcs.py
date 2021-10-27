# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python [conda env:nbarrier] *
#     language: python
#     name: conda-env-nbarrier-py
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

eof = xr.open_dataset('%s/full_eof_pacific_OOPE_latmax_%d_lonmin_%d_lonmax_%d.nc' %(dirin, latmax, lonmin, lonmax))
eof = eof.rename({'w': 'l'})
eof['l'] = const['l']

var = eof['eofvar']
var

pc = eof['eofpcs']
pc

letters = list(string.ascii_lowercase)
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = 120
lattext = 30

pc.shape

nino.shape

time = np.arange(nino.shape[0])

# +
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

fig = plt.figure(figsize=(12, 8))
axgr = AxesGrid(fig, 111,  aspect=False, nrows_ncols=(3, 2), axes_pad=(1.23, 0.4), cbar_pad=0.05)
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
        ax.fill_between(time, 0, nino, where=(nino > 0), facecolor='firebrick', interpolate=True)
        ax.fill_between(time, 0, nino, where=(nino < 0), facecolor='steelblue', interpolate=True)
        cs = ax.plot(time, pctemp, color='k')
        ax.plot(time, pciroll, color='cyan', linewidth=2)
        ax.plot(pci, color='Gold', linewidth=2)
        title = 'L=%.fcm, EOF %d (%.f' %(l, e + 1, vartemp) + '%' + ')'
        ax.set_title(title)
        ax.set_xlim(0, time.max())
        cpt += 1
plt.savefig('pcs_latmax_%d_lonmin_%d_lonmax_%d.png' %(latmax, lonmin, lonmax))
# -


