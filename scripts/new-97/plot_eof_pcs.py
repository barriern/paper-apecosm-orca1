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

latmax = 10
lonmin = 150
lonmax = -100
# -

dnino, nino = read_index(filename='../data/index/oni.data')

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred' : 'l'})
const['l'] = const['length'].values * 100
const

wstep = const['weight_step']
wstep

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'

eof = xr.open_dataset('%s/full_eof_pacific_OOPE_latmax_10_lonmin_150_lonmax_-100.nc' %dirin)
eof = eof.rename({'w': 'l'})
eof['l'] = const['l']

var = eof['eofvar']
var

pc = eof['eofpcs']
pc

cov = xr.open_dataset('%s/covariance_OOPE_anomalies_pcs_latmax_10_lonmin_150_lonmax_-100_all_sizes.nc' %dirin).isel(y=ilat)
cov = cov['__xarray_dataarray_variable__']
cov = cov.rename({'w': 'l'})
cov['l'] = const['l']
cov.name = 'cov'
cov

letters = list(string.ascii_lowercase)
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = 120
lattext = 30

# +
dictgrid = {'crs':ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()

fig = plt.figure(figsize=(12, 8))
axes_class = (GeoAxes, dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 2), axes_pad=(1.23, 0.4), label_mode='', cbar_mode='each', cbar_pad=0.05)
cbar_axes = axgr.cbar_axes

cpt = 0
for l in [3, 20, 90]:
    for e in range(2):
        
        temp = cov.sel(l=l, method='nearest').isel(eof=e) * wstep.sel(l=l, method='nearest')
        temp = temp.to_masked_array()
        vartemp = var.sel(l=l, method='nearest').isel(eof=e)
        pctemp = pc.sel(l=l, method='nearest').isel(eof=e).to_masked_array()

        perc = np.percentile(np.abs(np.ravel(temp[temp.mask == False])), 99)
        
        corr = np.corrcoef(pctemp, nino)[0, 1]
        if(corr < 0):
            temp *= -1
            
        ax = axgr[cpt]
        cs = ax.pcolormesh(lonf, latf, temp[1:, 1:], transform=proj2, cmap=plt.cm.RdBu_r)
        plt.colorbar(cs, cbar_axes[cpt])
        ax.add_feature(cfeature.COASTLINE, zorder=1001)
        ax.add_feature(cfeature.LAND, zorder=1000)
        cs.set_clim(-perc, perc)
        title = 'L=%.fcm, EOF %d (%.f' %(l, e + 1, vartemp) + '%' + ')'
        ax.set_title(title)
        ax.set_xlim(-60, 130)
        ax.set_ylim(-50, 50)
        
        cpt += 1
        
lonf.shape
