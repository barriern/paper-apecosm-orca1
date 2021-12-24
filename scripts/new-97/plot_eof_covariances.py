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
import matplotlib.ticker as mticker

latmax = 20
lonmin = 150
lonmax = -120

ilat = slice(None, -3)
# -

dnino, nino = read_index(filename='../data/index/oni.data')

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(y=ilat, z=0)
mesh

lonf = mesh['glamf']
latf = mesh['gphif']

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred' : 'l'})
const['l'] = const['length'].values * 100
const

wstep = const['weight_step']
wstep

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'
dirin = 'data/'

filename = '%s/full_eof_pacific_OOPE_latmax_%d_lonmin_%d_lonmax_%d.nc' %(dirin, latmax, lonmin, lonmax)
filename

eof = xr.open_dataset(filename)
eof = eof.rename({'w': 'l'})
eof['l'] = const['l']
eof

var = eof['eofvar']
var

pc = eof['eofpcs']
pc

cov = xr.open_dataset('%s/covariance_OOPE_anomalies_pcs_latmax_%d_lonmin_%d_lonmax_%d_all_sizes.nc' %(dirin, latmax, lonmin, lonmax)).isel(y=ilat)
cov = cov['__xarray_dataarray_variable__']
cov = cov.rename({'w': 'l'})
cov['l'] = const['l']
cov.name = 'cov'
cov

# +
letters = list(string.ascii_lowercase)
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = 110
lattext = 30

dictgrid = {'crs':ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
plt.rcParams['font.size'] = 15

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()

fig = plt.figure(figsize=(12, 8), facecolor='white')
axes_class = (GeoAxes, dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 2), axes_pad=(0.6, 0.4), label_mode='', cbar_mode='each', cbar_pad=0.05)
cbar_axes = axgr.cbar_axes

clim = {}
clim[3] = 70
clim[20] = 6
clim[90] = 3

cpt = 0
for l in [3, 20, 90]:
    for e in range(2):
        
        temp = cov.sel(l=l, method='nearest').isel(eof=e) * wstep.sel(l=l, method='nearest')
        temp = temp.where(temp != 0)
        temp = temp.to_masked_array().T
        vartemp = var.sel(l=l, method='nearest').isel(eof=e)
        pctemp = pc.sel(l=l, method='nearest').isel(eof=e).to_masked_array()

        perc = np.percentile(np.abs(np.ravel(temp[temp.mask == False])), 99)
        
        corr = np.corrcoef(pctemp, nino)[0, 1]
        if(corr < 0):
            temp *= -1
               
        perc = clim[l]
        
        ax = axgr[cpt]
        cs = ax.pcolormesh(lonf, latf, temp[1:, 1:], transform=proj2, cmap=plt.cm.RdBu_r)
        cb = plt.colorbar(cs, cbar_axes[cpt])
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND)
        cs.set_clim(-perc, perc)
        title = 'L=%.fcm, EOF %d (%.f' %(l, e + 1, vartemp) + '\%' + ')'
        ax.set_title(title)
        ax.set_extent([130, -60 + 360, -40, 40], crs=proj2)
        cb.set_label('J/m2')
        gl = ax.gridlines(**dictgrid)
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = False
        gl.left_labels = False
        if(cpt in [0, 2, 4]):
            gl.left_labels = True
        if(cpt > 3):
            gl.bottom_labels = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlocator = mticker.FixedLocator([150, 180, -150, -120, -90, -60])
        ax.text(lontext, lattext, letters[cpt] + ')', ha='center', va='center', bbox=dicttext)
        
        cpt += 1

plt.savefig('covariance_OOPE_pcs_latmax_%d_lonmin_%d_lonmax_%d_all_sizes.png' %(latmax, lonmin, lonmax), bbox_inches='tight')
# -
ax = plt.axes(projection=proj)
perc = np.percentile(np.abs(np.ravel(toto[toto.mask == False])), 99)
cs = ax.pcolormesh(lonf, latf, toto[1:, 1:], transform=proj2, cmap=plt.cm.RdBu_r)
cs.set_clim(-perc, perc)


