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
plt.rcParams['font.size'] = 15

ilat = slice(None, -3)
# -

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

filename = '%s/covariance_OOPE_anomalies_oni_sizes.nc' %(dirin)
filename

eof = xr.open_dataset(filename).isel(y=ilat)
eof = eof.rename({'w': 'l'})
eof['l'] = const['l']
eof

cov = eof['OOPE']
cov

# +
letters = list(string.ascii_lowercase)
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = 110
lattext = 27

dictgrid = {'crs':ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
plt.rcParams['font.size'] = 15

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()

fig = plt.figure(figsize=(12, 8), facecolor='white')
axes_class = (GeoAxes, dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 1), axes_pad=(0.6, 0.4), label_mode='', cbar_mode='each', cbar_pad=0.05)
cbar_axes = axgr.cbar_axes

cpt = 0
for l in [3, 20, 90]:

    temp = cov.sel(l=l, lags=0, method='nearest') * wstep.sel(l=l, method='nearest')
    temp = temp.where(temp != 0)
    temp = temp.to_masked_array()
    perc = np.percentile(np.abs(np.ravel(temp[temp.mask == False])), 99)
 
    ax = axgr[cpt]
    cs = ax.pcolormesh(lonf, latf, temp[1:, 1:], transform=proj2, cmap=plt.cm.RdBu_r)
    cb = plt.colorbar(cs, cbar_axes[cpt])
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND)
    cs.set_clim(-perc, perc)
    title = 'L=%.fcm' %(l)
    ax.set_title(title)
    cb.set_label('J/m2')
    gl = ax.gridlines(**dictgrid)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.left_labels = True
    if(cpt == 2):
        gl.bottom_labels = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    #gl.xlocator = mticker.FixedLocator([150, 180, -150, -120, -90])
    ax.text(lontext, lattext, letters[cpt] + ')', bbox=dicttext, ha='center', va='center')
    #ax.set_extent([130, -60, -40, 40], crs=proj2)
    ax.set_extent([130, -60 + 360, -40, 40], crs=proj2)

    cpt += 1

plt.savefig('covariance_OOPE_oni.png', bbox_inches='tight')
# -

