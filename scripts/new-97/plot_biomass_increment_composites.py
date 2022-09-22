# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python [conda env:nbarrier2]
#     language: python
#     name: conda-env-nbarrier2-py
# ---

# +
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.pyplot as plt
import numpy as np

ilat = slice(None, -3)

projout = ccrs.PlateCarree(central_longitude=180)
projin = ccrs.PlateCarree(central_longitude=0)

dictpbox = {'transform': projin, 'linestyle': '--', 'linewidth': 1, 'color':'k'}
gridparams = {'crs': ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}
axes_class = (GeoAxes, dict(map_projection=projout))
# -

const = xr.open_dataset('data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'w'})
length = const['length'] * 100
length = length.rename({'w': 'l'})

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(y=ilat, z=0)
mesh

lonf = mesh['glamf'].values
latf = mesh['gphif'].values
lonf.shape

seasons = ['JFM', 'AMJ', 'JAS', 'OND']
years = [0, 1, 2]
labels = ['%s, year %d' %(s, y) for y in years for s in seasons ]
labels

clim = [400000, 300, 3]

for varname in ['zadv_trend', 'madv_trend', 'zdiff_trend', 'mdiff_trend', 'predationTrend', 'growthTrend']:

    print('+++++++++++++++++ Processing variable ', varname)
    data = xr.open_dataset('data/seasonal_composite_trends_%s.nc' %varname).isel(y=ilat)[varname]
    data['l'] = length
    data

    cpt = 0
    for l in [3, 20, 90]:
        
        print('Processing l ', l, 'cm')
    
        datatemp = data.sel(l=l, method='nearest').values
        ccc = np.percentile(np.ravel(np.abs(datatemp[~np.isnan(datatemp)])), 95)
        ccc = clim[cpt]
        cpt += 1

        fig = plt.figure(figsize=(10, 10))
        cmap = plt.cm.get_cmap('RdBu_r')

        axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(4, 3), axes_pad=(0.1, 0.45), 
                        label_mode='', cbar_mode='single', cbar_size=0.1, cbar_pad=0.01, cbar_location="bottom")

        for i in data['season'].values:

            toplot = datatemp[i, :, :]
            cs = axgr[i].pcolormesh(lonf, latf, toplot[1:, 1:], transform=projin, cmap=cmap) 
            cbax = axgr.cbar_axes[i]
            cb = cbax.colorbar(cs)
            cs.set_clim(-ccc, ccc)
            axgr[i].set_title(labels[i])
            axgr[i].add_feature(cfeature.LAND)
            axgr[i].add_feature(cfeature.COASTLINE)

        plt.suptitle('%s, L = %d cm' %(varname, l), y=0.87)
        figname = 'seasonal_compo_%s_l_%dcm.png' %(varname, l)
        plt.savefig(figname, bbox_inches='tight')
        plt.close(fig)

