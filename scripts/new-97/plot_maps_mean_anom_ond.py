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

# # Plotting maps for mean OOPE and Nino anoms

# ## Imports

# +
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import string

yslice = slice(None, -3)
# -

# ## Loading constant fields

const = xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const

# The `l` coordinate is assigned with length in cm.

const.coords['l'] = const['length'] * 100

# Now we extract weight step and length

wstep = const['weight_step']
wstep

length = const['length'] * 100
length

# ## Loading the covariances

filename = 'data/covariance_OOPE_anomalies_oni_sizes.nc'
filename

cov = xr.open_dataset(filename).isel(y=yslice)
cov = cov.rename({'w': 'l'})
cov['l'] = const['l']
cov = cov['OOPE'].sel(lags=0) * wstep
cov = cov.sel(l=[3, 20, 90], method='nearest')
cov

# ## Loading the climatology

clim = xr.open_dataset('data/pacific_clim_OOPE.nc').rename({'w': 'l'}).isel(y=yslice)
clim

# We assign the `l` coordinate the values of length in cm.

clim.coords['l'] = length

# Now, the mean is computed over the 12 months, to have one map per size class.

varclim = (clim['OOPE'] * wstep).mean(dim='month').sel(l=[3, 20, 90], method='nearest')
varclim.name = 'OOPE_clim'
varclim

# ## Extracting the anomalies
#
# Now we extract the anomalies for the OND months of 97 (peak of El Nino). First, we read the raw OOPE data:

anom = xr.open_dataset('data/composite_OOPE_map.nc').isel(y=yslice)
anom

# Again, we rename the `w` dimension and change the length coordinate into `cm`

anom = anom.rename({'w': 'l'})
anom.coords['l'] = length

# Then, we extract the proper sizes:

anom = anom.sel(l=[3, 20, 90], method='nearest')
anom

# Now, we move from `J/m2/kg` into `J/m2` for the raw fields.

varanom = anom['OOPE'] * wstep
varanom.name = 'OOPE'
varanom

# ## Loading the mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0, y=yslice)
mesh

lonf = mesh['glamf'].values
lonf

latf = mesh['gphif'].values
latf

# ## Drawing the map
#
# First, define some variables

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

tmask_noatl = xr.open_dataset('data/corrected_mask_noatl.nc')
tmask_noatl = tmask_noatl['tmask'].values
tmask_noatl = (tmask_noatl == tmask_noatl.min()).astype(int)
plt.imshow(tmask_noatl)

# +
lont = mesh['glamt'].values
latt = mesh['gphit'].values
imask = np.nonzero(tmask_noatl == 1)
lon1d = np.ravel(lont[imask])
lat1d = np.ravel(latt[imask])

output = proj.transform_points(proj2, lon1d, lat1d)
lonout = output[..., 0]
latout = output[..., 1]
# -

gridparams = {'crs': ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}

lontext = 120
lattext = 30
lontext2 = 110
lattext2 = -30
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)


def manage_axes(gl, cpt):
    
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.left_labels = False
    if cpt in [0, 2, 4]:
        gl.left_labels = True
    if cpt in [4, 5]:
        gl.bottom_labels = True


# +
plt.rcParams['font.size'] = 15
textprops = {}
textprops.update(ha='center', va='center', transform=proj, bbox=dicttext, zorder=2000)
plt.rcParams['contour.negative_linestyle'] = 'solid'

step = 0.5
levels = np.arange(-1, 3 + step, step)
step = 1
levels = np.arange(-0.75, -0.75 + 6 * 0.75, 0.75)

fig = plt.figure(figsize=(16, 16))
plt.subplots_adjust(top=0.95)
axes_class = (GeoAxes, dict(map_projection=proj))

letters = list(string.ascii_lowercase)

axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 2), axes_pad=(0.5, 0.5), label_mode='', cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="bottom")

cptlet = 0

# ## Plotting the mean
# for l in range(3):
#     cpt = (l * 3)
#     ax = axgr[cpt]
    
#     temp = varclim.isel(l=l)
#     ltemp = float(temp.l)
#     temp = np.log10(temp)
    
#     cs = ax.pcolormesh(lonf, latf, temp[1:, 1:], transform=proj2, shading='auto', cmap=plt.cm.viridis)
#     #cs.set_clim(-1, 3)
#     ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
#     ax.add_feature(cfeature.COASTLINE, zorder=1001)

#     gl = ax.gridlines(**gridparams)

#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90])
#     ax.set_ylim(-40, 40)
#     ax.set_xlim(-60, 130)
#     ax.text(lontext, lattext, letters[cpt] + ")", **textprops)
#     ax.text(lontext2, lattext2, '%.f cm' %ltemp, **textprops, zorder=2000)
    
#     cbax = axgr.cbar_axes[cpt]
#     cb = cbax.colorbar(cs)
#     if(l == 0):
#         #cb.set_label('Mean biomass dens. (Log(J/m2))')
#         ax.set_title('Mean biomass dens.\nLog(J/m2)')
#     manage_axes(gl, cpt)
    
#     cpt += 1
    
## Plotting the anomalies
for l in range(3):
    
    cpt = (l * 2)
    ax = axgr[cpt]
    
    temp = varclim.isel(l=l)
    ltemp = float(temp.l)
    tpclim = np.log10(temp.to_masked_array())
    
    cl = ax.tricontour(lonout, latout, tpclim[imask], colors='k', linewidths=2, levels=levels)
    plt.clabel(cl)
    
    temp = varanom.isel(l=l)
    ltemp = float(temp.l)
    temp = temp.to_masked_array()
    temp = np.ma.masked_where(tmask_noatl==0, temp)
    
    ccc = np.percentile(np.abs(np.ravel(temp[temp.mask == False])), 99)
    print('ccc = ', ccc)

    if(l == 2):
        ccc = 7
    
    cs = ax.pcolormesh(lonf, latf, temp[1:, 1:], transform=proj2, shading='auto', cmap=plt.cm.RdBu_r)
    cs.set_clim(-ccc, ccc)
    ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
    ax.add_feature(cfeature.COASTLINE, zorder=1001)

    gl = ax.gridlines(**gridparams)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90])
    ax.set_ylim(-40, 40)
    ax.set_xlim(-60, 130)
    ax.text(lontext, lattext, letters[cptlet] + ")", **textprops)
    ax.text(lontext2, lattext2, '%.f cm' %ltemp, **textprops)
    
    cbax = axgr.cbar_axes[cpt]
    cb = cbax.colorbar(cs)
    if l == 0:
        #cb.set_label('97-OND biomass anoms (J/m2)')
        ax.set_title('El Nino biomass anoms (J/m2)')
    manage_axes(gl, cpt)
    
    cpt += 1
    cptlet += 1
    
## Plotting the anomalies
for l in range(3):
    
    cpt = (l * 2 + 1)
    ax = axgr[cpt]
    
    temp = varclim.isel(l=l)
    ltemp = float(temp.l)
    tpclim = np.log10(temp.to_masked_array())
    
    cl = ax.tricontour(lonout, latout, tpclim[imask], colors='k', linewidths=2, levels=levels)
    plt.clabel(cl)
    
    temp = cov.isel(l=l)
    ltemp = float(temp.l)
    temp = temp.to_masked_array()
    temp = np.ma.masked_where(tmask_noatl==0, temp)
    
    ccc = np.percentile(np.abs(np.ravel(temp[temp.mask == False])), 98)
    print('ccc = ', ccc)
    
    cs = ax.pcolormesh(lonf, latf, temp[1:, 1:], transform=proj2, shading='auto', cmap=plt.cm.RdBu_r)
    cs.set_clim(-ccc, ccc)
    ax.add_feature(cfeature.LAND, zorder=1000, color='lightgray')
    ax.add_feature(cfeature.COASTLINE, zorder=1001)

    gl = ax.gridlines(**gridparams)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator([150, 180, -180, -150, -120, -90])
    ax.set_ylim(-40, 40)
    ax.set_xlim(-60, 130)
    ax.text(lontext, lattext, letters[cptlet] + ")", **textprops)
    ax.text(lontext2, lattext2, '%.f cm' %ltemp, **textprops)
    
    cbax = axgr.cbar_axes[cpt]
    cb = cbax.colorbar(cs)
    if(l == 0):
        #cb.set_label('Cov. ONI/biomass anoms (J/m2)')
        ax.set_title('Cov. ONI/biomass anoms (J/m2)')
    manage_axes(gl, cpt)
    
    cptlet += 1
    cpt += 1    
    
    
plt.savefig('map_mean_anom_OND_97.png', bbox_inches='tight', facecolor='white')
plt.savefig('gr10.jpg', bbox_inches='tight', facecolor='white')
# -

