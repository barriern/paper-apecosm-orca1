# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
plt.rcParams['text.usetex'] = False

grid = varname = 'gamma1'
ilat = slice(None, -4);
# -

# ## Loading Apecosm length file

const  = xr.open_dataset('data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'w'})
lengths = const['length'] * 100
wstep = const['weight_step']
wstep

# ## Reading mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0, y=ilat)
mesh

lat = mesh['gphif']
lat

lon = mesh['glamf']
lon

# ## Loading the climatology

clim = xr.open_dataset('data/pacific_clim_%s.nc' %(grid)).isel(y=ilat)
clim

varclim = clim[varname]
varclim

# ## Loading the field

data = xr.open_dataset('data/pacific_nino97_%s.nc' %(grid)).isel(y=ilat)
data

nweights = data.dims['w']
nweights

ntime = data.dims['time']
ntime

var = data[varname]
var

# ## Computing the anomalies

anom = var.groupby('time.month') - varclim
anom

if varname == 'OOPE':
    print('Weighting with wstep')
    anom = anom * wstep
    anom.name = 'OOPE'
anom

anom = anom.rename({'w': 'length'})
anom

anom = anom.where(abs(lat) <= 20)
anom

anom['length'] = lengths.values
anom

# ## Plotting the anomalies

mapin = ccrs.PlateCarree()
mapout = ccrs.PlateCarree(central_longitude=180)
for w in [14, 45, 80]:
    print('@@@@@@@@@@@@@@@@@@@@@@@@ Processing w ', w)
    figname = 'pacific_anom_%s_lclass_%.3d.pdf' %(varname, w)
    temp = anom.isel(length=w)
    cmax = abs(temp).quantile(98 / 100).values
    cmax = float(cmax)
    with PdfPages(figname) as pdf:
        for t in range(ntime):
            print(t, '/', ntime)
            fig = plt.figure()
            toplot = temp.isel(time=t)
            ax = plt.axes(projection=mapout)
            cs = ax.pcolormesh(lon, lat, toplot.to_masked_array()[1:, 1:], cmap=plt.cm.RdBu_r, transform=mapin)
            cs.set_clim(-cmax, cmax)
            ax.add_feature(cfeature.LAND, zorder=100)
            ax.add_feature(cfeature.COASTLINE, zorder=101)
            cb = plt.colorbar(cs, shrink=0.5)
            cb.set_label(varname)
            gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=180), draw_labels=True,
                  linewidth=0.5, color='k', alpha=0.5, linestyle='--')
            gl.top_labels = False
            gl.left_labels = True
            gl.right_labels = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            title = '%4d-%.2d, size=%2.f cm' %(toplot['time.year'].values, toplot['time.month'].values, toplot['length'].values)
            ax.set_title(title)
            pdf.savefig(bbox_inches='tight')
            plt.close(fig)


