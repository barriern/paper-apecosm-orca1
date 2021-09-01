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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

grid = 'speed_V'
varname = 'vo'
ilat = slice(None, -4)
# -

# ## Reading mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0, y=ilat)
mesh

lat = mesh['gphif']
lat

lon = mesh['glamf']
lon

# ## Loading the climatology

clim = xr.open_dataset('data/pacific_clim_%s.nc' %(grid)).isel(olevel=0, y=ilat)
clim

varclim = clim[varname]
varclim

# ## Loading the field

data = xr.open_dataset('data/pacific_nino97_%s.nc' %(grid)).isel(olevel=0, y=ilat)
data

ntime = data.dims['time_counter']
ntime

var = data[varname]
var

# ## Computing the anomalies

anom = var.groupby('time_counter.month') - varclim
anom

# ## Plotting the anomalies

mapin = ccrs.PlateCarree()
mapout = ccrs.PlateCarree(central_longitude=180)
figname = 'pacific_surface_anom_%s.pdf' %(varname)
figname

cmax = abs(anom).quantile(99.5 / 100).values
cmax = float(cmax)
cmax

with PdfPages(figname) as pdf:
    for t in range(ntime):
        print(t, '/', ntime)
        fig = plt.figure()
        toplot = anom.isel(time_counter=t)
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
        title = '%4d-%.2d' %(toplot['time_counter.year'].values, toplot['time_counter.month'].values)
        ax.set_title(title)
        pdf.savefig(bbox_inches='tight')
        plt.close(fig)




