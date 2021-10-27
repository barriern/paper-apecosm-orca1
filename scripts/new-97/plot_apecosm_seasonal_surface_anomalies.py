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
plt.rcParams['text.usetex'] = False
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes

grid = varname = 'OOPE'
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

clim = xr.open_mfdataset('data/pacific_clim_*%s.nc' %(grid)).isel(y=ilat)
clim

if varname == 'adv_trend':
    print('summing trend')
    varclim = clim['madv_trend'] + clim['zadv_trend']
    varclim.name = varname
else:
    varclim = clim[varname]
varclim

# ## Loading the field

data = xr.open_mfdataset('data/pacific_nino97_*%s.nc' %(grid)).isel(y=ilat)
data

nweights = data.dims['w']
nweights

ntime = data.dims['time']
ntime

if varname == 'adv_trend':
    print('Summing trend')
    var = data['madv_trend'] + data['zadv_trend']
    var.name = varname
else:
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

# +
#anom = anom.where(abs(lat) <= 20)
#anom
# -

anom['length'] = lengths.values
anom

anomyear = anom.groupby('time.year').mean(dim='time')
anomyear

cmax_ar = abs(anomyear).chunk({'year': -1}).quantile(0.98, dim=['x', 'y', 'year']).compute()

# ## Plotting the anomalies

# +
mapin = ccrs.PlateCarree()
mapout = ccrs.PlateCarree(central_longitude=180)

fig = plt.figure(figsize=(14, 10))

axes_class = (GeoAxes, dict(map_projection=mapout))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 3), axes_pad=(0.7, 0.5), label_mode='', cbar_mode='each', cbar_size=0.1, cbar_pad=0.3, cbar_location="right")
cbar_axes = axgr.cbar_axes
wpred = [14, 45, 80]

cpt = 0
for y in range(1996, 1999):
    temp = anomyear.sel(year=y)
    for s in wpred:
        print(s)
        cmax = cmax_ar.isel(length=s).values
        tempfinal = temp.isel(length=s)
        ax = axgr[cpt]
        cs = ax.pcolormesh(lon, lat, tempfinal.to_masked_array()[1:, 1:], cmap=plt.cm.RdBu_r, transform=mapin)
        cs.set_clim(-cmax, cmax)
        ax.add_feature(cfeature.LAND, zorder=100)
        ax.add_feature(cfeature.COASTLINE, zorder=101)
        cb = plt.colorbar(cs, cbar_axes[cpt], shrink=0.5, )
        cb.set_label('J/m2')
        ax.set_title('%.4d, L=%.f cm' %(y, lengths[s]))
        cpt += 1
# -


