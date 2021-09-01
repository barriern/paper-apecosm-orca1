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
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import matplotlib.dates as mdates
plt.rcParams['text.usetex'] = False

grid = varname = 'mort_day'
latmax = 1
# -

# ## Loading Apecosm length file

const  = xr.open_dataset('data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'w'})
lengths = const['length'] * 100
wstep = const['weight_step']
wstep

# ## Reading mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
mesh

lat = mesh['gphit']
lat
lat.plot()

mesh = mesh.where(abs(lat) <= latmax)

lon = mesh['glamt'].mean(dim='y')
lon

lon = (lon + 360) % 360
lon

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmask']
surf

# ## Loading the climatology

clim = xr.open_dataset('data/pacific_clim_%s.nc' %(grid))
clim

varclim = clim[varname]
varclim

# ## Loading the field

data = xr.open_dataset('data/pacific_nino97_%s.nc' %(grid))
data

nweights = data.dims['w']
nweights

var = data[varname]
var

# ## Computing the anomalies

anom = var.groupby('time.month') - varclim
anom

# ## Computing the meridional anomalies

anomout = (anom * surf).sum(dim='y') / surf.sum(dim='y')
anomout

if(varname == 'OOPE'):
    print('Converting OOPE to density')
    anomout = anomout * wstep
    anomout.name = 'OOPE'
anomout

anomout['x'] = lon

anomout['w'] = lengths.values
anomout

anomout = anomout.rename({'w': 'length'})
anomout

# +
#with PdfPages('pacific_anom_hov_%s.pdf' %varname) as pdf:
#    for w in range(nweights):
#        print(w, '/', nweights)
#        fig = plt.figure()
#        cs = anomout.isel(length=w).plot(robust=True)
#        pdf.savefig(bbox_inches='tight')
#        plt.grid(True, color='k', linestyle='--')
#        plt.close()

# +
fig = plt.figure(figsize=(12, 10))
grid = ImageGrid(fig, 111,  # similar to subplot(111)
                 nrows_ncols=(2, 2),  # creates 2x2 grid of axes
                 axes_pad=[0.7, 0.3],  # pad between axes in inch.
                 cbar_mode='each', aspect=False, cbar_pad=0.1)
x = anomout['x']
y = data['time']
cbar_axes = grid.cbar_axes

cpt = 0
for w in [0, 14, 45, 90]:
    ax = grid[cpt]
    temp = anomout.isel(length=w)
    cmax = float(abs(temp).quantile(0.995))
    #cmax = float(abs(temp).max())
    cs = grid[cpt].pcolormesh(x.values, y.values, temp, shading='auto')
    cl = grid[cpt].contour(x.values, y.values, temp, levels=np.linspace(-cmax, cmax, 11), colors='k', linewidths=0.5)
    cl0 = grid[cpt].contour(x.values, y.values, temp, levels=0, colors='k', linewidths=1)
    cs.set_clim(-cmax, cmax)
    cb = plt.colorbar(cs, cbar_axes[cpt])
    ax.set_title('Length = %.2e cm' %lengths[w])
    p = grid[cpt].grid(True, color='k', linestyle='--', linewidth=.5)
    cpt += 1
sp = plt.suptitle(varname)
plt.savefig('pacific_anom_hov_%s.pdf' %varname, bbox_inches='tight')
# -




