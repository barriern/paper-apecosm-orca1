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

grid = varname = 'repfonct_day'
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

with PdfPages('pacific_anom_hov_%s.pdf' %varname) as pdf:
    for w in range(nweights):
        print(w, '/', nweights)
        fig = plt.figure()
        cs = anomout.isel(length=w).plot(robust=True)
        pdf.savefig(bbox_inches='tight')
        plt.grid(True, color='k', linestyle='--')
        plt.close()


