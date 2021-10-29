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
import os
from dask.diagnostics import ProgressBar

lonmin = 160
lonmax = -90
latmax = 5
mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(z=0)
mesh
# -

lon = mesh['glamt']
lat = mesh['gphit']

domain = ((lon > lonmin) | (lon < lonmax))
domain = domain & (abs(lat) < latmax)
domain.plot()

mask = mesh['tmask'].where(domain)
mask.plot()

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmask']
surf = surf.where(domain)
surf.plot()

scratch = os.getenv('SCRATCH')
data = xr.open_mfdataset(os.path.join(scratch, 'scratch', 'pacific*nc'))
data

oope = data['OOPE']
oope

ts = (oope * surf).sum(dim=('y', 'x'))
ts

dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'
fileout = 'apecosm_sardara_latmax_%d_lonmin_%d_lonmax_%d.nc' %(latmax, lonmin, lonmax)
fileout = os.path.join(dirout, fileout)
fileout

delayed = ts.to_netcdf(fileout, unlimited_dims='time', compute=False)

with ProgressBar():
    delayed.compute()


