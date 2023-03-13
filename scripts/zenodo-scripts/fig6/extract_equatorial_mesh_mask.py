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
from glob import glob
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False

ilat = 61
# -

mesh = xr.open_dataset('../data/static/pacific_mesh_mask.nc').isel(y=ilat)
mesh

lon0 = mesh['glamt'].isel(z=0)
lon0

mesh['x'] = (lon0 + 360) % 360
mesh = mesh.sortby(mesh['x'])
mesh

mesh['tmask'].plot()

delayed = mesh.to_netcdf('equatorial_mesh_mask.nc', compute=False)
with ProgressBar():
    delayed.compute()


