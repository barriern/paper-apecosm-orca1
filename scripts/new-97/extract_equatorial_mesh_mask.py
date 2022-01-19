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
from glob import glob
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False

ilat = 61
# -

mesh = xr.open_dataset('data/pacific_mesh_mask.nc').isel(y=ilat)
mesh

lon0 = mesh['glamt'].isel(z=0)
lon0

mesh['x'] = (lon0 + 360) % 360
mesh = mesh.sortby(mesh['x'])
mesh

mesh['tmask'].plot()

delayed = mesh.to_netcdf('data/equatorial_mesh_mask.nc', compute=False)
with ProgressBar():
    delayed.compute()


