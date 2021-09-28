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

import xarray as xr
import matplotlib.pyplot as plt


# ## Loading mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc')
mesh

volume = mesh['e2t'] * mesh['e1t'] * mesh['e3t_0'] * mesh['tmaskutil'] * mesh['tmask']
volume

lat = mesh['gphit']
lat

volume = volume.where(abs(lat) <= 5)
volume.isel(z=0).plot()

clim = xr.open_dataset('data/pacific_clim_FORAGE.nc').isel(dn=0)
clim

varclim = clim['FORAGE']
varclim

anoms = xr.open_dataset('data/pacific_nino97_FORAGE.nc').sel(time=slice('1997-10-01', '1997-12-31')).isel(dn=0)
anoms                                                             

varanoms = anoms['FORAGE']
varanoms

varanoms = varanoms.groupby('time.month') - varclim
varanoms

# ## Compute the mean climatology

varclim = (varclim * volume).sum(dim=['y', 'x']) / (volume.sum(dim=['y', 'x']))
varclim


