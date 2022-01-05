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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

chunk = {'y' : 126 //4, 'x': 126 // 2}
# -


# ## Loading mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc')
mesh

volume = mesh['e2t'] * mesh['e1t'] * mesh['e3t_0'] * mesh['tmaskutil'] * mesh['tmask']
volume

lat = mesh['gphit']
lat

volume = volume.where(abs(lat) == 0)
volume.isel(z=0).plot()

volume = volume.rename({'z': 'depth'})

clim = xr.open_dataset('data/pacific_clim_FORAGE.nc').isel(dn=0)
clim

varclim = clim['FORAGE']
varclim

varclim = varclim.chunk(chunk)
varclim

# Computation of anomalies:

anoms = xr.open_dataset('data/pacific_nino97_FORAGE.nc').sel(time=slice('1997-10-01', '1997-12-31')).isel(dn=0)
anoms                                                             

varanoms = anoms['FORAGE']
varanoms

varanoms = varanoms.chunk(chunk)
varanoms

varanoms = varanoms.groupby('time.month') - varclim
varanoms

# ## Compute the mean climatology

tsclim = ((varclim * volume).sum(dim=['y']) / (volume.sum(dim=['y']))).mean(dim='month')
tsclim

delayed_clim = tsclim.to_netcdf('mean_forage.nc', compute=False)

tsanoms = ((varanoms * volume).sum(dim=['y']) / (volume.sum(dim=['y']))).mean(dim='time')
tsanoms

delayed_anoms = tsanoms.to_netcdf('mean_forage_anomalies_ond_97.nc', compute=False)

with ProgressBar():
    delayed_clim.compute()

with ProgressBar():
    delayed_anoms.compute()
