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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar
from glob import glob

chunk = {'y' : 126 //4, 'x': 126 // 2}

filelist

varname = 'PLK'
filename = 'ptrc_T'

filelist = glob('data/*%s*nc' %filename)
# -


filelist
xr.open_dataset(filelist[0])

glob('*nc')

# ## Loading mesh mask

mesh = xr.open_dataset('data/pacific_mesh_mask.nc')
mesh

volume = mesh['e2t'] * mesh['e1t'] * mesh['e3t_0'] * mesh['tmaskutil'] * mesh['tmask']
volume

lat = mesh['gphit']
lat

volume = volume.where(abs(lat) <= 5)
volume.isel(z=0).plot()

volume = volume.rename({'z': 'olevel'})
volume

clim = xr.open_dataset('data/pacific_clim_%s.nc' %filename)
clim

if(varname == 'PLK'):
    varclim = clim['PHY2'] + clim['GOC'] + clim['ZOO'] + clim['ZOO2']
    varclim.name = 'PLK'
else:
    varclim = clim[varname]
varclim

# Computation of anomalies:

toread = 'data/pacific_nino97_%s.nc' %filename
toread
anoms = xr.open_dataset(toread).sel(time_counter=slice('1997-10-01', '1997-12-31'))
anoms                                                             

if(varname == 'PLK'):
    varanoms = anoms['PHY2'] + anoms['GOC'] + anoms['ZOO'] + anoms['ZOO2']
    varanoms.name = 'PLK'
else:
    varanoms = anoms[varname]
varanoms


varanoms = varanoms.groupby('time_counter.month') - varclim
varanoms

# ## Compute the mean climatology

tsclim = ((varclim * volume).sum(dim=['y']) / (volume.sum(dim=['y']))).mean(dim='month')
tsclim

delayed_clim = tsclim.to_netcdf('mean_%s.nc' %varname, compute=False)

tsanoms = ((varanoms * volume).sum(dim=['y']) / (volume.sum(dim=['y']))).mean(dim='time_counter')
tsanoms

delayed_anoms = tsclim.to_netcdf('mean_%s_anomalies_ond_97.nc' %varname, compute=False)

with ProgressBar():
    delayed_clim.compute()

with ProgressBar():
    delayed_anoms.compute()


