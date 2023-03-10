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

import xarray as xr
from dask.diagnostics import ProgressBar

mesh = xr.open_dataset('../../data/static/pacific_mesh_mask.nc')
mesh = mesh.isel(z=0)
lat = mesh['gphit']
mesh = mesh.where(abs(lat) <= 10)
surf = mesh['e1t'] * mesh['e2t'] * mesh['tmaskutil']
surf.plot()

const = xr.open_dataset('../../data/static/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

data = xr.open_mfdataset('../../data/apecosm/*_OOPE_Y*nc')
data

data = data.rename({'w': 'l'})
data['l'] = const['l']
data

output = (data['OOPE'] * const['weight_step'] * surf).sel(l=slice(30, 70)).sum(dim=['l', 'y'])
output.name = 'OOPE'
output['x'] = mesh['glamt'].mean(dim='y')
output

delayed = output.to_netcdf('integrated_biomass_30-70cm-10N-10S.nc', compute=False)
with ProgressBar():
    delayed.compute()
