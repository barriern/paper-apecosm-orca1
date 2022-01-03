# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import xarray as xr

dirin = '../new-97/data/'

const = xr.open_dataset('%s/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc' %dirin)
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

data = xr.open_dataset('%s/pacific_nino97_OOPE.nc' %dirin)
data = data.rename({'w': 'l'})
data['l'] = const['l']
data['OOPE']

mesh = xr.open_dataset('%s/pacific_mesh_mask.nc' %dirin).isel(z=0)
mesh

lat = mesh['gphit']

mesh = mesh.where(abs(lat) <= 10)

surf = mesh['e1t'] * mesh['e2t'] * mesh['tmask']
surf.plot()

hovoope = (data['OOPE'] * surf).sum(dim='y')
hovoope

hovoope = hovoope * const['weight_step']
hovoope

temp = hovoope.sel(l=slice(30, 70)).sum(dim='l')
temp['x'] = mesh['glamt'].mean(dim='y')
temp

temp.to_netcdf('hovmoller_apecosm.nc')
