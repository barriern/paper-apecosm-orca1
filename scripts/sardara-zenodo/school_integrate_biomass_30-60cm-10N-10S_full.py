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

import xarray as xr
from dask.diagnostics import ProgressBar

mesh = xr.open_dataset('/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0, z=0)
lat = mesh['gphit']
mesh = mesh.where(abs(lat) <= 10)
surf = mesh['e1t'] * mesh['e2t'] * mesh['tmaskutil']
surf.plot()

const = xr.open_dataset('/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
const = const.rename({'wpred': 'l'})
const['l'] = const['length'] * 100
const

data = xr.open_mfdataset('/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/*school_day*Y*nc')
data = data.isel(community=0)
data = data.rename({'w': 'l'})
data['l'] = const['l']
school = data['school_day']
school

data = xr.open_mfdataset('/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/*OOPE*Y*nc')
data = data.isel(community=0)
data = data.rename({'w': 'l'})
data['l'] = const['l']
data

output = (data['OOPE'] * const['weight_step'] * surf * school).sel(l=slice(30, 70)).sum(dim=['l', 'y'])
output.name = 'OOPE'
output['x'] = mesh['glamt'].mean(dim='y')
output

delayed = output.to_netcdf('school_integrated_biomass_30-70cm-10N-10S.nc', compute=False)
with ProgressBar():
    delayed.compute()
