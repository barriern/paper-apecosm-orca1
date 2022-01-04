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

import xarray as xr

data = xr.open_dataset('integrated_biomass_30-70cm-10N-10S.nc')
data

data['x'] = (data['x'] + 360) % 360
data = data.sortby(data['x'])
data

data['OOPE'].plot()

boolean = (data['x'] <= -120 + 360) & (data['x'] >= 120)
data = data.where(boolean, drop=True)
data

barycenter = (data['x'] * data['OOPE']).sum(dim='x') / (data['OOPE'].sum(dim='x'))
barycenter

barycenter.plot()

barycenter.to_netcdf('barycenter_apecosm.nc')
