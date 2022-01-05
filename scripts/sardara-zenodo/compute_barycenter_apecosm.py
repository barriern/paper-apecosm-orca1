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

prefix = ''

data = xr.open_dataset('%sintegrated_biomass_30-70cm-10N-10S.nc' %prefix)
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

import matplotlib.pyplot as plt
import numpy as np
plt.figure(figsize=(8, 12))
time = np.arange(barycenter.shape[0])
plt.plot(barycenter.values, time)
#plt.gca().set_xticks(np.arange(140, -120 + 360, 20))

barycenter.to_netcdf('%sbarycenter_apecosm.nc' %prefix)

barycenter.to_dataframe(name='lon').to_csv('%sbarycenter_apecosm.csv' %prefix)


