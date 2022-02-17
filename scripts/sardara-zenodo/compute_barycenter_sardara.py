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

# +
import xarray as xr

data = xr.open_dataset('regridded_catch_gear_PS.nc')
# -

data

catch = data['catch'].sum(dim='species')
catch

catch = catch.where(abs(data['lat']) <= 10, drop=True)
catch

catch = catch.where((data['lon'] >= 120) & (data['lon'] <= -120+360))
catch

barycenter = (catch['lon'] * catch).sum(dim=['lon', 'lat']) / catch.sum(dim=['lon', 'lat'])
barycenter.name = 'catch'
barycenter

barycenter.plot()

barycenter.to_netcdf('barycenter_sardara.nc')

barycenter.to_dataframe(name='lon').to_csv('barycenter_sardara.csv')

import matplotlib.pyplot as plt
import numpy as np
time = barycenter['time'].values
iok = np.nonzero((time >= 199001) & (time <= 201612))[0]
barycenter = barycenter.isel(time=iok)

test = xr.open_dataset('barycenter_old_sardara.nc')
time = test['time'].values
iok = np.nonzero((time >= 199001) & (time <= 201612))[0]
test = test['catch'].isel(time=iok)


# +
plt.figure(figsize=(12, 8))


plt.plot(test.values, label='new')
plt.plot(barycenter.values, label='zeno')
plt.legend(fontsize=15)
# -

boolean = np.isnan(test.values) | np.isnan(barycenter.values)
itest = np.nonzero(boolean == False)
print(np.corrcoef(test.values[itest], barycenter.values[itest]))


