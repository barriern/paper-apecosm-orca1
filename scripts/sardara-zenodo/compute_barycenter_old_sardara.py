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

data1 = xr.open_dataset('../sardara/data/regridded_catch_gear_PS_species_SKJ_1x1.nc')
data2 = xr.open_dataset('../sardara/data/regridded_catch_gear_PS_species_YFT_1x1.nc')
catch = data1['catch'] + data2['catch']
catch
# -

catch = catch.where(abs(data1['lat']) <= 10, drop=True)
catch

catch = catch.where((data1['lon'] >= 120) & (data1['lon'] <= -120+360))

barycenter = (catch['lon'] * catch).sum(dim=['lon', 'lat']) / catch.sum(dim=['lon', 'lat'])
barycenter.name = 'catch'
barycenter

barycenter.plot()

barycenter.to_netcdf('barycenter_old_sardara.nc')


