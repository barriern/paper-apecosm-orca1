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


