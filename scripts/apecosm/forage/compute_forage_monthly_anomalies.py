import xarray as xr
import numpy as np
import apecosm.ts as ts

data = xr.open_mfdataset("data/hov*nc", combine='by_coords')

forage = np.squeeze(data['FORAGE'].to_masked_array())

clim, anom = ts.get_monthly_clim(forage)
print(clim.shape, anom.shape)

dsout = xr.Dataset()
dsout['depth'] = data['depth']
dsout['x'] = data['x']
dsout['time'] = data['time']
dsout['clim'] = (['months', 'dn', 'x', 'depth', 'size_group'], clim)
dsout['anom'] = (['time', 'dn', 'x', 'depth', 'size_group'], anom)
dsout.to_netcdf('forage_anom.nc', unlimited_dims='time')

