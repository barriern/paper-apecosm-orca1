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

# +
import xarray as xr
import numpy as np

varname = 'predationTrend'
offset = 1
ylist = [1982, 1997, 2015]
# -

data = xr.open_dataset('data/equatorial_full_%s.nc' %varname)[varname]
data = data.groupby('time.month') - data.sel(time=slice('1971-01-01', '2000-12-31')).groupby('time.month').mean(dim='time')

for y in ylist:
    datestart = '%.4d-01-01' %y
    dateend = '%.4d-12-31' %(y + offset)
    temp = data.sel(time=slice(datestart, dateend))
    if y == ylist[0]:
        output = temp.values
    else:
        output += temp.values
output /= len(ylist)

dsout = xr.Dataset()
dsout['x'] = data['x']
dsout['time'] = np.arange((offset + 1) * 12)
dsout[varname] = (['time', 'x', 'l'], output)
dsout[varname].attrs['description'] = 'Composites of monthly anomalies for ' + ','.join(np.array(ylist, dtype=str))
dsout.to_netcdf('data/nino_equatorial_composites_%s.nc' %varname)
