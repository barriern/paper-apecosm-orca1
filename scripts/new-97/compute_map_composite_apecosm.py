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
#     display_name: Python [conda env:nbarrier]
#     language: python
#     name: conda-env-nbarrier-py
# ---

# +
import xarray as xr
from glob import glob
import numpy as np

for varname in ['madv_trend', 'zadv_trend', 'zdiff_trend', 'mdiff_trend']:

    filelist = []
    for y in [1982, 1997, 2015]:
        pattern = '/home1/scratch/nbarrier/pacific_ORCA1*%s*%d*.nc' %(varname, y)
        temp = glob(pattern)[0]
        filelist.append(temp)
    print(filelist)

    data = xr.open_mfdataset(filelist)[varname]
    print(data)

    clim = xr.open_dataset('data/pacific_clim_%s.nc' %varname)[varname]
    print(clim)

    anom = data.groupby('time.month') - clim
    print(anom)

    months = anom['time.month'].values
    iok = np.nonzero(months >=10)[0]
    print(iok)

    anom = anom.isel(time=iok)
    print(anom['time'])

    compo = anom.mean(dim='time')
    print(compo)

    fileout = 'data/composite_%s_map.nc' %varname
    print(fileout)

    compo.to_netcdf(fileout)
# -


