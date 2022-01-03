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

import pandas as pd
import numpy as np
import xarray as xr

data = pd.read_csv('processed_sardara_3.csv', index_col=0)
data

years = np.unique(data['year'].values)
years

m2d, y2d = np.meshgrid(np.arange(1, 13), years)
time = np.ravel(100 * y2d + m2d)
time[:24]

lonmin = data['lonmin'].values
lonmax = data['lonmax'].values
lonmin[lonmin < 0] += 360
lonmax[lonmax < 0] += 360

latmin = data['latmin'].values
latmax = data['latmax'].values

lat = np.arange(-30, 30 + 1, 1)
lat

lon = np.arange(130, -60 + 360 + 1, 1)
lon

ntime = len(time)
nlat = len(lat)
nlon = len(lon)
output = np.zeros((ntime, 2, nlat - 1, nlon - 1), dtype=float)

for i in range(data.shape[0]):
    
    if(i % 10000 == 0):
        print('%d / %d' %(i, data.shape[0]))
    
    temp = data.iloc[i, :]  # extract one row of data
    temptime = temp['year'] * 100 + temp['month']  # recover the time
    species = temp['species']
    if(species == 'SKJ'):
        ispecies = 0
    elif species == 'YFT':
        ispecies = 1
    itime = np.nonzero(time == temptime)[0][0]  # get the time index

    # recover the areas associated with the data and its coordinate
    lonamin = lonmin[i]
    lonamax = lonmax[i]

    # recover the areas associated with the data and its coordinate
    latamin = latmin[i]
    latamax = latmax[i]

    # we loop over the regular grid
    for y in range(nlat - 1):
        latcmin = lat[y]
        latcmax = lat[y + 1]
        for x in range(nlon - 1):
            loncmin = lon[x]
            loncmax = lon[x + 1]

            # if the area is out of the new cell, nothing is done
            if(lonamax <= loncmin):
                continue

            if(lonamin >= loncmax):
                continue

            if(latamax <= latcmin):
                continue

            if(latamin >= latcmax):
                continue

            recouvx = (min(loncmax, lonamax) - max(loncmin, lonamin)) / (lonamax - lonamin)
            recouvy = (min(latcmax, latamax) - max(latcmin, latamin)) / (latamax - latamin)

            output[itime, ispecies, y, x] += (recouvx * recouvy) * temp['value']
    

# +
dsout = xr.Dataset()
dsout

dsout['lon'] = (['lon'], 0.5 * (lon[:-1] + lon[1:]))
dsout['lat'] = (['lat'], 0.5 * (lat[:-1] + lat[1:]))
dsout['time'] = (['time'], time)
dsout['catch'] = (['time', 'species', 'lat', 'lon'], output)
dsout['catch'].attrs['species'] = '0=SKJ, 1=YFT'
dsout['catch'].attrs['gear'] = 'PS'
dsout

# +
fileout = 'regridded_catch_gear_PS.nc'
fileout

dsout.to_netcdf(fileout)
