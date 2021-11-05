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
import pandas as pd
import numpy as np
import re
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import re
from mpl_toolkits.axes_grid1 import ImageGrid
plt.rcParams['text.usetex'] = False
import xarray as xr

pacific_latmax = 40
pacific_lonmin = 120
pacific_lonmax = -60

gear = 'PS'
species = 'YFT'
# -

# ## Extracting the code of the gear

# Here, we extract the code that corresponds to the gear.

data = pd.read_csv('data/CL_FI_GEAR_LEVEL2.csv')
code = data[data['Abbreviation'] == gear].loc[:, 'Code']
code = code.values[0]
code

# ## Extracting the areas

areas = pd.read_csv('data/processed_sardara_areas_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax), index_col=0)
areas

# Now we convert the longitudes into a Pacific System:

loncen = areas.loc[:, 'loncen'].values
loncen[loncen < 0] += 360
loncen

lonmin = areas.loc[:, 'lonmin'].values
lonmin[lonmin < 0] += 360
lonmin

lonmax = areas.loc[:, 'lonmax'].values
lonmax[lonmax < 0] += 360
lonmax

areas = areas.assign(lonmin=lonmin, lonmax=lonmax, loncen=loncen)
areas

# Now we create the output grid, considering a 5x5 grid

lonmin = areas.loc[:, 'lonmin'].min()
lonmax = areas.loc[:, 'lonmax'].max()
lonout = np.arange(lonmin, lonmax + 5, 5)
lonout

latmin = areas.loc[:, 'latmin'].min()
latmax = areas.loc[:, 'latmax'].max()
latout = np.arange(latmin, latmax + 5, 5)
latout

lon2d, lat2d = np.meshgrid(lonout, latout)

# ## Data processing

# Now we process the data:

data = pd.read_csv('data/processed_sardara_data_L2_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax))
data

# We extract the time as an array

dates = [d.split('-') for d in data.loc[:, 'time_start']]

years = np.array([d[0] for d in dates], dtype=int)
years

months = np.array([d[1] for d in dates], dtype=int)
months

data = data.assign(years=years, months=months)
data

years = np.unique(years)
years

# We extract the data that corresponds to the species:

data = data[data['species'] == species]

data = data[data['gear'] == code]
data

# ## Regridding the data
#

lat2d.shape
lon2d.shape
ntime = len(years)

m2d, y2d = np.meshgrid(np.arange(1, 13), years)
time = np.ravel(100 * y2d + m2d)

nlat = latout.shape[0]
nlon = lonout.shape[0]
ntime = len(time)
output = np.zeros((ntime, nlat - 1, nlon - 1))
output.shape

for i in range(data.shape[0]):
    
    if(i % 1000 == 0):
        print('%d / %d' %(i, data.shape[0]))
    
    temp = data.iloc[i, :]  # extract one row of data
    temptime = temp['years'] * 100 + temp['months']  # recover the time
    itime = np.nonzero(time == temptime)[0][0]  # get the time index
    
    # recover the areas associated with the data and its coordinate
    temparea = areas[areas['code'] == temp['geographic_identifier']]
    lonamin = temparea['lonmin'].values[0]
    lonamax = temparea['lonmax'].values[0]
    latamin = temparea['latmin'].values[0]
    latamax = temparea['latmax'].values[0]
    
    # we loop over the regular grid
    for y in range(nlat - 1):
        latcmin = latout[y]
        latcmax = latout[y + 1]
        for x in range(nlon - 1):
            loncmin = lonout[x]
            loncmax = lonout[x + 1]
            
            # 
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
        
            output[itime, y, x] += (recouvx * recouvy) * temp['value']

dsout = xr.Dataset()
dsout

dsout['lon'] = (['lon'], 0.5 * (lonout[:-1] + lonout[1:]))
dsout['lat'] = (['lat'], 0.5 * (latout[:-1] + latout[1:]))
dsout['time'] = (['time'], time)
dsout['catch'] = (['time', 'lat', 'lon'], output)
dsout['catch'].attrs['species'] = species
dsout['catch'].attrs['gear'] = gear

fileout = 'data/regridded_catch_gear_%s_species_%s.nc' %(gear, species)
fileout

dsout.to_netcdf(fileout)


