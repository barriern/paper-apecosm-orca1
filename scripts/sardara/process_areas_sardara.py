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

# # Processing Sardara area data
#
# Data provided by Jonathan Rault (Datarmor: `/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1_REF/analyses/sardara/data-analysis/sardara/julien-dataset/data`)
#
# ## Imports

# +
import pandas as pd
import numpy as np
import re
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

pacific_latmax = 40
pacific_lonmin = 140
pacific_lonmax = -60
# -

# ## Reading the CSV

datatemp = pd.read_csv('data/geographic_id_wkt.csv', sep=',')
datatemp

# ## Working on the polygons
#
# ### Removing of wrong polygons

polygons_list_temp = datatemp.loc[:, 'st_astext'].values
polygons_list_temp

# We extract the polygons which are really string, not int (i.e we discard the `NaN` value of `UNK`)

test = np.array([isinstance(a, str) for a in polygons_list_temp])
test

iok = np.nonzero(test == True)[0]
iok

data = datatemp.iloc[iok, :]
data

# ### Extract polygon coordinates

# Now we extract the longitudes based on the strig value, using regular expressions

polygons_list = data.loc[:, 'st_astext'].values
polygons_list

pattern = 'MULTIPOLYGON\(\(\((-?[0-9]+)\ (-?[0-9]+),(-?[0-9]+)\ (-?[0-9]+),(-?[0-9]+)\ (-?[0-9]+),(-?[0-9]+)\ (-?[0-9]+).*\)\)\)'
regex = re.compile(pattern, re.VERBOSE)
regex

# +
lonmin = []
lonmax = []
latmin = []
latmax = []
for temp in polygons_list[:]:
    match = regex.match(str(temp))
    if(match is None):
        print(temp, ' not processed')
    else:
        coords = np.array(match.groups()).astype(float)
    lon = coords[::2]
    lat =  coords[1::2]
    latmin.append(lat.min())
    latmax.append(lat.max())
    lonmin.append(lon.min())
    lonmax.append(lon.max())

lonmin = np.array(lonmin)
latmin = np.array(latmin)
lonmax = np.array(lonmax)
latmax = np.array(latmax)
# -

# We now extract the barycenter of the polygons

loncen = 0.5 * (lonmin + lonmax)
loncen
latcen = 0.5 * (latmin + latmax)

# We convert them into Atlantic coordinates (-180/180)

loncen[loncen > 180] -= 360
loncen[loncen < -180] += 360
loncen.min(), loncen.max()

# Now we extract the Pacific data:

testlat = np.abs(latcen) < pacific_latmax
testlon1 = (loncen <= pacific_lonmax)
testlon2 = (loncen >= pacific_lonmin)
test = testlat & ((testlon2 | testlon1))
iok = np.nonzero(test == True)[0]

# Now, we plot the resulting coordinates for the area

ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
plt.scatter(loncen[iok], latcen[iok], transform=ccrs.PlateCarree(), marker='.', color='r')
ax.coastlines()
ax.set_global()
plt.savefig('areas.png')

data = data.iloc[iok, :]
data

# ## Creating output area dataset

# Now, we write in the dataset the coordinates of the box. But first, we add the data coordinates that we computed.

data = data.assign(lonmin=lonmin[iok])
data

data = data.assign(lonmax=lonmax[iok])
data

data = data.assign(latmax=latmax[iok])
data

data = data.assign(latmin=latmin[iok])
data

data = data.assign(latcen=latcen[iok], loncen=loncen[iok])
data

# Finally, we write the new dataframe

output_ds = data.loc[:, ['code', 'lonmin', 'lonmax', 'latmin', 'latmax', 'loncen', 'latcen']]
output_ds

output_ds.to_csv('data/processed_sardara_areas_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax))
