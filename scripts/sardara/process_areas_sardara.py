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

pacific_latmax = 5
pacific_lonmin = 160
pacific_lonmax = -90
# -

# ## Reading the CSV

data = pd.read_csv('data/geographic_id_wkt.csv', sep=',')
data

# ## Working on the polygons
#
# ### Removing wrong polygon

polygons_list = data.loc[:, 'st_astext'].values
polygons_list

# We extract the polygons which are really string, not int (i.e we discard the `NaN` value of `UNK`)

test = np.array([isinstance(a, str) for a in polygons_list])
test

iok = np.nonzero(test == True)[0]
iok

data = data.iloc[iok, :]
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
    match = regex.match(temp)
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

# We extract the codes for Pacific data

lonmin[lonmin > 180] -= 360
lonmax[lonmax > 180] -= 360
lonmin[lonmin < -180] += 360
lonmax[lonmax < -180] += 360

testlat = ((np.abs(latmax) < pacific_latmax)) & ((np.abs(latmin) < pacific_latmax))
testlon1 = (lonmax <= pacific_lonmax) & (lonmin >= 0)
testlon2 = (lonmin >= pacific_lonmin) & (lonmin <= 180)
test = testlat & (testlon1 | testlon2)
iok = np.nonzero(test == True)[0]
len(iok)

lonmin = lonmin[iok]
lonmax = lonmax[iok]
latmin = latmin[iok]
latmax = latmax[iok]

ax = plt.axes(projection=ccrs.PlateCarree())
for i in range(20):
    x = [lonmin[i], lonmax[i], lonmax[i], lonmin[i], lonmin[i]]
    y = [latmin[i], latmin[i], latmax[i], latmax[i], latmin[i]]
    plt.plot(x, y, transform=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()

data = data.iloc[iok, :]
data

# ## Creating output area dataset

# Now, we write in the dataset the coordinates of the box.

data = data.assign(lonmin=lonmin)
data

data = data.assign(lonmax=lonmax)
data

data = data.assign(latmax=latmax)
data

data = data.assign(latmin=latmin)
data

# Finally, we write the new dataframe

output_ds = data.loc[:, ['code', 'lonmin', 'lonmax', 'latmin', 'latmax']]
output_ds

output_ds.to_csv('data/processed_sardara_areas_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax))
