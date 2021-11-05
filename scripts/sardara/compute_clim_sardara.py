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

pacific_latmax = 40
pacific_lonmin = 140
pacific_lonmax = -60

datestart = 2006
dateend = 2015
years = list(range(datestart, dateend + 1))
ntime = len(years) * 12
ntime
# -

# First, the pre-processed file is loaded:

data = pd.read_csv('data/processed_sardara_data_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax), index_col=0)
data

# Then, we convert the time-start column into a dates using regular expressions

dates = data.loc[:, 'time_start']
dates

regex = re.compile('([0-9][0-9][0-9][0-9])\-([0-9][0-9])\-([0-9][0-9])')
yymmdd = np.array([regex.match(d).groups() for d in dates]).astype(int)
yymmdd.shape

# The years and months are added to the file:

data = data.assign(year=yymmdd[:, 0])
data = data.assign(month=yymmdd[:, 1])
data

# Now we extract the years that correspond to the climatology:

iok = np.nonzero((yymmdd[:, 0] >= datestart) & (yymmdd[:, 0] <= dateend))[0]
data.iloc[iok, :]

# Now, we sum the catches that have the same dates, geographic area and species (this is to take into account that different countries may fish on the same areas at the same dates).

dataint = data.iloc[iok, :].groupby(['time_start', 'geographic_identifier', 'species']).sum().loc[:, 'value']
dataint

# Now, we can compute the temporal mean:

dataclim = dataint.groupby(['geographic_identifier', 'species'])
dataclim


datamean = dataclim.mean().rename('mean')
datamean

nvalues = dataclim.count().rename('count') / ntime
nvalues.max()
ntime

dsout = pd.DataFrame([datamean, nvalues]).T
dsout

dsout.to_csv('data/clim_sardara_data_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax))


