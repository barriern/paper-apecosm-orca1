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

# # Second processing of the Sardara file.
#
# This script aims at providing a new processing to the Sardara file. Starting from the data file extracted for the Pacific, this script aims at:
# - Extract monthly data and discard yearly and 3-monthly data
# - Aggregate the catches that share the same area, starting date, gear and species.

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
pacific_lonmin = 120
pacific_lonmax = -60
# -

# First, the pre-processed file is loaded:

data = pd.read_csv('data/processed_sardara_data_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax), index_col=0)
data


# Then, we convert the time-start column into a dates using regular expressions. First, we create a function that does that:

def make_dates(data, arg):
    dates = data.loc[:, 'time_' + arg]
    regex = re.compile('([0-9][0-9][0-9][0-9])\-([0-9][0-9])\-([0-9][0-9])')
    yymmdd = np.array([regex.match(d).groups() for d in dates]).astype(int)
    yy = yymmdd[:, 0]
    mm = yymmdd[:, 1]
    dd = yymmdd[:, 2]
    dictout = {'year_' + arg:yy, 'month_' + arg: mm, 'day_' + arg: dd, 'date_' + arg: 10000*yy + mm * 100 + dd}
    data = data.assign(**dictout)
    return data


# The years and months are added to the file:

data = make_dates(data, 'start')
data = make_dates(data, 'end')
data

# Now we remove data which are not monthly ones:

test = data.loc[:, 'month_end'] - data.loc[:, 'month_start']
np.unique(test)

data = data[test == 0]
data

# Now, we sum the catches that have the same dates, geographic area and species (this is to take into account that different countries may fish on the same areas at the same dates).

dataint = data.groupby(['time_start', 'geographic_identifier', 'species', 'gear']).sum().loc[:, 'value']
dataint

dataint.to_csv('data/processed_sardara_data_L2_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax))
