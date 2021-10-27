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
import pandas as pd
import numpy as np
import re
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from pandas.api.types import union_categoricals
plt.rcParams['text.usetex'] = False

pacific_latmax = 5
pacific_lonmin = 160
pacific_lonmax = -90
# -

data = pd.read_csv('data/processed_sardara_data_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax), index_col=0)
data

output_data = data.groupby(['time_start', 'species',  'catchtype', 'schooltype']).sum().loc[:, ['value']]
output_data

output_data.to_csv('data/analysed_sardara_data_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax))


