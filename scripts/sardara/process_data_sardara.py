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

# # Processing Sardara data
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

# ## Reading processed areas

areas = pd.read_csv('data/processed_sardara_areas_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax))
areas

areacodes = areas.loc[:, 'code'].values
areacodes

# ## Reading the CSV

data = pd.read_csv('data/global_catch_ird_level2.csv', sep=',')
data

# ## Removing wrong geographic identifier

# Now we remove areas which are not identified from integers

area = data.iloc[:, 4].values.astype(str)
area

pattern = '[0-9]+'
regex = re.compile(pattern)

test = np.array([(regex.match(a) is not None) for a in area])

datanew = data.iloc[test, :]
datanew

datanew = datanew.assign(geographic_identifier=datanew.loc[:, 'geographic_identifier'].astype(int))
datanew
