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
# The aim of this script is to:
# - Remove the data with bad geographic identifiers (countries instead of codes)
# - Remove the data which are not in the Pacific Ocean, based on the comparison with geographic identifiers
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
pacific_lonmin = 120
pacific_lonmax = -60
# -

# ## Reading processed areas

areas = pd.read_csv('data/processed_sardara_areas_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax), index_col=0)
areas

areacodes = areas.loc[:, 'code'].values
areacodes

# ## Reading the CSV

data = pd.read_csv('data/global_catch_ird_level2.csv', sep=',')
data

# ## Removing wrong geographic identifier

# Now we remove areas which are not identified from integers, i.e when the identifier is a string and not a number.

area = data.loc[:, 'geographic_identifier'].values.astype(str)
area

pattern = '[0-9]+'
regex = re.compile(pattern)

test = np.array([regex.match(a) is not None for a in area])

datanew = data.iloc[test, :]
datanew

# Now we convert the identifier from string to integer.

datanew = datanew.assign(geographic_identifier=datanew.loc[:, 'geographic_identifier'].astype(int))
datanew

# Now we extract the data whose codes are in the Pacific domain. To do so, we compare the input geographic identifier with the ones of the Pacific area files.

test_index = np.array([v in areacodes for v in datanew['geographic_identifier'].values])
iok = np.nonzero(test_index == True)[0]
len(iok)

# Finally, we extract the data for these locations and we save the new file.

dataout = datanew.iloc[iok, :]
dataout

dataout.to_csv('data/processed_sardara_data_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax))
