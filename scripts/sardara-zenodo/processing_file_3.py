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

# # Processing Sardara #3
#
# The aim here is to the polygon columns into lonmin/lonmax/latmin/latmax columnsm in order to facilitate the regridding. Besides, some polygon informations are missing and are discarded.

# +
import pandas as pd
import numpy as np
from datetime import datetime

data = pd.read_csv('processed_sardara_2.csv', index_col=0)
data
# -

polygon = data['geom_wkt'].values.astype(str)
iok = np.nonzero(polygon != 'nan')[0]

# +
lonmin = []
lonmax = []
latmin = []
latmax = []

for polstr in polygon[iok]:

    try:
        temp = polstr.replace('POLYGON((', '').replace('))', '')
        temp
    except AttributeError:
        print('Error with')
        print(polstr)
        

    pol = temp.split(',')
    pol

    output = []
    for p in pol:
        output.append(p.split(' '))
    output = np.array(output, dtype=float)
    output
    lonmin.append(output[:, 0].min())
    lonmax.append(output[:, 0].max())
    latmin.append(output[:, 1].min())
    latmax.append(output[:, 1].max())
# -

dataout = data.iloc[iok, :]
dataout = dataout.assign(lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax)
dataout

dataout.to_csv('processed_sardara_3.csv')


