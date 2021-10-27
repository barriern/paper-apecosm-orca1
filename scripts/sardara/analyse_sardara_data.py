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
import re
from mpl_toolkits.axes_grid1 import ImageGrid
plt.rcParams['text.usetex'] = False

pacific_latmax = 5
pacific_lonmin = 160
pacific_lonmax = -90
# -

data = pd.read_csv('data/processed_sardara_data_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax), index_col=0)
data

output_data = data.groupby(['time_start', 'species']).sum().loc[:, ['value']]
output_data

date = output_data.index.to_frame().loc[:, 'time_start']
date = np.unique(date.values)
date

regex = re.compile('([0-9][0-9][0-9][0-9])\-([0-9][0-9])\-([0-9][0-9])')
yymmdd = np.array([regex.match(d).groups() for d in date]).astype(int)
yymmdd.shape

ymin = yymmdd[:, 0].min()
ymax = yymmdd[:, 0].max()
ymin, ymax

datestr = ['%.4d-%.2d' %(y,m) for y in range(ymin, ymax+1) for m in range(1, 13)]
datestr

nyears = len(range(ymin, ymax + 1))
nyears

skj = output_data.loc[(slice(None, None), 'SKJ'), 'value']
skj

yft = output_data.loc[(slice(None, None), 'YFT'), 'value']
yft


def build_array(df):
    ddd = [l[0] for l in df.index.to_list()]
    output = np.zeros(nyears * 12, dtype=float)
    cpt = 0
    for y in range(ymin, ymax + 1):
        for m in range(1, 13):
            key = '%.4d-%.2d-%.2d' %(y, m, 1)
            if key in ddd:
                print(key)
                output[cpt] = df.loc[key].values[0]
            cpt += 1 
    output = np.ma.masked_where(output==0, output)
    return output


array_yft = build_array(yft)

array_skj = build_array(skj)

# +
time = np.arange(nyears * 12)
fig = plt.figure(figsize=(12, 12))

ax1 = plt.subplot(3, 1, 1)
ax1.plot(time, array_skj)
ax1.set_title('SKJ')
ax1.grid(True, 'both')
ax1.set_ylabel('MT')

ax2 = plt.subplot(3, 1, 2, sharex=ax1)
ax2.plot(time, array_yft)
ax2.set_title('YFT')
#ax2.xaxis.set_visible(False)
ax2.grid(True, 'both')
ax2.set_ylabel('MT')

ax = plt.subplot(3, 1, 3, sharex=ax1)
ax.plot(time, array_yft + array_skj)
stride = 12 * 5
t = ax.set_xticks(time[::stride])
l = ax.set_xticklabels(datestr[::stride])
ax.grid(True)
ax.set_ylabel('MT')

plt.savefig('captures_sardara_latmax_%d_lonmin_%d_lonmax_%d.png' %(pacific_latmax, pacific_lonmin, pacific_lonmax), bbox_inches='tight')
# -


