# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import matplotlib.pyplot as plt
import sys
sys.path.append('../nino')
from extract_nino import read_index
import numpy as np

datestart = 199601
dateend = 200112
# -

date, nino = read_index('../data/index/oni.data')

iok = np.nonzero((date >= datestart) & (date <= dateend))[0]
date[iok]

date = date[iok]
nino = nino[iok]

years = date // 100
months = (date - 100 * years).astype(int)
datestr = np.array(['%.4d-%.2d' %(y, m) for y, m in zip(years, months)])
datestr

ntime = len(datestr)
time = np.arange(ntime)
time

stride = 3
al = 0.5
fig = plt.figure(figsize=(14, 8))
ax = plt.axes()
ax.plot(time, nino, color='k', marker='.')
ax.fill_between(time, 0, nino, where=(nino > 0), interpolate=True, facecolor='FireBrick', alpha=al)
ax.fill_between(time, 0, nino, where=(nino < 0), interpolate=True, facecolor='SteelBlue', alpha=al)
ticks = ax.set_xticks(time[::stride])
ticklabels = ax.set_xticklabels(datestr[::stride], rotation=45, ha='right')
ax.grid(True)
ax.set_xlim(time.min(), time.max())
ax.set_ylim(-2.5, 2.5)
plt.savefig('nino.pdf', bbox_inches='tight')
#plt.close(fig)


