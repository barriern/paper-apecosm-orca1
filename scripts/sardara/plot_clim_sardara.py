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

# # Plotting of Sardara climatology

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
import matplotlib as mpl

pacific_latmax = 40
pacific_lonmin = 140
pacific_lonmax = -60

species = 'SKJ'

cmap = plt.cm.jet
# -

# ## Area extraction
#
# The area files are extracted:

areas = pd.read_csv('data/processed_sardara_areas_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax), index_col=0)
areas

# The cell resolution is also extracted and added to the area dataframe

areas = areas.assign(dx=abs(areas.loc[:, 'lonmin'] - areas.loc[:, 'lonmax']))
areas

# ## Loading the data file

# Then, the pre-processed file is loaded:

data = pd.read_csv('data/clim_sardara_data_latmax_%d_lonmin_%d_lonmax_%d.csv' %(pacific_latmax, pacific_lonmin, pacific_lonmax), index_col=[0, 1])
data

# Now we extract the data for all the areas and one species:

temp = data.loc[(slice(None, None), species), :]
temp

# Now we recover the area indexes for all the input data that we process:

indexes = [l[0] for l in temp.index]

# ## Sorting areas by resolution

# We extract the resolution of the area:

dxareas = areas['dx'].values
dxareas

# We loop over the different resolutions and store the dataframe in a dictionnary:

# +
output = {}
for r in np.unique(dxareas):
    print('Processing resolution ', r)
    # index of the areas with the given resolution
    iareas = np.nonzero(dxareas == r)[0]
    # we extract the codes of the areas with the given resolution
    temparea = areas.iloc[iareas, :].loc[:, 'code'].values
    test = np.array([t in temparea for t in indexes])
    if(np.all(test == False)):
        print('No areas with resolution ', r)
        continue
    
    iok = np.nonzero(test == True)[0]
    
    output[r] = temp.iloc[iok, :]

output.keys()
# -

# ## Normalisation by the surface

res = list(output.keys())
res.sort()
res

# +
for r in res:
    
    print('Processing res = ', r)
    
    dy = Rt * np.deg2rad(r)
    tempbis = output[r]
    tempbis

    iii = [l[0] for l in tempbis.index]
    
    norm = [] 
    
    for i in iii[:]:
        temparea = areas[areas['code'] == i]
        latcen = np.deg2rad(temparea['latcen'].values[0])
        dx = Rt * np.cos(np.deg2rad(latcen)) * np.deg2rad(r)
        norm.append(output[r].loc[i, 'mean'].values[0] / (dx * dy))
        
    output[r] = output[r].assign(norm_mean=norm)
    
output[r]
# -

# Now we search for minimum and maximum values to plot on the colorbar:

# +
mmm = []
MMM = []
for r in res:
    print(r)
    mval = output[r].loc[:, 'norm_mean'].min()
    MVAL = output[r].loc[:, 'norm_mean'].max()
    mmm.append(mval)
    MMM.append(MVAL)

mmm = np.min(mmm)
MMM = np.max(MMM)
MMM - mmm
# -

# Now we process loop over all the resolutions in reverse order:

# +
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude = 180))

Rt = 6371.009 # km

for r in res[::-1]:
    
    print(r)
    
    tempbis = output[r]
    tempbis

    indexes = [l[0] for l in tempbis.index]
    len(indexes)

    vmax = tempbis.max()
    vmax

    #Now, we draw the figure:

    projin = ccrs.PlateCarree()
    projout = ccrs.PlateCarree(central_longitude=180)
    toplot = 'norm_mean'

    cpt = 0
    nindex = len(indexes)
    for i in indexes[:]:

        if(cpt % 20 == 0):
            print(cpt, ' / ', nindex)

        cpt += 1

        temparea = areas[areas['code'] == i]

        xxx = np.array([temparea['lonmin'].values[0], temparea['lonmax'].values[0], temparea['lonmax'].values[0], temparea['lonmin'].values[0]])
        yyy = np.array([temparea['latmin'].values[0], temparea['latmin'].values[0], temparea['latmax'].values[0], temparea['latmax'].values[0]])
        
        tp = tempbis.loc[(i, species), toplot]
        #tp = (tp - mmm) / (MMM - mmm)
        tp = (np.log10(tp) - np.log10(mmm)) / (np.log10(MMM) - np.log10(mmm))
        
        col = cmap(tp)
        
        points = projout.transform_points(projin, xxx, yyy)
        xout = points[:, 0]
        yout = points[:, 1]
        xy = np.array([xout, yout]).T

        poly = mpl.patches.Polygon(xy, closed=True, facecolor=col)
        ax.add_artist(poly)

ax.coastlines()
ax.set_global()

ax = plt.axes([0.95, 0.2, 0.02, 0.6])
norm = mpl.colors.Normalize(vmin=np.log10(mmm), vmax=np.log10(MMM))
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
cb.set_label('%s catch (MT)' %species)
plt.savefig('%s_catch_resolution_%f.png' %(species, r), bbox_inches='tight')
# -


