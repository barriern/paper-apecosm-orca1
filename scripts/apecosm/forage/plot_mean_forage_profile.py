import string
import xarray as xr
import matplotlib.pyplot as plt
import apecosm
import numpy as np
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.mpl.gridliner as gridliner

_DEGREE_SYMBOL = u'\u00B0'

def _east_west_formatted(longitude, num_format='g'):
    fmt_string = u'{longitude:{num_format}}{degree}{hemisphere}'
    if(longitude == 180):
        longitude = 0
    output = fmt_string.format(longitude=abs(longitude), num_format=num_format,
                             hemisphere=gridliner._lon_hemisphere(longitude),
                             degree=_DEGREE_SYMBOL)
    return output

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
letters = list(string.ascii_lowercase)
lontext = 150
depthtext = -175

lontext2 = 265
depthtext2 = -175

classes = ['0-3cm', '3cm-20cm', '20cm-90cm', '90cm-200cm']
comm = ['Epi.', 'Mig.', 'Meso.']
dn = ['Day', 'Night']

mesh = xr.open_dataset('/home/barrier/Work/apecosm/ORCA1/DATA_APE_ORCA1/pisces/data/corrected_mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0)
depth = mesh['gdept_1d'].values
e3t = mesh['e3t_1d'].values
lon = mesh['glamt'].values
lat = mesh['gphit'].values
latmax = 2

idepth = np.nonzero(depth <= 1000)[0]
depth = depth[idepth] 
e3t = e3t[idepth]

ilat, ilon = np.nonzero(np.abs(lat) <= latmax)
jmin, jmax = ilat.min(), ilat.max() + 1
ilat = slice(jmin, jmax)

lon = lon[ilat, :]
lon = np.mean(lon, axis=0)
lon0 = lon.copy()

lon[lon < 0] += 360

data = xr.open_dataset("data/yearly_mean_forage_year.nc")
data = data.isel(community=0, depth=idepth)
data['x'] = lon
data = data.sortby(data['x'])
lon = data['x'].values
forage = data['FORAGE'].to_masked_array()
forage = np.squeeze(forage)

ilon = np.nonzero((lon <= 300) & (lon >= 130))[0]

lon = lon[ilon]
forage = forage[:, ilon, :, :]

fig = plt.figure(figsize=(10, 14))

axgrid = AxesGrid(fig, 111, nrows_ncols=(4, 2), axes_pad=(0.5, 0.4), cbar_mode="each", aspect=False, cbar_pad=0.05)
cax = axgrid.cbar_axes

darray = [0, 1, 0, 1, 0,1, 0, 1]
sarray = [0, 0, 1, 1, 2, 2, 3, 3]

for i, ax in enumerate(axgrid):

    d = darray[i]
    s = sarray[i]

    temp = forage[d, :, :, s].T
    temp = np.ma.masked_where(temp == 0, temp)
    itemp = np.nonzero((temp.mask == False) & (temp != 0))
    cmin = np.percentile(temp[itemp], 1)
    cmax = np.percentile(temp[itemp], 99)

    cmin = temp[itemp].min()
    cmax = temp[itemp].max()

    cs = ax.pcolormesh(lon, -depth, temp, cmap=plt.cm.jet)
    cs.set_clim(cmin, cmax)
    cb = cax[i].colorbar(cs)
    #cl = ax.contour(lon, -depth, temp, 11, colors='k', linewidths=0.5)
    #cb.add_lines(cl)
    cb.set_label_text('$J.m^{-3}$')
    ax.set_title(classes[s])
    ax.set_ylim(-200, 0)
    ax.text(lontext, depthtext, letters[i] + ")", bbox=dicttext, ha='center', va='center')
    ax.text(lontext2, depthtext2, dn[d], bbox=dicttext, ha='center', va='center')
    
    if i in [0, 2, 4, 6]:
        ax.set_ylabel('Depth (m)')
    
    if i in [6, 7]:
        ax.set_xlabel('Lon')
    
    ax.set_xlim(130, 300)
    labels = ['150', '180', '-150', '-120', '-90', '-60']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    labels = [_east_west_formatted(float(l)) for l in labels]

    axgrid[i].set_xticklabels(labels)
    axgrid[i].set_xticks(xticks)
    ax.grid()
        

plt.savefig('forage_mean.png', bbox_inches='tight')

plt.close(fig)
        
"""
test = np.sum(forage * e3t[np.newaxis, np.newaxis, :, np.newaxis], axis=2)

plt.figure()

for s in range(4):
    plt.subplot(2, 2, s + 1)

    plt.plot(test[0, :, s], label='0')
    plt.plot(test[1, :, s], label='1')

plt.legend()

plt.savefig('test.png')
"""
