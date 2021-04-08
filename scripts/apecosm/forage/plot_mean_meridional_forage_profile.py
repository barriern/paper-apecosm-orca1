import string
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.mpl.gridliner as gridliner

lonmax = -150
latmax = 20
off = 5

_DEGREE_SYMBOL = u'\u00B0'

def _north_south_formatted(latitude, num_format='g'):
    fmt_string = u'{latitude:{num_format}}{degree}{hemisphere}'
    return fmt_string.format(latitude=abs(latitude), num_format=num_format,
                             hemisphere=gridliner._lat_heimisphere(latitude),
                             degree=_DEGREE_SYMBOL)

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
letters = list(string.ascii_lowercase)
lontext = -latmax + off
depthtext = -175

lontext2 = latmax - off
depthtext2 = -175

classes = ['0-3cm', '3cm-20cm', '20cm-90cm', '90cm-200cm']
comm = ['Epi.', 'Mig.', 'Meso.']
dn = ['Day', 'Night']

#mesh = xr.open_dataset('/home/barrier/Work/apecosm/ORCA1/DATA_APE_ORCA1/pisces/data/corrected_mesh_mask_eORCA1_v2.2.nc')
mesh = xr.open_dataset('/Users/Nicolas/Work/sent/apecosm/ORCA1/corrected_mesh_mask_eORCA1_v2.2.nc')
mesh = mesh.isel(t=0)
depth = mesh['gdept_1d'].values
e3t = mesh['e3t_1d'].values
lon = mesh['glamt'].values
lat = mesh['gphit'].values

idepth = np.nonzero(depth <= 1000)[0]
depth = depth[idepth] 
e3t = e3t[idepth]

data = xr.open_dataset("data/meridional_%.f_forage.nc" %lonmax)
data = data.isel(depth=idepth)
lat = data['y'].values
forage = data['FORAGE'].to_masked_array()
forage = np.squeeze(forage)

ilat = np.nonzero(np.abs(lat) <= latmax)[0]

lat = lat[ilat]
forage = forage[:, ilat, :, :]  # dn, lat, depth, w

fig = plt.figure(figsize=(10, 14))

axgrid = AxesGrid(fig, 111, nrows_ncols=(4, 2), axes_pad=(0.9, 0.4), cbar_mode="each", aspect=False, cbar_pad=0.05)
cax = axgrid.cbar_axes

darray = [0, 1, 0, 1, 0,1, 0, 1]
sarray = [0, 0, 1, 1, 2, 2, 3, 3]

ccc = [60, 25, 10, 1.5]

for i, ax in enumerate(axgrid):

    d = darray[i]
    s = sarray[i]

    temp = forage[d, :, :, s].T
    temp = np.ma.masked_where(temp == 0, temp)
    itemp = np.nonzero((temp.mask == False) & (temp != 0))
    perc = 0.5
    cmin = np.percentile(temp[itemp], 0.5)
    cmax = np.percentile(temp[itemp], 100 - perc)
    cmin = 0 
    cmax = ccc[s]

    #cmin = temp[itemp].min()
    #cmax = temp[itemp].max()

    cs = ax.pcolormesh(lat, -depth, temp, cmap=plt.cm.jet, shading='auto')
    cs.set_clim(cmin, cmax)
    cb = cax[i].colorbar(cs)
    cl = ax.contour(lat, -depth, temp, 11, colors='k', linewidths=0.5)
    #cb.add_lines(cl)
    cb.set_label_text('$J.m^{-3}$')
    ax.set_title(classes[s])
    ax.set_ylim(-200, 0)
    ax.text(lontext, depthtext, letters[i] + ")", bbox=dicttext, ha='center', va='center')
    ax.text(lontext2, depthtext2, dn[d], bbox=dicttext, ha='center', va='center')
    
    if i in [0, 2, 4, 6]:
        ax.set_ylabel('Depth (m)')
    
    if i in [6, 7]:
        ax.set_xlabel('Lat')
    
    ax.set_xlim(-latmax, latmax)
    step = 5
    labels = np.arange(-latmax, latmax + step, step)
    xticks = np.array([float(l) for l in labels])
    labels = [_north_south_formatted(float(l)) for l in labels]

    axgrid[i].set_xticklabels(labels)
    axgrid[i].set_xticks(xticks)
    ax.grid(linestyle='--', linewidth=0.5)
        

plt.savefig('forage_meridional_mean_%.f.png' %lonmax, bbox_inches='tight')

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
