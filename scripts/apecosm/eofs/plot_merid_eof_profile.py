import xarray as xr
import sys
sys.path.append('../../nino')
from extract_nino import read_index
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import ImageGrid
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

def _north_south_formatted(latitude, num_format='g'):
    fmt_string = u'{latitude:{num_format}}{degree}{hemisphere}'
    return fmt_string.format(latitude=abs(latitude), num_format=num_format,
                             hemisphere=gridliner._lat_heimisphere(latitude),
                             degree=_DEGREE_SYMBOL)

dnino, nino = read_index(filename='data/index/oni.data')
ynino = dnino // 100
iok = np.nonzero((ynino <= 2018) & (ynino >= 1958))
nino = nino[iok]
dnino = dnino[iok]

latmax = 20
lonmax = -150
step = 5

data = xr.open_dataset('data/eof_merid_profile_%.f.nc' %lonmax)
lat = data['y'].values
depth = data['z'].values
idepth = np.nonzero(depth <= 200)[0]
data = data.isel(z=idepth)
depth = depth[idepth]

eof = data['eofmap'].values
pc = data['eofpc'].values
eofvar = data['eofvar'].values * 100
nd, ncom, nbins, neof, nx, nz = eof.shape
nbins = 4
ncom = 1

classes = ['0-3cm', '3cm-20cm', '20cm-90cm', '90cm-200cm']
comm = ['Epi.', 'Mig.', 'Meso.']
dn = ['Day', 'Night']

for nd in range(2):

    fig = plt.figure(figsize=(10, 14))
    axgr = ImageGrid(fig, 111,  nrows_ncols=(4, 2), axes_pad=(0.5, 0.4), label_mode='L', cbar_mode='each', cbar_pad=0.05, aspect=False)
    cbar_axes = axgr.cbar_axes
    axout = list(enumerate(axgr))
    axout = [p[1] for p in axout]

    cpt = 0
    for c in range(ncom):
        for s in range(nbins):
            temp = eof[nd, c, s]
            cmax = np.percentile(np.abs(temp), 99.5)
            for e in range(2):
                ax = axout[cpt]
                temp = eof[nd, c, s, e]

                corrcoef = np.corrcoef(pc[nd, c, s, e], nino)[0, 1]
                if(corrcoef < 0):
                    temp *= -1

                cs = ax.pcolormesh(lat, -depth, temp.T, cmap=plt.cm.RdBu_r)
                cs.set_clim(-cmax, cmax)
                cbar_axes[cpt].colorbar(cs)
                #cs = ax.contour(lat, -depth, temp.T, 11, colors='k', linewidths=0.5)

                ax.set_xlim(-latmax, latmax)
                labels = np.arange(-latmax, latmax + step, step)
                xticks = np.array([float(l) for l in labels])
                labels = [_north_south_formatted(float(l)) for l in labels]
                ax.set_xticks(xticks)
                ax.set_xticklabels(labels)
                ax.grid(True, linestyle='--', linewidth=0.5)
                ax.set_xlabel('Latitude')
                ax.set_ylabel('Depth')

                title = 'Size=%s, EOF %d (%.2f' %(classes[s], e + 1, eofvar[nd, c, s, e]) + r'\%)'
                ax.set_title(title)

                cpt += 1

    
    plt.suptitle(dn[nd], y=0.94)
    plt.savefig('eof_merid_profile_dn_%d_%.f' %(nd, lonmax), bbox_inches='tight')
