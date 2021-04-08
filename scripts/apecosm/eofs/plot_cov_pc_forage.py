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
lonmax = -150
latmax = 30
step = 10

mesh = xr.open_dataset('../../data/mesh_mask_eORCA1_v2.2.nc')
depth = mesh['gdept_1d'].values[0]

lat = mesh['gphit'].values[0]
lon = mesh['glamt'].values[0]

ilat, ilon = np.nonzero(np.abs(lat) <= latmax)
ilat = slice(ilat.min(), ilat.max() + 1)
lon0 = lon[ilat, :].mean(axis=0)
print(lon0.shape)

ilon = np.nonzero((lon0 >= -150 - 2) & (lon0 <= -150 + 2))[0]
ilon = slice(ilon.min(), ilon.max() + 1)
lat = np.mean(lat[:, ilon], axis=-1)
print(lat.shape)
print(depth.shape)

idepth = np.nonzero(depth <= 200)[0]
ilat = np.nonzero(np.abs(lat) <= latmax)[0]

classes = ['0-3cm', '3cm-20cm', '20cm-90cm', '90cm-200cm']
comm = ['Epi.', 'Mig.', 'Meso.']
dn = ['Day', 'Night']

toplot = xr.open_dataset("data/cov_forage_pc.nc")
toplot = toplot.isel(y=ilat, z=idepth)
toplot = toplot['cov'].values # com, size, eof, dn, y, depth
ncom, nbins, neof, nd, ny, nz = toplot.shape 
print(toplot.shape)

depth = depth[idepth]
lat = lat[ilat]

ncom = 1

for nd in range(2):

    fig = plt.figure(figsize=(10, 14))
    axgr = ImageGrid(fig, 111,  nrows_ncols=(3, 2), axes_pad=(0.5, 0.4), label_mode='L', cbar_mode='each', cbar_pad=0.05, aspect=False)
    cbar_axes = axgr.cbar_axes
    axout = list(enumerate(axgr))
    axout = [p[1] for p in axout]

    cpt = 0
    for c in range(ncom):
        for s in range(nbins):
            temp2 = toplot[c, s, :, nd, :, :]  # eof, y, x
            cmax = np.percentile(np.abs(temp2), 99.5)
            for e in range(2):
                ax = axout[cpt]
                temp = temp2[e, :, :]

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

                title = 'Size=%s, EOF %d' %(classes[s], e + 1) 
                ax.set_title(title)

                cpt += 1

    
    plt.suptitle(dn[nd], y=0.94)
    plt.savefig('cov_pc_merid_profile_dn_%d_%.f' %(nd, lonmax), bbox_inches='tight')
