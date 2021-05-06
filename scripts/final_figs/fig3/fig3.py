import xarray as xr
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
import matplotlib.pyplot as plt
import numpy as np
import envtoolkit.map as amap
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

mesh = xr.open_dataset('../../data/mesh_mask_eORCA1_v2.2.nc')
lat = mesh['gphit'].values[0]
lon = mesh['glamt'].values[0]

ilat, ilon = np.nonzero(np.abs(lat) <= 2)
ilat = slice(ilat.min(), ilat.max() + 1)

lon = np.mean(lon[ilat, :], axis=0)
z1d = mesh['gdept_1d'].values[0]

print(lon.shape)
print(z1d.shape)
zmax = 600
iz = np.nonzero(np.abs(z1d) <= zmax)[0]

#axgr = ImageGrid(fig, 111, nrows_ncols=(3, 1), ngrids=None, direction='row', axes_pad=(0.05, 0.3), share_all=False, aspect=True, label_mode='L', cbar_mode='each', cbar_location='right', cbar_pad='5%', cbar_size='5%', cbar_set_cax=True, axes_class=None)

nlevels = 5

for v in ['thetao', 'O2']:

    print(v)

    data = xr.open_dataset('data/equatorial_%s_mean.nc' %v)
    data['x'] = lon
    data = amap.lonflip('x', data)
    datam = data['%s_mean' %v].values  # x, z

    data = xr.open_dataset('data/zonal_monthly_covariance_monthly_%s_enso.nc' %v) 
    data['x'] = lon
    data = amap.lonflip('x', data)
    cov = data['covariance'].values  # x, z

    lonf = data['x'].values
    print(lon)

    cov = np.ma.masked_where(cov == 0, cov)

    plt.figure()
    ax = plt.gca()
    cs = plt.pcolormesh(lonf, -z1d[iz], cov.T[iz, :])
    cl = plt.contour(lonf, -z1d[iz], datam.T[iz, :], nlevels, linewidths=0.5, colors='k')
    ax.set_ylim(-500, 0)
    plt.colorbar(cs)

    ax.set_xlim(130, 300)
    labels = ['150', '180', '-150', '-120', '-90', '-60']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    labels = [_east_west_formatted(float(l)) for l in labels]

    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)

    plt.savefig('toto')

