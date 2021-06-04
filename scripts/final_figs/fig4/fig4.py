import xarray as xr
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
import matplotlib.pyplot as plt
import numpy as np
import envtoolkit.map as amap
import cartopy.mpl.gridliner as gridliner
import string

dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)

letters = list(string.ascii_lowercase)

plt.rcParams['image.cmap'] = 'RdBu_r'
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['font.size'] = 12

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
zmax = 300
iz = np.nonzero(np.abs(z1d) <= zmax)[0]

z1d = z1d/ 1000.

titles = \
["Biomass dens. [0-3cm]", 
"Biomass dens. [3cm-20cm]", 
"Biomass dens. [20cm-90cm]", 
"Temperature", 
"Oxygen",
"Plankton"]

units = 3 * ['J/m3'] + ['Celcius'] + 2 * ['mmol/m3']

######### Reading Mean Pisces fields
dataglob = xr.open_mfdataset('data/equatorial_*_mean.nc')
dataglob['PLK_mean'] = dataglob['GOC_mean'] + dataglob['PHY2_mean'] + dataglob['ZOO_mean'] + dataglob['ZOO2_mean']
dataglob['CHL_mean'] = dataglob['NCHL_mean'] + dataglob['DCHL_mean']

# converting from atl to pac systems
dataglob['x'] = lon
dataglob = amap.lonflip('x', dataglob)

########## Reading mean Forage and convert from Pac. to Atl
forage_mean = xr.open_dataset('data/equatorial_forage_mean.nc')
forage_mean['x'] = lon
forage_mean = amap.lonflip('x', forage_mean)
forage_mean = forage_mean['forage_mean'].values # size, com, z, x, dn
forage_mean = forage_mean[:3, 0, :, :, 0]

########## Reading cov Forage and convert from Pac. to Atl
forage_cov = xr.open_dataset('data/zonal_monthly_covariance_yearly_enso.nc')
forage_cov['x'] = lon
forage_cov = amap.lonflip('x', forage_cov)
lonf = forage_cov['x'].values
forage_cov = forage_cov['covariance'].values  # sizes, com, z, x, dn
forage_cov = forage_cov[:3, 0, :, :, 0]
forage_cov = np.ma.masked_where(forage_cov == 0, forage_cov)
forage_mean = np.ma.masked_where(forage_mean == 0, forage_mean)
    
fig = plt.figure(figsize=(12, 12))

axgr = ImageGrid(fig, 111, nrows_ncols=(3, 2), ngrids=None, direction='row', axes_pad=(0.8, 0.4), share_all=False, aspect=False, label_mode='L', cbar_mode='each', cbar_location='right', cbar_pad='5%', cbar_size='5%', cbar_set_cax=True, axes_class=None)
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

nlevels = 5

cpt = 0

ccc = [15, 6, 1, 1.5, 15, 0.45]

axorder = [0, 2, 4, 1, 3, 5]

lontext = 140
lattext = -0.25

for s in range(3):

    ilon = np.nonzero((lonf >= 130) & (lonf <= 300) & (np.ma.getmaskarray(forage_cov[s, 0, :]) == False))[0]

    temp = forage_cov[s, iz, :][:, ilon]
    datam = forage_mean[s, iz, :][:, ilon]
    print(datam.shape)
    axi = axorder[cpt]
    ax = axout[axi]
    cs = ax.pcolormesh(lonf[ilon], -z1d[iz], temp)
    cl = ax.contour(lonf[ilon], -z1d[iz], datam, nlevels, linewidths=1, colors='k')
    ax.clabel(cl, cl.levels, fmt="%.f", manual=False)
    cs.set_clim(-ccc[cpt], ccc[cpt])
    ax.text(lontext, lattext, letters[cpt] + ")", ha='center', va='center', bbox=dicttext)

    ax.set_xlim(130, 300)
    labels = ['150', '180', '-150', '-120', '-90', '-60']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    labels = [_east_west_formatted(float(l)) for l in labels]

    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)
    ax.set_title(titles[cpt])

    cbar = axgr.cbar_axes[axi].colorbar(cs)
    cbar.set_label_text(units[cpt])
    ax.grid(True, linestyle='--')
    ax.set_ylabel('Depth (km)')

    cpt += 1
    

for v in ['thetao', 'O2', 'PLK']:

    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ", v)

    datam = dataglob['%s_mean' %v].values  # x, z

    data = xr.open_dataset('data/zonal_monthly_covariance_monthly_%s_enso.nc' %v) 
    data['x'] = lon
    data = amap.lonflip('x', data)
    cov = data['covariance'].values  # x, z
    lonf = data['x'].values

    print(lonf.shape, cov.shape)

    cov = np.ma.masked_where(cov == 0, cov)  # lon, depth
    ilon = np.nonzero((lonf >= 130) & (lonf <= 300) & (np.isnan(cov[:, 0]) == False))[0]

    temp = cov[ilon, :][:, iz].T

    axi = axorder[cpt]
    ax = axout[axi]
    cs = ax.pcolormesh(lonf[ilon], -z1d[iz], temp)
    cl = ax.contour(lonf[ilon], -z1d[iz], datam.T[iz, :][:, ilon], nlevels, linewidths=1, colors='k')
    ax.clabel(cl, cl.levels, fmt="%.1f", manual=False)
    #ax.set_ylim(-z1d[iz].min(), 0)
    cs.set_clim(-ccc[cpt], ccc[cpt])
    ax.text(lontext, lattext, letters[cpt] + ")", ha='center', va='center', bbox=dicttext)

    ax.set_xlim(130, 300)
    labels = ['150', '180', '-150', '-120', '-90', '-60']
    xticks = np.array([float(l) for l in labels])
    xticks[xticks < 0] += 360
    labels = [_east_west_formatted(float(l)) for l in labels]

    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)
    ax.set_title(titles[cpt])
    ax.grid(True, linestyle='--')

    cbar = axgr.cbar_axes[axi].colorbar(cs)
    cbar.set_label_text(units[cpt])
    ax.set_ylabel('Depth (km)')

    cpt += 1

plt.savefig('fig4.png', bbox_inches='tight')

