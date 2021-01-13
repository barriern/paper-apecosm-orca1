import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import string
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes

# Define the settings of the grid 
def set_grid():

    gl = ax.gridlines(**dictgrid)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

# converts to recover the 
def get_cpt(cpt):
    conversion = [1, 3, 5, 2, 4, 6]
    return conversion[cpt - 1] - 1

letters = list(string.ascii_lowercase)

dirin = '/home/barrier/Work/apecosm/ORCA1/final_figs/data'

ilength = [14, 45, 80]

dictgrid = {'crs':ccrs.PlateCarree(central_longitude=180), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
dictcbar = {}
lontext = 120
lattext = 30

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

# reading mesh mask
mesh = xr.open_dataset('%s/mesh_mask_eORCA1_v2.2.nc' %dirin)
mesh = mesh.isel(t=0, z=0)
tmask = mesh['tmask'].to_masked_array()
lon = mesh['glamt'].values
lat = mesh['gphit'].values

# reading length / weight step
const = xr.open_dataset('%s/ORCA1_JRA_CO2_CYC4_ConstantFields.nc' %dirin)
wstep = const['weight_step'].values
length = const['length'].values * 100
nlength = len(length)

# init figures
fig = plt.figure(figsize=(12, 8))
axes_class = (GeoAxes,
              dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 2), axes_pad=(1, 0.4), label_mode='', cbar_mode='each', cbar_pad=0.05)
cbar_axes = axgr.cbar_axes
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]

cpt = 1

length = length[ilength]

# processing mean maps
dirin = '/home/barrier/Desktop/correlation_maps'

data = xr.open_dataset("%s/debugged_corr_mask_yearly_mean_OOPE.nc" %dirin)
cov = data['OOPE'].to_masked_array()  # lat, lon, w
cov = np.squeeze(cov)

cov = np.ma.masked_where(cov == 0, cov)
cov = cov * wstep
cov = np.log10(cov)

cov = cov[:, :, ilength]


for p in range(3):

    temp = cov[:, :, p]

    iok = np.nonzero(temp.mask == False)
    cmax = np.percentile(temp[iok], 99)
    cmin = np.percentile(temp[iok], 1)

    ax = axout[get_cpt(cpt)]
    cs = ax.pcolormesh(lon, lat, temp, transform=proj2, cmap=plt.cm.jet)
    ax.add_feature(cfeature.LAND, zorder=100, color='lightgray')
    ax.add_feature(cfeature.COASTLINE, zorder=101)
    cb = cbar_axes[get_cpt(cpt)].colorbar(cs, **dictcbar)
    cb.set_label_text('Log10 Biom. dens. ($J.m^{-2}$)')
    ax.set_ylim(-40, 40)
    ax.set_xlim(-60, 130)
    cs.set_clim(cmin, cmax)
    ax.set_title('Length = %.2e cm' %length[p])
    ax.text(lontext, lattext, letters[cpt - 1] + ")", ha='center', va='center', transform=proj, bbox=dicttext)
    set_grid()

    cpt += 1

# processing covariance
data = xr.open_dataset("%s/debugged_corr_mask_covariance_yearly_enso_epis_OOPE.nc" %dirin)
cov = data['covariance'].to_masked_array()  # lat, lon, w
cov = np.ma.masked_where(cov == 0, cov)
cov = cov * wstep
cov = cov.T  # w, lon, lat

cov = cov[ilength]

for p in range(3):

    temp = cov[p].T

    iok = np.nonzero(temp.mask == False)
    perc = np.percentile(np.abs(temp[iok]), 99)

    #ax = plt.subplot(sp0, sp1, get_cpt(cpt), projection=proj)
    ax = axout[get_cpt(cpt)]
    cs = ax.pcolormesh(lon, lat, temp, transform=proj2, cmap=plt.get_cmap('RdBu_r'))

    ax.add_feature(cfeature.LAND, zorder=100, color='lightgray')
    ax.add_feature(cfeature.COASTLINE, zorder=101)
    cb = cbar_axes[get_cpt(cpt)].colorbar(cs, **dictcbar)
    cs.set_clim(-perc, perc)
    cb.set_label_text('Biom. Nino. cov. ($J.m^{-2}$)')
    ax.set_ylim(-40, 40)
    ax.set_xlim(-60, 130)
    ax.set_title('Length = %.2e cm' %length[p])
    set_grid()
    
    ax.text(lontext, lattext, letters[cpt - 1] + ")", ha='center', va='center', transform=proj, bbox=dicttext)

    cpt += 1


plt.savefig('covariance_mean_epi.png', bbox_inches='tight')
