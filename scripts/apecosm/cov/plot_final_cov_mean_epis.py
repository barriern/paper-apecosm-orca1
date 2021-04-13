import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import string
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.ticker as mticker
from cycler import cycler

colcyc =  cycler(color=["cyan", 'magenta', 'gold', 'plum', 'green', 'orange'])
plt.rc('axes', prop_cycle=colcyc)
plt.rcParams['lines.markersize'] = 10

def plot_domain(ax, cpt):
    if(False):
    #if cpt == 1:
        _plot_domain(ax, lonmin0, latmin0, label='Dom 1')
        _plot_domain(ax, lonmin1, latmin1, label='Dom 2')
        #ax.legend()
    #elif cpt == 3:
        _plot_domain(ax, lonmin2, latmin2, label='Dom 3')
        _plot_domain(ax, lonmin3, latmin3, label='Dom 4')
        #ax.legend()
    #elif cpt == 5:
        _plot_domain(ax, lonmin4, latmin4, label='Dom 5')
        _plot_domain(ax, lonmin5, latmin5, label='Dom 6')
        #ax.legend()

def _plot_domain(axxx, lonmin, latmin, label=''):

    lonout = [lonmin]
    latout = [latmin]

    proj = ccrs.PlateCarree(central_longitude=0)
    proj2 = ccrs.PlateCarree(central_longitude=180)
    axxx.plot(lonout, latout, transform=proj, zorder=1000, marker='.', label=label)

lonmin0, latmin0 = 172.5, 4.5
lonmin1, latmin1 = -150.5, -4
lonmin2, latmin2 = 165.5, 8
lonmin3, latmin3 = -175.5, -4

lonmin4, latmin4 = -178.5, 3
lonmin5, latmin5 = 145.5, 5

# Define the settings of the grid 
def set_grid(cpt):

    gl = ax.gridlines(**dictgrid)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator([150, 180, -150, -120, -90, -60])

# converts to recover the value of the ax to used
def get_cpt(cpt):
    conversion = [1, 3, 5, 2, 4, 6]
    #return cpt - 1
    return conversion[cpt - 1] - 1

# recovers the letters in lower case
letters = list(string.ascii_lowercase)

########################################################## shared settings

# parameters to display the grid for lon/klat
dictgrid = {'crs':ccrs.PlateCarree(central_longitude=0), 'draw_labels':True, 'linewidth':0.5, 'color':'gray', 'alpha':0.5, 'linestyle':'--'}

# parameters to display the text
dicttext = dict(boxstyle='round', facecolor='lightgray', alpha=1)
lontext = 120
lattext = 30

# length to display
ilength = [14, 45, 80]

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree(central_longitude=0)

########################################################### reading mesh mask
dirin = '/home/barrier/Work/apecosm/ORCA1/final_figs/data'
dirin = '/home/barrier/Work/scientific_communication/articles/ongoing/paper-apecosm-orca1/scripts/data'
dirin = '../../data'

mesh = xr.open_dataset('%s/mesh_mask_eORCA1_v2.2.nc' %dirin)
mesh = mesh.isel(t=0, z=0)
tmask = mesh['tmask'].to_masked_array()
lon = mesh['glamt'].values
lat = mesh['gphit'].values
lonf = mesh['glamf'].values
latf = mesh['gphif'].values

lon[lon < 0] += 360

# reading length / weight step
const = xr.open_dataset('%s/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc' %dirin)
wstep = const['weight_step'].values
length = const['length'].values * 100
nlength = len(length)

########################################################### init figures
fig = plt.figure(figsize=(12, 8))
axes_class = (GeoAxes, dict(map_projection=proj))
axgr = AxesGrid(fig, 111,  axes_class=axes_class, nrows_ncols=(3, 2), axes_pad=(1, 0.4), label_mode='', cbar_mode='each', cbar_pad=0.05)
cbar_axes = axgr.cbar_axes
axout = list(enumerate(axgr))
axout = [p[1] for p in axout]
cpt = 1


########################################################### processing mean maps
dirin = '../../correlation_maps'
dirin = 'data/'

data = xr.open_dataset("%s/final-runs_yearly_mean_OOPE.nc" %dirin)
data = data.isel(community=0)
cov = data['OOPE'].to_masked_array()  # lat, lon, w
cov = np.squeeze(cov)

print(cov.shape)  # lat, lon, com, w

cov = np.ma.masked_where(cov == 0, cov)
cov = cov * wstep

'''
covout = np.zeros((332, 362, 3))

imin = 0
for p in range(3):
    print("++++++++++++++++++++++++++ ", p)
    imax = ilength[p]
    sss = slice(imin, imax)
    covout[:, :, p] = np.sum(cov[:, :, sss], axis=-1)
    if(p == 2): break
    imin = ilength[p]
    imax = ilength[p + 1]

covout = np.ma.masked_where(cov.mask, covout)
'''

cov = cov[:, :, ilength]
cov = np.log10(cov)
length = length[ilength]


latmask = (np.abs(lat) <= 40)
lonmask = (lon >= 130) & (lon <= 300)
dommask = latmask & lonmask
dommask = dommask[1:, 1:]

for p in range(3):

    temp = cov[1:, 1:, p]

    iok = np.nonzero((temp.mask == False) & (dommask == True))
    cmax = np.percentile(temp[iok], 99)
    cmin = np.percentile(temp[iok], 1)

    ax = axout[get_cpt(cpt)]
    cs = ax.pcolormesh(lonf, latf, temp, transform=proj2, cmap=plt.cm.jet, shading='auto')
    ax.add_feature(cfeature.LAND, zorder=100, color='lightgray')
    ax.add_feature(cfeature.COASTLINE, zorder=101)
    cb = cbar_axes[get_cpt(cpt)].colorbar(cs)
    cb.set_label('Log10 Biom. dens. ($J.m^{-2}$)')
    cs.set_clim(cmin, cmax)
    ax.set_title('Length = %.2e cm' %length[p])
    ax.text(lontext, lattext, letters[cpt - 1] + ")", ha='center', va='center', transform=proj, bbox=dicttext)
    set_grid(cpt)
    plot_domain(ax, get_cpt(cpt))
    ax.set_ylim(-40, 40)
    ax.set_xlim(-60, 130)
    #ax.vlines(-150 + 360, -30, 30)
    #ax.vlines(-100, -30, 30)
    #ax.vlines(-100, -30, 30)

    cpt += 1

########################################################### processing covariance
data = xr.open_dataset("%s/final-runs_covariance_yearly_enso_epis_OOPE.nc" %dirin)
cov = data['covariance'].to_masked_array()  # lat, lon, w
cov = np.ma.masked_where(cov == 0, cov)
cov = cov * wstep
cov = cov.T  # w, lon, lat

cov = cov[ilength]

for p in range(3):

    temp = cov[p].T
    temp = temp[1:, 1:]

    iok = np.nonzero((temp.mask == False) & (dommask == True))
    perc = np.percentile(np.abs(temp[iok]), 99)

    ax = axout[get_cpt(cpt)]
    cs = ax.pcolormesh(lonf, latf, temp, transform=proj2, cmap=plt.get_cmap('RdBu_r'))

    ax.add_feature(cfeature.LAND, zorder=100, color='lightgray')
    ax.add_feature(cfeature.COASTLINE, zorder=101)
    cb = cbar_axes[get_cpt(cpt)].colorbar(cs)
    cs.set_clim(-perc, perc)
    cb.set_label('Biom. Nino. cov. ($J.m^{-2}$)')
    ax.set_ylim(-40, 40)
    ax.set_xlim(-60, 130)
    ax.set_title('Length = %.2e cm' %length[p])
    set_grid(cpt)
    
    ax.text(lontext, lattext, letters[cpt - 1] + ")", ha='center', va='center', transform=proj, bbox=dicttext)
    plot_domain(ax, get_cpt(cpt))

    cpt += 1

plt.savefig('covariance_mean_epi.png', bbox_inches='tight')
