import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

proj = ccrs.PlateCarree(central_longitude=180)
proj2 = ccrs.PlateCarree()

mesh = xr.open_dataset("/Users/Nicolas/Work/sent/apecosm/ORCA1/mesh_mask_eORCA1_v2.2.nc")
lonf = mesh['glamf'].values[0]
latf = mesh['gphif'].values[0]

data = xr.open_dataset("data/annual_eof_oni_annual_density.nc")
eof = data['eofmap'].to_masked_array()
pc = data['eofpc'].values
var = data['eofvar'].values * 100
time = data['time'].values

ncom, nbins, neof = var.shape
coms = ['Epi.']
sizes = ['3cm', '20cm', '90cm']

for c in range(1):
    for s in range(nbins):
        for e in range(2):

            plt.figure()
            ax = plt.subplot(2, 1, 1, projection=proj)
            temp = np.abs(eof[c, s, e, 1:, 1:])
            clim = temp.max()
            
            cs = ax.pcolormesh(lonf, latf, eof[c, s, e, 1:, 1:], transform=proj2, cmap=plt.cm.RdBu_r)
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.COASTLINE)
            ax.set_ylim(-40, 40)
            ax.set_xlim(-60, 130)
            cs.set_clim(-clim, clim)
            ax.set_title('%s, %s, EOF %d (%.2f' %(coms[c], sizes[s], e + 1, var[c, s, e]) + '%)')
            plt.colorbar(cs, orientation='horizontal')

            ax = plt.subplot(2, 1, 2)
            plt.plot(time, pc[c, s, e, :], color='black')
            ax.set_xlim(time.min(), time.max())
            ax.set_ylim(-3, 3)

            plt.savefig('annual_eof_com_%d_size_%s_eof_%d' %(c, s, e), bbox_inches='tight')
