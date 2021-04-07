import numpy as np
from datetime import date
import os.path
import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt
import cartopy.crs as crs
from scipy import stats
import xarray as xr
import apecosm.ts as ts
import scipy.signal as sig

lonmax = -150

# Load the mesh mask
mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc")
mesh = mesh.isel(t=0)
e3t = mesh['e3t_1d'].values
z = mesh['gdept_1d'].values

idepth = np.nonzero(z <= 1000)[0]
z = z[idepth]
e3t = e3t[idepth]
e3t /= np.sum(e3t)
e3t = np.sqrt(e3t)

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data/'
data = xr.open_mfdataset('%s/meridional_%.f_forage_year_**.nc' %(dirin, lonmax))
time = data['time']
lat = data['y'].values
ilat = np.nonzero(np.abs(lat) <= 20)[0]
data = data.isel(y=ilat, depth=idepth)
lat = data['y'].values
dens = data['FORAGE'].to_masked_array()  # time, dn, y, depth, com, size

mask = np.ma.getmaskarray(dens)
clim, dens = ts.get_monthly_clim(dens) # time, dn, x, depth, comm, size
del(clim)

dens = np.ma.masked_where(mask == True, dens)
dens = np.transpose(dens, (1, 4, 5, 0, 2, 3))   # dn, com, size, time, x, z
ndn, ncom, nbins, ntime, nx, nz = dens.shape

neofs = 2

eofmap = np.zeros((ndn, ncom, nbins, neofs, nx, nz))
eofts = np.zeros((ndn, ncom, nbins, neofs, ntime))
eofvar = np.zeros((ndn, ncom, nbins, neofs))

for d in range(ndn):
    for c in range(ncom):
        for s in range(nbins):

            print(dens.shape)
            temp = dens[d, c, s, :, :, :]  # ntime, nx, nz
            print(temp.shape)

            temp = sig.detrend(temp.T).T
            print(temp.shape)
            print(e3t.shape)

            solver = Eof(temp, weights=e3t)
            # EOFs are multiplied by the square - root of
            # their eigenvalues ( units are in Pa )
            #maps = solver.eofs(eofscaling=2, neofs=neofs)
            maps = solver.eofsAsCovariance(neofs=neofs)  # neof, x, z
            print(maps.shape)

            eofmap[d, c, s, :, :, :] = maps
            eofts[d, c, s, :, :] = solver.pcs(pcscaling=1, npcs=neofs).T
            eofvar[d, c, s, :] = solver.varianceFraction(neigs=neofs)

dirout = 'data/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

dsout = xr.Dataset({'eofmap':(['nd', 'com', 'bins', 'eof', 'x', 'z'], eofmap), 
                    'eofpc':(['nd', 'com', 'bins', 'eof', 'time'], eofts),
                    'eofvar':(['nd', 'com', 'bins', 'eof'], eofvar)})
dsout['time'] = (['time'], time)
dsout['y'] = (['y'], lat)
dsout['z'] = (['z'], z)
#dsout['sizes'] = (['sizes'], bins)
dsout.attrs['creation_date'] = str(date.today())
dsout.attrs['script'] = os.path.realpath(__file__)
dsout.attrs['description'] = 'EOF for densities by class'
dsout.to_netcdf('%s/eof_merid_profile_%.f.nc' %(dirout, lonmax))
