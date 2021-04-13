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

dirout = './'
ilon_glob = slice(58, 228)
ilat_glob = slice(140, 233)

# Load the mesh mask
mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc")
mesh = mesh.isel(y=ilat_glob, x=ilon_glob)
tmask = np.squeeze(mesh['tmask'].values[0, 0])
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values
surf = e1t * e2t
surf = np.squeeze(surf)
weights = surf

latmax = 20

tmask = xr.open_dataset('eof_mask_%d.nc' %latmax)
tmask = tmask.isel(y=ilat_glob, x=ilon_glob)
tmask = tmask['mask'].values
print(tmask.shape)
print(weights.shape)
ilat, ilon = np.nonzero(tmask == 1)

weights_tot = np.sum(weights[ilat, ilon])
weights[ilat, ilon] /= weights_tot
print(np.sum(weights[ilat, ilon]))
weights[ilat, ilon] = np.sqrt(weights[ilat, ilon])

data = xr.open_mfdataset('/home1/scratch/nbarrier/*OOPE*nc')
dens = data['OOPE'].to_masked_array() 
print(dens.shape)
mask = dens.mask
time = data['time']
clim, dens = ts.get_monthly_clim(dens)  # time, y, x, w
del(clim)
dens = np.ma.masked_where(mask == True, dens)
del(mask)

dens = np.transpose(dens, (3, 0, 1, 2))   # w, time, y, x
nbins, ntime, nlat, nlon = dens.shape

neofs = 2

eofmap = np.zeros((nbins, neofs, nlat, nlon))
eofts = np.zeros((nbins, neofs, ntime))
eofvar = np.zeros((nbins, neofs))

for s in range(nbins):

    print(dens.shape)
    print(ilat.shape)
    temp = dens[s, :, ilat, ilon].T  # ntime, ndims
    print('temp ', temp.shape)
    iok = np.nonzero(temp[0].mask == False)
    temp[:, iok] = sig.detrend(temp[:, iok].T).T

    print(weights[ilat, ilon].shape)
    solver = Eof(temp, weights=weights[ilat, ilon])
    # EOFs are multiplied by the square - root of
    # their eigenvalues ( units are in Pa )
    #maps = solver.eofs(eofscaling=2, neofs=neofs)
    maps = solver.eofsAsCovariance(neofs=neofs)

    eofmap[s, :, ilat, ilon] = maps.T
    eofts[s, :, :] = solver.pcs(pcscaling=1, npcs=neofs).T
    eofvar[s, :] = solver.varianceFraction(neigs=neofs)

eofmap = eofmap * tmask[np.newaxis, np.newaxis, :, :]
eofmap = np.ma.masked_where(eofmap == 0, eofmap)

dirout = 'data/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

dsout = xr.Dataset({'eofmap':(['bins', 'eof', 'y', 'x'], eofmap), 
                    'eofpc':(['bins', 'eof', 'time'], eofts),
                    'eofvar':(['bins', 'eof'], eofvar)})
dsout['time'] = (['time'], time)
#dsout['sizes'] = (['sizes'], bins)
dsout.attrs['creation_date'] = str(date.today())
dsout.attrs['script'] = os.path.realpath(__file__)
dsout.attrs['ilon'] = str(ilon_glob)
dsout.attrs['ilat'] = str(ilat_glob)
dsout.attrs['description'] = 'EOF for densities by class'
dsout.to_netcdf('%s/eof_full_density_%d.nc' %(dirout, latmax))
