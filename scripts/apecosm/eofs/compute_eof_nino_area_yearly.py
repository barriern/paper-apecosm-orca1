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

# Load the mesh mask
mesh = xr.open_dataset("/home/datawork-marbec-pmod/forcings/APECOSM/ORCA1_HINDCAST/corrected_mesh_mask_eORCA1_v2.2.nc")
tmask = np.squeeze(mesh['tmask'].values[0, 0])
e1t = mesh['e1t'].values
e2t = mesh['e2t'].values
surf = e1t * e2t
surf = np.squeeze(surf)
weights = np.sqrt(surf)
nlat, nlon = surf.shape
#lon = mesh['glamt'].values
#lat = mesh['gphit'].values
#lon = np.squeeze(lon)
#lat = np.squeeze(lat)

tmask = xr.open_dataset('eof_mask.nc')
tmask = tmask['mask'].values
ilat, ilon = np.nonzero(tmask == 1)

data = xr.open_mfdataset('/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/yearly/*OOPE*nc')
data = data.isel(w=[14, 45, 80])
print(data)
dens = data['OOPE'].to_masked_array() # time, y, x, c, w
ilat, ilon, ic, iw = np.nonzero(dens[0].mask == False)
print(dens[:, ilat, ilon, ic, iw].shape)
dens[:, ilat, ilon, ic, iw] = sig.detrend(dens[:, ilat, ilon, ic, iw].T).T
time = data['time']
dens = dens - np.mean(dens, axis=0, keepdims=True)

dens = np.transpose(dens, (3, 4, 0, 1, 2))   # com, w, time, y, x
print(dens.shape)
ncom, nbins, ntime, nlat, nlon = dens.shape

neofs = 3

eofmap = np.zeros((ncom, nbins, neofs, nlat, nlon))
eofts = np.zeros((ncom, nbins, neofs, ntime))
eofvar = np.zeros((ncom, neofs, nbins))

for c in range(ncom):
    for s in range(nbins):

        print(dens.shape)
        print(ilat.shape)
        temp = dens[c, s, :, ilat, ilon].T
        print(temp.shape)

        print(weights[ilat, ilon].shape)
        solver = Eof(temp, weights=weights[ilat, ilon])
        # EOFs are multiplied by the square - root of
        # their eigenvalues ( units are in Pa )
        #maps = solver.eofs(eofscaling=2, neofs=neofs)
        maps = solver.eofsAsCovariance(neofs=neofs)

        eofmap[c, s, :, ilat, ilon] = maps.T
        eofts[c, s, :, :] = solver.pcs(pcscaling=1, npcs=neofs).T
        eofvar[c, s, :] = solver.varianceFraction(neigs=neofs)

eofmap = eofmap * tmask[np.newaxis, np.newaxis, np.newaxis, :, :]
eofmap = np.ma.masked_where(eofmap == 0, eofmap)

dirout = 'data/'
dirout = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

dsout = xr.Dataset({'eofmap':(['com', 'bins', 'eof', 'y', 'x'], eofmap), 
                    'eofpc':(['com', 'bins', 'eof', 'time'], eofts),
                    'eofvar':(['com', 'bins', 'eof'], eofvar)})
dsout['time'] = (['time'], time)
#dsout['sizes'] = (['sizes'], bins)
dsout.attrs['creation_date'] = str(date.today())
dsout.attrs['script'] = os.path.realpath(__file__)
dsout.attrs['description'] = 'EOF for densities by class'
dsout.to_netcdf('%s/annual_eof_oni_annual_density.nc' %(dirout))
