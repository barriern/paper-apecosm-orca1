import xarray as xr
import apecosm.ts as ts 
import scipy.signal as sig
import numpy as np

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

data = xr.open_dataset('%s/eof_oni_annual_density_20.nc' %dirin)
pc = data['eofpc'].values  #  com, bins, eof, time 
ncom, nbins, neof, ntime = pc.shape
print('pc.shape = ', pc.shape)

forage = xr.open_mfdataset("%s/meridional_-150_forage_year_*.nc" %dirin)
lat = forage['y'].values
forage = forage['FORAGE'].to_masked_array()   # time, dn, y, depth, comm, size
print('Forage 1 ', forage.shape)
clim, forage = ts.get_monthly_clim(forage)   # get monthly anomalies
del(clim)
print('Forage 2 ', forage.shape)
ntime, ndn, ny, ndepth, ncom, ns2 = forage.shape
forage = np.transpose(forage, (4, 5, 1, 2, 3, 0))   # comm, size, dn, y, depth, time
forage = sig.detrend(forage)
print('Forage 3 ', forage.shape)

covout = np.zeros((ncom, nbins, neof, ndn, ny, ndepth))
print('Covout.shape ', covout.shape)

for c in range(ncom): 
    for s in range(nbins):
        for e in range(neof):
            pctemp = pc[c, s, e]   # time
            for d in range(ndn):
                for y in range(ny):
                    for z in range(ndepth):
                        ftemp = forage[c, s, d, y, z]
                        covout[c, s, e, d, y, z] = np.cov(pctemp, ftemp)[0, 1]
#                        break
#                    break
#                break
#            break
#        break
#    break

dsout = xr.Dataset()
dsout['cov'] = (['c', 's', 'e', 'd', 'y', 'z'], covout)
dsout['c'] = (['y'], lat)
dsout.to_netcdf('%s/cov_forage_pc.nc' %dirin)

