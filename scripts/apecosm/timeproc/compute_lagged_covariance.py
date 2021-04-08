from glob import glob
import xarray as xr
import numpy as np
import scipy.signal as sig
import sys
import os
sys.path.append('/home1/datahome/nbarrier/apecosm/configs/apecosm_orca1/diags/paper-apecosm-orca1/scripts/nino')
from extract_nino import read_index
import apecosm.ts as ts

dmin, dmax = 195801, 201812
daten, nino = read_index(filename='/home1/datahome/nbarrier/apecosm/configs/apecosm_orca1/diags/nino/index/oni.data')
iok = np.nonzero((daten >= 195801) & (daten <= 201812))
daten, nino = daten[iok], nino[iok]
ntime = len(nino)

lags = sig.correlation_lags(ntime, ntime)
ilags = np.nonzero(np.abs(lags) <= 3*12)[0]
lags = lags[ilags]
nlags = len(lags)

dirout = '/home/datawork-marbec-pmod/outputs/APECOSM/ORCA1/final-runs/output/tsanalysis'
dirin = dirout

def compute_lagged_covariance(varname):
    
    data = xr.open_dataset('%s/flipped_%s.nc' %(dirin, varname))
    nw = data.dims['w']
    nx = data.dims['x']
    ny = data.dims['y']
    data = data[varname]
    
    cov = np.zeros((nw, nx, ny, nlags))
    for w in range(nw):
        for i in range(nx):
            for j in range(ny):
                temp = data.isel(w=w, x=i, y=j).values
                if(np.isnan(temp[0])):
                    continue
                clim, temp = ts.get_monthly_clim(temp)
                temp = sig.detrend(temp)
                cov[w, i, j, :] = sig.correlate(temp, nino)[ilags] / ntime
                
    dsout = xr.Dataset()
    dsout['cov'] = (['w', 'x', 'y', 'lags'], cov)
    dsout['lags'] = (['lags'], lags)
    dsout.attrs['file'] = os.path.realpath(__file__)
    dsout.to_netcdf('%s/covariance_%s.nc' %(dirin, varname))



if __name__ == '__main__':

    compute_lagged_covariance('OOPE')
    #compute_lagged_covariance('repfonct_day')
