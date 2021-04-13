import numpy as np
import xarray as xr
import sys
sys.path.append('../../nino')
from extract_nino import read_index
import scipy.signal as sig
import os.path

const =  xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values
print(length)

dnino, nino = read_index('../../data/index/oni.data')
inino = np.nonzero((dnino >= 195801) & (dnino <= 201812))
nino = nino[inino]
dnino = dnino[inino]

data = xr.open_dataset('data/eof_full_density_20.nc')
pc = data['eofpc'].values  # bins, eof, time
pc = pc[:, 0, :]   # bins, time
nbins = pc.shape[0]

tpi = xr.open_dataset('../../nino/filt_tpi.nc')
tpi = tpi['tpi_filt'].to_masked_array()
iok = np.nonzero(tpi.mask == False)
tpi = tpi[iok]
tpi = (tpi - tpi.mean()) / tpi.std()
ntime2 = tpi.shape[0]


ntime = dnino.shape[0]

lags = sig.correlation_lags(ntime, ntime)
ilags = np.nonzero((lags >= 0) & (lags <= 5 * 12))[0]
ilags = np.nonzero(np.abs(lags) <= 5 * 12)[0]
lags = lags[ilags]
nlags = len(lags)

lags2 = sig.correlation_lags(ntime2, ntime2)
ilags2 = np.nonzero(np.abs(lags2) <= 5 * 12)[0]

corrtpi = np.zeros((nbins, nlags))
corroni = np.zeros((nbins, nlags))

for l in range(nbins):

    pctemp = pc[l]

    corrtemp = sig.correlate(pctemp, nino)[ilags] / ntime
    corroni[l] = corrtemp
    
    print(tpi.shape)
    print(pctemp[iok].shape)
    print(ntime2)
    corrtemp = sig.correlate(pctemp[iok], tpi)[ilags2] / ntime2
    corrtpi[l] = corrtemp

print(corrtpi[45])

dsout = xr.Dataset()
dsout['corrtpi'] = (['l', 'lags'], corrtpi)
dsout['corroni'] = (['l', 'lags'], corroni)
dsout['lags'] = (['lags'], lags)
dsout['l'] = (['l'], length * 100)
dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('full_eof_tpi_oni_corr.nc')

