# +
import numpy as np
import xarray as xr
import sys
sys.path.append('../nino')
from extract_nino import read_index
import scipy.signal as sig
import os.path

latmax = 20
lonmin = 150
lonmax = -120
eof = 0
# -

signs = -np.ones((100))
index = list(range(84, 93 + 1)) + [99]
index = np.array(index)
signs[index] *= -1

const =  xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values
print(length)

dnino, nino = read_index('../data/index/oni.data')
inino = np.nonzero((dnino >= 195801) & (dnino <= 201812))
nino = nino[inino]
dnino = dnino[inino]

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'
dirin = 'data/'

filename = '%s/full_eof_pacific_OOPE_latmax_%d_lonmin_%d_lonmax_%d.nc' %(dirin, latmax, lonmin, lonmax)
filename

data = xr.open_dataset(filename)
data

pc = data['eofpcs'].values  # bins, eof, time
pc = pc[:, :, eof]   # bins, time
pc.shape

#pc = (pc.T * signs).T
nbins = pc.shape[0]

tpi = xr.open_dataset('../data/filt_tpi.nc')
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

    corrtemp = sig.correlate(pctemp[iok], tpi)[ilags2] / ntime2
    corrtpi[l] = corrtemp

print(corrtpi[45])

dsout = xr.Dataset()
dsout['corrtpi'] = (['l', 'lags'], corrtpi)
dsout['corroni'] = (['l', 'lags'], corroni)
dsout['lags'] = (['lags'], lags)
dsout['l'] = (['l'], length * 100)
#dsout.attrs['file'] = os.path.realpath(__file__)
dsout.to_netcdf('full_eof_tpi_oni_corr_eof_%d.nc' %(eof + 1))

