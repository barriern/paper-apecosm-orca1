# +
import numpy as np
import xarray as xr
import sys
sys.path.append('../nino')
from extract_nino import read_index
import scipy.signal as sig
import os.path
import matplotlib.pyplot as plt
import scipy.stats as stats

latmax = 20
lonmin = 150
lonmax = -120
eof = 0
ddof = 12
# -

const =  xr.open_dataset('../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values

dnino, nino = read_index('../data/index/oni.data')
inino = np.nonzero((dnino >= 195801) & (dnino <= 201812))
nino = nino[inino]
dnino = dnino[inino]
nino = (nino - nino.mean()) / nino.std()

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/processed_pacific'
dirin = 'data/'

filename = '%s/full_eof_pacific_OOPE_latmax_%d_lonmin_%d_lonmax_%d.nc' %(dirin, latmax, lonmin, lonmax)
filename

data = xr.open_dataset(filename)
data

pc = data['eofpcs'].values  # bins, eof, time
pc = pc[:, :, eof]   # bins, time
# normalisation to insure that autocorr at lag 0 is 1.0
pc = (pc - np.mean(pc, axis=-1, keepdims=True)) / np.std(pc, axis=-1, keepdims=True)

nbins = pc.shape[0]

tpi = xr.open_dataset('../data/filt_tpi.nc')
tpi = tpi['tpi_filt'].to_masked_array()
iok = np.nonzero(tpi.mask == False)
tpi = tpi[iok]
tpi = (tpi - tpi.mean()) / tpi.std()
ntime2 = tpi.shape[0]
ntime = dnino.shape[0]

# Computation of the lag for the ONI time-series (complete) and the filtered TPI time-series (incomplete)

lags = sig.correlation_lags(ntime, ntime)
ilags = np.nonzero((lags >= 0) & (lags <= 5 * 12))[0]
ilags = np.nonzero(np.abs(lags) <= 5 * 12)[0]
lags = lags[ilags]
nlags = len(lags)

lags2 = sig.correlation_lags(ntime2, ntime2)
ilags2 = np.nonzero(np.abs(lags2) <= 5 * 12)[0]
lags2 = lags2[ilags2]


# Computation of the autocorrelation at lag 1 for TPI and ONI indexes.

def autocorr(pc):
    output = np.corrcoef(pc[1:], pc[:-1])[0, 1]
    return output


autocorr_oni = autocorr(nino)
autocorr_oni

autocorr_tpi = autocorr(tpi)
autocorr_tpi

corrtpi = np.zeros((nbins, nlags))
corroni = np.zeros((nbins, nlags))

sigtpi = np.zeros((nbins, nlags))
sigoni = np.zeros((nbins, nlags))

nptcom = ntime - np.abs(lags)
nptcom2 = ntime2 - np.abs(lags2)


def bretherton(a, b):
    res = (1 - a * b) / (1 + a*b)
    return res


tdis = stats.t
alpha = 0.95

for l in range(nbins):

    pctemp = pc[l]
    auto_pc = autocorr(pctemp)
    
    breth = bretherton(auto_pc, autocorr_oni)
    breth2 = bretherton(auto_pc, autocorr_tpi)
    
    dof = (nptcom - ddof) * breth
    dof2 = (nptcom2 - ddof) * breth2
    
    tlim = tdis.interval(alpha=alpha, df = dof - 2)[0]
    tlim2 = tdis.interval(alpha=alpha, df = dof2 - 2)[0]
    
    rlim = tlim / np.sqrt(dof - 2 + tlim**2)
    rlim2 = tlim2 / np.sqrt(dof2 - 2 + tlim2**2)
    
    sigoni[l] = rlim
    sigtpi[l] = rlim2
    
    corrtemp = sig.correlate(pctemp, nino)[ilags] / nptcom
    corroni[l] = corrtemp
    
    corrtemp = sig.correlate(pctemp[iok], tpi)[ilags2] / nptcom2
    corrtpi[l] = corrtemp

dsout = xr.Dataset()
dsout['corrtpi'] = (['l', 'lags'], corrtpi)
dsout['corroni'] = (['l', 'lags'], corroni)
dsout['sigtpi'] = (['l', 'lags'], sigtpi)
dsout['sigoni'] = (['l', 'lags'], sigoni)
dsout['lags'] = (['lags'], lags)
dsout['l'] = (['l'], length * 100)
dsout.to_netcdf('full_eof_tpi_oni_corr_eof_%d.nc' %(eof + 1))

