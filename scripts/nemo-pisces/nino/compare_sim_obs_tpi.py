import xarray as xr
import matplotlib.pyplot as plt
from envtoolkit.ts import Lanczos 
import numpy as np

nWgt      = 157                          
fca       = 1./156                      
fca = 1 / fca

lanc = Lanczos('lp', pca=fca, nwts=nWgt)

datamod = xr.open_dataset('model_tpi_index.nc')
mod = datamod['tpi'].values
dmod = datamod['time'].values

dataobs = xr.open_dataset('../../data/filt_tpi.nc')
obs = dataobs['tpi_raw'].values
fobs = dataobs['tpi_filt'].values
dobs = dataobs['time'].values

test_filt = lanc.wgt_runave(obs)

ntime = len(dobs)

print(np.corrcoef(obs, mod)[0, 1])

plt.figure()
#plt.plot(mod, label='raw mod')
#plt.plot(obs, label='raw obs')
plt.plot(test_filt, label='my filt')
plt.plot(fobs, label='their filt')
plt.legend()
plt.savefig('compare_tpi.png')
