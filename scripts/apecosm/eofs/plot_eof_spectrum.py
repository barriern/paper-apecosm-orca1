import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import scipy.signal as sig
from scipy import interpolate
import envtoolkit.spectral as spec
from scipy.stats import lognorm, beta, norm
from scipy.optimize import curve_fit

def lnorm(x, mu, sigma, scale):

    return scale / (x * sigma * np.sqrt(2*np.pi)) * np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))

def norm(x, mu, sigma, scale):
    return scale / (sigma * np.sqrt(2*np.pi)) * np.exp(-0.5 * ((x - mu) / (sigma))**2)

const = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values * 100 # length = cm

data = xr.open_dataset("data/eof_full_density_20.nc", decode_times=False)
pc = data['eofpc'].values  # com, bins, eof, time
pc = pc[:, 0, :]  # extraction of first EOF for first 
isize = [0, 20, 30, 40, 50, 60, 70, 80]
isize = [0, 30, 70]
pc = pc[:, :]
t = np.arange(data.dims['time'])

sizes = np.array(['%.1ecm' %l for l in length])
print(sizes)

nbandw = 11
deltat = 1.

# init. of the output spectra (raw and smoothed)
peri = []
smoothed = []

# loop over all the EOFs.
for p in range(pc.shape[0]):

    pctemp = pc[p]
    pctemp -= pctemp.mean()
    pctemp /= pctemp.std()
    [spectrum2, freq2, error2] = spec.multitaper(pctemp, nbandw=nbandw, deltat=deltat)
    
    # removing 0 frequency
    spectrum2 = spectrum2[1:]
    freq2 = freq2[1:]
    
    # saving of normalized frequency
    y = spectrum2 * freq2
    peri.append(y)

    # fitting with a scaled normalized lnorm function 
    opt, cov = curve_fit(lnorm, freq2, y)

    ytemp = lnorm(freq2, opt[0], opt[1], opt[2])

    smoothed.append(ytemp)

smoothed = np.array(smoothed)
peri = np.array(peri)
x = np.log(freq2)

plt.figure()
x = freq2
x = np.log(freq2)
#x = freq2

'''
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.title('PCs')
plt.plot(t/12, pc.T)
plt.xlabel('Years')
plt.ylim(-3, 3)
plt.ylabel('PC')
plt.legend(sizes, fontsize=7, ncol=3)
'''

x0 = np.arange(1, 30 + 4, 4) # x tick labels
print(x0)
x0 = np.array([1, 5, 20, 40])
x0 *= 12
xticks = np.log(1/x0)

ax = plt.gca()
plt.title('Spectrum')

l = plt.plot(x, (peri).T[:, isize])
colors = [lt.get_color() for lt in l]
print(colors)

for p in range(len(isize)):
    print(isize[p])
    print(smoothed)
    lbis = plt.plot(x, smoothed[isize[p]])
    #lbis = plt.plot(xnew, smoothed[isize[p]], linestyle='--')

plt.legend(sizes[isize], fontsize=7)

plt.xlabel('Period (years)')
plt.ylabel('PSD * f')
ax.set_xticks(xticks)
ax.set_xticklabels(x0 / 12, rotation=45, ha='right')
ax.set_xlim(x.min(), x.max())
#ax.set_xlim(ax.get_xlim()[::-1])
plt.savefig('spectrum.pdf', bbox_inches='tight')
