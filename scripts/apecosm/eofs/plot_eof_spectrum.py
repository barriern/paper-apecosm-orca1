import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import scipy.signal as sig
from scipy import interpolate
import envtoolkit.spectral as spec
from scipy.stats import lognorm, beta, norm
from scipy.optimize import curve_fit

''' Log normal PDF '''
def lnorm(x, mu, sigma, scale):
    return scale / (x * sigma * np.sqrt(2*np.pi)) * np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))

''' Normal PDF ''' 
def norm(x, mu, sigma, scale):
    return scale / (sigma * np.sqrt(2*np.pi)) * np.exp(-0.5 * ((x - mu) / (sigma))**2)

const = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values * 100 # length = cm

data = xr.open_dataset("data/eof_full_density_20.nc", decode_times=False)
pc = data['eofpc'].values  # bins, eof, time
pc = pc[:, 0, :]  # extraction of first EOF for first 
pc = pc[:, :]
t = np.arange(data.dims['time'])

sizes = np.array(['%.1ecm' %l for l in length])

nbandw = 11
deltat = 1.

# init. of the output spectra (raw and smoothed)
peri = []
smoothed = []

# loop over all the size classes.
for p in range(pc.shape[0]):

    pctemp = pc[p]
    #pctemp -= pctemp.mean()
    #pctemp /= pctemp.std()
    [spectrum2, freq2, error2] = spec.multitaper(pctemp, nbandw=nbandw, deltat=deltat)
    
    # removing 0 frequency
    spectrum2 = spectrum2[1:]
    freq2 = freq2[1:]
    
    # saving of normalized frequency
    y = spectrum2 * freq2
    peri.append(y)

    # fitting with a scaled normalized lnorm function 
    opt, cov = curve_fit(lnorm, freq2, y)

    # computes and save the fitting curve
    ytemp = lnorm(freq2, opt[0], opt[1], opt[2])
    smoothed.append(ytemp)

# conversion of the outputs from list to array
smoothed = np.array(smoothed)
peri = np.array(peri)
x = np.log(freq2)

# selection of some years to draw
x0 = np.array([0.5, 1, 5, 20, 40]) # years
x0 = np.arange(1, 10)
print(x0)
x0 *= 12   # months
xticks = np.log(1/x0)

plt.figure()

isize = [0, 20, 30, 40, 50, 60, 70, 80]
isize = [1, 7]

ax = plt.gca()
plt.title('Spectrum')

l = plt.plot(x, (peri).T[:, isize])
colors = [lt.get_color() for lt in l]

if True:
    for p in range(len(isize)):
        lbis = plt.plot(x, smoothed[isize[p]])

plt.legend(sizes[isize], fontsize=7)

plt.xlabel('Period (years)')
plt.ylabel('PSD * f')
ax.set_xticks(xticks)
ax.set_xticklabels(x0 / 12, rotation=45, ha='right')
ax.set_xlim(x.min(), x.max())
#ax.set_xlim(ax.get_xlim()[::-1])
#ax.set_xscale('log')
#ax.set_yscale('log')
plt.savefig('spectrum.pdf', bbox_inches='tight')

imax = np.argmax(smoothed, axis=1)
freqmax = freq2[imax]  
tmax = 1 / freqmax
tmax = tmax / 12.
tmax = freqmax
tmax = np.log(freqmax)

x0 = np.array([1, 2, 3, 4, 5, 7, 9])
x0 *= 12   # months
xticks = np.log(1/x0)

plt.figure()
ax = plt.subplot(1, 2, 1)
plt.plot(tmax, length, marker='.')
plt.ylim(length.min(), length.max())
ax.set_yscale('log')
plt.ylim(length.min(), 100)
plt.ylabel('Length (cm)')
plt.xlabel('Period (years)')
ax.set_xticks(xticks)
ax.set_xticklabels(x0 / 12, rotation=45, ha='right')
ax.set_xlim(-5, -2)
ax.grid(True)

step = 0.025
levels = np.arange(0, 0.3 + step, step)
ax2 = plt.subplot(1, 2,  2, sharey=ax, sharex=ax)
cs = plt.contourf(np.log(freq2), length, peri, cmap=plt.cm.jet, levels=levels, extend='both')
plt.yscale('log')
plt.ylim(length.min(), 100)
plt.colorbar(cs)
ax2.set_xticks(xticks)
ax2.set_xticklabels(x0 / 12, rotation=45, ha='right')
plt.xlabel('Period (years)')
ax2.grid(True)
plt.savefig('pcolor_power.pdf', bbox_inches='tight')

