import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import scipy.signal as sig
import envtoolkit.spectral as spec

const = xr.open_dataset('../../data/ORCA1_JRAC02_CORMSK_CYC3_FINAL_ConstantFields.nc')
length = const['length'].values * 100 # length = cm

data = xr.open_dataset("data/eof_full_density_20.nc", decode_times=False)
pc = data['eofpc'].values  # com, bins, eof, time
pc = pc[:, 0, :]  # extraction of first EOF for first 
size = [0, 10, 20, 30, 40, 50, 60, 70]
pc = pc[size, :]
t = np.arange(data.dims['time'])

length = length[size]

sizes = ['%.fcm' %l for l in length]
print(sizes)

nbandw = 5
deltat = 1.

# init. of the output spectra (raw and smoothed)
peri = []
smoothed = []

# loop over all the EOFs.
for p in range(pc.shape[0]):

    pctemp = pc[p]
    [spectrum2, freq2, error2] = spec.multitaper(pctemp, nbandw=nbandw, deltat=deltat)
    
    #freq2, spectrum2 = sig.welch(pctemp, nperseg=51)
    
    # removing 0 frequency
    spectrum2 = spectrum2[1:]
    freq2 = freq2[1:]


    peri.append(spectrum2 * freq2)
    x = np.log(freq2)
    
    y = spectrum2 * freq2
    pol = np.polyfit(x, y, 4)
    ytemp = np.polyval(pol, x)

    smoothed.append(ytemp)

smoothed = np.array(smoothed)
peri = np.array(peri)
x = np.log(freq2)

print(x.shape)
print(peri.shape)
print(smoothed.shape)

plt.figure()
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

#ax = plt.subplot(212)
ax = plt.gca()
plt.title('Spectrum')
l = plt.plot(x, (peri).T)
colors = [lt.get_color() for lt in l]
print(colors)
plt.legend(sizes, fontsize=7)
#for p in range(peri.shape[0]):
#    lbis = plt.plot(x, smoothed[p], linestyle='--', color=colors[p])
plt.xlabel('ln(f)')
plt.ylabel('PSD * f')
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xticks(xticks)
ax.set_xticklabels(x0 / 12, rotation=45, ha='right')
ax.set_xlim(ax.get_xlim()[::-1])
plt.savefig('spectrum.pdf', bbox_inches='tight')
