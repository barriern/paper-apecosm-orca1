import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import apecosm.misc as misc

plt.rcParams['text.usetex'] = False
zmax = 500
dirin = 'data/'

convert = {'PHY':'Nano. phy.', 'PHY2': 'Diat.', 'ZOO': 'Micro. zoo.', 'ZOO2': 'Meso. zoo.', 'thetao': 'Temp.', 'O2': 'Oxy.'}
convert['PAR'] = 'PAR'
convert['GOC'] = 'Big part. matter'
units = {}
units['O2'] = 'mmol/m3'
units['PHY'] = 'mmol/m3'
units['PHY2'] = 'mmol/m3'
units['ZOO'] = 'mmol/m3'
units['ZOO2'] = 'mmol/m3'
units['thetao'] = 'C'
units['PAR'] = 'W/m2'
units['GOC'] = 'mmol/m3'

props = dict(boxstyle='round', facecolor='lightgray', alpha=1.)

def plot_var(varname, letter, cmax=None):

    print('Plotting ', varname)

    data2 = xr.open_dataset('data/mean_profile_yearly_%s.nc' %(varname))
    datamean = data2[varname].values

    data = xr.open_dataset('%s/covariance_yearly_enso_profile_%s.nc' %(dirin, varname))
    cov = data['covariance'].values
    lon = data['lon'].values
    depth = -data['depth'].values
    ndepth, nlon = cov.shape
    lon[lon < 0] += 360

    lonmin = 130
    lonmax = 300
    iok = np.nonzero((lon <= lonmax) & (lon >= lonmin))[0]
    lon = lon[iok]
    cov = cov[:, iok]
    print(cov.shape)

    datamean = datamean[:, iok]

    iz = np.nonzero(np.abs(depth) <= zmax)[0]
    depth = depth[iz]
    cov = cov[iz, :]
    datamean = datamean[iz, :]
    
    depth /= 1000

    cov = np.ma.masked_where(cov == 0, cov)
    if cmax is None:
        cmin, cmax = misc.find_percentile(np.abs(cov))
    cmin = -cmax

    cs = plt.pcolormesh(lon, depth, cov, cmap='RdBu_r')
    cs.set_clim(cmin, cmax)
    levels = np.linspace(cmin, cmax, 11)
    #cl = plt.contour(lon, depth, cov, colors='k', levels=levels, linewidths=0.5)
    #cl2 = plt.contour(lon, depth, cov, colors='k', levels=0, linewidths=1)
    cl = plt.contour(lon, depth, datamean, colors='k', linewidths=1)
    plt.clabel(cl, fmt='%.1f', fontsize=6)
    cb = plt.colorbar(cs)
    #cb.add_lines(cl)
    cb.set_label(units[varname])
    #cb.set_ticks(levels)
    plt.ylabel('Depth (km)')
    #plt.xlabel('Longitude')
    plt.title(convert[varname])
    plt.grid(True)

    xticks = np.arange(150, 150 + 6 * 30, 30)
    print(xticks)

    xticks1 = [150]
    labels1 = ['%dE' %(x) for x in xticks1]
    
    xticks2 = [0]
    labels2 = ['%d$^\circ$' %(x) for x in xticks2]
    
    xticks3 = [210, 240, 270, 300]
    labels3 = ['%dW' %(np.abs(360-x)) for x in xticks3]

    print(xticks)

    yticks = np.arange(-400, 0, 100)/1000
    plt.gca().set_yticks(yticks)

    labels = labels1 + labels2 + labels3

    plt.gca().set_xticks(xticks)
    plt.gca().set_xticklabels(labels)

    plt.xlim(lonmin, lonmax)

    plt.text(lon.max()-10, -400/1000, letter, ha='right', va='center', bbox=props, fontsize=10)

def hidex(ax1):
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_xlabel('')
def hidey(ax2):
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.set_ylabel('')

plt.figure(figsize=(10, 8))
plt.subplots_adjust(hspace=0.2, wspace=0.1)
cpt = 1

let = ['', 'a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)']

ax1 = plt.subplot(4, 2, cpt)
plot_var('thetao', let[cpt])
hidex(ax1)
cpt += 1

ax1 = plt.subplot(4, 2, cpt)
plot_var('O2', let[cpt])
hidex(ax1)
hidey(ax1)
cpt += 1

ax1 = plt.subplot(4, 2, cpt)
plot_var('PHY', let[cpt])
hidex(ax1)
cpt += 1

ax2 = plt.subplot(4, 2, cpt)
plot_var('PHY2', let[cpt])
hidey(ax2)
hidex(ax2)
cpt += 1

ax1 = plt.subplot(4, 2, cpt)
plot_var('ZOO', let[cpt])
hidex(ax1)
cpt += 1

ax2 = plt.subplot(4, 2, cpt)
plot_var('ZOO2', let[cpt])
hidey(ax2)
hidex(ax2)
cpt += 1

ax1 = plt.subplot(4, 2, cpt)
plot_var('GOC', let[cpt])
cpt += 1

ax2 = plt.subplot(4, 2, cpt)
plot_var('PAR', let[cpt])
hidey(ax2)
cpt += 1

plt.savefig('covariance_profiles.pdf', bbox_inches='tight')
