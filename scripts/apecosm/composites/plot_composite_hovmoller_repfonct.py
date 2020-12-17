import xarray as xr
import numpy as np
import apecosm.misc as misc
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 13

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'
dirin = 'data/'

constant = xr.open_dataset('/home/barrier/Work/apecosm/ORCA1/figures_orca1/plot_hov/ORCA1_JRA_CO2_CYC4_ConstantFields.nc')
constant = constant.isel(w=[14, 45, 80])
wstep = constant['weight_step'].values

pref = 'corr_mask'
pref = 'debugged_corr_mask'
#pref = 'climTemp'
#pref = 'climPlk'

for varname in ['repfonct_day']:

    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@ ", varname)

    data = xr.open_dataset('%s/%s_hovmoller_composite_%s.nc' %(dirin, pref, varname))
    lon = data['lon'].values
    lon[lon < 0] += 360
    lonmin = 130
    lonmax = 300
    iok = np.nonzero((lon <= lonmax) & (lon >= lonmin))[0]

    time = data['time'].values + 1
    lon[lon < 0] += 360

    xticks1 = np.arange(180 + 30, lon.max() + 30, 30)
    labels1 = ['%dE' %(x - 180) for x in xticks1]

    xticks2 = np.array([180 - 2*30, 180 - 30])
    labels2 = ['%dW' %(-x + 180) for x in xticks2]

    xticks3 = np.array([180])
    labels3 = [0]

    xticks = np.concatenate((xticks2, xticks3, xticks1))
    labels = labels2 + labels3 + labels1

    sizes = np.array([2, 20, 90])
    comm = ['Meso.', 'Mig.', 'Epi.'][::-1]

    oope = data['compo_mean'].values[:, iok, :, :] * 100
    lon = lon[iok]

    #xticks = time
    monthticks = np.array(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
    monthticks = np.ravel(np.tile(monthticks, (4, 1)))

    cpt = 0
       
    plt.figure(figsize=(12, 12))
    plt.subplots_adjust(left=0.04, right=0.95, bottom=0.05, top=0.98, wspace=0.2)

    cmaxlist = [0, 250, 30, 6, 80, 40, 40, 250, 50, 15]
    if(pref == 'debugged_corr_mask'): 
        cmaxlist = [0, 
                    80, 8, 3, 
                    40, 10, 20, 
                    30, 3, 4]

    for c in range(0, 3):

        for s in range(3):

            cpt += 1
            print('+++++++++++++++++++++++++++++++++++++++ c=%d, s=%d' %(c, s))

            ax = plt.subplot(3, 3, cpt)
            temp = oope[:, :, c, s] 
            temp = np.ma.masked_where(temp == 0, temp)
            #print(temp)
            #ax.set_facecolor('black')
     
            cmin, cmax = misc.find_percentile(np.abs(temp), 1)
            cmax = cmaxlist[cpt] 

            cmin = -cmax
            levels = np.linspace(cmin, cmax, 11)

            title = '%s, %s cm' %(comm[c], sizes[s])
            plt.title(title)
            cs = plt.pcolormesh(lon, time, temp, cmap='RdBu_r')
            #cl = plt.contour(lon, time, temp, colors='k', levels=levels, linewidths=1)
            #cs.set_clim(cmin, cmax)
            cb = plt.colorbar(cs)
            #cb.add_lines(cl)

            if cpt in [3, 6, 9]:
                cb.set_label('')

            if cpt in [1, 1 + 3, 1 + 6]:
                plt.ylabel('Time (month)')
            if cpt in [7, 8, 9]:
                plt.xlabel('Longitude')
            plt.grid(True)

            ax.set_yticks(np.arange(12, 48, 12))

            xticks = np.arange(150, 150 + 6 * 30, 30)

            xticks1 = [150]
            labels1 = ['%dE' %(x) for x in xticks1]

            xticks2 = [0]
            labels2 = ['%d$^\circ$' %(x) for x in xticks2]

            xticks3 = [210, 240, 270, 300]
            labels3 = ['%dW' %(np.abs(360-x)) for x in xticks3]

            labels = labels1 + labels2 + labels3

            plt.gca().set_xticks(xticks)
            plt.gca().set_xticklabels(labels)

            ax.set_xticks(xticks)
            ax.set_xticklabels(labels)
            ax.set_xlim(lonmin, lonmax)

    plt.savefig('%s_hovmoller_composites_%s.pdf' %(pref, varname), bbox_inches='tight')
