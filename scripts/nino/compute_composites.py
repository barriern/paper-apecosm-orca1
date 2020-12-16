import xarray as xr
import numpy as np
import apecosm.misc as misc
import matplotlib.pyplot as plt
from extract_nino import read_index
import os.path
from datetime import datetime

dirin = '/home1/datawork/nbarrier/apecosm/apecosm_orca1/diags/data'

def get_compo_year():

    daten, nino = read_index()
    nino = np.ma.masked_where(np.isnan(nino), nino)

    thres = 1

    inino = np.nonzero(nino > 1.75)[0]
    inina = np.nonzero(nino < -1)[0]
    ineutral = np.nonzero(np.abs(nino) < 0.5)[0]

    ynino = daten // 100
    mnino = (daten - 100 * ynino).astype(np.int)

    yfinal = []

    for y in np.unique(ynino):

        test = (ynino == y)
        test = test & (mnino >= 10) & (mnino <= 12)
        iok = np.nonzero(test)[0]
        temp = np.mean(nino[iok])
        if(temp > 1):
            yfinal.append(y)
    yfinal.remove(1986)

    return yfinal

if __name__ == '__main__':

    yfinal = get_compo_year()

    daten, nino = read_index()
    ynino = daten // 100
    mnino = (daten - 100 * ynino).astype(np.int)
    year = ynino
    month = mnino

    output = []

    for y in yfinal:
        istart = np.nonzero((year == y) & (month == 1))[0][0]
        iend = np.nonzero((year == y + 3) & (month == 12))[0][0]
        iend += 1
        output.append(nino[istart:iend])

    output = np.array(output)

    print(output.shape)

    plt.figure()
    month = np.arange(1, output.shape[1] + 1)
    print(month.shape)
    for p in range(output.shape[0]):
        plt.plot(month, output[p], label=yfinal[p])
        plt.fill_betweenx([-3, 3], [10], [12], color='lightgray')

    plt.xlim(month.min(), month.max())
    plt.ylim(-3, 3)
    plt.grid(True)

    plt.legend(ncol=3)

    plt.savefig('nino_index_composites.png', bbox_inches='tight')
