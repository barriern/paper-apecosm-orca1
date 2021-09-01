# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['text.usetex'] = False

varname = 'PHY2'
# -

# ## Loading the file

data = xr.open_dataset('data/pacific_profile_anoms_%s.nc' %varname).drop('month')
data

data = data.rename({'__xarray_dataarray_variable__': varname})
data

ntime = data.dims['time_counter']
ntime

anom = data[varname]
anom

# ## Plotting the anomalies

figname = 'pacific_profile_anom_%s.pdf' %(varname)
figname

cmax = abs(anom).quantile(99.5 / 100).values
cmax = float(cmax)
cmax

with PdfPages(figname) as pdf:
    for t in range(ntime):
        #print(t, '/', ntime - 1)
        fig = plt.figure()
        toplot = anom.isel(time_counter=t)
        ax = plt.axes()
        cs = toplot.plot()
        cs.set_clim(-cmax, cmax)
        ax.set_ylim(ax.get_ylim()[::-1])
        pdf.savefig(bbox_inches='tight')
        plt.close(fig)


