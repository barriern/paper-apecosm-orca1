import Ngl
import Nio
import xarray as xr
import numpy as np
import copy
from misc import find_percentile

prefix = 'debugged_corr_mask'

constant = xr.open_dataset('/home/barrier/Work/apecosm/ORCA1/figures_orca1/plot_hov/ORCA1_JRA_CO2_CYC4_ConstantFields.nc')
constant = constant.isel(w=[14, 45, 80])
wstep = constant['weight_step'].values
print(wstep)

datamean = xr.open_dataset('data/%s_yearly_mean_OOPE.nc' %prefix, engine='pynio')
datamean = datamean.isel(time=0)
datamean = datamean['OOPE'].to_masked_array()
for c in range(3):
    datamean[:, :, :, c] *= wstep[c]
datamean = np.log10(datamean, where=((datamean.mask == False) & (datamean > 0)))


mesh = xr.open_dataset('../data/mesh_mask_eORCA1_v2.2.nc', engine='pynio')
mesh = mesh.isel(t=0, z=0)
lon = mesh['glamt'].values
lat = mesh['gphit'].values
tmask = mesh['tmask'].values

lon2 = lon.copy()
lon2[lon2<0] += 360
tttt = (lon2 >= 130) & (lon2 <= 300)
tttt = tttt & (np.abs(lat) <= 40)


data = xr.open_dataset('data/%s_covariance_yearly_enso_OOPE.nc' %prefix, engine='pynio')
cov = data['covariance'].to_masked_array()

resngl = Ngl.Resources()
resngl.wkColorMap = "BlueDarkRed18"
wks = Ngl.open_wks("pdf", "%s_covariance_maps_OOPE" %prefix , resngl)

res = Ngl.Resources()

res.nglDraw = False
res.nglFrame = False

# Set map resources.
res.mpLimitMode = "LatLon"     # limit map via lat/lon
res.mpMinLatF = -40         # map area
res.mpMaxLatF = 40    # latitudes
res.mpMinLonF = 130     # and
res.mpMaxLonF = 300     # longitudes
res.mpCenterLonF = 180

print(res.mpMinLonF, res.mpMaxLonF)

res.mpFillOn = True
res.mpLandFillColor = "LightGray"
res.mpOceanFillColor = -1
res.mpInlandWaterFillColor = "LightBlue"
res.mpGeophysicalLineThicknessF = 1.5
res.mpOutlineOn   = True
res.mpGridAndLimbOn = False
res.mpOutlineBoundarySets = "National"
res.mpFillOn                    = True

res.lbTitleOffsetF = 0.5

res.lbLabelStride = 4
res.lbOrientation = "Vertical" # vertical colorbar 
res.lbTitlePosition = "Right" # cbar title position
res.lbTitleAngleF = 90
res.lbTitleDirection = "Across"
res.lbTitleFontHeightF = 0.02
res.lbLabelFontHeightF = 0.02

# Scalar field resources
res.sfXArray        = lon
res.sfYArray        = lat

res.cnFillOn                    = True
res.cnLinesOn                   = False
res.cnLineLabelsOn              = False
res.cnFillMode = "CellFill"
res.nglSpreadColorEnd = 2 # index of first color for contourf
res.nglSpreadColorStart = 27 # index of last color for contourf

res.nglSpreadColorEnd, res.nglSpreadColorStart = res.nglSpreadColorStart, res.nglSpreadColorEnd

res.cnLevelSelectionMode="ExplicitLevels"
res.pmLabelBarWidthF = 0.06 # cbar width

# map tick managements
res.tmYROn       = True
res.tmYRLabelsOn = False
res.tmYLOn       = True
res.tmYLLabelsOn = True

res.tmXTOn       = True
res.tmXBOn       = True
res.tmXTLabelsOn = False
res.tmXBLabelsOn = False
res.tmYLLabelsOn = False
res.tmYTLabelsOn = False

res.tmXTLabelFontHeightF = 0.02
res.tmXBLabelFontHeightF = 0.02
res.tmYLLabelFontHeightF = 0.02
res.tmYRLabelFontHeightF = 0.02

plot = []

comm = ['Epi.', 'Mig.', 'Meso.']
size = [3, 20, 90]

cmaxlist = [260, 20, 3, 30, 10, 10, 150, 10, 5]
cmaxlist = [60, 5, 1, 15, 4, 8, 20, 3, 2]

res2 = copy.deepcopy(res)
res2.cnFillOn                    = False
res2.cnLinesOn                   = True
res2.cnLineLabelsOn              = True
res2.cnInfoLabelOn  = False
res2.cnLineLabelDensityF  = 2.0
res2.cnLineLabelPlacementMode = "Constant"
res2.cnFillMode = "AreaFill"
res2.cnLevelSelectionMode = 'AutomaticLevels'
res2.cnMaxLevelCount = 11

cpt = 0 

for c in range(3):
    for w in range(3):
    
        if  cpt < 6:
            res.tmXBLabelsOn       = True
        else:
            res.tmXBLabelsOn = True
        
        if  cpt < 3:
            res.tmXTLabelsOn       = False
        else:
            res.tmXTLabelsOn = False
        
        if  cpt in [0, 3, 6]:
            res.tmYLLabelsOn       = True
        else:
            res.tmYLLabelsOn = True

        if cpt in [2, 2 + 3, 2 + 6]:
            res.lbTitleString = "Density (J/m2)" # cbar title string
        else:
            res.lbTitleString = "" # cbar title string
            #res.lbTitleString = "Density (J/m2)" # cbar title string

        res.tiMainString = '%s, size = %d cm' %(comm[c], size[w])

        temp =  cov[:,:, c, w] * wstep[w]
        temp = np.ma.masked_where(tmask==0, temp)
        temp = np.ma.masked_where(tttt == False, temp)

        temp2 = datamean[:, :, c, w]
        temp2 = np.ma.masked_where(tttt == False, temp2)

        cmax = np.max(np.abs(temp))
        #cmin, cmax = find_percentile(np.abs(temp))
        cmax = cmaxlist[cpt]
        res.cnLevels = np.linspace(-cmax, cmax, 21)

        temp = Ngl.contour_map(wks, temp, res)
        temp2 = Ngl.contour(wks, temp2, res2)
        #Ngl.overlay(temp, temp2)

        plot.append(temp)

        cpt += 1

panelres = Ngl.Resources()
#panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
#panelres.nglPanelLabelBarLabelFontHeightF = 0.015    # Labelbar font height
#panelres.nglPanelLabelBarHeightF          = 0.1750   # Height of labelbar
#panelres.nglPanelLabelBarWidthF           = 0.700    # Width of labelbar
#panelres.lbLabelFont                      = "helvetica-bold" # Labelbar font
#panelres.nglPanelTop                      = 0.935
panelres.nglPanelFigureStrings            = ["a)", "b)", "c)", "d)", "e)", "f)", 'g)', 'h)', 'i)']
panelres.nglPanelFigureStringsJust        = "BottomLeft"
panelres.nglPanelFigureStringsFontHeightF  = 0.01
panelres.nglPanelRight = 0.95
panelres.nglPanelXWhiteSpacePercent = 1.5
panelres.nglPanelYWhiteSpacePercent = 3

Ngl.panel(wks,plot, [3, 3], panelres)  

#Ngl.draw_colormap(wks)

Ngl.end()
