import copy
import Ngl
import Nio
import xarray as xr
import numpy as np

mesh = xr.open_dataset('../data/mesh_mask_eORCA1_v2.2.nc', engine='pynio')
mesh = mesh.isel(t=0, z=0)
lon = mesh['glamt'].values
lat = mesh['gphit'].values
tmask = mesh['tmask'].values

hadley = xr.open_dataset('data//covariance_yearly_enso_surface_sst.nc', engine='pynio')
lonhad = hadley['longitude'].values
lathad = hadley['latitude'].values
hadley = hadley['covariance'].to_masked_array()
hadley = np.ma.masked_where(hadley == 0, hadley)

model = xr.open_dataset('data//covariance_yearly_enso_surface_thetao.nc', engine='pynio')
model = model['covariance'].to_masked_array()
model = np.ma.masked_where(model == 0, model)

print(model.shape, lon.shape, lat.shape, hadley.shape)

resngl = Ngl.Resources()
resngl.wkColorMap = "BlueDarkRed18"
wks = Ngl.open_wks("pdf", "covariance_maps_hadley_model" , resngl)

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

res.lbTitleOffsetF = 0.1

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

res.pmLabelBarDisplayMode = "Never"

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
res.tmXBLabelsOn = True
res.tmYLLabelsOn = True
res.tmYTLabelsOn = False

val = 0.012
res.tmXTLabelFontHeightF = val
res.tmXBLabelFontHeightF = val
res.tmYLLabelFontHeightF = val
res.tmYRLabelFontHeightF = val

res.tiMainFontHeightF = 0.015

plot = []

cpt = 0 

res.lbTitleString = "Temperature (C)" # cbar title string
#res.lbTitleString = "" # cbar title string
#res.lbTitleString = "Density (J/m2)" # cbar title string

#res.tiMainString = '%s, size = %d cm' %(comm[c], size[w])
cmax = 1.6
step = 0.1
res.cnLevels = np.arange(-cmax, cmax + step, step)

res2 = copy.deepcopy(res)
res2.cnFillOn                    = False
res2.cnLinesOn                   = True
res2.cnLineLabelsOn              = False
res2.cnFillMode = "AreaFill"

res.sfXArray        = lonhad
res.sfYArray        = lathad
res.tiMainString = 'Hadley SST'
temp = Ngl.contour_map(wks, hadley, res)
plot.append(temp)


res.tiMainString = 'Model SST'
res.sfXArray        = lon
res.sfYArray        = lat

temp = Ngl.contour_map(wks, model, res)
temp2 = Ngl.contour(wks, model, res2)
#Ngl.overlay(temp, temp2)
plot.append(temp)

cpt += 1

panelres = Ngl.Resources()
panelres.nglPanelLabelBar                 = True     # Turn on panel labelbar
#panelres.nglPanelLabelBarLabelFontHeightF = 0.015    # Labelbar font height
#panelres.nglPanelLabelBarHeightF          = 0.1750   # Height of labelbar
panelres.nglPanelLabelBarWidthF           = 0.700    # Width of labelbar

#panelres.lbLabelFont                      = "helvetica-bold" # Labelbar font
#panelres.nglPanelTop                      = 0.935
panelres.nglPanelFigureStrings            = ["a)", "b)"]
panelres.nglPanelFigureStringsJust        = "BottomLeft"
panelres.nglPanelFigureStringsFontHeightF  = 0.02
panelres.nglPanelBottom = 0.05
#panelres.nglPanelXWhiteSpacePercent = 1.5
#panelres.nglPanelYWhiteSpacePercent = 3

txres               = Ngl.Resources()
txres.txFontHeightF = 0.017
Ngl.text_ndc(wks, "Temperature anomalies (C)", 0.52, 0.03, txres)


Ngl.panel(wks,plot, [2, 1], panelres)  

#Ngl.draw_colormap(wks)

Ngl.end()
