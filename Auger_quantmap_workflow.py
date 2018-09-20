# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 17:10:19 2016

@author: tkc
"""
import pandas as pd
import numpy as np
import sys, glob
import scipy.stats
import matplotlib.pyplot as plt
import os
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules')
import Auger_smdifquant_functions as AESsmquant
import Auger_quantmap_functions as QM
from Auger_utility_functions import pickelemsGUI
import Auger_utility_functions as AESutils
from scipy.signal import medfilt
os.chdir('C:\\Temp\\AugerQM') 

#%% CREATE PHI FILES FOR AUTOTOOL, SPATIAL AREAS, MULTIPLEX CONDITIONS
# A few of these are also stored in Auger import main (to allow QM data combination prior to quant)
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\AESquantparams.csv', encoding='utf-8')
AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437')

# create pixarray file (correlating spe file w/ pixel position) and Autotool + spatial area files
QMpixarray=QM.QMarray_setup()
QMpixarray=QMarray_setup()
# Save and name pixarray
QMpixarray.to_csv('QMpixarray_50x50scr20m.csv', index=False)

# instead create a rectangular array (i.e. scan over FIB section area)
QMpixarray=QMrectarray_setup()
QMpixarray.to_csv('QMpixarray_rectangle.csv', index=False)

# Choose element for quant map
Elements=AESutils.pickelemsGUI(AESquantparams, Smdifpeakslog, Integquantlog)
Elements=['C','O','Fe','Mg','Si']

# Create custom quantmap multiplex file
QM.QMmultiplex_setup(Elements, AESquantparams)
multiplex_setup = QMmultiplex_setup(Elements, AESquantparams)

# Make annotation of image with overlaid mapped region (last arg is margin )
QM.showQMregion('1Sep17.247.jpg', AugerParamLog, 0.2)
# Interactive make annotation from pixarray (separate cropped jpg)
superimposearr(QMpixarray, allregs=True, crop=True)

# Reload pixarray definition file 
QMpixarray=QM.loadQMpixarray()
npfiles=glob.glob('*.npy')
basename='Acfer094map'

# Associate multiplex spe data files with correct pixels from quantmap (after data collection)
QMpixarray=QM.linkfilename(QMpixarray, 'GC1_20Oct17', startnum=115)

QMpixarray.to_csv('GC_rectangle_QMpixarray', index=False) # requires manual save
QMpixarray=QM.loadQMpixarray() # reload after pixel positions are linked with spe files

# Get spectral region details and multiplex energy values from quantmap multiplex file
spectralregs, energy = QM.get_spectral_regs(QMpixarray)

# Make 3D numpy array (specimage) with all spectral data; energy x values in energy list
specimage, energy =QM.makespecimage(QMpixarray, spectralregs)

specimage=QM.loadspecimage(os.getcwd()) # load existing spectral image 

np.save('GC3_specimage.npy', specimage) # save numpy array to file
specimage=np.load('Acfer094_Omap_area8_pixscrambled.npy')
specimage2=np.load('Acfer094map_area8.npy') # Reload array from disk

# Generate full element data associated with spectral images 
# (element list, peak/low/hiback index #s dand eV ranges)
Elemdata=QM.getelemdata(spectralregs, AESquantparams)

# Use for odd setups (like single continuous peak-background regions)
kwargs={'elems':['O']}
Elemdata=QM.getelemdata(spectralregs, AESquantparams, **kwargs)

# Find charge across mapped region using O map/scan
chargemap, peakamplmap=QM.findnegpeaks(specimage, Elemdata, 'O') # for wider scan w/ s7d7
chargemap, peakamplmap=findnegpeaks(specimage, Elemdata, 'O')

# new combined method w/ deriv and integcounts
amplmaps, shiftmaps, integmaps=QM.findallpeaks(specimage, Elemdata)

# Save list of maps in standard np stack 
QM.savemaps(amplmaps, shiftmaps, integmaps, 'test') # save as uniquestr +'_amplmaps.npy'
# Reload of already saved stacks of maps
amplmaps, shiftmaps, integmaps, elemmaps=QM.loadmaps(os.getcwd())

# Quick look at spatial maps of charging and underlying peak amplitudes (raw version) 
QM.plotcharging(shiftmaps[1], amplmaps[1])

# Compare deriv based shift with integ based shift
# element 1 in list, for shiftmap 0 is deriv-based, 1 is integ-based
QM.plot_2_maps(shiftmaps[1][:,:,0], amplmaps[1][:,:,0]) 

# quick compare of deriv based and integ based peak shifts (wrong label on plot 2)
QM.plot_2_maps(shiftmaps[2][:,:,0], shiftmaps[2][:,:,1])

# Make histogram of peak shift/ charging values
QM.plothisto(chargemap, 15) # number of bins
partial=specimage[:,:,0]

# Summary statistics describing charging behavior
scipy.stats.describe(shiftmaps[1][:,:,1], axis=None)
scipy.stats.describe(amplmaps[2][:,:,3], axis=None) # element 2 sm-diff ampl.
scipy.stats.describe(integmaps[3][:,:,4], axis=None) # integ-based ampl.
scipy.stats.describe(newmap, axis=None)

# Mask weakest or strongest signal subsets and look at spatial distribution
weird=np.ma.masked_where(chargemap<=np.percentile(chargemap, 5), chargemap)
lowvals=np.ma.masked_where(chargemap<=135, chargemap)
weird=np.ma.masked_where(np.logical_or(chargemap<=135, chargemap>=165), chargemap)
highvals=np.ma.masked_where(chargemap>=10, chargemap)
highvals=np.ma.masked_where(chargemap>=np.percentile(chargemap, 95), chargemap)
realvals=np.ma.masked_where(chargemap==np.nan, chargemap)
weak=np.ma.masked_where(peakamplmap <=np.percentile(peakamplmap, 10), peakamplmap)
strong=np.ma.masked_where(peakamplmap >=np.percentile(peakamplmap, 90), peakamplmap)
np.ma.count_masked(lowvals) # counts number of masked values

# Apply median filter (or uniform filter) to raw chargemap
charge_medfilt=medfilt(chargemap, kernel_size=3)
peakampl_medfilt=medfilt(peakamplmap, kernel_size=3)
smoothspec=QM.uniformfilter(specimage, size=3) # odd window/kernel size (1 is no transform)

# Spatial plotting of various masked subsets
fig, axes = plt.subplots(nrows=1, ncols=1, squeeze=False)
plt.imshow(chargemap)
plt.imshow(peakamplmap)
plt.imshow(weak)
plt.imshow(strong)
plt.imshow(highvals)
plt.imshow(weird)
plt.imshow(newamplmap)

# Histogram plot of charging values (sometimes reveals erroneous ones)
QM.plothisto(chargemap, 15)
QM.plothisto(newmap, 15)

# Interactive plot of single chosen pixel (tk calling plotpixels)
kwargs={}
kwargs=plotpix_tk(specimage, energy, Elemdata, spectralregs, amplmaps, integmaps, 
           shiftmaps, AESquantparams, **kwargs)

# Look at subset of masked pixels (normally checking underlying spectra from extrema)
pixlist=QM.makepixlist(highvals) # get masked pixels back
pixlist=makepixlist(highvals)
pixlist=makepixlist(weird)
pixlist=QM.pickrandompixels(specimage, 5)
pixlist=[[0,13]]

# plot report of counts, deriv or both for subset of chosen pixels
QM.pixelreport_tk(specimage, pixlist, energy, Elemdata, spectralregs, amplmaps, 
                  integmaps, shiftmaps, AESquantparams, **kwargs)

# Replace any masked pixels with median filtered values (what about edges?) 
newmap=QM.replacemaskpix(lowvals, charge_medfilt)
# Can also filter peak amplitude map with bad shift values from charge map (replace again w/ median)
chargemap=QM.replacemaskpix2(peakampl, peakampl_medfilt, lowvals) # filtering w/ 3rd array
chargemap=replacemaskpix2(chargemap, charge_medfilt, highvals) 

# if bad values at map's edge, one can mask again (by value) and replace w/ h
lowvals=np.ma.masked_where(newmap<=135, newmap)
newmap=np.ma.filled(lowvals, 150)

plt.imshow(newmap)
plt.imshow(chargemap)
plt.imshow(charge_medfilt)
# Save of modified chargemap (O peak shift w/ bad pix modifications
np.save('Acfer094_Omap_area8_chargemap.npy', newmap) 
np.save('Acfer094_Omap_area8_peakamplmap.npy', peakamplmap)

# Load each element's s7d7 amplitude into amplmapdict (usually after remove/replace of defects)
amplmapdict={}
amplmapdict.update({'O': peakamplmap})
amplmapdict.update({'Si': peakamplmap})

# Show peak amplitudes (s7d7 style ) for all elements
showpeakampls(amplmapdict, Elements)

ratio=ratioimages(amplmapdict, 'O','Fe')
fig, axes = plt.subplots(nrows=1, ncols=1, squeeze=False)
plt.imshow(ratio)

# Save charging for each pixel, reorder areas by charge magnitude and create new spatial area files
QMpixarray=reorder_spatial_array(QMpixarray, chargemap, '50x50_scrambled')

# Create series of shifted multiplex files and associated Autotool
QMpixarray, Autotool=QM.multi_shift_setup(Elements, AESquantparams, QMpixarray)

QMpixarrayshift=QM.multi_shift_setup(Elements, AESquantparams, QMpixarray)
QMpixarrayshift=multi_shift_setup(Elements, AESquantparams, QMpixarray)

QMpixarrayshift.to_csv('QMpixarray_shifted.csv',index=False)

'''create autotool to map in order of observed charge (separate shifted multiplex files)
instead of normal mapping in X, Y pixel order
'''
# tk interface for making shifted multiplexes
shiftQMarray=make_shifted_multiplexes(QMpixarray)

# Find peak shift from ideal position (in AESquantparams) across entire spectral image
# shiftdict holds peak shift array (50x50x1) for each element
AESqp=QM.shiftAESparams(AESquantparams, 150)
shiftdict, peakstats=QM.calcshifts(specimage, Elemdata, AESquantparams)
shiftdict, peakstats=calcshifts(specimage, Elemdata, AESqp)

shiftdict, peakstats=QM.calcshifts(smoothspec, Elemdata, AESquantparams)

shiftarr=shiftdict.get('O',[]) # 50x50 array of values with each peak shift

# Plotting of shift images (should reveal spatially dependent charging)
# Use single value for shift?  apply global value?  use local value but with smoothed image
# Note: shift can be negative but image display vals all pos 
# trends correct but display numbers are different
QM.showshiftimages(shiftdict, peakstats)

test=shiftdict.get('O')
# Perform background fits over all elemental regions (returns array of slopes & intercepts)
backarray=QM.calcbackgrounds(specimage, energy, Elemdata)
backarray2=QM.calcbackgrounds(smoothspec, energy, Elemdata) # run same on background averaged
np.save('Acfer094map_backfits.npy', backarray) 
np.save('Acfer094map_backfits_smoothed.npy', backarray) 

# Interface for plotting multiplex spectra from chosen list of pixels 
kwargs={}
shiftdict={}
backarray=pd.DataFrame()

# retrieve and plot single element numpy array 
Ocnts=sumdict.get('O',[])[:,:,2]
Ocnts=sumdict2.get('O',[])[:,:,2]
scipy.stats.describe(Ocnts)

# plot all intensity maps using integrated, subtracted peak data
QM.showintegimages(sumdict, Elemdata)
QM.showintegimages(sumdict2, Elemdata)

# After peak amplitude maps are available for each peak (w/ bad vals mask/replaced) can run calccomposition 

#%% Combine separate spe files to single quantmap (QM) file
# Combine separate spe/csv files into single quant map file containing all areas (AugerParamLog autosaved)
AugerParamLog=QM.combineQMdata(AugerParamLog,'139-143',QMname='')

# Renumber/rename spatial areas for quant map files
SpatialAreas=pd.read_csv('spatialareaslog.csv') 
SpatialAreas=QM.renumberQMareas(SpatialAreas, '148-167',QMname='') # copy spatial areas for combined QM file, autosaved 
#%%
includemask=AugerParamLog['Comments'].str.contains('quantmap', case=False, na=False) # for selecting files if doing any batch reports
qmlist=AugerParamLog.loc[includemask]
qmlist=qmlist[2:]

Elements=['Si', 'Mg','Fe'] # list of multiplex elements (same as used for quant)

Mgmap=elementmaps[0][1] # extra
Mgmap.mean()
Mgmap.max()
Mgmap.histogram()
hist=np.histogram(Mgmap, bins=15)
plt.imshow(hist)

# Create amplitude maps for given elements list of quant map for selected filenumber
elementmaps=QM.createampmaps(221222, Elements, Smdifpeakslog, SpatialAreas) # returns separate 512x512 intensity map 
elementmaps=createampmaps(148167, Elements, Smdifpeakslog, SpatialAreas) 

# Plot the above set of element maps in single figure (optional save)
QM.plotmaps(elementmaps, Elements, savename='Mg_map_221222.jpg') # separate heatmaps of each element's sm-diff amplitude
QM.plotmaps(elementmaps, Elements, savename='') # plot only 

# Make and save image with elemental hotspot positions indicated, saved as savestr+combined filenumber
QM.showsubareas(Mgrich, SpatialAreas, image='tr184.187.jpg', savestr='Mgrich',label=True)

# Overlay selected element map on SE image with variable opacity (alpha) and optional direct save
QM.overlayimage(elementmaps, 'Mg', 'tr184.220.jpg', savename='Mgrich_area331.jpg', alphaval=0.3)

# Probably easiest to just have quant from these pixels in the same master smdifpeak logs
# Also could easily have optional separate peak log for quantmaps

# Perform sm-diff quant on qmlist in usual way
Elements=['Mg','In','Si'] # same set from quantmap multiplex
Backregs=[121,200,405,800,1475,1850] 
# LOAD AESQUANTPARAMS (peak positions, kfactors etc.)
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Params\\AESquantparams.csv', encoding='utf-8')
QMpeakslog=AESsmquant.smdifbatchquant(qmlist, Elements, Backregs, AESquantparams)

# Slice image and find/display subareas rich in selected element

