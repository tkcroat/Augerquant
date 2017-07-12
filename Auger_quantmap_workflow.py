# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 17:10:19 2016

@author: tkc
"""
import pandas as pd
import numpy as np
import sys

if 'C:\\Users\\tkc\\Documents\\Python_Scripts' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts')
import Auger_smdifquant_functions as AESsmquant
import Auger_quantmap_functions as QM

#%% CREATE PHI FILES FOR AUTOTOOL, SPATIAL AREAS, MULTIPLEX CONDITIONS
# A few of these are also stored in Auger import main (to allow QM data combination prior to quant)

AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437')
# Generate square pixel arrays of specified dimension and create associated .phi files (margin as percentage of field)
kwargs={'imageregint':2} # image registration interval (value of 1 means register every 20 areas)
QMpixarray, autotool=QM.makesquarearray(margin=0.2, arraysize=100,basename='array100x100file', **kwargs)
QMpixarray, autotool=makesquarearray(margin=0.2, arraysize=100,basename='array100x100file', **kwargs)

QM.writemultiplex('QM_multiplex.phi', dwelltime=20, numcycles=3, reginterval=1, regmode='Areas')
QM.writemultiplex('QM_multiplex_100x100.phi', dwelltime=5, numcycles=1, reginterval=50, regmode='Areas')
QM.writeautotool(autotool, 'QM_autotool20x20.phi') # write phi autotool file with correct spatial area names

QMpixarray.to_csv('QMpixarray_acfer094area14.csv', index=False)

# Make annotation of image with overlaid mapped region (last arg is margin )
QM.showQMregion('Acfer094.157.jpg', AugerParamLog, 0.2)

# Reload pixarray definition file 
QMpixarray=loadQMpixarray()
# Associate multiplex spe data files with correct pixels from quantmap (after data collection)
QMpixarray=linkfilename(QMpixarray, 'Acfer094map', startnum=101)
# Make 3D numpy array (specimage) with all spectral data; energy x values in energy list
specimage, energy =makeQMarray(QMpixarray, 'Acfer094map')

# 
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

