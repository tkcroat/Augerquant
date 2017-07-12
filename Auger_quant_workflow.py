# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:49:58 2016
Auger peak finding and quantitative routines ... batch processing
@author: tkc
First get it working for single file.
"""
#%%
import pandas as pd
import os, sys
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules')

import Auger_smdifquant_functions as AESsmquant
import Auger_integquant_functions as AESintquant
import Auger_utility_functions as AESutils
import Auger_plot_functions as AESplot
''' AESsmquant contains functions related to peak finding in smooth-differentiated spectra
whereas AESquant contains background fitting and integration over peaks for direct from counts '''
# import Auger_integquant_functions as AESquant

#%%  SET DATA DIRECTORY and load main files (and sub files if desired)
os.chdir('C:\Temp\Auger')

# Reloading of key log files (required after alterations to csv files or restart)
AugerParamLog, spelist, Smdifpeakslog, Integquantlog, Backfitlog, AESquantparams= AESutils.loadmainfiles()
AugerParamLog, spelist, Smdifpeakslog, Integquantlog, Backfitlog, AESquantparams= loadmainfiles()
subspelist, Smdifpeakslogsubs, Backfitlogsubs, Integquantlogsubs=AESutils.loadsubfiles()
#%% Optional filtering of datasets (according to any custom criteria)
spelist=spelist[(spelist['Filenumber']>=10000)]
spelist=spelist[spelist['Filename'].str.contains('54-9')]

# Interactive choice of spe files by filenumber
myfiles=AESutils.pickspectraGUI(spelist)
#%% MAIN BATCH SMDIFF (DERIVATIVE) QUANT LOOP from smooth-differentiated peaks  (via Multipak S7D7 algorithm) of all selected files in dataframe

#Element and background region setup
Elements=AESutils.pickelemsGUI(AESquantparams) # interactive element selection
Elements=AESutils.pickelemsGUI(AESquantparams) # interactive element selection

# MAIN BATCH QUANT LOOP function which returns Smdifpeakslog dataframe with amplitude of every peak in above list
# reprocess=False ... loads prior smdifpeakslog and skips files already processed
Smdifpeakslog=AESsmquant.smdifbatchquant(spelist, Elements, AESquantparams, reprocess=False)
Smdifpeakslog=smdifbatchquant(spelist, Elements, AESquantparams, reprocess=False)
Smdifpeakslogsubs=AESsmquant.smdifbatchquant(subspelist, Elements, AESquantparams, reprocess=True)

# Save smdifpeaks log result run below (this is not auto-saved!)
Smdifpeakslog.to_csv('Smdifpeakslog.csv', index=False) # saves above peaks file
Smdifpeakslogsubs.to_csv('Smdifpeakslog_subs.csv', index=False)

# get peak statistics
peakstats=AESsmquant.getpeakstats(Smdifpeakslog, Elements)
peakstatssubs=AESsmquant.getpeakstats(Smdifpeakslogsubs, Elements)
peakstats.to_csv('peakstats.csv', index=False)

# histogram plots
AESplot.plothist(Smdifpeakslog, spelist, Elements, col='Amplitude')
AESplot.plothist(Integquantlog, spelist, Elements, col='Adjcnts')

#%% DIRECT INTEGRAL QUANT METHOD 
# Now for integral method directly on counts (but guided by peak positions as determined above with smooth-diff spectrum)
Elements=AESutils.pickelemsGUI(AESquantparams) # interactive element selection

# MAIN BATCH INTEGRAL METHOD fitting and direct integration loops (returns background fitting params and integrated quant results)
# WARNING ... overwrite=True will overwrite background fit and peaks columns in source csv file (counts obviously unchanged)
Backfitlog, Integquantlog = AESintquant.integbatchquant(spelist, Smdifpeakslog, AESquantparams, Elements, reprocess=True, overwrite=True)
# set overwrite to false if redoing quant to add a group of new elements
Backfitlog, Integquantlog = AESintquant.integbatchquant(spelist, Smdifpeakslog, AESquantparams, Elements, reprocess=False, overwrite=False)

# perform this explicitly on sub spe files
Backfitlogsubs, Integquantlogsubs = AESintquant.integbatchquant(subspelist, Smdifpeakslogsubs, AESquantparams, Elements, reprocess=True, overwrite=True)
# the background fit results and integration results must be manually saved (not-autosaved)
Backfitlog.to_csv('Backfitlog.csv', index=False)
Integquantlog.to_csv('Integquantlog.csv', index=False)

Backfitlogsubs.to_csv('Backfitlog_subs.csv', index=False)
Integquantlogsubs.to_csv('Integquantlog_subs.csv', index=False)

# Calculated (or recalculate after AESquantparams change) adjusted counts based on k-factor2 for direct integration method
Integquantlog=AESintquant.calcadjcounts(Integquantlog, AESquantparams, sig=2, kerrors=True) 
#%% CREATE PLOT REPORTS TO VIEW DERIVATIVE AND INTEGRAL METHOD RESULTS
# Reloading of key log files (required after alterations to csv files or restart)
AugerParamLog, spelist, Smdifpeakslog, Integquantlog, Backfitlog, AESquantparams= AESutils.loadmainfiles()
subspelist, Smdifpeakslogsubs, Backfitlogsubs, Integquantlogsubs=AESutils.loadsubfiles()

# Choose elements list to plot and label
plotelems=AESutils.pickelemsGUI(AESquantparams) # interactive element selection
plotelems=['30-2130'] # can also just enter a range
plotelems.append('100-500') # can also include eV range along with above list

spelist=AESutils.pickspectraGUI(spelist)# interactively grab subset of filenumbers
# Grab compositional subsets for separate plot reports (example)
Mgrich=Smdifpeakslog[(Smdifpeakslog['PeakID']=='Mg') & (Smdifpeakslog['Amplitude']>5000)]

# General smooth-differentiated plot report for selected elements
AESplot.reportSD(spelist, Smdifpeakslog, plotelems, AESquantparams) # default name is SDplots_report_09Jun17.pdf
AESplot.reportSD(Mgrich, Smdifpeakslog, plotelems, AESquantparams, PDFname='Mgrich_areas.pdf')
AESplot.reportSD(subspelist, Smdifpeakslogsubs, plotelems, AESquantparams, PDFname='SDreportsubs_3Jan17.pdf')

# Plot report for direct peaks (option to include direct background fits from integral quant method)

# Create a counts/background plot report over selected files 
kwargs={} # kwargs fed back to set defaults in case of rerun
kwargs=AESplot.countsbackreport_tk(spelist, Elements, Backfitlog, AESquantparams, **kwargs) 
kwargs=countsbackreport_tk(spelist, Elements, Backfitlog, AESquantparams, **kwargs) 
kwargs=AESplot.countsbackreport_tk(subspelist, Elements, Backfitlogsubs, AESquantparams, **kwargs) 

# Plot report of derivative and counts on same page (top and bottom)
kwargs={}
kwargs=AESplot.countsderivreport_tk(spelist, Elements, Smdifpeakslog, Backfitlog, AESquantparams, **kwargs)
kwargs=countsderivreport_tk(spelist, Elements, Smdifpeakslog, Backfitlog, AESquantparams, **kwargs)

# Plot report of counts - background for selected peaks (all areas )
myfiles=AESutils.pickspectraGUI(spelist)
AESplot.reportpeaksall(myfiles, plotelems, AESquantparams, PDFname='Peaks_report.pdf')
# If erroneous spectra are found, add phrase exclude to comments in AugerParamLog for this filenumber... reload and then can be excluded from future processing
# For refitting problematic subset and partial redo of those quant results ... see Auger_quant_alt.py

# Single interactive plot of selected filenumber(s) .. calls AESplot1
kwargs={}
kwargs=AESplot.AESplot_tk(spelist, plotelems, Smdifpeakslog, Backfitlog, AESquantparams, **kwargs)
kwargs=AESplot_tk(spelist, plotelems, Smdifpeakslog, Backfitlog, AESquantparams, **kwargs)

#%% COMPUTING ELEMENTAL COMPOSITIONS based on derivative quant method or integral quant method 

# Calculate (or  recalculate after AESquantparams change the adjusted amplitude for deriv method
Smdifpeakslog=AESsmquant.calcamplitude(Smdifpeakslog, AESquantparams)
Smdifpeakslogsubs=AESsmquant.calcamplitude(Smdifpeakslogsubs, AESquantparams)
Smdifpeakslog.to_csv('Smdifpeakslog.csv', index=False) # save above result (not auto-saved)
Smdifpeakslogsubs.to_csv('sub\\Smdifpeakslog_subs.csv', index=False) 

# Now calculate composition (at %) for selected elements
Elements=AESutils.pickelemsGUI(AESquantparams) 

# Calculate composition based on derivative method (smooth-diff comp)
# Non-zero threshold excludes weak lines from quant results (threshold=1 means excluded if peak amplitude < noise amplitude
Smdifcomp=AESsmquant.calccomposition(spelist, Smdifpeakslog, Elements, threshold=0.0)
Smdifcompsubs=AESsmquant.calccomposition(subspelist, Smdifpeakslogsubs, Elements, threshold=0.0)

Smdifcompsub=Smdifcompsubs[Smdifcompsubs['AESbasis']>0]
# saves compositions and underlying elements/k-factors to auto-named xls file
AESutils.writecomps(Smdifcomp, AESquantparams, Elements)
# read back previously saved compositional determinations (smdiffcomp or integcomp)
Smdifcomp=pd.read_excel('Smdiffcomp_03Jan17.xlsx', sheetname='smdiffcomp')
Smdifcompsubs=pd.read_excel('sub\\smdiffcomp_04Jan17.xlsx', sheetname='smdiffcomp')

Integcomp=pd.read_excel('integcomp_02Jan17.xlsx', sheetname='integcomp')

Smdifcompsubs=pd.read_csv('Smdifcompsubs.csv')
Mgoutlierquant=calccomposition(Mgoutliers, Smdifpeakslog, Elements, threshold=0.0)
Mgoutlierquant.to_csv('Mgoutlier_quant.csv',index=False)

# Compare compositions from combine-averaged spe file and the underlying component spe files 
Compresult, Fullcompresult=AESsmquant.compareavg_subs(Smdifcomp, Smdifcompsubs)
Compresult, Fullcompresult=compareavg_subs(Smdifcomp, Smdifcompsubs)
Compresult.to_csv('smdiff_avgcomp_vs_subs_5Jan17.csv',index=False)
Fullcompresult.to_csv('Compare_comp_full.csv',index=False)
mycols=['Feampl','Feampl1','Feampl2','Feampl3']
Fesub=Fullcompresult[mycols]
Fesub.to_csv('temp.csv', index=False)

Smdifcomp.to_csv('Smdifcomp.csv',index=False)
Smdifcompsubs.to_csv('Smdifcompsubs.csv',index=False)
Ferich=Smdifcomp[Smdifcomp['Feampl']>1500]
Sirich=Sirich.sort_values(['sigSi'], ascending=False)
Sirich=Sirich[0:20]

# Generic plotting of above compositional comparisons with outliers returned
Compresult=pd.read_csv('Smdifcompavg_vs_subs.csv')
xcol='Feampl'
ycol='Feamplavg'
xcol='Mg'
ycol='Mgavg'
outliers=AESplot.scatplot(Compresult, xcol, ycol, thresh=0.1)

# CALCULATE COMPOSITIONS from direct integration quant (normally includes at.% error)
# First calculate or recalculated elemental peak adjusted by sensitivity k-factor
Integquantlog=AESintquant.calcadjcounts(Integquantlog, AESquantparams, sig=2, kerrors=True) 
Integquantlogsubs=AESintquant.calcadjcounts(Integquantlogsubs, AESquantparams, sig=2, kerrors=True) 
# kerrors true propogates estimated error in k-factor defined in AESquantparams into all errors 
# sig refers to incl 2 sig errors for kfactor
Integquantlog.to_csv('Integquantlog.csv', index=False)
Integquantlogsubs.to_csv('sub\\Integquantlog_subs.csv', index=False)

# Calculate compositions in at. % based on integral method
Integcomp2=calccomp(spelist, Integquantlog, Elements, AESquantparams)
Integcomp=AESintquant.calccomp(spelist, Integquantlog, Elements, AESquantparams)
Integcompsubs=AESintquant.calccomp(subspelist, Integquantlogsubs, Elements, AESquantparams)
Integcomp=calccomp(spelist, Integquantlog, Elements, AESquantparams)

# Write compositions to xls (along with k-factor details)
AESutils.writecomps(Integcompsubs, AESquantparams, Elements)

# save above result (not autosaved)
Integcomp.to_csv('Integcomp.csv', index=False)
Integcompsubs.to_csv('Integcomp_subs.csv', index=False)

Compresult, Fullcompresult=AESintquant.compareavg_subs(Integcomp, Integcompsubs)
Compresult, Fullcompresult=compareavg_subs(Integcomp, Integcompsubs)
Compresult.to_csv('integ_avgcomp_vs_subs_5Jan17.csv', index=False)

#%% COMPARISON PLOTS OF DERIVATIVE VS INTEGRAL COMPOSITIONS
# Reload composition files for sm diff (deriv) and integral methods 
Smdifcomp=pd.read_csv('Smdifcomp.csv', encoding='cp437')
Smdifcomp=pd.read_excel('Smdiffcomp_08Jan17.xlsx', sheetname='smdiffcomp')
Smdifcompsubs=pd.read_excel('sub\\smdiffcomp_04Jan17.xlsx', sheetname='smdiffcomp')
# Read back calculated and saved compositions
Integcomp=pd.read_excel('integcomp_02Jan17.xlsx', sheetname='integcomp')
Integcompsubs=pd.read_excel('sub\\integcomp_02Jan17.xlsx', sheetname='integcomp')

Smdifcomp=AESutils.dropexcluded(Smdifcomp,spelist) # refilter any df by filenumber based on filtered spe list
Integcomp=AESutils.dropexcluded(Integcomp,spelist) 
Smdifcomp.to_csv('temp.csv', index=False)

# Scatter compositional comparison plots
# Directly compare deriv and integ on exact same data files with multi-element scatter plots
compdata, outliers=AESplot.scattercompplot_tk(Smdifcomp,Integcomp, Elements)
compdata, outliers=AESplot.scattercompplot_tk(Smdifcomp, Smdifcomp2, Elements)

compdata.to_csv('compdata_full.csv', index=False)
compdata=pd.read_csv('compdata_full.csv') # read back temp saved data

# Various ways of plotting/saving subsets from above comparative compositional data
SiN=outliers[outliers['%Sib']>0.4] # subset of outliers from sulfur plot
Feoutliers=outliers[(outliers['Feb']>100000) & (outliers['Element']=='Fe')]
Sioutliers=outliers[outliers['Element']=='Si']

Mgoutliers=Mgoutliers.sort_values(['Pval'])
temp=Mgoutliers[0:11]

SiN.to_csv('SiN.csv', index=False)
WeakFe3_5=compdata[(compdata['sigFeb']>3)&(compdata['sigFeb']<5)] # filtefsring based on elemental line significance level (# sig deviation)
WeakFe3_5=WeakFe3_5[WeakFe3_5['Element']=='S'] # eliminate 5x duplicate

compset.to_csv('compset_full.csv', index=False)
outliers.to_csv('outliers.csv', index=False)

# excluding weaker spectra based on AESbasis (which sums strength of all peaks present)
compdata=compdata[compdata['Element']=='Si'] 
compset=compdata[compdata['AESbasis']>200] # remove weak spectra
compset=compset[compset['Si']<2000]
compset=compset[compset['Ti']>800]

outliers=compdata[compdata['AESbasis']>200]
thismask=(compdata['Fe']<20)&(compdata['Feb']<20)
compdata=compdata.loc[~thismask] # knock out 15 more with overlall weak Fe

Smdifcomp=Smdifcomp[Smdifcomp['Filenumber']==106]

# make a standard ternary plot
ternelems=['Fe', 'Si', 'Mg+Ca']
ternelems=['Ga', 'Pt', 'Mg+Si+In']
AESplot.plotternary(Smdifcomp, ternelems)

#%% COMPARE COMPOSITIONAL CONSISTENCY BETWEEN DUPLICATED SPECTRA FROM SAME SAMPLE
integcomp=pd.read_csv('integcomp.csv', encoding='cp437')

duplicates=integcomp.duplicated(['Sample'], keep=False)
duplicates=integcomp.loc[duplicates]
duplicates=duplicates.sort_values(['Sample'])
duplicates=duplicates[duplicates['AESbasis']>3000] # remove weak data (w/ integ AES basis <3000)
duplicates=duplicates.reset_index()  # all possible duplicates incl. multi-area ones

duplicates=AESutils.findduplicates(integcomp) # returns duplicated but removes those with only multiple areas in single filenumber (maybe heterogeneous)

duplicate_dataset=AESutils.compareduplicatecomps(duplicates, Elements) # returns avg comp in cases w/ multiple determinations
duplicate_dataset.to_csv('Duplicated_compositions.csv', index=False)

# Plot report with ratios from duplicates 
ternelems=['Fe', 'Si', 'Mg']
ternelems=['Fe', 'S', 'Mg']
ternelems=['Mg+Ca','Fe', 'Si']

plotduplicatesternary(duplicates, ternelems)

# Constructing plots of various outliers (found in above duplicate dataset)
thissample=spelist[spelist['Sample']=='f2-2c12']

AESplot.reportderivcntall(Tirich, Smdifpeakslog, plotelems, AESquantparams, backfitdf=Backfitlog, PDFname='Tirich.pdf')
#%% COMPARE AUGER DERIVED COMPOSITIONS WITH SEM-EDX OR TEM_EDX compositions of same craters
SEMcomp=pd.read_csv('C2010W_SEMquant.csv', encoding='cp437')
TEMcomp=pd.read_csv('C2010W_FIBTEMcomp.csv', encoding='cp437')

samplelist=['f2-1c5', 'f2-1c6', 'f2-1c7', 'f2-2c5', 'f7-2c2', 'f7-2c3', 'f7-3c2']

#%%  OTHER DATA FILTERING AND PLOTTING 
# SELECTING SUBSETS BASED ON ELEMENTAL COMPOSITION
# Get spectra with strongest peaks (top ten) for any element
Carich=CompC2010WselectavgFe.sort('Ca',ascending=False)
Carich=Carich.head(10)
Mgrich=CompC2010WselectavgFe.sort('Mg',ascending=False)
Mgrich=Mgrich.head(10)

compFeredo=C2010WAESmajor2[C2010WAESmajor2['Fe']>0]
compFeredo2=C2010WAESmajor2Fe2[C2010WAESmajor2Fe2['Fe']>0]

C2010sel_integ.to_csv('C2010sel_integ2.csv', index=False)
C2010sel_smdif.to_csv('C2010sel_smdif.csv', index=False)
compdataset.to_csv('Quantcomp_set.csv', index=False)

# Reload of previously calculated Auger compositions  `
C2010sel_integ=pd.read_csv('C2010sel_integ.csv', encoding='utf-8')
C2010sel_smdif=pd.read_csv('C2010sel_smdif.csv', encoding='utf-8')

# Scatterplot comparison of two sets of compositions 
# uses element names (not peak names like Fe2); basis True compares basis not at.%

compdataset=scatterratioplot(C2010sel_smdif, C2010sel_integ, 'Si', 'Mg', basis=False)

# Save altered log files
C2010sel_integ.to_csv('C2010sel_integ2.csv', index=False)

compdataset.to_csv('SEM_AESinteg_comparison.csv', index=False)
Smdifpeakslog.to_csv('Smdifpeakslog.csv', index=False)

# if multiple spectra exist, choose one with largest AESbasis for given Elements list
spelist=AESsmquant.getbestduplicate(spelist, Smdifpeakslog, Elements) 
AESutils.outputduplicates(spelist,'Sample') # simple utility to output samples with duplicate spectra into console

AESutils.copyselectfiles(spelist,'best') # moves all underlying data into new subfolder (named by string argument)

# Copying ten best spectra (based on AESbasis for chosen elements into another folder
spelist=spelist.sort('AESbasis',ascending=False)
best=spelist.head(10)
AESutils.copyselectfiles(best,'best') #

# Comparing effects of different noise thresholds in sm-diff calculation
integcompwithCa=integcomp[integcomp['Ca']>0]
integcompwithCa=integcompwithCa.sort_values(['Ca'], ascending=True)
temp=integcompwithCa.head(10)
integcompwithCa.to_csv('Ca_containing.csv')

# Subset of higher mag images (cropped)
AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
images=AugerParamLog[AugerParamLog['Filename'].str.contains('.sem')]
images=images[images['FieldofView']<25]
images['Filename']=images['Filename'].replace('.sem','.jpg', regex=True)
images=images.drop_duplicates('Sample')
# only single images 
images=images[~images['Sample'].str.contains('and')]
images=images[~images['Sample'].str.contains(',')]
images=images[~images['Sample'].str.contains('/')]
AESutils.copyselectfiles(images,'images')
# Drop multi-image craters
singleim=images[~images['Sample'].str.contains('and')]

cratersizes=pd.read_excel('C:\\Users\\tkc\Desktop\\Research\\Stardust_craters\\C2010W_crater_logbook.xlsx', sheetname='Crater_positions')
# Autocrop each image based on max diam measurement
 