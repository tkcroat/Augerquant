# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:49:58 2016
Auger peak finding and quantitative routines ... batch processing
@author: tkc
First get it working for single file.
"""
#%%
import pandas as pd
import numpy as np
import os, sys, shutil, glob, re
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules')

import Auger_smdifquant_functions as AESsmquant
import Auger_integquant_functions as AESintquant
import Auger_utility_functions as AESutils
import Auger_plot_functions as AESplot
''' AESsmquant contains functions related to peak finding in smooth-differentiated spectra
whereas AESquant contains background fitting and integration over peaks for direct from counts '''
# import Auger_integquant_functions as AESquant

#%%  SET DATA DIRECTORY 
os.chdir('C:\Temp\Auger')
# load master AugerParamLog for this data directory (created during batch import)
AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437') # use Excel compatible encoding cp437 instead of utf-8 
AugerParamLogsubs=pd.read_csv('sub\\Augerparamlog_subs.csv', encoding='cp437')
AugerParamLog=AugerParamLog.drop_duplicates(['Filenumber'])
AugerParamLog.to_csv('Augerparamlog.csv',index=False)
AugerParamLog=updatepaths(AugerParamLog)

# Reloading of key log files (required after alterations to csv files or restart)
AugerParamLog, spelist, Smdifpeakslog, Integquantlog, Backfitlog, AESquantparams= AESutils.loadmainfiles()
Smdifpeakslogsubs, Backfitlogsubs, Integquantlogsubs=AESutils.loadsubfiles()
#%%
# FILTERING DATASET (according to any custom criteria)
#select all spectral files for quant from log 
spelist=AugerParamLog[(AugerParamLog['Areas']>=1)]  # remove image & map files 
spelist=AugerParamLog[(AugerParamLog['Filenumber']>=10000)]
spelist=spelist[(spelist['Filenumber']>=10000)]

subspelist=AugerParamLogsubs[(AugerParamLogsubs['Areas']>=1)]
subspelist=subspelist.loc[~subspelist['Comments'].str.contains('exclude', case=False, na=False)]

# remove files excluded by log
spelist=spelist.loc[~spelist['Comments'].str.contains('exclude', case=False, na=False)]
temp=spelist.loc[[85]]
tempspe=subspelist[0:51]

excllist=np.ndarray.tolist(Excluded.Filenumber.unique())
spelist=AugerParamLog.loc[~AugerParamLog.Filenumber.isin(excllist)]

# Remove files excluded for various reasons (has exclude in comments fields)
Excluded=AugerParamLog2[AugerParamLog2['Comments'].str.contains('exclude', case=False, na=False)]
# Exclude from spelist any spe files with "exclude" in comments
spelist=spelist.loc[~AugerParamLog['Comments'].str.contains('exclude', case=False, na=False)]
includemask=spelist['Comments'].str.contains('avg', case=False, na=False)
spelist=spelist.loc[includemask]
includemask=spelist['Comments'].str.contains('quantmap', case=False, na=False) # quantmap is combined file, subqm is sub array files
# select on sub spe files (combined via avg into other files)
subspelist=spelist.loc[spelist['FilePath'].str.contains('sub', case=False, na=False)]

spelist=spelist.loc[~spelist['FilePath'].str.contains('sub', case=False, na=False)]
spelist=spelist.sort_values(['Filenumber'], ascending=True)

# Select sublist of files by number for quant reprocess quant
filelist=[114]
spelist=spelist[spelist.Filenumber.isin(filelist)]

# check for missing data files
filelist=glob.glob('*.csv')
AESutils.checklog(filelist, AugerParamLog) 

tempspe=spelist[0:21]  # different slices from Auger parameters log
combinelist=combinelist.loc[320:330] 
#%% MAIN BATCH SMDIFF (DERIVATIVE) QUANT LOOP from smooth-differentiated peaks  (via Multipak S7D7 algorithm) of all selected files in dataframe

#Element and background region setup
Elements=['S','C','Ca','Ti','O','Fe1','Fe2','Fe','Na','Mg','Al','Si'] # using dominant Auger peak for each element... see AESquantparams.py for associated energies
Elements=['S','C','Ca','Ti','Ti2','O','Fe1','Fe2','Fe','Na','Mg','Al','Si','In','F'] 
Elements=['S','C','Ca','Ti','N','O','Fe1','Fe2','Fe','Ti2','Mg','Al','Si']
Elements=['S','C','Ca','Ti','O','Fe1','Fe2','Fe','Ti2','Mg','Al','Si','Cs','Cs2']
Elements=['Mg']
Elements=['C','O','Si', 'W']
Elements=['Pt','Ga','Si', 'Mg','In','O']
Elements=['Si', 'Mg','In','Fe']


# list of eV values where background is searched for smooth-diff peaks (gives method of determining significance of peak amplitudes)

# LOAD AESQUANTPARAMS (peak positions, kfactors etc.)
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Params\\AESquantparams.csv', encoding='utf-8')
# Can also load a local version (if modifications are needed for this sample's quant
AESquantparams=pd.read_csv('AESquantparams.csv', encoding='utf-8')

# MAIN BATCH QUANT LOOP function which returns Smdifpeakslog dataframe with amplitude of every peak in above list
# reprocess=False ... loads prior smdifpeakslog and skips files already processed
Smdifpeakslog=AESsmquant.smdifbatchquant(spelist, Elements, AESquantparams, reprocess=True)
Smdifpeakslogsubs=AESsmquant.smdifbatchquant(subspelist, Elements, Backregs, AESquantparams, reprocess=True)
Smdifpeakslogsubs=smdifbatchquant(subspelist, Elements, Backregs, AESquantparams, reprocess=True)

Smdifpeakslogsubs=smdifbatchquant(subspelist, Elements, Backregs, AESquantparams, reprocess=True)

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
# now for integral method directly on counts (but guided by peak positions as determined above with smooth-diff spectrum)
Elements=['S','Ca','Fe2','Fe','Mg','Si'] # choose major elements
Elements=['Fe2','Fe','Mg','Si','Au']
Elements=['Ca','N', 'Mg','Si','Ti','Fe','Fe2']
Elements=['S', 'Ca','Mg','Si','Ti','Fe','Fe2']
Elements=['S','Ca','Fe2','Fe','Mg','Si','Ti','Ti2']
Elements=['Pt','Ga','Si', 'Mg','In']
Elements=['S']

# reload smooth-diff peak positions (if not already open)
Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437') 
Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')

# load peak positions, kfactors, integration params, etc.
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Params\\AESquantparams.csv', encoding='utf-8')


# MAIN BATCH INTEGRAL METHOD fitting and direct integration loops (returns background fitting params and integrated quant results)
# WARNING ... overwrite=True will overwrite background fit and peaks columns in source csv file (counts obviously unchanged)
Backfitlog, Integquantlog = AESintquant.integbatchquant(spelist, Smdifpeakslog, AESquantparams, Elements, reprocess=True, overwrite=True)
Backfitlog, Integquantlog = integbatchquant(spelist, Smdifpeakslog, AESquantparams, Elements, reprocess=True, overwrite=True)
# set overwrite to false if redoing quant to add a group of new elements
Backfitlog, Integquantlog = AESintquant.integbatchquant(spelist, Smdifpeakslog, AESquantparams, Elements, reprocess=False, overwrite=False)

# perform this explicitly on sub spe files
Backfitlogsubs, Integquantlogsubs = AESintquant.integbatchquant(subspelist, Smdifpeakslogsubs, AESquantparams, Elements, reprocess=True, overwrite=True)
# the background fit results and integration results must be manually saved (not-autosaved)
Backfitlog.to_csv('Backfitlog.csv', index=False)
Integquantlog.to_csv('Integquantlog.csv', index=False)

Backfitlogsubs.to_csv('Backfitlog_subs.csv', index=False)
Integquantlogsubs.to_csv('Integquantlog_subs.csv', index=False)

integbatchquant(spelist, Smdifpeakslog, AESquantparams, Elements, reprocess=False, overwrite=True)

# Calculated (or recalculate after AESquantparams change) adjusted counts based on k-factor2 for direct integration method
Integquantlog=AESintquant.calcadjcounts(Integquantlog, AESquantparams, sig=2, kerrors=True) 
#%% CREATE PLOT REPORTS TO VIEW DERIVATIVE AND INTEGRAL METHOD RESULTS
# Reloading of key log files (required after alterations to csv files or restart)
AugerParamLog, spelist, Smdifpeakslog, Integquantlog, Backfitlog, AESquantparams= AESutils.loadmainfiles()
Smdifpeakslogsubs, Backfitlogsubs, Integquantlogsubs=AESutils.loadsubfiles()

AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437')
Smdifpeakslogsubs=pd.read_csv('sub\\Smdifpeakslog_subs.csv', encoding='cp437')
Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
Integquantlogsubs=pd.read_csv('sub\\Integquantlog_subs.csv', encoding='cp437')
Backfitlog=pd.read_csv('Backfitlog.csv', encoding='cp437')
Backfitlogsubs=pd.read_csv('sub\\Backfitlog_subs.csv', encoding='cp437')
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Params\\AESquantparams.csv', encoding='utf-8') # global version
AESquantparams=pd.read_csv('AESquantparams.csv', encoding='utf-8') # load local version instead

# DATA SLICE AND DICE
# Select subset of quality Auger spectra
spelist=AugerParamLog[(AugerParamLog['Areas']>=1)] # selects only spe files
excludemask=spelist['Comments'].str.contains('exclude', case=False, na=False)
spelist=spelist.loc[~excludemask]
spelist=spelist.sort_values(['Filenumber'], ascending=True)
# Select only subfiles before combination via averaging
AugerParamLog=updatepaths(AugerParamLog) # ensures that log properly reflects which csv files are in subdir
AugerParamLog.to_csv('AugerParamLog.csv', index=False) # optional save after path update
submask=spelist['FilePath'].str.contains('sub', case=False, na=False)
subspelist=spelist.loc[submask]
spelist=spelist.loc[~submask] # yields only unique or combined files in base dir (not sub-files)
tempspe=spelist[0:10]
Tirich=Smdifpeakslog[Smdifpeakslog['PeakID']=='Ti']
Tirich=Tirich[Tirich['Amplitude']>5000]

# Choose elements list to plot and label
plotelems=['S', 'C', 'Ca','O', 'Fe', 'Fe2','Mg','Si']
plotelems=['S','C','Ca','O', 'Fe', 'Fe2','Fe1','Mg','Si'] # list of elements/lines to label on plots
plotelems=['C','O', 'Fe', 'Fe2','Fe1', 'Na', 'Mg','Si','Au', 'Au2']
plotelems=['Ti', 'Ti2', 'N', 'Ca', 'Fe2','Fe', 'Mg','Si']
plotelems=['C','Ca','Fe', 'Fe2','Mg','Si','O','Cs','Cs2']
plotelems=['C','O','Fe','Mg','Si']
plotelems=['S','Ca','O','Fe','Fe2','Mg','Si']
plotelems=['S','Ca','Fe','Fe2','Mg','Si','O']
plotelems=['Ti', 'N', 'Ti2','Si', 'Fe', 'Ca','S','Mg']
plotelems=['100-2000'] # can also just enter a range
plotelems=['100-500','O', 'Fe', 'Fe2'] # mixture of ranges and elements is allowed
plotelems=['Fe','Fe2','Mg','Si', 'Au','Au2']
plotelems=['In','O','Mg','Si', 'Pt','Ga']
plotelems=['In','O','Mg','Si', 'Pt','Ga','Fe','Fe2','S','Ca','Al']
plotelems=['Mg','Si','Fe']
plotelems=['C','O', 'Si']
plotelems=['30-2130']

testlist=spelist[0:1] # select subset from above list of spectra
# General smooth-differentiated plot report for selected elements

AESplot.reportSD(spelist, Smdifpeakslog, plotelems, AESquantparams, PDFname='SDreport_2Feb17.pdf')
AESplot.reportSD(Mgrich, Smdifpeakslog, plotelems, AESquantparams, PDFname='Mgrich_areas.pdf')
AESplot.reportSD(subspelist, Smdifpeakslogsubs, plotelems, AESquantparams, PDFname='SDreportsubs_3Jan17.pdf')
AESplot.reportSD(tempspe, Smdifpeakslog, plotelems, AESquantparams, PDFname='SDreportsubs_test.pdf')

# plot report for direct peaks (option to include direct background fits from integral quant method)
AESplot.reportcountsback(spelist, plotelems, AESquantparams, backfitdf=Backfitlog, PDFname='Countsback_report_30Dec16.pdf')
AESplot.reportcountsback(tempspe, plotelems, AESquantparams, backfitdf=Backfitlogsubs, PDFname='Countsback_subs_test.pdf')
AESplot.reportcountsback(spelist, plotelems, AESquantparams, backfitdf=False, PDFname='countsback_report7Dec16.pdf') # just counts without backgrounds

AESplot.reportderivcntall(Mgrich, plotelems, AESquantparams, Smdifdf=Smdifpeakslog, backfitdf=False, PDFname='Derivcnts_MgrichQM188207_12Dec16.pdf')
AESplot.reportderivcnt(Ferich, plotelems, AESquantparams, Smdifdf=Smdifpeakslog, backfitdf=False, PDFname='Derivcnts_FerichQM188207_12Dec16.pdf')

AESplot.reportpeaksall(spelist, plotelems, AESquantparams, PDFname='Peaks_report.pdf')
# If erroneous spectra are found, add phrase exclude to comments in AugerParamLog for this filenumber... reload and then can be excluded from future processing
# For refitting problematic subset and partial redo of those quant results ... see Auger_quant_alt.py

#%% COMPUTING ELEMENTAL COMPOSITIONS based on derivative quant method or integral quant method 

# Calculate (or  recalculate after AESquantparams change the adjusted amplitude for deriv method
Smdifpeakslog=AESsmquant.calcamplitude(Smdifpeakslog, AESquantparams)
Smdifpeakslogsubs=AESsmquant.calcamplitude(Smdifpeakslogsubs, AESquantparams)
Smdifpeakslog.to_csv('Smdifpeakslog.csv', index=False) # save above result (not auto-saved)
Smdifpeakslogsubs.to_csv('sub\\Smdifpeakslog_subs.csv', index=False) 
# Now calculate composition (at %) for selected elements
Elements=['S','Ca','Fe', 'Mg','Si','Fe2']
Elements=['S','Ca','Fe', 'Mg','Si']
Elements=['S', 'Fe','Mg','Si']
Elements=['S','Fe','Mg','Si']
Elements=['S','Mg','Fe','Ca','Si','C','O','Fe2']
Elements=['Au']
Elements=['S','Mg','Fe','Ca','Si','Ti', 'N'] # for refractory mix craters
Elements=['Mg','Si','Fe1', 'O']
Elements=['S','Mg','Ca','Si','Ti']
Elements=['S','Mg','Ca','Si','Ti','Fe']
Elements=['S','Ti','Ca','Si','Fe2']
Elements=['Mg','Fe','Au','Si','Fe2', 'Au2']
Elements=['Mg','Fe','Au','Si','O']
Elements=['Mg','Cs','Cs2','Si','O','C']
Elements=['In','Mg','Si', 'Pt','Ga']
Elements=['In','Mg','Si']
# Calculate composition based on derivative method (smooth-diff comp)
# Non-zero threshold excludes weak lines from quant results (threshold=1 means excluded if peak amplitude < noise amplitude
Smdifcomp=AESsmquant.calccomposition(spelist, Smdifpeakslog, Elements, threshold=0.0)
Smdifcompsubs=AESsmquant.calccomposition(subspelist, Smdifpeakslogsubs, Elements, threshold=0.0)


Smdifcompsub=Smdifcompsubs[Smdifcompsubs['AESbasis']>0]
# saves compositions and underlying elements/k-factors to auto-named xls file
AESutils.writecomps(Smdifcomp, AESquantparams, Elements)

Smdifcomp=pd.read_excel('Smdiffcomp_03Jan17.xlsx', sheetname='smdiffcomp')
Smdifcompsubs=pd.read_excel('sub\\smdiffcomp_04Jan17.xlsx', sheetname='smdiffcomp')

Integcomp=pd.read_excel('integcomp_02Jan17.xlsx', sheetname='integcomp')

Smdifcompsubs=pd.read_csv('Smdifcompsubs.csv')
Mgoutlierquant=calccomposition(Mgoutliers, Smdifpeakslog, Elements, threshold=0.0)
Mgoutlierquant.to_csv('Mgoutlier_quant.csv',index=False)

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

# Set keyword arguments (dict) for plotting functions (some are optional others have built-in defaults)
kwargs={}
kwargs.update({'basis':False}) # can plot Auger basis (before at% conversion) or atomic percent
kwargs.update({'errbars':'xy'}) # Optional error bars on x, y or xy
# Join options for comparison of comp1 and comp2 
# For comparing smdiff and integ comps on exact same file (default value)
kwargs.update({'joinlist':['Filenumber','Areanumber']}) 
kwargs.update({'joinlist':['Sample']}) # Compare all spectra from same sample (more points); 
# Sample join can be erroneous (i.e. if some areanumbers in spe file are background not sample)
kwargs.update({'thresh':0.1}) # threshold for defining/returning point as outlier (higher gives more outliers)

# Directly compare deriv and integ on exact same data files with multi-element scatter plots
compdata, outliers=AESplot.scattercompplot(Smdifcomp,Integcomp, Elements, joinlist=['Filenumber','Areanumber'], thresh=0.1, basis=True, errbars='y')
compdata, outliers=AESplot.scattercompplot(Smdifcomp,Integcomp, Elements, joinlist=['Filename','Areanumber'], thresh=0.05, basis=True, errbars='y')
compdata, outliers=AESplot.scattercompplot(Smdifcomp, Smdifcomp2, Elements, joinlist=['Filename','Areanumber'], thresh=0.05, basis=True, errbars='y')

compdata, outliers=scattercompplot(Smdifcomp,Integcomp, Elements, joinlist=['Filenumber','Areanumber'], thresh=0.1, basis=True, errbars='y')
# basis=true plots elemental basis (strength of line) whereas false plots at.%
# Joinlist -- can also merge with ['Sample', 'Areanumber' 
# thresh --- larger number returns more outliers 
# errbars -- includes errorbars (which are generally only available for integral method )
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

# Plotting counts/backfit and derivatives for outliers (from Filenumber and Areanumber, not loop through areas)
plotelems=['S','Ca','Fe', 'Mg','Si','Ti','Ti2']
plotelems=['S','Ca','Fe', 'Fe2', 'Mg','Si']
plotelems=['S','Fe', 'Ca', 'Fe2','Fe1','Mg','Si']
plotelems=['S','Ti','Ti2','Ca', 'Fe2','Fe','Si']
plotelems=['S','Ca','Fe', 'Fe2', 'Mg','Si']
plotelems=['30-1730'] # can also just pass a single ev range instead of discrete element ranges

# plot/report for all in spe list (optional pass and plot of logs with smdif and background determination points)
AESplot.reportderivcntall(spelist, plotelems, AESquantparams, Smdifdf=Smdifpeakslog, backfitdf=Backfitlog, PDFname='cnts_deriv_report.pdf')
AESplot.reportderivcntall(spelist, plotelems, AESquantparams, Smdifdf=False, backfitdf=False, PDFname='cnts_deriv_report.pdf')

# Plot/report for selected filenum/areanums and selected elements
AESplot.reportderivcnt(temp,  plotelems, AESquantparams, Smdifdf=Smdifpeakslog, backfitdf=Backfitlog, PDFname='Mgoutlier_plots2.pdf')
AESplot.reportderivcnt(Sioutliers,  plotelems, AESquantparams, Smdifdf=Smdifpeakslog, backfitdf=False, PDFname='Srich_areas.pdf')

# make a standard ternary plot
ternelems=['Fe', 'Si', 'Mg+Ca']
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


