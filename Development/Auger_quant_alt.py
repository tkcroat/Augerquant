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
if 'C:\\Users\\tkc\\Documents\\Python_Scripts' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts')

import Auger_smdifquant_functions as AESsmquant
import Auger_integquant_functions as AESintquant
import Auger_utility_functions as AESutils
import Auger_plot_functions as AESplot
''' AESsmquant contains functions related to peak finding in smooth-differentiated spectra
whereas AESquant contains background fitting and integration over peaks for direct from counts '''
# import Auger_integquant_functions as AESquant

#%% REFIT of problematic peaks
# Manual refitting of failed fits on single peaks (usually Ca)
# filter with SPE list above 

AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437')
Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
Backfitlog=pd.read_csv('Backfitlog.csv', encoding='cp437') 
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\AESquantparams.csv', encoding='utf-8') # global version
# ALTERNATE QUANT and FIT PARAMS (which are sometimes used if problems arise)
AESquantparams=pd.read_csv('AESquantparams.csv', encoding='utf-8') # load local version instead
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\AESquantparams_Ca_refit.csv', encoding='utf-8')
# Change type of Ca fit to 'Carefit'... linear 

# Pick out the spectra for processing
spelist=AugerParamLog[(AugerParamLog['Areas']>=1)] # selects only spe files
excludemask=spelist['Comments'].str.contains('exclude', case=False, na=False)
spelist=spelist.loc[~excludemask]

Elements=['Ca'] 
# Create smoothed column for all csvs who lack them (now added during import process)
AESutils.addsmoothloop(spelist) # autosaves smcounts column for each area (7pt adjacent averaging)

# Comparing Ti and Ti2 peak magnitudes
Tidata=Integquantlog[Integquantlog['Element'].str.contains('Ti')]
Tidata=Integquantlog[Integquantlog['Element']=='Ti']
Ti2data=Integquantlog[Integquantlog['Element']=='Ti2']
Tidata=Tidata.dropna(subset=['Erradjcnts'])
Ti2data=Ti2data.dropna(subset=['Erradjcnts'])
Ticomp=pd.merge(Tidata, Ti2data, how='inner',on=['Filenumber','Area'], suffixes=('','_2'))
Ticomp.to_csv('Tidata.csv',index=False)

# Are peak shifts correlated in strong spectra?  Stdev in shift col after strength filter
Fe=Smdifpeakslog[Smdifpeakslog['PeakID']=='Fe']
Fe['Amplitude'].hist(bins=20)
Fe['Amplitude'].describe()
Smdifpeakslog.columns

# For peak shift determinations, use smdiff negpeak location or direct peak
Shiftlog.plot.scatter(x='Amplitude',y='Avgshift')

# Older/ alt plotting schemes (now replaced by single plot/report interface)
# General smooth-differentiated plot report for selected elements
AESplot.reportSD(spelist, Smdifpeakslog, plotelems, AESquantparams) # default name is SDplots_report_09Jun17.pdf
AESplot.reportSD(Mgrich, Smdifpeakslog, plotelems, AESquantparams, PDFname='Mgrich_areas.pdf')
AESplot.reportSD(subspelist, Smdifpeakslogsubs, plotelems, AESquantparams, PDFname='SDreportsubs_3Jan17.pdf')

# Create a counts/background plot report over selected files 
kwargs={} # kwargs fed back to set defaults in case of rerun
kwargs=AESplot.countsbackreport_tk(spelist, Elements, Backfitlog, AESquantparams, **kwargs) 
kwargs=countsbackreport_tk(spelist, Elements, Backfitlog, AESquantparams, **kwargs) 
kwargs=AESplot.countsbackreport_tk(subspelist, Elements, Backfitlogsubs, AESquantparams, **kwargs) 

# Plot report of derivative and counts on same page (top and bottom)
kwargs={}
kwargs=AESplot.countsderivreport_tk(myfiles, Elements, Smdifpeakslog, Backfitlog, AESquantparams, **kwargs)
kwargs=countsderivreport_tk(myfiles, Elements, Smdifpeakslog, Backfitlog, AESquantparams, **kwargs)

# Plot report of counts - background for selected peaks (all areas )
myfiles=AESutils.pickspectraGUI(spelist)
AESplot.reportpeaksall(spelist, plotelems, AESquantparams, PDFname='Peaks_report.pdf')
