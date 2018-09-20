# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 10:14:59 2018

@author: tkc
"""

#%% Load modules
import os, glob,sys # already run with functions 
import pandas as pd
import numpy as np
#import re, struct (only used in sub-functions)
# import csv, fileinput
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules')
import Auger_batch_import_functions as Auger
import Auger_utility_functions as AESutils
import Auger_quantmap_functions as QM
import Auger_plot_functions as AESplot
import Auger_smdifquant_functions as AESsmquant
import Auger_integquant_functions as AESintquant

#%% Set data directory with Auger data files (via cd in Ipython console, using tinker or using os.chdir)
os.chdir('C:\\Temp\\AugerQM\\28Sep17')
# load filenames of interest
filelist=glob.glob('*.sem')+glob.glob('*.map')+glob.glob('*.spe')
filelist=glob.glob('*.map')
filelist=glob.glob('*.sem')
filelist=glob.glob('*.spe')

# Find and read Auger logbook (xls file that contains phrase "Auger_logbook")
# this file is to associate sample/project names with Auger file numbers and can be used to combine/average multiple spe files
Augerlogbook=Auger.openorcreatelogbook(filelist)

kwargs={}
kwargs={'move':False} # option to not move files to /sub
# Extract params and process SEM, SPE and MAP files
AugerParamLog = Auger.Augerbatchimport(filelist, Augerlogbook, **kwargs) 

# log files are overwritten if they exist by default (but data binaries are not)
AugerParamLog.to_csv('Augerparamlog.csv',index=False) # save all params to new more complete log file (not autosaved)

# Create jpg images annotated with spatial areas for spe files (assumes existence of .sem image taken just before .spe file)
# easiest to do this before moving combined files (into separate log)
SpatialAreas=pd.read_csv('spatialareaslog.csv') # open automatically created spatial areas log for spe files
Auger.makeannotatedjpg(AugerParamLog, SpatialAreas) # makes annotated jpgs for all spe files with prior sem image

# Reloading of key log files (required after alterations to csv files or restart)
AugerParamLog, spelist, Smdifpeakslog, Integquantlog, Backfitlog, AESquantparams= AESutils.loadmainfiles()

# Element and background region setup
Peaks=AESutils.pickelemsGUI(AESquantparams, Smdifpeakslog, Integquantlog) # interactive element selection

# MAIN BATCH QUANT LOOP function which returns Smdifpeakslog dataframe with amplitude of every peak in above list
Smdifpeakslog=AESsmquant.smdiffquant_gui(spelist, Peaks, AESquantparams, Smdifpeakslog)

# MAIN BATCH INTEGRAL METHOD fitting and direct integration loops (returns background fitting params and integrated quant results)
# WARNING ... overwrite=True will overwrite background fit and peaks columns in source csv file (counts obviously unchanged)
Backfitlog, Integquantlog, Shiftlog = AESintquant.integbatchquant(spelist, Smdifpeakslog, 
                AESquantparams, Peaks, reprocess=True, overwrite=True)
Backfitlog.to_csv('Backfitlog.csv', index=False)
Integquantlog.to_csv('Integquantlog.csv', index=False)
Shiftlog.to_csv('Shiftlog.csv', index=False) # Compares smdiff quant shift and integquant shift

# Usually batch plot report to PDF (many si)
kwargs={}
kwargs=AESplot.AESreport(spelist, Peaks, Smdifpeakslog, Backfitlog, Integquantlog, AESquantparams, **kwargs)
kwargs=AESreport(spelist, Peaks, Smdifpeakslog, Backfitlog, Integquantlog, AESquantparams, **kwargs)

#%% COMPUTING ELEMENTAL COMPOSITIONS based on derivative quant method or integral quant method 
# Choose full set of elements for quant
Peaks=AESutils.pickelemsGUI(AESquantparams, Smdifpeakslog, Integquantlog)
# Reloads elemental peaks used previously
Peaks=np.ndarray.tolist(Smdifpeakslog.PeakID.unique())
# Optional recalc (if different k-factor set is desired)
Smdifpeakslog=AESsmquant.calcamplitude(Smdifpeakslog, AESquantparams)
Smdifpeakslog.to_csv('Smdifpeakslog.csv', index=False)
# Calculate composition based on derivative method (smooth-diff comp) 
# Non-zero threshold excludes weak lines from quant results (threshold=1 means excluded if peak amplitude < noise amplitude
Smdifcomp=AESsmquant.calccomposition(spelist, Smdifpeakslog, Peaks, threshold=0.0)

# Save smdiff-based compositions in xls file (along w/ kfactors) 
AESutils.writecomps(Smdifcomp, AESquantparams, Peaks)
# Abbreviated compositional summary w/ selected element subset
Peaksexcl=['C']
Peaksexcl=['C','O'] # Fe2 already not present

# Compositional summary
Smdifsummary=AESutils.compsummary(Smdifcomp, Peaks, Peaksexcl)
Smdifsummary.to_csv('Smdiff_compsummary.csv', index=False)
Smdifsummary.to_csv('Smdiff_compsummary_metals.csv', index=False)

# CALCULATE COMPOSITIONS from direct integration quant (normally includes at.% error)
# First calculate or recalculated elemental peak adjusted by sensitivity k-factor
Integquantlog=AESintquant.calcadjcounts(Integquantlog, AESquantparams, sig=2, kerrors=True) 
# kerrors true propogates estimated error in k-factor defined in AESquantparams into all errors 
# sig refers to incl 2 sig errors for kfactor
Integquantlog.to_csv('Integquantlog.csv', index=False)

# Calculate compositions in at. % based on integral method
Integcomp=AESintquant.calccomp(spelist, Integquantlog, Peaks, AESquantparams)

# Write integration-based compositions to xls (along with k-factor details)
AESutils.writecomps(Integcomp, AESquantparams, Peaks)

''' Scatter compositional comparison plots
Directly compare deriv and integ on exact same data files with multi-element scatter plots
'''
compdata, outliers=AESplot.scattercompplot_tk(Smdifcomp,Integcomp, Peaks)
compdata, outliers=scattercompplot_tk(Smdifcomp,Integcomp, Peaks)

# Interactive quantplotter 
# TODO finish me
launch_plotter(os.getcwd())
