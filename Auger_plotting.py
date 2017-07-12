# -*- coding: utf-8 -*-
"""
Created on Wed May 11 08:09:38 2016
Smdiff quant visual checking using plots: need to plot 
1) spectra 
2) Exact points from which peak amplitude was determined
3) Nearest background regions and amplitude
- must be able to pass and select subregion for the plot (xrange - either string of numbers or an element)
Sequential Also would be good to create a loop through all peaks 

@author: tkc
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
if 'C:\\Users\\tkc\\Documents\\Python_Scripts' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts')
import Auger_plot_functions as AESplot
from matplotlib.backends.backend_pdf import PdfPages
#%% Import parameters log and peaks log 
# change to data directory of interest (use console or )
# set to interactive plots
AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp1252') # if altered by Excel (otherwise utf-8)
Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437') # load peak params file
Backfitlog=pd.read_csv('Backfitlog.csv', encoding='cp437') 
Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437') 
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\AESquantparams.csv', encoding='utf-8')
#AugerParamLog,Smdifpeakslog = truncatenumber(AugerParamLog,Smdifpeakslog) # truncate combined file numbers i.e. from 138140 to 138
# AugerParamLog=AugerParamLog.reset_index(drop=True)

#%% Filter files within AugerParamLog to plot only desired those .spe files 

spelist=AugerParamLog[(AugerParamLog['Areas']>=1)] # selects only spe files

excludemask=spelist['Comments'].str.contains('exclude', case=False, na=False)
spelist=spelist.loc[~excludemask]
spelist=spelist.sort_values(['Filenumber'], ascending=True)

spelist=spelist[0:10]  # different slices from Auger parameters log
spelist=spelist.loc[320:330] 

spelist=spelist.reset_index(drop=True)

#%% Batch PDF report of major element plots (either smoothdiff or raw counts) for each spectrum in 

plotelems=['S','C','Ca','O', 'Fe', 'Fe2','Fe1','Mg','Si'] # list of elements/lines to label on plots
plotelems=['C','O', 'Fe', 'Fe2','Fe1', 'Na', 'Mg','Si','Au', 'Au2']
plotelems=['S', 'O', 'C', 'Ca', 'Fe', 'Fe1', 'Fe2','Mg','Si']
plotelems=['C','Ca','O', 'Fe', 'Fe2','Fe1','Mg','Al','Si']
plotelems=['C','Ca','O', 'Fe', 'Fe2', 'F','In','Mg','Al','Si']
plotelems=['S','Ca','Fe', 'Fe2','Mg','Si']

# General smooth-differentiated plot report for selected elements
AESplot.reportSD(spelist, Smdifpeakslog, plotelems, AESquantparams, PDFname='SDreport.pdf')

# plot report for direct peaks (option to include direct background fits from integral quant method)
AESplot.reportcountsback(spelist, plotelems, AESquantparams, backfitdf=False, PDFname='countsback_report.pdf')
AESplot.reportcountsback(spelist, plotelems, AESquantparams, backfitdf=Backfitlog, PDFname='countsback_report.pdf')
reportcountsback(spelist, plotelems, AESquantparams, backfitdf=Backfitlog, PDFname='countsback_report.pdf')

reportpeaks(spelist, plotelems, AESquantparams, PDFname='Peaks_report.pdf')
#%% 
Smdifpeakslog=copyproblemcomments(Integquantlog, Smdifpeakslog)
Smdifpeakslog.to_csv('smdifpeakslog.csv',index=False)
Integquantlog.to_csv('integquantlog.csv',index=False)
#%% Comparison of compositions calculated in different ways
Elements=['S','Ca','Fe', 'Mg','Si']
datacomp=AESplot.scattercompplot(SmdifcompFe,SmdifcompavgFe, Elements, joinlist=['Sample','Areanumber'], basis=False)

# For internal comparison (ie. Fe from Fe alone vs Fe+Fe2 average)
outliers=scattercompplot(comp1,comp2, Elements, joinlist=['Filenumber','Areanumber'], basis=False)
outliers.to_csv('outliers.csv')

plotthese=outliers[outliers['Element']=='%Ca']
AESplot.reportSD(plotthese, Smdifpeakslog, plotelems, AESquantparams, PDFname='SDoutliers.pdf')
#%% Other plotting functions
AESplot.reportcountsmajor(spelist, Smdifpeakslog,PDFname='countsplot_report.pdf') # legacy funct (replaced by above)

AESplot.reportcountsback(spelist,Linearfitlog) # report with major peaks from all spectra incl background region fits
AESplot.reportpeaks(spelist, addgauss=True) # report with major peaks from all spectra incl background region fits

# PDF report for all available files (can be pretty slow)
AESplot.reportSDmajor(spelist, Smdifpeakslog, PDFname='SDplots_report.pdf') # Pass all files and produce PDF report
AESplot.reportcountsmajor(spelist, Smdifpeakslog) # Pass all files and produce PDF report

#%% Ways to select single data file of interest
# Slice by filenumber
filenumber=282284
Params=AugerParamLog[(AugerParamLog['Filenumber']==filenumber)]
Params=Params.squeeze() # converts df to series
Peaks=Smdifpeakslog[(Smdifpeakslog['Filenumber']==filenumber)]

# Slice by row # in logfile
rownumber=1
Params=AugerParamLog.iloc[rownumber] # this file's params as Series
filenumber=Params.Filenumber # retrieve filenumber
Peaks=Smdifpeakslog[(Smdifpeakslog['Filenumber']==filenumber)] # retrieve assoc. subset of peaks data
#%% # work on ways to select multiple data files for simultaneous plotting (under development)

filelist='138, 143, 198a2' # select list of files by number (area 1 if not specified)
# construct slices from main dataset based on above filestring... fileareas is tuple with filenum and associated area
Params, Peaks, fileareas=fileslice(filelist, AugerParamLog,Smdifpeakslog) # slice based on above
# TODO pass fileslice to version of plotmajor 

#%% 
#%% Set plot range by element, by evrange or all (max range); choose Auger peaks to add
plotrange='Fe'
plotelems='Fe'

plotrange='100-500'
plotrange='all'
plotrange='Fe'

areas='all' # select which spatial areas to plot
areas='1-4,8'
plotelems='C O Mg S Fe Si Al N' # these ideal peak positions will appear as vert lines on plots
#%%
# Plot single selected file by filenumber (defaults to area 1) .. choose plot range and elem lines to be displayed
filenumber=282284
plotrange='all' # set this to desired range 
plotelems='C O Mg S Fe Si Al N' 
Params=AugerParamLog[(AugerParamLog['Filenumber']==filenumber)]
Params=Params.squeeze() # flatten to Series
Peaks=Smdifpeakslog[(Smdifpeakslog['Filenumber']==filenumber)]
AESplot.AESplot1(Params,Peaks, plotrange, plotelems=plotelems)

#%%
# Plot current Auger filenumber -- counts and S7D7 stack (all areas) 
AESplotstack(Params,Peaks, plotrange, plotelems=plotelems)

# Single plot of direct counts and background of Auger filenumber incl. range for S, C/Ca, O, Fe, Mg, Si
plotcntsmajor(Params, areanum=1)
AESplot.plotcntsmajor(Params, areanum=1)

reportcountsback(AugerParamLog)  # pdf with major peaks direct counts and fitted background

# default call uses above defined areas, ranges, but if string is null function uses defaults

# plotting of dataframe counts along with numpy ndarrays (derivs of above data)

#%% MAJOR ELEMENTS PLOT -- Manual loop through stack of files ... 6-plot of all major peaks (S, C/Ca, O, Fe, Mg, Si)
rownumber=0  #screwed up rows: 34
areanum=1

plt.close('all') # close all open plot windows
areanum=areanum+1
rownumber=rownumber+1

Params=AugerParamLog.iloc[rownumber] # this file's params as Series
filenumber=Params.Filenumber # retrieve filenumber
print('Numareas in file ', filenumber,' from row ', rownumber,' is ', Params.Areas)
Peaks=Smdifpeakslog[(Smdifpeakslog['Filenumber']==filenumber)] # retrieve assoc. subset of peaks data
Peaks=Peaks[(Peaks['Areanum']==areanum)] # then only this areanum
AESplot.plotSDmajor(Params, Peaks, areanum)
# temp troubleshooting version
plotSDmajor(Params, Peaks, areanum) 