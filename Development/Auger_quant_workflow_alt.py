# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 18:08:48 2017

@author: tkc
"""
import pandas as pd

# load master AugerParamLog for this data directory (created during batch import)
AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437') # use Excel compatible encoding cp437 instead of utf-8 
AugerParamLogsubs=pd.read_csv('sub\\Augerparamlog_subs.csv', encoding='cp437')
AugerParamLog.to_csv('Augerparamlog.csv',index=False)
# reload smooth-diff peak positions (if not already open)
Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437') 
Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
Backfitlog=pd.read_csv('Backfitlog.csv', encoding='cp437')

Smdifpeakslogsubs=pd.read_csv('sub\\Smdifpeakslog_subs.csv', encoding='cp437')
Integquantlogsubs=pd.read_csv('sub\\Integquantlog_subs.csv', encoding='cp437')
Backfitlogsubs=pd.read_csv('sub\\Backfitlog_subs.csv', encoding='cp437')

# load peak positions, kfactors, integration params, etc.
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Params\\AESquantparams.csv', encoding='utf-8')
AESquantparams=pd.read_csv('AESquantparams.csv', encoding='utf-8') # load local version instead


# Various filtering operations
AugerParamLog=AugerParamLog.drop_duplicates(['Filenumber']) # drop accidental filenumber duplicates

subspelist=AugerParamLogsubs[(AugerParamLogsubs['Areas']>=1)]
subspelist=subspelist.loc[~subspelist['Comments'].str.contains('exclude', case=False, na=False)]

# remove files excluded by log
spelist=spelist.loc[~spelist['Comments'].str.contains('exclude', case=False, na=False)]
spelist=spelist.loc[[66]]
tempspe=subspelist[0:51]

excllist=np.ndarray.tolist(Excluded.Filenumber.unique())
spelist=AugerParamLog.loc[~AugerParamLog.Filenumber.isin(excllist)]

# Exclude from spelist any spe files with "exclude" in comments
spelist=spelist.loc[~AugerParamLog['Comments'].str.contains('exclude', case=False, na=False)]
includemask=spelist['Comments'].str.contains('avg', case=False, na=False)
spelist=spelist.loc[includemask]
includemask=spelist['Comments'].str.contains('quantmap', case=False, na=False) # quantmap is combined file, subqm is sub array files
# select on sub spe files (combined via avg into other files)
subspelist=spelist.loc[spelist['FilePath'].str.contains('sub', case=False, na=False)]

spelist=spelist.loc[~spelist['FilePath'].str.contains('sub', case=False, na=False)]
spelist=spelist.sort_values(['Filenumber'], ascending=True)

AugerParamLog=AESutils.updatepaths(AugerParamLog) # update paths 

# Check for missing data files
filelist=glob.glob('*.csv')
AESutils.checklog(filelist, AugerParamLog) 

# Manual ways of picking file subsets
tempspe=spelist[0:21]  # different slices from Auger parameters log
combinelist=combinelist.loc[320:330] 

# LOAD AESQUANTPARAMS (peak positions, kfactors etc.)
AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\AESquantparams.csv', encoding='utf-8')
# Can also load a local version (if modifications are needed for this sample's quant
AESquantparams=pd.read_csv('AESquantparams.csv', encoding='utf-8')

# Manually choosing sets of elements

Elements=['S','C','Ca','Ti','O','Fe1','Fe2','Fe','Na','Mg','Al','Si'] # using dominant Auger peak for each element... see AESquantparams.py for associated energies
Elements=['S','C','Ca','Ti','Ti2','O','Fe1','Fe2','Fe','Na','Mg','Al','Si','In','F'] 
Elements=['S','C','Ca','Ti','O','Fe1','Fe2','Fe','Ti2','Mg','Al','Si','Cs','Cs2']
Elements=['Mg']
Elements=['C','O','Si', 'W']
Elements=['Pt','Ga','Si', 'Mg','In','O']
Elements=['Si', 'Mg','In','Fe']
Elements=['S','Ca','Fe2','Fe','Mg','Si'] # choose major elements
Elements=['Fe2','Fe','Mg','Si','Au']
Elements=['S','Ca','Fe2','Fe','Mg','Si','Ti','Ti2']
Elements=['S']
Elements=['S','Mg','Fe','Ca','Si','Ti', 'N'] # for refractory mix craters

plotelems=['S', 'C', 'Ca','O', 'Fe', 'Fe2','Mg','Si']
plotelems=['S','C','Ca','O', 'Fe', 'Fe2','Fe1','Mg','Si'] # list of elements/lines to label on plots
plotelems=['C','O', 'Fe', 'Fe2','Fe1', 'Na', 'Mg','Si','Au', 'Au2']
plotelems=['C','Ca','Fe', 'Fe2','Mg','Si','O','Cs','Cs2']
plotelems=['C','O','Fe','Mg','Si']
plotelems=['S','Ca','Fe','Fe2','Mg','Si','O']
plotelems=['Ti', 'N', 'Ti2','Si', 'Fe', 'Ca','S','Mg'] 
plotelems=['100-500','O', 'Fe', 'Fe2'] # mixture of ranges and elements is allowed
plotelems=['Fe','Fe2','Mg','Si', 'Au','Au2']
plotelems=['In','O','Mg','Si', 'Pt','Ga','Fe','Fe2','S','Ca','Al']
plotelems=['C','O', 'Si']

# Non-gui interactive plot (old way)
kwargs={}
kwargs.update({'plotelems':Elements})
kwargs.update({'xrange':[100, 800]}) # auger range in ev
kwargs.update({'smdifpeaks':Smdifpeakslog}) # optional plotting of found peak params
kwargs.update({'backfitdf':Backfitlog}) # optional plots of background fits
kwargs.update({'areas':'1-4'}) # areas for inclusion (csv and/or hyphenated range)
kwargs.update({'plotcol':'S7D7'})  # either 'Counts' or 'S7D7' (deriv)
AESplot.AESplot1('101103, 106108', spelist, AESquantparams, **kwargs) # filenumber(s) as string

# older method plot/report for all in spe list (optional pass and plot of logs with smdif and background determination points)
AESplot.reportderivcntall(spelist, plotelems, AESquantparams, Smdifdf=Smdifpeakslog, backfitdf=Backfitlog, PDFname='cnts_deriv_report.pdf')
AESplot.reportderivcntall(spelist, plotelems, AESquantparams, Smdifdf=False, backfitdf=False, PDFname='cnts_deriv_report.pdf')

# Plot/report for selected filenum/areanums and selected elements
AESplot.reportderivcnt(temp,  plotelems, AESquantparams, Smdifdf=Smdifpeakslog, backfitdf=Backfitlog, PDFname='Mgoutlier_plots2.pdf')
AESplot.reportderivcnt(Sioutliers,  plotelems, AESquantparams, Smdifdf=Smdifpeakslog, backfitdf=False, PDFname='Srich_areas.pdf')

kwargs={}
kwargs.update({'PDFname':'Derivcnts_report.pdf'}) # optional custom PDF report name 
kwargs.update({'bf':Backfitlog}) # optional plot of points used for integ background fits
kwargs.update({'smd':Smdifpeakslog}) # optional plot of points in smooth-diff used for amplitude calcs
AESplot.reportderivcnt(Ferich, plotelems, AESquantparams, **kwargs)

# old kwarg construction for scatter comp plots