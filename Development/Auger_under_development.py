# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 10:16:40 2016

@author: tkc
"""
#%% Load modules
import matplotlib.pyplot as plt
import re, os, glob, sys, csv # already run with functions 
import pandas as pd
import numpy as np
from pylab import *
import scipy
import scipy.stats
from matplotlib.backends.backend_pdf import PdfPages # Only needed for plotting
from collections import defaultdict
from math import factorial # used by Savgol matrix
from PIL import Image, ImageDraw, ImageFont # needed for jpg creation
import datetime

import statsmodels.formula.api as smf # has OLS module with convenient outliers tests
#%%
                
def packagequant():
    ''' Write backfit, integ, smdiff, AESquantparams to single Excel file '''
    
# Use statsmodels package for outlier detection
comp1=pd.read_csv('Smdifcomp.csv', encoding='cp437')
comp2=pd.read_csv('SmdifcompavgFe.csv', encoding='cp437')

# Find a way to do linear fits with least squares but with fixed intercept (through origin)
# web suggestion sez linalg.lstsq but that didn't work initially
# slope,intercept=np.polyfit(data1, data2, 1)  numpy version
A=np.vstack([xdata, np.ones(len(xdata))]).T  # x needs to be a column vector fo linalb.lstsqNot sure why these arrays need stacking    
slope, _, _, _ = np.linalg.lstsq(A, ydata) # fixes intercept to zero? 


fig, axes = plt.subplots(nrows=1, ncols=1)
Augerfile.plot(x='Energy',y='Backfit1', color='r',ax=axes)
Augerfile.plot(x='Energy',y='Counts1', ax=axes)

print(slope1, errslope1, slope2, errslope2)

diff=comparelinfits(lowfitparams,upperfitparams, thresh=0.4)

# TESTING of linear peak fits
# params are Augerfile, lowfitparams, upperfitparams, fitbounds (4 points)

# setting up data range for fits 
ptlist = [i for i in range(fitbounds[0],fitbounds[1])]
ptlist2 = [i for i in range(fitbounds[2],fitbounds[3])]
fullpts=ptlist+ptlist2

# slice the dataframe
Augerslice=Augerfile[Augerfile.index.isin(ptlist)]
Augerslice=Augerfile[Augerfile.index.isin(fullpts)]
# plot dataframe slice
fig, axes = plt.subplots(nrows=1, ncols=1)
Augerslice.plot.scatter(x='Energy',y='Counts1', ax=axes)
Augerslice.plot.scatter(x='Energy',y='Backfit1', ax=axes)

# plot dataframe slice
fig, axes = plt.subplots(nrows=1, ncols=1)
Compresult.plot.scatter(x='Sampl',y='Samplavg', ax=axes)

# Plot lower line
xdata=Augerslice.Energy
xdata=xdata.loc[fitbounds[0]:fitbounds[1]+1]
ydata=lowfitparams[0]*xdata+lowfitparams[2]
plt.plot(xdata, ydata, color='b')
# Plot w/ 2 stdev different slope 
ymax=(lowfitparams[0]+lowfitparams[1])*xdata+lowfitparams[2]
plt.plot(xdata, ymax, color='r')
ymin=(lowfitparams[0]-lowfitparams[1])*xdata+lowfitparams[2]
plt.plot(xdata, ymin, color='r')

# Test fits with linregress
xdata=Augerslice.Energy
ydata=Augerslice.Counts1

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xdata, ydata) 

# NUMPY LINEAR FITS
z, residuals, rank, singular_values, rcond = np.polyfit(xdata, ydata, 1, full=True)
(slope,intercept), cov = np.polyfit(xdata, ydata, 1, cov=True)
sd_slope=np.sqrt(cov[0,0]) # variance of slope is element 0,0 in covariance matrix
sd_intercept=np.sqrt(cov[1,1]) # same for standard deviation of intercept (sqrt of variance)


# Plot upper line
xdata=Augerslice.Energy
xdata2=xdata.loc[fitbounds[0]:fitbounds[1]+1]
ydata2=lowfitparams[0]*xdata+lowfitparams[2]
plt.plot(xdata, ydata2, color='b', ax=axes)
ymax2=(lowfitparams[0]+lowfitparams[1])*xdata+lowfitparams[2]
plt.plot(xdata, ymax, color='r', ax=axes)
ymin2=(lowfitparams[0]+lowfitparams[1])*xdata+lowfitparams[2]
plt.plot(xdata, ymin2, color='r', ax=axes)
# saving slice of data (pass to origin)
Augerslice.to_csv('testfits.csv',index=False)


smooth_width=7
x1 = np.linspace(-3,3,smooth_width)
norm = np.sum(np.exp(-x1**2)) * (x1[1]-x1[0])  # ad hoc normalization
y1 = (4*x1**2 - 2) * np.exp(-x1**2) / smooth_width *8 #norm*(x1[1]-x1[0])

xdata=lowerfitpeak['Energy']
ydata=lowerfitpeak['Counts1'] # Counts1, Counts2 or whatever

regression= smf.ols("data ~ x", data=dict(data=ydata, x=xdata)).fit()
outliers=regression.outlier_test()
intercept,slope=regression.params
colnames=['Resid','Pval','Bonf']
outliers.columns=colnames # rename columns

xvals=outliers.index.values
outliers.plot.scatter(x=xvals, y='Resid')   # same indices as input data
plt.scatter(outliers.index, outliers.Resid)

def findcombofiles(spelist):
    '''Find repeated spe files to combine via averaging.. 
    also ensure x,y,z stage match  ''' 
    # get list of basenames that have >1 associated spe
    spelist['Basename']='' # add basename column
    for index, row in spelist.iterrows():
        fname=spelist.loc[index]['Filename']
        basename=fname.split('.')[0]
        spelist=spelist.set_value(index,'Basename',basename)
    basenamelist=spelist.Basename.unique()
    basenamelist=np.ndarray.tolist(basenamelist) # list of unique base names
    combinelist=[]
    for i, bname in enumerate(basenamelist):
        match=spelist[spelist['Basename']==bname]
        if len(match)>1: # only combinable if multiple matching spes with same basename
            # ensure they're at same x, y,z positions 
            xvals=match.X.unique()
            if len(xvals)!=1:
                print(len(match),' spe files with basename ', bname,' but ', len(xvals), ' unique X stage positions!... not combined')
                continue
            firstfile=match.Filenumber.min() # assumes consecutive file numbering w/ Autotool
            lastfile=match.Filenumber.max()
            combinelist.append([bname, firstfile, lastfile])
    return combinelist # list with basename, firstfile, lastfile


from operator import itemgetter
from itertools import groupby

              
def autocombinespe(AugerParamLog, movefiles=False):
    ''' automatic method of finding combinable measurements (same X,Y, areas, basename) and combine all areas separately
    via averaging 
    naming convention of basename.1.csv, basename.2.csv, basename.3.csv etc to basename.13.csv
                if str(AugerParamLog.loc[index]['Comments'])=='nan':
                AugerParamLog=AugerParamLog.set_value(index,'Comments','')
    '''
    if 'Basename' not in AugerParamLog: # in case basename not already added
        AugerParamLog['Basename']=''
        for index, row in AugerParamLog.iterrows():
            fname=AugerParamLog.loc[index]['Filename']
            basename=fname.split('.')[0]
            AugerParamLog=AugerParamLog.set_value(index,'Basename',basename)
    for index, row in AugerParamLog.iterrows(): # avoids nan problems with avg string filter 
        if str(AugerParamLog.loc[index]['Comments'])=='nan':
            AugerParamLog=AugerParamLog.set_value(index,'Comments','')
    # generate list of combinable files (same basename)
    basenamelist=AugerParamLog.Basename.unique()
    basenamelist=np.ndarray.tolist(basenamelist) # list of unique base names
    for i, bname in enumerate(basenamelist): # loop through each combinable basename
        logmatches=AugerParamLog[AugerParamLog['Basename']==bname]
        excludemask=logmatches['Comments'].str.contains('avg', case=False, na=False)
        logmatches=logmatches.loc[~excludemask] # fails if all are nan  
        # Exclude files already produced by averaging (avg in comments)
        
        if len(logmatches)>1: # only combinable if multiple matching spes with same basename
            xvals=logmatches.X.unique() # ensure they're all at same x and y positions 
            xvals=np.ndarray.tolist(xvals) # most likely len=1[ but possibly different
            for j, xval in enumerate(xvals):
                xmatches=logmatches[logmatches['X']==xval]
                if len(xmatches)==1: # single file and thus not combinable
                    continue
                # ensure all have same number of areas (true for combinable spectra made with Autotool loops)
                xmatches=xmatches[(xmatches['Areas']>=1)] # convenient way to drop .sem and .map 
                if len(xmatches.Areas.unique())!=1:
                    print('Different # of areas for spectra', bname, str(xmatches.Filenumber.min()),' to ', str(xmatches.Filenumber.min()))
                    continue
                # TODO fix this check if filenumbers are consecutive (true for Autotool)
                '''
                filenums=xmatches.Filenumber.unique()
                filenums=np.ndarray.tolist(filenums)
                for key, group in groupby(enumerate(filenums), lambda x: x[0]-x[1]):
                    print(group)
                    mygroup = map(itemgetter(1), group)
                    print (map(itemgetter(1),group))
                '''
                # now ready to combine this set of files
                firstnum=xmatches.Filenumber.min() # assumes consecutive file numbering w/ Autotool
                lastnum=xmatches.Filenumber.max()
                csvname=bname+str(firstnum)+str(lastnum)+'.csv' # string for combined file based on first-last filenumber
                if not os.path.isfile(csvname):  # skip combined file creation if it already exists (however log is still regenerated)    
                    avgcombinespe(xmatches,csvname) # combines above files and makes new averaged csv 
                    print('Average-combined Auger spectra ', csvname, ' created.')
                    # make new entry for Augerparamlog  
                    avgentry = avglogentry(xmatches) # create new logbook entry (series) for averaged spectrum (mostly from firstfile's info)
                    # append new Series entry to end of AugerParamLog
                    AugerParamLog=AugerParamLog.append(avgentry, ignore_index=True)
                else:
                    print(csvname,' already exists')
                if movefiles==True: # Moves all csv subfiles (those combined) into /sub directory
                    #TODO change this so that single AugerParamLog at top level remains... subs can be filtered out by "sub' in path
                    AugerParamLog=movespes(filelist, AugerParamLog) # shuffle csv sub files to subfolder and split Paramslog accordingly
    return AugerParamLog
 
logmatches.to_csv('logmatches.csv',index=False) 
    
def integpeaks(Augerfile, Backfitparams, areanum, Elemdata, Shifts, logmatch):
    ''' Background fit for each direct peak, shift is list of energy shifts of negpeak (same order as Eledata (opens source spectrum as Augerfile, 
    fits peak backgrounds above and below using Elemdata, saves background to source csv (overwrites existing fits), also saves linear fit params to logdataframe with position/amplitude/etc;
    desired elements out of data range are skipped (in prior findindices function)
    '''
    #create Smdifpeaks dataframe for temp storage of each peak's params
    Backfitparams=Backfitparams.dropna(subset=['Rval1']) # skip Gaussian fit if background fit fails
        
    AugerFileName=logmatch.Filename # 
    # Create temp df to hold and pass linear fit data
    mycols=['Filenumber', 'Filename', 'Filepath', 'Sample', 'Comments', 'Areanumber', 'Element', 'Integcounts', 
    'Backcounts', 'Significance', 'Xc', 'Width', 'Peakarea', 'Y0','Rsquared','Numchannels']
    Integresults=pd.DataFrame(columns=mycols) # empty df for all integ results for elems in this spe file
    peakname='Peaks'+str(areanum) # this is counts - background (but only calculated in vicinity of known elemental peaks)
    backfitname='Backfit'+str(areanum)
        # all fit regions modify fit region boundaries for this spectrum based on smooth-differentiated peak (2nd deriv, Savgol (poly=2, pts=11))
        # global shifts from smdifpeaks and local shift based on smoothed 2nd derivative 
        # already incorporated into Elemdata values (lower1,2 and upper1,2 fully adjusted)
                
        # loop through and fit all peaks for each element in this spatial area            
    for i, (elem, fittype, integpeak, lower1, lower2, upper1, upper2, kfact, integwidth, siglevel) in enumerate(Elemdata):
        # linear fit below this elem's peak (shifts and adjustments already made)
        Integresult=pd.DataFrame(index=np.arange(0,1),columns=mycols) # blank df row for this element
        if i in Backfitparams.index: # skips integ calc if backfit is n/a (but does put n/a entries in datalog)
            fitregion=Augerfile[lower1:upper2+1]
            if fitregion.empty==True: # skip if no data present (already should be skipped in Elemdata)
                continue
            # addgauss if save of gaussian peak fit in Augerfile is desired
            # Probably could skip Gaussian fitting entirely if peak is weak (check smdiff)
            fitregion, fitparams, rsquared, ier = fitgauss(fitregion, areanum, integwidth, elem, AugerFileName, addgauss=True)
            addgauss=True # maybe pass this arg from elsewhere
            if addgauss==True and ier in [1,2,3,4]: # copy gaussian fit over to csv file if successful
                gaussname="Gauss"+str(areanum)
                if gaussname not in Augerfile.dtypes.index: # add col if not already present                
                    Augerfile[gaussname]='' # add col for gaussian fit
                # Copy gaussian fit to Augerfile... fitregion only modified in new Gauss peak fit column
                Augerfile.loc[fitregion.index,fitregion.columns]=fitregion 
            # if gaussian fit is successful set center integration channel to index nearest xc
            # ier flag of 1,2,3,4 if fit succeeds but rsquared threshold is better
            if rsquared!='n/a': # skip integcounts calc but do put 'n/a' entries in df
                if rsquared>0.4:
                    xc=fitparams[0] # center of gaussian fit
                    center=int(round(xc,0))
                    tempdf=fitregion[fitregion['Energy']==center]
                    try:
                        centerindex=tempdf[peakname].idxmax() # corresponding index # of peak maximum
                    except:
                        print('Gaussian fit center out of data range for ', elem, ' in ', AugerFileName)                
                        centerindex=integpeak # backup method of finding center of integration region
                else: # indication of poor Gaussian fit R2<0.4 (use prior knowledge of peak position)
                    print('Failed gaussian fit for ', elem, ' in ', AugerFileName)
                    # set center integration channel to value passed by integpeak 
                    # this is ideal energy value but adjusted by shift found using smooth-diff quant method
                    centerindex=integpeak # already stores index number of central peak (ideal - sm-diff shift value)
                    # Still do the counts integration for poor gaussian fits 
                # perform integration over peak center channel + integwidth on either side 
                Augerpeak=Augerfile[centerindex-integwidth:centerindex+integwidth+1]
                integcounts=Augerpeak[peakname].sum() # get counts sum 
                backgroundcnts=Augerpeak[backfitname].sum() # sum counts over identical width in background fit
                # Used for peak significance i.e. typically 2 sigma of background integration over identical width
                # full integ width is 1.2*FWHM but integwidth here is closest integer half-width
                # Write fit params from tuple over to Integresult df            
                Integresult.iloc[0]['Integcounts']=integcounts
                Integresult.iloc[0]['Backcounts']=backgroundcnts
                Integresult.iloc[0]['Significance']=round(integcounts/(np.sqrt(backgroundcnts)),3)
        # TODO add 2/sqrt(n) calc of associated percent error (also can calculate later)
        Integresult.iloc[0]['Numchannels']=integwidth
        Integresult.iloc[0]['Rsquared']=rsquared
        Integresult.iloc[0]['Element']=elem         
        # These will be n/a if fit fails
        Integresult.iloc[0]['Xc']=fitparams[0]
        Integresult.iloc[0]['Width']=fitparams[1]
        Integresult.iloc[0]['Peakarea']=fitparams[2]
        Integresult.iloc[0]['Y0']=fitparams[3]
        Integresults=Integresults.append(Integresult, ignore_index=True) # add row to list with valid 
        # end of loop through each element 

    # assign params that are common to all areas/all peaks into rows of df (copied from original log)
    for index,row in Integresults.iterrows():
        Integresults.loc[index]['Filenumber']=logmatch.Filenumber   
        Integresults.iloc[index]['Filename']=logmatch.Filename
        Integresults.iloc[index]['Filepath']=logmatch.FilePath
        Integresults.iloc[index]['Sample']=logmatch.Sample
        Integresults.iloc[index]['Comments']=logmatch.Comments
        Integresults.loc[index]['Areanumber']=areanum 
    Integresults=Integresults[mycols] # put back in original order
    return Augerfile, Integresults # df with direct peak fitting info for all areas/ all elements

Augerfile, smcountname, evbreaks

def packagequant():
    ''' Write backfit, integ, smdiff, AESquantparams to single Excel file '''
    
# Use statsmodels package for outlier detection
comp1=pd.read_csv('Smdifcomp.csv', encoding='cp437')
comp2=pd.read_csv('SmdifcompavgFe.csv', encoding='cp437')

# Find a way to do linear fits with least squares but with fixed intercept (through origin)
# web suggestion sez linalg.lstsq but that didn't work initially
# slope,intercept=np.polyfit(data1, data2, 1)  numpy version
A=np.vstack([xdata, np.ones(len(xdata))]).T  # x needs to be a column vector fo linalb.lstsqNot sure why these arrays need stacking    
slope, _, _, _ = np.linalg.lstsq(A, ydata) # fixes intercept to zero? 
            
def checklog(filelist, AugerParamLog):
    ''' Checks the user Auger parameter log matches the actual data file list from directory
    prints out filenumbers that have a problem to console''' 
    spe=[] # list of filenumbers of spe files in directory
    semmap=[] # list of filenumbers for sem or map files in directory
    for i, name in enumerate(filelist): # deals with 3 cases (spe, sem or map)
        match=re.search(r'\d+.csv',name) # spe files converted to csv
        if match:
            tempstring=match.group(0)
            num=int(tempstring.split('.')[0])
            spe.append(num) # add to spe files in dir list    
    spelist=AugerParamLog[(AugerParamLog['Areas']>=1)]
	loglist=spelist.Filename.unique()
	loglist=np.ndarray.tolist(loglist)
    missingdata=[i for i in loglist if i not in spe] # list comprehension for missing data file (but present in excel logbook)
    missingentry=[i for i in spe if i not in loglist] # list comprehension for data file with missing log entry (could also convert to sets and compare)
    for i, val in enumerate(missingdata):
        print ('Data file number ', val, ' present in Auger params log but missing from directory')
    for i, val in enumerate(missingentry):
        print ('Data file number ', val, ' present in directory but missing from Auger params log')
    # check for duplicate entries in logbook
    myset=set([x for x in loglist if loglist.count(x) > 1])
    for i in myset:
        print('Duplicate entry for file number', i, ' in Auger params log.')
    return 

def plotdffits(df, colname, nparray):
    ''' Plot a dataframe  slice on top and a numpy array (i.e its derivative) on the bottom '''
    fig, axes = plt.subplots(nrows=2, ncols=1) # axes is array
    df.plot.scatter(x='Energy', y='Counts1', ax=axes[0,0], color='r')
    df.plot.scatter(x='Energy', y='Counts1', color='r')
    
    df.plot(x='Energy', y='Counts1', ax=axes[0,0], color='r')
    
    xvals=np.arange(df.Energy.min(),df.Energy.max(),1)
    plt.plot(xvals, )
      

# MORE PLOTTING STUFF UNDER DEVELOPMENT
def addternary(df, elemlist):
    ''' Working from df with AES basis and at.%
    S, Mg+Ca, Fe '''
    # TODO get ternary calculation working from df with at.%
    

