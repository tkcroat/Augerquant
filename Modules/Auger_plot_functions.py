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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re, os, sys
from matplotlib.backends.backend_pdf import PdfPages
import scipy
import scipy.stats # load in this sequence to get linregress working
import math
from statsmodels.formula.api import ols # ordinary least squares
import ternary # python-ternary plotting module
import datetime
import tkinter as tk
from sympy import Interval, Union # for overlapping plot range removal
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules')
from Auger_utility_functions import pickelemsGUI

plt.rcParams.update({'font.size':24})
plt.rcParams.update({'legend.fontsize':18})

def openspefile(AugerFileName):
    ''' Open and return  '''
    if os.path.isfile(AugerFileName):
        Augerfile=pd.read_csv(AugerFileName)
    elif os.path.isfile('sub//'+AugerFileName):
        Augerfile=pd.read_csv('sub//'+AugerFileName)
    else: # If not found skip to next file
        print(AugerFileName, ' not found.')
        Augerfile=pd.DataFrame() # Return empty frame
    return Augerfile

def getplotboundaries(Augerfile, plotrange, AESquantparams, colname='Energy'):
    ''' Gets typical boundary of plots of given line from quantparams from plotrange but remove duplicates
    or remove if no data in range
    do not return range for duplicates (ie. Fe2 is in Fe plot range, Ca in C plot range so return 0s 
    defaults to checking for vals in energy column but can check peaks or other cols
    can also just pass string with ev range 
	passed plotrange is a list that can contain elements, ev ranges or both '''
    plotranges=[] # returns list of length 2 lists for valid elements
    
    for i, thisrange in enumerate(plotrange):
        if '-' in thisrange: # direct ev range specified in length 1 list (not list of elements)
            plotranges.append([int(thisrange.split('-')[0]),int(thisrange.split('-')[1])])
        else: # if not hyphenated it's a named Auger peak 
            thiselemdata=AESquantparams[AESquantparams['element']==thisrange]
            if len(thiselemdata)==1:
                thisrange=thiselemdata.iloc[0]['plotrange']
                try:
                    match= re.finditer(r'\d+', thisrange)
                    if match:
                        thisrange=[m.group(0) for m in match] # Parse range into lower upper
                        thisrange=[int(i) for i in thisrange]
                        Augerslice=Augerfile[(Augerfile[colname]>thisrange[0]) & (Augerfile[colname]<thisrange[1])]
                        if not Augerslice.empty:
                            plotranges.append(thisrange)
                except:
                    print ('Problem finding plot range for element ', thisrange)
                    pass # combining multiple lines in single window isn't an error
            # knock out duplicated plot ranges (i.e. just plot Fe and Fe2 together)
    # now knock out any duplicate or overlapping plot ranges, such as C and Ca, Fe and Fe2
    # this is controllable via plotrange specified in AESquantparams
    if len(plotranges)>1:
        plotranges=combineranges(plotranges)
    return plotranges

def combineranges(mylist):
    ''' Remove exact duplicates and combine ranges if overlapping
    uses sympy functions Interval (real intervals) and Union'''
    intervals=[Interval(begin,end) for [begin, end] in mylist]
    u=Union(*intervals)
    # Convert all the FiniteSets (distinct elements in unions) into nested list with min, max values for subplots
    for i in range(0, len(u.args)):
        newrange=[]
        for i in range(0, len(u.args)):
            pair=[]
            for i in iter(u.args[i].boundary):
                pair.append(int(i))
            newrange.append(pair)
    return newrange

def parsefilenums(filestr): 
    ''' Takes string containing filenumbers or ranges of filenums and returns as int list''' 
    myfiles= set()
    invalid = set()
    # tokens are comma seperated values
    tokens = [x.strip() for x in filestr.split(',')]
    for i, tok in enumerate(tokens):
        try:
            # typical tokens are ints
            myfiles.add(int(tok))
        except: # if not, then it might be a range
            try:
                token = [int(k.strip()) for k in tok.split('-')]
                if len(token) > 1:
                    token.sort()
                    # try to build a valid range from hyphenated items
                    first = int(token[0])
                    last = int(token[len(token)-1])
                    for x in range(first, last+1):
                        myfiles.add(x)
            except:
                # not an int and not a range...
                invalid.add(tok)
    # Report invalid tokens before returning valid selection
    if invalid: # test if any invalid 
        print ('Invalid set of files: ' + str(invalid))
    myfiles=list(myfiles)
    return myfiles
   
def parseareas(areas, numareas): 
    ''' takes string containing areas and returns set of ints for chosen areas, defaults to all areas if not specified''' 
    myareas= set()
    invalid = set()
    # tokens are comma seperated values
    tokens = [x.strip() for x in areas.split(',')]
    if areas=='' or areas=='all': # defaults to all areas
        for i in range(1,numareas+1):
            myareas.add(int(i))
        return myareas
    for i in tokens:
        try:
            # typical tokens are ints
            if int(i)<=numareas: # skips if value inadvertently exceeds # of areas
                myareas.add(int(i))
        except: # if not, then it might be a range
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # try to build a valid range from hyphenated items
                    first = token[0]
                    last = token[len(token)-1]
                    for x in range(first, last+1):
                        if x<=numareas:
                            myareas.add(x)
            except:
                # not an int and not a range...
                invalid.add(i)
    # Report invalid tokens before returning valid selection
    if invalid: # test if any invalid 
        print ('Invalid set of areas: ' + str(invalid))
    myareas=list(myareas) # convert from set to list
    return myareas

def getelemenergy(plotelems, plotrange, AESquantparams, deriv=False):
    ''' Pass plotted data range as tuple or string and list of desired elements for plotting, return paired 
    list of elemental and ideal energies only returns element if in range of the given plot; can return either 
    direct peak ideal position or deriv 
    peak ideal position (using deriv switch)
    negpeak is energy of peak deriv peak whereas direct peak is negpeak- integpeak (expressed as eV shift'''

    elemlines=[] # energies list
    # Return empty list of plotelems only if eV range is passed
    if '-' in plotelems[0]:
        return elemlines # return empty list if plotelems is just an eV range
    if '-' in plotrange: # determine range for plot if passed as hyphenated string
        plotrange=(float(plotrange.split('-')[0]),float(plotrange.split('-')[1])) # convert to 2 element tuple
    if len(plotrange)!=2: # passed as tuple
        print('Problem parsing plotrange lower and upper limits.')
        return
    AESquantparams=AESquantparams[AESquantparams['element'].isin(plotelems)] # select rows in element list
    AESquantparams=AESquantparams[(AESquantparams['negpeak']>plotrange[0]) & (AESquantparams['negpeak']<plotrange[1])]
    for index,row in AESquantparams.iterrows():
        try:
            if deriv==False: # return direct peak positions ()
                elemlines.append([int(AESquantparams.loc[index]['negpeak']+AESquantparams.loc[index]['integpeak']),AESquantparams.loc[index]['element']])
            else: # return ideal deriv negative peak position
                elemlines.append([AESquantparams.loc[index]['negpeak'],AESquantparams.loc[index]['element']])
        except: # some error with this element's params
            print('Problem finding ideal peak position for ', AESquantparams.loc[index]['element'])
    return elemlines # list with [energy, element name] pairs for plotting

def findevbreaks(Params, Augerfile):
    ''' Pull index # of multiplex evbreaks and return coords of necessary vertical lines for matplotlib
    Augerslice as df will have same indices as original
    logmatch should already contain info for this spectrum (evbreak index #s needed)'''   
    evbreaks=Params.Evbreaks # needed to avoid endpoints or energy breaks during peak searches
    if type(evbreaks)==str: # create list of ints from loaded string of format [0, 115, 230, 361]
        tempstring=evbreaks.split('[')[1]
        tempstring=tempstring.split(']')[0] # remove brackets commonly around written list
        evbreaks=[int(s) for s in tempstring.split(',')] # int conversion
    # find energy x values from known evbreak indices (Augerslice.index matches original index #s)
    energyvals=[] # new list for corresponding energy x vals within range
    for i in range(0,len(Augerfile)):
        index=Augerfile.index[i]
        if index in evbreaks:
            energyvals.append(Augerfile.iloc[i]['Energy']) # find energy vals assoc with index # 
    return energyvals
    
def findelemlines(plotelems, xrange, AESquantparams, peaktype='Deriv'):
    ''' Pass list of element peaks and plot range; return energies if in range
    use negpeak values if derivative plot or negpeak+integpeak for direct counts plots
    ''' 
    # plotelems=plotelems.replace('Fe','Fe Fe1 Fe2') # choose all Fe lines
    # elems=[str(s) for s in plotelems.split(' ')]
    
    elemlines=[] # energies list    
    peaks=AESquantparams[AESquantparams['element'].isin(plotelems)] # select rows in element list
    peaks=peaks[(peaks['negpeak']>xrange[0]) & (peaks['negpeak']<xrange[1])]
    if peaktype=='Deriv':
        elemlines=peaks['negpeak'].tolist()
    elif peaktype=='Counts':
        # Integpeak energy is relative to negpeak so convert to actual ev value
        peaks['integpeak']=peaks['integpeak']+peaks['negpeak']
        elemlines=peaks['integpeak'].tolist()
    return elemlines

def maketitlestring(plotpts):
    ''' Extract label for plot containing element peak amplitude, noise amplitudes (at lower and higher eV), and significance (elem amplitude/noise amplitude)
    add as title string in subplots for each peak; plotpts dataframe slice usually has 1 (but sometimes 2) elements'''
    ampl=plotpts.iloc[0]['Amplitude']
    noisel=plotpts.iloc[0]['Lowbackamplitude']
    noiseh=plotpts.iloc[0]['Highbackamplitude']
    sigl=ampl/noisel
    sigh=ampl/noiseh
    amplstr=str(int(ampl))
    noisestr=str(int(noisel))+','+str(int(noiseh))
    sigstr='%.1f' %sigl+','+'%.1f' %sigh
    if len(plotpts)>1:
        ampl2=plotpts.iloc[1]['Amplitude'] # secondary elem peak amplitude if present
        amplstr=amplstr+','+str(int(ampl2))    
        sig2l=ampl2/noisel
        sig2h=ampl2/noiseh
        sigstr=sigstr + ':'+ '%.1f' %sig2l+','+'%.1f' %sig2h
    titlestring='Peak: ' + amplstr + ' low/hi: ' + noisestr + ' Sig: ' + sigstr
    return titlestring
    
def makenoisebars(plotpts):
    ''' Create vertical line/bars showing noise amplitude determination at lower and higher energy 
    TODO do this in a bit more systematic way
    '''
    lnoiseamp=int(plotpts.iloc[0]['Lowbackamplitude'])
    hnoiseamp=int(plotpts.iloc[0]['Highbackamplitude'])
    lowenergy=int(plotpts.iloc[0]['Lowback'])
    highenergy=int(plotpts.iloc[0]['Highback'])
    tup1=(lowenergy,-lnoiseamp/2, lnoiseamp/2)
    tup2=(highenergy,-hnoiseamp/2, hnoiseamp/2)
    noisebars=[tup1,tup2]
    return noisebars
    
def organizecomp(df):
    '''Get rid of common duplicated columns  '''
    removelist=['Projectb','Filenameb','FilePathb','Sampleb','Commentsb']
    singleelemlist=['Ca','Mg','Si','S']
    for i, val in enumerate(removelist):
        if val in df:
            df=df.drop(val,axis=1)
    for i, val in enumerate(singleelemlist): # don't drop basis as these will differ if comparing smdif and integ compositions
        if val+'amplb' in df:
            df=df.drop(val+'amplb',axis=1)
    return df
            
def returnoutliers(xdata, ydata):
    '''pass xcol and ycol as lists, makes plot and return outliers with pvals below specified threshold'''    
    # convert pandas series to lists
    regression= ols("data ~ x", data=dict(data=ydata, x=xdata)).fit()
    outliers=regression.outlier_test()    
    # df with cols as student_resid, unadj_p and bonf (bonferroni)
    colnames=['Resid','Pval','Bonf']
    outliers.columns=colnames # rename columns
    return outliers


def plothist(df, spelist, Elements, col='Amplitude'):
    '''Histogram plot of amplitudes (for smdiff quant data) or adjusted counts (for integ data)
    for all passed elements ''' 
    # Set up figure of proper size
    filelist=np.ndarray.tolist(spelist.Filename.unique())
    df=df[df['Filename'].isin(filelist)] # Keep only file subset in passed list
    # Set up histogram plots
    numcols=min(len(Elements),2) # 1 or 2 columns
    numrows=math.ceil(len(Elements)/2)
    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
    for i, elem in enumerate(Elements):
        thisrow=i%numrows
        thiscol=i//numrows
        axindex=thisrow,thiscol
        if col=='Amplitude':
            tempdf=df[df['PeakID']==elem] # select only this element subset
        elif col=='Adjcnts':
            tempdf=df[df['Element']==elem]
        else:
            print('No option for histograms of column', col)
            return
        # tempdf=tempdf.rename(columns={col:elem+' '+col})
        # histogram plot
        pd.DataFrame.hist(tempdf, column=[col], bins=15, ax=axes[axindex])
        axes[axindex].set_title(elem+' '+col, fontsize=20) # change axis label to show element
    return
    
def scatplot(df, xcol='Energy', ycol='Counts', thresh=0.1):
    '''Single plot of explicitly listed cols from composition dataset (after merging and usually on compresult df) '''
    fig, axes = plt.subplots(nrows=1, ncols=1)
    df.plot.scatter(x=xcol,y=ycol, ax=axes)
    xdata=df[xcol].as_matrix()
    ydata=df[ycol].as_matrix()
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xdata, ydata) 
    # set x range for linear plot 
    text1=str(round(slope,3))+' *x +' + str(round(intercept,3))
    text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,4))
    xmax=max(xdata)
    # adding linear fit to data
    x=np.linspace(0,xmax,100) 
    axes.text(0.025,0.9, text1, fontsize=12, transform=axes.transAxes)
    axes.text(0.025,0.8, text2, fontsize=12, transform=axes.transAxes)
    # TODO set axes as active??
    plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
    # Now look for outliers and plot in separate color
    theseoutliers=returnoutliers(xdata.tolist(), ydata.tolist()) # params to find outliers
    df= pd.concat([df, theseoutliers], axis=1, join='inner')
    outliers=df[df['Pval']<thresh]
    if not outliers.empty : # plot outliers in red
        outliers.plot.scatter(x=xcol, y=ycol, ax=axes, color='r')    
    return outliers  

def getbackfitpts(backfitdf, AugerFileName, areanum):
    ''' Find points over which background fit was performed (for scatter plots)
    called by various plot and reporting functions'''
    if backfitdf.empty:
        return []
    thisfilebackpts=backfitdf[(backfitdf['Filename']==AugerFileName) & (backfitdf['Areanumber']==areanum)]
    indexptslist=[] # this gets all the lower1, lower2, upper1, upper2 index point boundaries
    thisarr=thisfilebackpts.Lower1.unique()
    thislist=np.ndarray.tolist(thisarr)
    indexptslist.extend(thislist)
    thisarr=thisfilebackpts.Lower2.unique()
    thislist=np.ndarray.tolist(thisarr)
    indexptslist.extend(thislist)
    thisarr=thisfilebackpts.Upper1.unique()
    thislist=np.ndarray.tolist(thisarr)
    indexptslist.extend(thislist)
    thisarr=thisfilebackpts.Upper2.unique()
    thislist=np.ndarray.tolist(thisarr)
    indexptslist.extend(thislist)
    indexptslist=[int(i) for i in indexptslist] 
    indexptslist.sort()
    return indexptslist

def get_plotkwargs(mykwargs, elem):
    ''' Gets subset of keyword arguments that are needed as plot arguments
    pass elem to construct x and y err col names which may be needed 
    some unpacking of compressed kwargs is needed
    
    '''
    kwargdict={} # new dict for pandas plot options only
    if 'errbars' in mykwargs:
        val=mykwargs.get('errbars')
        hasxerr=mykwargs.get('errx', False) # true if x error col exists
        hasyerr=mykwargs.get('erry', False) 
        if 'x' in val and hasxerr:
            # col name after merge depends on both hasxerr and hasyerr
            if hasyerr:
                colname='err'+elem+'1'
            else:
                colname='err'+elem
            kwargdict.update({'xerr':colname}) # name of x error column 
        if 'y' in val and hasyerr: # either errSi or errSib depending on merge
            if hasxerr:
                colname='err'+elem+'2' # suffix only appended for duplicate cols
            else:
                colname='err'+elem
            kwargdict.update({'yerr':colname}) # name of x error column 
    if 'log' in mykwargs: # Also can use log scale on plots 
        val=mykwargs.get('log')
        if val=='xy':
            kwargdict.update({'loglog':True})
        elif val=='x':
            kwargdict.update({'logx':True})
        elif val=='y':
            kwargdict.update({'logx':True})
    return kwargdict

''' Plot unit testing
comp1=Smdifcomp
comp2=Integcomp
kwargs1=get_plotkwargs(kwargs, elem)
kwargs.update({'errx':False})
compdata, outliers=AESplot.scattercompplot(Smdifcomp,Integcomp, Elements, **kwargs)
'''

def scattercompplot(comp1, comp2, elemlist, **kwargs):
    '''Pass two versions of composition calculation (using different lines or whatever) and compare 
    major elements using scatter graphs .. single point for each sample
    uses inner merge to select only subset with values from each df
    elemlist: set of elements to compare in scatter plots (multiple subplots created) 
    Keyword args:
    joinlist: list of columns to use for merge of two composition files
        filenumber, areanumber is default and used to compare compositions from exact same spe files 
        but computed via sm-diff or integ methods 
        sample - returns all sample inner merge matches (some of which could be ), such as if a 
        background region were also measured as one of the areas in an spe file
    thresh - always used; distinguishes outlier points in scatter plot; defaults to 0.1
          higher gives more outliers which are returned and plotted separately
    basis - bool (true means plot basis, false is defaul at.% plot)
    errbars - optional plotting of x, y or xy error bars
    '''
    
    elemlist=[re.match('\D+',i).group(0) for i in elemlist] 
    # strip number from peaks like Fe2 if present; columns will be element names (Fe) not peak names (Fe2)
    basis=kwargs.get('basis', False) # boolean for plot of at.% (default) or direct basis
    if basis: # plot basis if true, default at.% if false
        elemlist=['%'+s for s in elemlist]            
    numregions=len(elemlist)
    # Check if error columns are present in passed comp dataset (test w/ elem1 but should be same for all)
    tempname='err'+elemlist[0]
    if tempname in comp1:
        kwargs.update({'errx':True})
    if tempname in comp2:
        kwargs.update({'erry':True})
    # set nrows and ncols for figure of proper size
    cols=divmod(numregions,2)[0]+ divmod(numregions,2)[1]
    if numregions>1:
        rows=2
    else:
        rows=1        
    fig, axes = plt.subplots(nrows=rows, ncols=cols, squeeze=False) # axes is array
    plt.tight_layout()  # normally tight layout needed
    # Get choice for how to merge two separate compositions 
    joinlist=kwargs.get('joinlist',['Filenumber','Areanumber']) # defaults to file# & area#
    # Merge dfs with comp1 and comp2 using inner join 
    compdata=pd.merge(comp1, comp2, how='inner', on=joinlist, suffixes=('','b'))
    mycols=compdata.dtypes.index.tolist() # list of columns from full dataset
    # Dataframes for full returned datasets
    outliers=pd.DataFrame(columns=mycols) # empty dataframe for outlying points only
    fulldata=pd.DataFrame(columns=mycols)
    newcols=['Resid','Pval','Bonf', 'Element'] # single column for residuals but separate value needed for each row per element
    mycols.extend(newcols)
    for i, cname in enumerate(newcols):
        outliers[cname]=''
        fulldata[cname]=''
    # now create scatter plot for each requested element
    for i,elem in enumerate(elemlist):
        # new version of base compositional data for each loop (otherwise reindexing problems)
        compdata=pd.merge(comp1, comp2, how='inner', on=joinlist, suffixes=('1','2')) 
         # determine which subplot to use
        if (i+1)%2==1:
            rownum=0
        else:
            rownum=1
        colnum=int((i+1)/2.1)       
        xcol=elem+'1'
        ycol=elem+'2' # same element from second dataset
        axindex=rownum, colnum # 
        # Check if error bars needed and if so check/find/return error bar columns
        kwargs1=get_plotkwargs(kwargs, elem)
        compdata.plot.scatter(x=xcol, y=ycol, ax=axes[axindex], **kwargs1) # single plot axes has no [#,#]
        # linear regression: fitting, plot and add labels
        xdata=compdata[xcol].as_matrix() # this data column as np array
        ydata=compdata[ycol].as_matrix()
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xdata, ydata) # imported from scipy.stats
        # set x range for linear plot 
        text1=str(round(slope,3))+' *x +' + str(round(intercept,3))
        text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,4))
        x=np.linspace(0,max(xdata)*1.05,100) # set range for lin regress plot (slightly larger than xmax)
        # Add labels to correct position in this subplot 
        axes[rownum,colnum].text(0.025,0.9, text1, fontsize=12, transform=axes[axindex].transAxes)
        axes[rownum,colnum].text(0.025,0.8, text2, fontsize=12, transform=axes[axindex].transAxes)
        plt.axes(axes[axindex]) # set correct axes as active
        plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        # Now test, plot and return outliers
        theseoutliers=returnoutliers(xdata.tolist(), ydata.tolist()) # index # 
        # Add residual and pval to each compositional comparison line
        compdata= pd.concat([compdata, theseoutliers], axis=1, join='inner') # same length so just join by index
        compdata['Element']=elem # set after concat operation
        thresh=kwargs.get('thresh',0.1) # 0.1 default value (always using some outlier threshold)
        theseoutliers=compdata[compdata['Pval']<thresh] # now filter by threshold for outliers (with all cols)         

        # Now set correct axes name which depends on # of regions
        if not theseoutliers.empty: # plot outliers on top in different color
            theseoutliers.plot.scatter(x=xcol, y=ycol, ax=axes[axindex], color='r', **kwargs1)
            # kwargs1 contains x and/or y error columns if requested
        # hide subplots with no data
        for i in range(0,rows*cols):
            if i>len(elemlist)-1: # hide this axis
                if (i+1)%2==1: # clunky way of finding proper 
                    thisrow=0
                else:
                    thisrow=1
                thiscol=int((i+1)/2.1)  
                axindex=thisrow, thiscol # tuple to index axes 
                axes[axindex].set_visible(False)     
        outliers=outliers.append(theseoutliers) # outliers from all elements
        fulldata=fulldata.append(compdata)
        # possibly could use ignore index but probably either is fine
    # plt.tight_layout()
    # TODO fix occasional error "Nonetype has no attribute is_bbox"
    fulldata=fulldata[mycols] # put back in original column order
    outliers=outliers[mycols]
    fulldata=organizecomp(fulldata) # duplicated elements for same filenum have element-specific fit results
    outliers=organizecomp(outliers)
    return fulldata, outliers

def scatterratioplot(comp1, comp2, elem1, elem2, basis=False):
    '''Pass two versions of composition calculation (using different lines or whatever) and compare 
    major elements using scatter graphs .. single point for each sample
    uses inner merge to select only subset with values from each df'''
    elemlist=[elem1, elem2] 
    # strip number from peaks like Fe2 if present; columns will be element names (Fe) not peak names (Fe2)
    if basis==False: # use atomic % (which is the default), not basis for each element
        elemlist=['%'+s for s in elemlist]
          
    fig, axes = plt.subplots(nrows=1, ncols=1) # axes is array
    # merge dfs with comp1 and comp2 using inner join
    comp1['Ratio']=comp1[elem1]/comp1[elem2] 
    # replace infinite vals with nan and drop them
    comp1=comp1.replace([np.inf, -np.inf], np.nan).dropna(subset=['Ratio'])
    comp2['Ratio']=comp2[elem1]/comp2[elem2] 
    # replace infinite vals with nan and drop them
    comp2=comp2.replace([np.inf, -np.inf], np.nan).dropna(subset=['Ratio'])
    
    df=pd.merge(comp1, comp2, how='inner', on=['Sample'], suffixes=('','b'))    
    
    xcol='Ratio'
    ycol='Ratiob' # same element from second dataset
    
    df.plot.scatter(x=xcol, y=ycol, ax=axes) # single plot axes has no [#,#]
    data1=df[xcol]
    data2=df[ycol]
    # slope,intercept=np.polyfit(data1, data2, 1)  numpy version
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(data1, data2) # imported from scipy.stats
    # set x range for linear plot 
    text1=str(round(slope,2))+' *x +' + str(round(intercept,2))
    text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,3))
    xmax=max(max(data1),max(data2))*1.1 # set to slightly larger than max of dataset
    x=np.linspace(0,xmax,100) # setting range for 
    axes.text(0.025,0.9, text1, fontsize=12, transform=axes.transAxes)
    axes.text(0.025,0.8, text2, fontsize=12, transform=axes.transAxes)
    plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
    return df

def plotternary(comp1, ternelems):
    ''' Take compositional data, compute as 3-tuples and plot on ternary diagram '''
    # Calculate ternary quantities for all in duplicates dataframe (allows each axis to be sum of elements)
    hyphensplit = re.compile('(\+[a-zA-Z]+)').split
    ternlist=[part for img in ternelems for part in hyphensplit(img) if part]
    ternlist=[x.replace('+','') for x in ternlist]
    try:
        comp1['Tbasis']=comp1[ternlist].sum(axis=1) # sum a list of columns 
    except:
        print('Failed sum.. missing data for given element?')
    # calculate T1, T2, T3 values for ternary plot (auto normalized to unity)
    for i, val in enumerate(ternelems):
        num=str(i+1)
        if '+' not in val:
            comp1['T'+num]=comp1[val]/comp1['Tbasis']
        else:
            elems=[str(s) for s in val.split('+')]
            comp1['T'+num]=comp1[elems].sum(axis=1)/comp1['Tbasis']
    # create list with 3 points as tuples
    plotpts=[]
    for index, row in comp1.iterrows():
        plotpts.append((comp1.loc[index]['T1'],comp1.loc[index]['T2'],comp1.loc[index]['T3']))
        
    # create ternary plot
    figure,tax = ternary.figure(scale=1.0)
    fontsize=20
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=0.1, linewidth=0.5) # denser lines with higher multiple
    tax.set_title("Ternary Composition Plot", fontsize=fontsize)
    tax.left_axis_label(ternelems[2], fontsize=fontsize)
    tax.right_axis_label(ternelems[1], fontsize=fontsize)
    tax.bottom_axis_label(ternelems[0], fontsize=fontsize)
    tax.ticks(axis='lbr', linewidth=1, multiple=0.1) # set ticks
    tax.clear_matplotlib_ticks() # remove default Matplotlib axes (non-existent x and y axes)
    tax.scatter(plotpts, marker='s', s=40, color='b')  # s is point size
    return # return with ternary compositional data

def AESplot1(files, spelist, AESquantparams, **kwargs): # single file already selected by number
    '''Plot a single frame over desired range with passed filenumbers/ areas 
    plotelems kwarg is list of elements who's ideal position will be shown as blue vert line on plot
    filenumbers - Single # or list of numbers
    TODO add yrange as kwarg?
    kwargs:  xrange, yrange 
        plotelems - set of elements to label 
        backfitdf -- optional background fitting  
        smdifpeaks -- optional plot of pts in S7D7 used for amplitude determination
        areas - text string with areas for inclusion in plot (csv or ranges with dash)
        plotcol - Counts or S7D7 (or Both)
        '''
    # convert passed files string or int into actual filenumber list 
    if isinstance(files, int):
        filenums=[files] # convert single filenumber to len 1 list
    elif isinstance(files, str):
        filenums= parsefilenums(files) # convert string of vaious and ranges to int list of filenumbers
    else:
        print('Files for plotting must be string w or w/o ranges or single integer')
        return
    myfiles=spelist[spelist.Filenumber.isin(filenums)] # filename if logmatch is series
    if len(myfiles)==0:
        print('No files of specified number selected.')
        return
    # convert pass4ed areas string 
    areas=kwargs.get('areas',1)
    numareas=int(myfiles.iloc[0]['Areas'])
    if isinstance(areas, str):
        areas = parseareas(areas, numareas)
    elif isinstance(areas, int):
        areas = [areas] # len 1 list
    else:
        print('Files for plotting must be string w or w/o ranges or single integer')
        return
    # Plot of all filenumbers in list and all areas in areas list
    xrange=kwargs.get('xrange', [100,1900]) # string or list 
    # Things added to plots (breaks in multiplex scans, points used for S7D7 amplitudes)
    
    AugerFileName=myfiles.iloc[0]['Filename']
    Augerfile=openspefile(AugerFileName)
    if Augerfile.empty: # file not found problem
        return
    energyvals=findevbreaks(myfiles.iloc[0], Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
           
    plotcol=kwargs.get('plotcol','Counts') # defaults to counts
    if plotcol=='Both':
        fig, axes = plt.subplots(nrows=2, ncols=1, squeeze=False, sharex=True) # axes is array
    else:
        fig, axes = plt.subplots(nrows=1, ncols=1, squeeze=False) # axes is array
    
    # Use color and linestyles to differentiate filenumber(s) and areanumber(s)
    mylegend=[]
    colorlist=['b','r','g','c','m','y','k', 'olive','pink','purple']
    linestyles=['solid','dashed','dashdot','dotted','--',':']
    '''
    # TODO does colorlist represent different areas or different filenumbers?
    # Default is to use color for areas (different sample regions) as filenumbers are often 
    temporal reruns (time slices of same underlying spectrum possibly w/ eneryg shifts)
    TODO switch to handling these using groupby
    '''      
    filecount=0
    
    for index, row in myfiles.iterrows():
        try:
            AugerFileName=myfiles.loc[index]['Filename']
            Augerfile=openspefile(AugerFileName)
            if Augerfile.empty: # file not found
                continue
            # Also check in sub directory
            filecount+=1
            thisfilenum=myfiles.loc[index]['Filenumber']
            Augerslice=Augerfile[(Augerfile['Energy']>xrange[0]) & (Augerfile['Energy']<xrange[1])]
            areacount=0 # differs from areanum in that always starting from zero
            for i, areanum in enumerate(areas):
                # default is color for areas and linestyle for 
                plkwargs={'color':colorlist[areacount], 'linestyle':linestyles[filecount%6]}
                if len(areas)<3: # switch this if using small # of areas
                      plkwargs={'color':colorlist[filecount], 'linestyle':linestyles[areacount%6]}
                areacount+=1
                colname=plotcol+str(areanum)
                mylegend.append(str(thisfilenum)+' a'+str(areanum))
                # Augerslice.plot(x='Energy', y=colname, ax=axes, linestyle=thisls, color=colorlist[areanum])
                if plotcol=='Both':
                    Augerslice.plot(x='Energy', y='Counts'+str(areanum), ax=axes[0,0], **plkwargs)
                    Augerslice.plot(x='Energy', y='S7D7'+str(areanum), ax=axes[1,0], **plkwargs)
                else:
                    Augerslice.plot(x='Energy', y=colname, ax=axes[0,0], **plkwargs)
                # Optional additional plotting of backgrounds (only for direct counts plotting)
                if 'backfitdf' in kwargs and plotcol!='S7D7': # do for counts and both selections
                    # plot background itself always in axes 0, 0
                    Augerslice.plot(x='Energy', y='Backfit'+str(areanum), ax=axes[0,0])
                    # scatter plot points over which background was fitted
                    backfitdf=kwargs.get('backfitdf', pd.DataFrame())
                    indexptslist=getbackfitpts(backfitdf, AugerFileName, areanum)
                    backpts=Augerslice[Augerslice.index.isin(indexptslist)] # gets background fitted pts but only from this data slice
                    if not backpts.empty: # show fitted pts from counts
                        backpts.plot.scatter(x='Energy', y='Counts'+str(areanum), ax=axes[0,0])
                # Optional plotting of smdiff points for amplitude calc (only for differentiated plots)
                if 'smdifpeaks' in kwargs and plotcol!='Counts':
                    smdiffpts=kwargs.get('smdifpeaks', pd.DataFrame())
                    # filter for this filenumber and areanumber
                    smdiffpts=smdiffpts[(smdiffpts['Filenumber']==thisfilenum) & (smdiffpts['Areanumber']==areanum)]
                    smdiffpts=smdiffpts[(smdiffpts['Peakenergy']>xrange[0]) & (smdiffpts['Peakenergy']<xrange[1])]
                    # plot background itself
                    if not smdiffpts.empty: # show fitted pts from counts
                        smdiffpts['Pospeak']=smdiffpts['Peakenergy']-smdiffpts['Peakwidth']
                        if plotcol=='Both':
                            axindex=axes[1,0] # plot in lower 
                        else:
                            axindex=axes[0,0] # main single window
                        smdiffpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axindex, color='r')
                        smdiffpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axindex, color='r')
        except:
            print('Error plotting file', str(thisfilenum))
    # Alter legend with filenumber+ areanumber
    axes[0,0].legend(mylegend, loc='best')
    for axis in ['top', 'bottom','left','right']:
        axes[0,0].spines[axis].set_linewidth(2.5) # change axes thickness
        try: # sometimes plotting both counts and deriv
            axes[1,0].spines[axis].set_linewidth(2.5) 
        except:
            pass
    # Energy breaks (for multiplex files and elemental peak locations same for all plots
    for i, val in enumerate(energyvals):
        axes[0,0].axvline(x=val, color='r') 
        try:
            axes[1,0].axvline(x=val, color='r') 
        except:
            pass
    # Add elemental peaks as vert lines (but slightly different position for deriv vs counts)
    # TODO add elem/peak label next to line
    plotelems=kwargs.get('plotelems',[])
    if len(plotelems)>0:
        if plotcol=='Counts':
            elemlines = findelemlines(plotelems, xrange, AESquantparams, peaktype='Counts')
            for i, val in enumerate(elemlines):
                axes[0,0].axvline(x=val, color='b')
        elif plotcol=='S7D7':
            elemlines = findelemlines(plotelems, xrange, AESquantparams, peaktype='Deriv')
            for i, val in enumerate(elemlines):
                axes[0,0].axvline(x=val, color='b')
        elif plotcol=='Both': # counts plot is on top
            elemlines = findelemlines(plotelems, xrange, AESquantparams, peaktype='Counts')
            for i, val in enumerate(elemlines):
                axes[0,0].axvline(x=val, color='b')
            elemlines = findelemlines(plotelems, xrange, AESquantparams, peaktype='Deriv')
            for i, val in enumerate(elemlines):
                axes[1,0].axvline(x=val, color='b')
    return

def reportSD(paramlog, Smdifpeakslog, plotelems, AESquantparams, PDFname='SDplot_report'):
    ''' generalized newer version of SDmajor, plots all spe files, associated quant points, and 
    labels selected elemental lines  
    '''
    plt.ioff()
    if PDFname=='SDplot_report': # using default name with date
        now=datetime.datetime.now()
        PDFname=PDFname+'_'+now.strftime('%d%b%y')+'.pdf'
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])
            Augerfile=openspefile(AugerFileName)
            if Augerfile.empty:
                continue
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            thisfilepeaks=Smdifpeakslog[(Smdifpeakslog['Filename']==AugerFileName)] # retrieve assoc. subset of peaks data
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # determine size of plot for this filenumber (same for all areas)
            plotranges=getplotboundaries(Augerfile, plotelems, AESquantparams) # returns plot ranges for all regions with data from plotelems

            for i in range(0,numareas): # create new plot for each area 
                areanum=i+1 #
                thisareapeaks=thisfilepeaks[(thisfilepeaks['Areanumber']==areanum)] # then only this areanum
                # set plot rows and columns
                numcols=min(len(plotranges),2) # 1 or 2 columns
                numrows=math.ceil(len(plotranges)/2)
                # new plot for each spatial area for each filenumber
                try:
                    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
                    colname='S7D7'+str(areanum)
                    mytitle=AugerFileName+' '+str(filenumber)+' area #'+str(areanum)
                    plt.suptitle(mytitle)
                    for j, bounds in enumerate(plotranges):
                        [lower, upper]=bounds                    
                        thisrow=j%numrows
                        thiscol=j//numrows
                        if len(plotranges)<7:
                            axindex=thisrow,thiscol
                        else:
                            axindex=thiscol,thisrow # flips to 2 row, n column style
                        Augerslice=Augerfile[(Augerfile['Energy']>lower) & (Augerfile['Energy']<upper)] # already known that this isn't empty
                        Augerslice.plot(x='Energy', y=colname, ax=axes[axindex])
                        # Section for labeling plotelements
                        elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=True) # can pass plot range as lower,upper tuple
                        # list of tuples with energy,elemname
                        for k, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                            try:
                                axes[axindex].axvline(x=elemtuple[0], color='b') # O line
                                axes[axindex].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                            except:
                                print('Problem labeling elements')
                        # Section for adding smooth-diff quant data
                        plotpts=thisareapeaks[(thisareapeaks['Peakenergy']>lower) & (thisareapeaks['Peakenergy']<upper)]
                        if not plotpts.empty:
                            try:
                                plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[axindex], color='r')
                                plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[axindex], color='r')                                        
                                titlestring=maketitlestring(plotpts)
                                axes[axindex].set_title(titlestring, fontsize=10)
                                noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                                axes[axindex].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                                axes[axindex].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
                            except:
                                print('Problem adding points from smdif quant calcs' )
                        # add red vert line at multiplex energy break if present                    
                        for l, val in enumerate(energyvals):
                            if val > lower and val < upper: 
                                axes[axindex].axvline(x=val, color='r') # on counts plot
                    # now hide empty subplots (common for odd number of elemental regions)
                    for subplt in range(0,numrows*numcols):
                        if subplt>len(plotranges)-1:
                            thisrow=subplt//numcols
                            thiscol=subplt%numcols
                            axindex=thisrow, thiscol # tuple to index axes 
                            axes[axindex].set_visible(False)
                    pdf.savefig(fig)
                    plt.close(fig) # closes recently displayed figure (since it's saved in PDF report)
                    print(AugerFileName,' area', areanum, 'plotted')
                except:
                    print('Unknown problem plotting area ', str(areanum),' for ', AugerFileName)
    plt.ion()
    return

def reportcountsback(spelist, plotrange, AESquantparams, Backfitlog, PDFname='Countsback_report', **kwargs):
    ''' Plot of list of passed elements 
    args:
        spelist - log entries for spe files (default all plotting but can be trimmed with kwarg)
        plotrange- list of elements, hyphenated eV ranges or both
        AESquantparams- peak params for known Auger peaks
        Backfitlog always passed (even if empty): contains pts used for backfits (optionally plotted)
        PDFname- optional name for report (or use default of Countsback_report_7May17.pdf)
    
    kwargs: 
        'filesubset' - filenumbers to plot (defaults to all)
        'addbackfit' - optional plot of background fit col (stored in each Augerfile)
        'backfitpts' - plot of points used for background fit (useful for troubleshooting)
		'plotelems; - list of elements/peaks for labelling
                ... dataframe passed as value (not a bool)
    Evbreaks from multiplex plotted as red lines (often different dwell times for different elements)
    '''
    plt.ioff()
    if 'filesubset' in kwargs: # only plot report for subset of filenumbers (defaults to all)
        filenumlist=kwargs.get('filesubset',[])
        spelist=spelist[spelist.Filenumber.isin(filenumlist)]

    with PdfPages(PDFname) as pdf:
        for index,row in spelist.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=spelist.loc[index]['Filename']
            numareas=int(spelist.loc[index]['Areas'])
            # find and open file (or skip if not found)
            Augerfile=openspefile(AugerFileName)
            if Augerfile.empty:
                continue
            Params=spelist.loc[index] # grab row for this spe file as Series
            # filenumber=Params.Filenumber # retrieve filenumber            
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # Determine plot ranges from string of elements, evranges or both 
            # Automatic knockout of regions with no data;  same range for all areas within given filenumber
            plotranges=getplotboundaries(Augerfile, plotrange, AESquantparams) # returns plot ranges for all regions with data from plotelems
            # TODO optional plot of all areas on same page (currently all on separate)
            for i in range(0,numareas): # Create separate plot page for each area 
                areanum=i+1
                # optional plotting of points used to fit backgrounds (often a useful diagnostic)
                if kwargs.get('backfitpts', False):
                    indexptslist=getbackfitpts(Backfitlog, AugerFileName, areanum)
                # Set plot rows and columns
                numcols=min(len(plotranges),2) # 1 or 2 columns
                numrows=math.ceil(len(plotranges)/2)
                # new plot for each spatial area
                try:
                    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by 3 axes array
                    colname='Counts'+str(areanum)
                    mytitle=AugerFileName +' area #'+str(areanum)
                    plt.suptitle(mytitle)
                    # Now loop over the elemental plot ranges
                    for j, bounds in enumerate(plotranges):
                        [lower, upper]=bounds                    
                        thisrow=j%numrows
                        thiscol=j//numrows
                        axindex=thisrow, thiscol
                        Augerslice=Augerfile[(Augerfile['Energy']>=lower) & (Augerfile['Energy']<=upper)] # already known that this isn't empty
                        Augerslice.plot(x='Energy', y=colname, ax=axes[axindex]) # plot counts
                        if kwargs.get('addbackfit', False)==True: # optional plot of background fits
                            backfitcol='Backfit'+str(areanum)
                            # shouldn't matter if no vals are in backfit col
                            Augerslice.plot(x='Energy', y=backfitcol, ax=axes[axindex])
                        if kwargs.get('backfitpts', False):
                            # Now add scatter plot points at fit region boundaries
                            backpts=Augerslice[Augerslice.index.isin(indexptslist)] # gets background fitted pts but only from this data slice
                            if not backpts.empty: # show fitted pts from counts
                                backpts.plot.scatter(x='Energy', y=colname, ax=axes[axindex])
                        # Section for labeling plotelements
                        if 'plotelems' in kwargs:
                            plotelems=kwargs.get('plotelems',[]) # list of element peaks for labeling
                            elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=False) # can pass plot range as lower,upper tuple
                            # list of tuples with energy,elemname
                            for k, elemtuple in enumerate(elemlines):
                                # elemtuple[0] is energy and [1] is element symbol
                                # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                try:
                                    axes[axindex].axvline(x=elemtuple[0], color='b') # O line
                                    yval=(Augerslice[colname].max()-Augerslice[colname].min())*0.9+Augerslice[colname].min()
                                    axes[thisrow,thiscol].text(elemtuple[0],yval, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                                except:
                                    print('Problem labeling elements')
                            # fold 
                        # Add red vertical lines at multiplex energy breaks 
                        for l, val in enumerate(energyvals):
                            if val > lower and val < upper: 
                                axes[thisrow,thiscol].axvline(x=val, color='r') # on counts plot 
                    for subplt in range(0,numrows*numcols): # hide any empty subplots
                        if subplt>len(plotranges)-1:
                            axes[subplt%numrows,subplt//numrows].set_visible(False)
                    pdf.savefig(fig)
                    plt.close(fig) # closes recently displayed figure (since it's saved in PDF report)
                    print(AugerFileName,' area', areanum, 'plotted') # end of long try plotting entire area
                except:
                    print('Problem plotting file ', AugerFileName, 'area', areanum, ' likely no data for specified elements.')
    plt.ion()
    return

def reportderivcnt(spelist, plotrange, AESquantparams, Backfitlog, Smdifpeakslog, PDFname='this_report.pdf', **kwargs):
    ''' Comparison plots for both derivative and counts itself (don't have to worry about axes w/o indexing)
    plots selected filenumber for all area (looped)
    plots all spe files, associated quant points, and labels selected elemental lines
    args:
        spelist - log entries for spe files (default all plotting but can be trimmed with kwarg)
        plotrange- list of elements, hyphenated eV ranges or both
        AESquantparams- peak params for known Auger peaks
        Backfitlog always passed (even if empty): contains pts used for backfits (optionally plotted)
        Smdifpeakslog -always passed (even if empty): contains pts used for amplitude calcs
        PDFname- optional name for report (or use default of Countsback_report_7May17.pdf)
    kwargs:
        'filesubset' - filenumbers to plot (defaults to all)
        'addbackfit' - optional plot of background fit col (stored in each Augerfile)
        'backfitpts' - plot of points used for background fit (useful for troubleshooting)
        'smdifpts' - plot of points used for background fit (useful for troubleshooting)
                ... dataframe passed as value (not a bool)
    '''
    plt.ioff()
    if 'filesubset' in kwargs: # only plot report for subset of filenumbers (defaults to all)
        filenumlist=kwargs.get('filesubset',[])
        spelist=spelist[spelist.Filenumber.isin(filenumlist)]
    with PdfPages(PDFname) as pdf:
        for index,row in spelist.iterrows():
            AugerFileName=spelist.loc[index]['Filename']
            numareas=int(spelist.loc[index]['Areas'])
            Augerfile=openspefile(AugerFileName)
            if Augerfile.empty:
                continue
            # Find evbreaks if multiplex (identical for each filenumber)
            energyvals=findevbreaks(spelist.loc[index], Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # find ranges for subplots (from element list or eV range or both)
            plotranges=getplotboundaries(Augerfile, plotrange, AESquantparams) # returns plot ranges for all regions with data from plotelems
            if 'smdifpts' in kwargs:
                thisfilepeaks=Smdifpeakslog[Smdifpeakslog['Filename']==AugerFileName] # retrieve assoc. subset of peaks data
            # SKIP energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # determine size of plot for this filenumber (same for all areas); plotranges combines some plotelems (i.e. C and Ca together)
            for i in range(0,numareas):
                areanum=i+1
                # set plot rows and columns
                numrows=2 # counts on top and S7D7 on bottom
                numcols=len(plotranges)
                if kwargs.get('backfitpts', False):
                    indexptslist=getbackfitpts(Backfitlog, AugerFileName, areanum)
                try:
                    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
                    mytitle=AugerFileName.replace('.csv','') +' area #'+str(areanum)                
                    if len(plotranges)>4:
                        plt.tight_layout() # shrinks to fit axis labels
                    plt.suptitle(mytitle)
                    for j, bounds in enumerate(plotranges):
                        [lower, upper]=bounds                    
                        thiscol=j
                        Augerslice=Augerfile[(Augerfile['Energy']>lower) & (Augerfile['Energy']<upper)] # already known that this isn't empty
                        Augerslice.plot(x='Energy', y='S7D7'+str(areanum), ax=axes[0,thiscol]) # deriv in upper plot
                        Augerslice.plot(x='Energy', y='Counts'+str(areanum), ax=axes[1,thiscol]) # counts in lower 
                        if kwargs.get('addbackfit', False)==True: # optional plot of background fits
                            backfitcol='Backfit'+str(areanum)
                            # plot backfit in lower along with counts
                            if backfitcol in Augerslice: # avoid error if no backfit yetperformed
                                Augerslice.plot(x='Energy', y=backfitcol, ax=axes[1,thiscol])
                        if kwargs.get('backfitpts', False):
                            # Now add scatter plot points at fit region boundaries
                            backpts=Augerslice[Augerslice.index.isin(indexptslist)] # gets background fitted pts but only from this data slice
                            if not backpts.empty: # show fitted pts from counts
                                backpts.plot.scatter(x='Energy', y='Counts'+str(areanum), ax=axes[1,thiscol])
                        # Section for labeling plot elements
                        if 'plotelems' in kwargs:
                            plotelems=kwargs.get('plotelems',[]) # list of element peaks for labeling
                            elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=True) # can pass plot range as lower,upper tuple
                            # list of tuples with energy,elemname
                            for k, elemtuple in enumerate(elemlines):
                                # elemtuple[0] is energy and [1] is element symbol
                                # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                try:
                                    axes[0,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                    axes[0,thiscol].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                                except:
                                    print('Problem labeling elements')
                            # same labeling of elements in upper plot (w/ peaks slightly shifted)
                            elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=False) # can pass plot range as lower,upper tuple
                            # list of tuples with energy,elemname
                            for k, elemtuple in enumerate(elemlines):
                                # elemtuple[0] is energy and [1] is element symbol
                                # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                try:
                                    axes[1,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                    axes[1,thiscol].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                                except:
                                    print('Problem labeling elements')

                        # Section for adding smooth-diff quant data
                        if 'smdifpts' in kwargs:
                            thisareapeaks=thisfilepeaks[thisfilepeaks['Areanumber']==areanum] 
                            plotpts=thisareapeaks[(thisfilepeaks['Peakenergy']>lower) & (thisareapeaks['Peakenergy']<upper)]
                            if not plotpts.empty:
                                try:
                                    plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,thiscol], color='r')
                                    plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,thiscol], color='r')                                        
                                    titlestring=maketitlestring(plotpts)
                                    axes[0,thiscol].set_title(titlestring, fontsize=10)
                                except:
                                    print('Problem adding points from smdif quant calcs for ', AugerFileName,'area # ', areanum )
                        # add red vert line at multiplex energy break if present
                        for l, val in enumerate(energyvals):
                            if val > lower and val < upper: 
                                axes[0,thiscol].axvline(x=val, color='r') # on deriv plot
                                axes[1,thiscol].axvline(x=val, color='r') # on counts plot 
                    print(AugerFileName, str(areanum), 'plotted.')
                    pdf.savefig(fig)
                    plt.close(fig)
                except:
                    print('Unknown problem plotting', AugerFileName,' area #', areanum)
    plt.ion()
    return
    
def reportpeaksall(spelist, plotelems, AESquantparams, PDFname='Peaks_report.pdf'):
    ''' Generalized version of reportpeaksmajor working from spelist (loop over spatial areas 
    plots subtracted data for all spe files
    '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in spelist.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=spelist.loc[index]['Filename']
            numareas=int(spelist.loc[index]['Areas'])  
            Augerfile=openspefile(AugerFileName)
            if Augerfile.empty:
                continue
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            filenumber=spelist.loc[index]['Filenumber'] # retrieve filenumber
            # skip evbreaks as surely no breaks are present in narrow region around peak
            # determine size of plot for this filenumber (same for all areas)
            plotranges=getplotboundaries(Augerfile, plotelems, AESquantparams, colname='Peaks1') # returns plot ranges for data-containing regions
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 #
                # Determine # rows and columns from len(plotranges)
                numcols=min(len(plotranges),2) # 1 or 2 columns
                numrows=math.ceil(len(plotranges)/2)
                try:
                    if len(plotranges)<7:                    
                        fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
                    else:
                        fig, axes = plt.subplots(nrows=numcols, ncols=numrows, figsize=(16,9), squeeze=False) # switch to 2 row style for >7 subplots
                    colname='Peaks'+str(areanum)
                    mytitle=str(filenumber)+' area #'+str(areanum)
                    plt.suptitle(mytitle)
                    for j, bounds in enumerate(plotranges):
                        [lower, upper]=bounds                    
                        thisrow=(j)//2 # correct zero based indexing
                        thiscol=(j)%2
                        if len(plotranges)<7:
                            axindex=thisrow,thiscol
                        else:
                            axindex=thiscol,thisrow # flips to 2 row, n column style
                        Augerslice=Augerfile[(Augerfile['Energy']>lower) & (Augerfile['Energy']<upper)] # already known that this isn't empty
                        Augerslice.plot(x='Energy', y=colname, ax=axes[axindex])
                        # Set tighter plot range for peaks only                        
                        Augerslice=Augerslice.dropna(subset=[colname]) # drop na values
                        myplotrange=[Augerslice['Energy'].min()*.95,Augerslice['Energy'].max()*1.05] # 
                        axes[axindex].set_xlim(myplotrange) 
                        # Section for labeling plotelements
                        elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=False) # can pass plot range as lower,upper tuple
                        axes[axindex].axhline(y=0)
                        # list of tuples with energy,elemname
                        for k, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            try:
                                axes[axindex].axvline(x=elemtuple[0], color='b') # O line
                                axes[axindex].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                            except:
                                print('Problem labeling elements')                            
                        # skip evbreaks (or could get them from AugerParamLog if really necessary)
                except:
                    print('Unknown problem with single element multiplex.')
                print(AugerFileName, 'area', str(areanum), 'plotted.')
                pdf.savefig(fig)
                plt.close(fig) # closes recently displayed figure (since it's saved in PDF report)
    plt.ion()
    return

def reportpeaks(dfselect, plotelems, AESquantparams, PDFname='Peaks_report.pdf'):
    ''' Peaks report for subset of files (i.e. outliers) and only plots listed filenum and areanum (no all areas loop
    plots subtracted data for all spe files... seldom used (normally reportpeaksall)
    '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in dfselect.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=dfselect.loc[index]['Filename']
            areanum=dfselect.loc[index]['Areanumber']  
            Augerfile=openspefile(AugerFileName)
            if Augerfile.empty:
                continue

            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            filenumber=dfselect.loc[index]['Filenumber'] # retrieve filenumber

            # determine size of plot for this filenumber (same for all areas)
            plotranges=getplotboundaries(Augerfile, plotelems, AESquantparams, colname='Peaks1') # returns plot ranges for data-containing regions

            # Determine # rows and columns from len(plotranges)
            numcols=min(len(plotranges),2) # 1 or 2 columns
            numrows=math.ceil(len(plotranges)/2)
            try:
                if len(plotranges)<7:                    
                    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
                else:
                    fig, axes = plt.subplots(nrows=numcols, ncols=numrows, figsize=(16,9), squeeze=False) # switch to 2 row style for >7 subplots
                colname='Peaks'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                for j, bounds in enumerate(plotranges):
                    [lower, upper]=bounds                    
                    thisrow=(j)//2 # correct zero based indexing
                    thiscol=(j)%2
                    if len(plotranges)<7:
                        axindex=thisrow,thiscol
                    else:
                        axindex=thiscol,thisrow # flips to 2 row, n column style
                    Augerslice=Augerfile[(Augerfile['Energy']>lower) & (Augerfile['Energy']<upper)] # already known that this isn't empty
                    Augerslice.plot(x='Energy', y=colname, ax=axes[axindex])
                    # Set tighter plot range for peaks only                        
                    Augerslice=Augerslice.dropna(subset=[colname]) # drop na values
                    myplotrange=[Augerslice['Energy'].min()*.95,Augerslice['Energy'].max()*1.05] # 
                    axes[axindex].set_xlim(myplotrange) 
                    # Section for labeling plotelements
                    elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=True) # can pass plot range as lower,upper tuple
                    axes[axindex].axhline(y=0)
                    # list of tuples with energy,elemname
                    for k, elemtuple in enumerate(elemlines):
                        # elemtuple[0] is energy and [1] is element symbol
                        try:
                            axes[axindex].axvline(x=elemtuple[0], color='b') # O line
                            axes[axindex].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                        except:
                            print('Problem labeling elements')                            
                    # skip evbreaks (or could get them from AugerParamLog if really necessary)
            except:
                print('Unknown problem with single element multiplex.')
            pdf.savefig(fig)
            plt.close(fig) # closes recently displayed figure (since it's saved in PDF report)
    plt.ion()
    return

def AESplot_tk(spelist, Elements, Smdifflog, Backfitlog, AESquantparams, **kwargs):
    ''' tk interface for args/kwargs of AESplot1 function 
    all args/dataframes must be passed through to plot functions 
    kwargs: fileset, areas, xrangestr -- usually just holds values entered during prior run  '''
    # first print out existing info in various lines
    root = tk.Tk()
    filestr=tk.StringVar() # comma separated or range of filenumbers for plot
    filestr.set(kwargs.get('fileset','')) # get entry from prior run
    areastr=tk.StringVar()  # area numbers for inclusion
    areastr.set(kwargs.get('areas','')) # use val from prior run if present
    xrangestr=tk.StringVar()  # energy range in eV 
    xrangestr.set(kwargs.get('xrangestr',''))
    plotcol=tk.StringVar() # column for plotting (Counts or deriv S7D7)
    smdifbool=tk.BooleanVar() # Bool for plotting of peak locations in smdifpeakslog
    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    plotelemstr=tk.StringVar()
    mytext='Peaks to be labeled:'+', '.join(Elements)
    plotelemstr.set(mytext)
    newelems=tk.StringVar() # string for new element choices if made
    newelems.set('')
    choice=tk.StringVar()  # plot or abort

    a=tk.Label(root, text='Enter filenumbers for plotting').grid(row=0, column=0)
    b=tk.Entry(root, textvariable=filestr).grid(row=0, column=1)
    a=tk.Label(root, text='Enter areanumbers for inclusion').grid(row=1, column=0)
    b=tk.Entry(root, textvariable=areastr).grid(row=1, column=1)
    a=tk.Label(root, text='Enter xrange in eV (default 100-1900eV)').grid(row=2, column=0)
    b=tk.Entry(root, textvariable=xrangestr).grid(row=2, column=1)
    
    radio1 = tk.Radiobutton(root, text='Counts', value='Counts', variable = plotcol).grid(row=3, column=0)
    radio1 = tk.Radiobutton(root, text='S7D7', value='S7D7', variable = plotcol).grid(row=4, column=0)
    
    c=tk.Checkbutton(root, variable=smdifbool, text='Plot smdiff peak locations?')
    c.grid(row=3, column=1)
    d=tk.Checkbutton(root, variable=backfitbool, text='Plot background fits?')
    d.grid(row=3, column=2)
    
    # option to reselect labeled elemental peaks 
    # TODO fix this nested tk interface ... chosen elements are not changing
    def changeelems(event):
        # pickelemsGUI imported from AESutils
        newelemlist=pickelemsGUI(AESquantparams) # get new elements/peaks list 
        newelems.set(', '.join(newelemlist))
        newtext='Peaks to be labeled: '+', '.join(newelemlist)
        plotelemstr.set(newtext)
        
    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    a=tk.Label(root, text=plotelemstr.get()).grid(row=5, column=0)
    d=tk.Button(root, text='Change labelled element peaks')
    d.bind('<Button-1>', changeelems)
    d.grid(row=6, column=0)

    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=6, column=2)

    d=tk.Button(root, text='Plot')
    d.bind('<Button-1>', plot)
    d.grid(row=6, column=1)

    root.mainloop()
        
    mychoice=choice.get()
    
    if mychoice=='plot':
        # Set up kwargs for plot 
        kwargs={}
        kwargs.update({'plotcol':plotcol.get()})        
        # TODO fix possible label element updating 
        if len(Elements)>0: # if not updated use passed list 
            kwargs.update({'plotelems':Elements})
        if xrangestr.get()!='' and '-' in xrangestr.get(): # should be hyphenated range if altered
            xmin=int(xrangestr.get().split('-')[0])
            xmax=int(xrangestr.get().split('-')[1])
            kwargs.update({'xrange':(xmin, xmax)})       
            kwargs.update({'xrangestr':xrangestr.get()})
        if areastr!='': # defaults to all areas...include in kwargs only if entry is made
            kwargs.update({'areas':areastr.get()})
        if smdifbool.get():
            kwargs.update({'smdifpeaks':Smdifflog}) 
        if backfitbool.get():
            kwargs.update({'backfitdf':Backfitlog})
        if filestr.get()!='': # also put into kwargs to initialize next run
            kwargs.update({'fileset':filestr.get()})
        print(kwargs)
        AESplot1(filestr.get(), spelist, AESquantparams, **kwargs)

    return kwargs

def countsbackreport_tk(spelist, Elements, Backfitlog, AESquantparams, **kwargs):
    ''' tk interface for args/kwargs of countsback plot reporting 
    all args/dataframes must be passed through to plot functions 
    Options: 
    1) filenums -- use all or choose subsets by filenumbers (spelist normally has all spe filenumbers in folder)
    2) xrange-- can be list of elements, list of ev ranges or combination (parsed by setplotboundaries)
    3) custom pdf name option
    4) plot background fits -bool
    TODO: subsetting of areas? 
    '''
    # TODO use prior kwargs to set defaults
    # first print out existing info in various lines
    root = tk.Tk()
    filestr=tk.StringVar() # comma separated or range of filenumbers for plot
    # Variable for entering eV range(s).. optionally used
    xrangestr=tk.StringVar()  # energy range in eV 
    xrangestr.set('')
    now=datetime.datetime.now()
    PDFname='Countsback_report'+'_'+now.strftime('%d%b%y')+'.pdf'
    PDFnamestr=tk.StringVar()  # energy range in eV 
    PDFnamestr.set(PDFname)
    # elements list only displayed for now... change with separate tk selector
    plotelemstr=tk.StringVar()
    mytext='Peaks to be labeled:'+', '.join(Elements)
    plotelemstr.set(mytext)
    rangechoice=tk.StringVar()

    # todo implement possible change of elements
    newelems=tk.StringVar() # string for new element choices if made
    newelems.set('')

    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    backfitptsbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    plotelemsbool=tk.BooleanVar()  # optional labelling of elements
    plotelemsbool.set(True) # default to true
    choice=tk.StringVar()  # plot or abort
    # set defaults based on prior run
    # set defaults based on prior run
    if 'addbackfit' in kwargs:
        backfitbool.set(True)
    if 'backfitpts' in kwargs:
        backfitptsbool.set(True)
    if 'plotelems' not in kwargs:
        backfitptsbool.set(False)
    a=tk.Label(root, text='Enter filenumbers for plot report (default all)').grid(row=0, column=0)
    b=tk.Entry(root, textvariable=filestr).grid(row=0, column=1)
    a=tk.Label(root, text='Optional xrange in eV (default 100-1900eV)').grid(row=1, column=0)
    b=tk.Entry(root, textvariable=xrangestr).grid(row=1, column=1)
    a=tk.Label(root, text='Enter PDF report name').grid(row=2, column=0)
    b=tk.Entry(root, textvariable=PDFnamestr).grid(row=2, column=1)
    # printout of peaks to be plotted (if chosen)
    a=tk.Label(root, text=plotelemstr.get()).grid(row=3, column=0)
    # Radio button for range (use element list, use ev strings or both)
    radio1 = tk.Radiobutton(root, text='Plot element ranges', value='elems', variable = rangechoice).grid(row=0, column=2)
    radio1 = tk.Radiobutton(root, text='Plot eV range', value='evrange', variable = rangechoice).grid(row=1, column=2)
    radio1 = tk.Radiobutton(root, text='Join elements and eV range', value='both', variable = rangechoice).grid(row=2, column=2)
    radio1 = tk.Radiobutton(root, text='Use full data range', value='full', variable = rangechoice).grid(row=3, column=2)
    d=tk.Checkbutton(root, variable=backfitbool, text='Plot background fits?')
    d.grid(row=4, column=0)
    d=tk.Checkbutton(root, variable=backfitptsbool, text='Plot points used to fit background?')
    d.grid(row=5, column=0)    
    d=tk.Checkbutton(root, variable=plotelemsbool, text='Label element peaks?')
    d.grid(row=6, column=0)    
    # option to reselect labeled elemental peaks 
    # TODO fix this nested tk interface ... chosen elements are not changing
    def changeelems(event):
        newelemlist=AESutils.pickelemsGUI(AESquantparams) # get new elements/peaks list 
        newelems.set(', '.join(newelemlist))
        newtext='Peaks to be labeled: '+', '.join(newelemlist)
        plotelemstr.set(newtext)
        
    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    d=tk.Button(root, text='Change labelled element peaks')
    d.bind('<Button-1>', changeelems)
    d.grid(row=7, column=0)

    d=tk.Button(root, text='Plot')
    d.bind('<Button-1>', plot)
    d.grid(row=7, column=1)
    
    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=7, column=2)

    root.mainloop()
        
    mychoice=choice.get()
    
    if mychoice=='plot':
        # Set up kwargs for plot 
        kwargs={}
        if filestr.get()!='': # optional subset of files based on entered number
            filenumlist=parsefilenums(filestr.get())
            kwargs.update({'filesubset':filenumlist})
        if backfitbool.get(): # plot backgrounds (stored in each Augerfile)
            kwargs.update({'addbackfit':True})
        if backfitptsbool.get(): # add points used for fitting background regions
            kwargs.update({'backfitpts':True})
        if plotelemsbool.get(): # add points used for fitting background regions
            kwargs.update({'plotelems':Elements})
        # Now create plotrange arg string (elems + evrange) based on radiobutton choice
        plotrange=[]
        if rangechoice.get()=='elems':
            plotrange.extend(Elements)
        elif rangechoice.get()=='evrange': # hyphenated ev range entered
            tempstr=xrangestr.get()
            plotrange.extend(tempstr.split(',')) # parse into strings if necessary (list must be passed)
        elif rangechoice.get()=='both':
            plotrange.extend(Elements) # start with elements list
            tempstr=xrangestr.get()
            plotrange.extend(tempstr.split(',')) # add ev ranges
        elif rangechoice.get()=='full': 
            # set to maximum possible range.. will be reduced if larger than data range
            plotrange.append('0-2500') # parse into strings if necessary (list must be passed)
        # get and pass PDF name string as arg (default name is prefilled)
        myPDFname=PDFnamestr.get()
        '''
        if newelems!='': # use newly assigned values
            kwargs.update({'plotelems':newelems.get()})       
        '''

        reportcountsback(spelist, plotrange, AESquantparams, Backfitlog, PDFname=myPDFname, **kwargs)
        
    return kwargs

def countsderivreport_tk(spelist, Elements, Smdifpeakslog, Backfitlog, AESquantparams, **kwargs):
    ''' tk interface for args/kwargs of countsderiv reports (deriv on top and counts below for all areas)
    all args/dataframes must be passed through to plot functions 
    Options: 
    1) filenums -- use all or choose subsets by filenumbers (spelist normally has all spe filenumbers in folder)
    2) xrange-- can be list of elements, list of ev ranges or combination (parsed by setplotboundaries)
    3) custom pdf name option
    4) plot background fits -bool
    5) plot smdifpeak amplitude points -bool (for deriv plot below)
    6) area subsets
    TODO: subsetting of areas? 
    '''
    # TODO use prior kwargs to set defaults
    # first print out existing info in various lines
    root = tk.Tk()
    filestr=tk.StringVar() # comma separated or range of filenumbers for plot
    # Variable for entering eV range(s).. optionally used
    xrangestr=tk.StringVar()  # energy range in eV 
    xrangestr.set('')
    now=datetime.datetime.now()
    PDFname='Countsderiv_report'+'_'+now.strftime('%d%b%y')+'.pdf'
    PDFnamestr=tk.StringVar()  # energy range in eV 
    PDFnamestr.set(PDFname)
    # elements list only displayed for now... change with separate tk selector
    plotelemstr=tk.StringVar()
    mytext='Peaks to be labeled:'+', '.join(Elements)
    plotelemstr.set(mytext)
    rangechoice=tk.StringVar()

    # todo implement possible change of elements
    newelems=tk.StringVar() # string for new element choices if made
    newelems.set('')

    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    backfitptsbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    smdifbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    plotelemsbool=tk.BooleanVar()  # optional labelling of elements
    plotelemsbool.set(True) # default to true
    choice=tk.StringVar()  # plot or abort
    # set defaults based on prior run
    if 'addbackfit' in kwargs:
        backfitbool.set(True)
    if 'backfitpts' in kwargs:
        backfitptsbool.set(True)
    if 'plotelems' not in kwargs:
        backfitptsbool.set(False)
    a=tk.Label(root, text='Enter filenumbers for plot report (default all)').grid(row=0, column=0)
    b=tk.Entry(root, textvariable=filestr).grid(row=0, column=1)
    a=tk.Label(root, text='Optional xrange in eV (default 100-1900eV)').grid(row=1, column=0)
    b=tk.Entry(root, textvariable=xrangestr).grid(row=1, column=1)
    a=tk.Label(root, text='Enter PDF report name').grid(row=2, column=0)
    b=tk.Entry(root, textvariable=PDFnamestr).grid(row=2, column=1)
    # printout of peaks to be plotted (if chosen)
    a=tk.Label(root, text=plotelemstr.get()).grid(row=3, column=0)
    # Radio button for range (use element list, use ev strings or both)
    radio1 = tk.Radiobutton(root, text='Plot element ranges', value='elems', variable = rangechoice).grid(row=0, column=2)
    radio1 = tk.Radiobutton(root, text='Plot eV range', value='evrange', variable = rangechoice).grid(row=1, column=2)
    radio1 = tk.Radiobutton(root, text='Join elements and eV range', value='both', variable = rangechoice).grid(row=2, column=2)
    radio1 = tk.Radiobutton(root, text='Use full data range', value='full', variable = rangechoice).grid(row=3, column=2)
    d=tk.Checkbutton(root, variable=backfitbool, text='Plot background fits?')
    d.grid(row=4, column=0)
    d=tk.Checkbutton(root, variable=backfitptsbool, text='Plot points used to fit background?')
    d.grid(row=5, column=0)    
    d=tk.Checkbutton(root, variable=plotelemsbool, text='Label element peaks?')
    d.grid(row=6, column=0)
    d=tk.Checkbutton(root, variable=smdifbool, text='Label amplitude pts in sm-diff?')
    d.grid(row=7, column=0)
    # option to reselect labeled elemental peaks 
    # TODO fix this nested tk interface ... chosen elements are not changing
    def changeelems(event):
        newelemlist=AESutils.pickelemsGUI(AESquantparams) # get new elements/peaks list 
        newelems.set(', '.join(newelemlist))
        newtext='Peaks to be labeled: '+', '.join(newelemlist)
        plotelemstr.set(newtext)
        
    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    d=tk.Button(root, text='Change labelled element peaks')
    d.bind('<Button-1>', changeelems)
    d.grid(row=8, column=0)

    d=tk.Button(root, text='Plot')
    d.bind('<Button-1>', plot)
    d.grid(row=8, column=1)
    
    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=8, column=2)

    root.mainloop()
        
    mychoice=choice.get()
    
    if mychoice=='plot':
        # Set up kwargs for plot 
        kwargs={}
        if filestr.get()!='': # optional subset of files based on entered number
            filenumlist=parsefilenums(filestr.get())
            kwargs.update({'filesubset':filenumlist})
        if backfitbool.get(): # plot backgrounds (stored in each Augerfile)
            kwargs.update({'addbackfit':True})
        if backfitptsbool.get(): # add points used for fitting background regions
            kwargs.update({'backfitpts':True})
        if plotelemsbool.get(): # add points used for fitting background regions
            kwargs.update({'plotelems':Elements})
        if smdifbool.get(): # add points used for fitting background regions
            kwargs.update({'smdifpts':True})
        # Now create plotrange arg string (elems + evrange) based on radiobutton choice
        plotrange=[]
        if rangechoice.get()=='elems':
            plotrange.extend(Elements)
        elif rangechoice.get()=='evrange':
            tempstr=xrangestr.get()
            plotrange.extend(tempstr.split(',')) # parse into strings if necessary (list must be passed)
        elif rangechoice.get()=='both':
            plotrange.extend(Elements) # start with elements list
            tempstr=xrangestr.get()
            plotrange.extend(tempstr.split(',')) # add ev ranges
        elif rangechoice.get()=='full': 
            # set to maximum possible range.. will be reduced if larger than data range
            plotrange.append('0-2500') # parse into strings if necessary (list must be passed)
        # get and pass PDF name string as arg (default name is prefilled)
        myPDFname=PDFnamestr.get()
        '''
        if newelems!='': # use newly assigned values
            kwargs.update({'plotelems':newelems.get()})       
        '''
        reportderivcnt(spelist, plotrange, AESquantparams, Backfitlog, Smdifpeakslog, PDFname=myPDFname, **kwargs)
    return kwargs

def scattercompplot_tk(comp1, comp2, elemlist, **kwargs):
    ''' tk interface for args/kwargs of scattercompplot function (compares compositions derived two different ways)
    all args/dataframes must be passed through to plot functions 
    Keyword args:
    joinlist: list of columns to use for merge of two composition files
        filenumber, areanumber is default and used to compare compositions from exact same spe files 
        but computed via sm-diff or integ methods 
        sample - returns all sample inner merge matches (some of which could be ), such as if a 
        background region were also measured as one of the areas in an spe file
    thresh - always used; distinguishes outlier points in scatter plot; defaults to 0.1
          higher gives more outliers which are returned and plotted separately
    basis - bool (true means plot basis, false is defaul at.% plot)
    errbars - optional plotting of x, y or xy error bars'''
    # first print out existing info in various lines
    root = tk.Tk()
    thresh=tk.StringVar()  # threshold for outliers returned (string not doublevar)
    thresh.set(0.1)
    errbars=tk.StringVar() # x, y, xy or none
    jointype=tk.StringVar()
    plottype=tk.StringVar() # for basis or at.% comparison
    # radiobuttons for how to join separate datasets

    choice=tk.StringVar()  # plot or abort
    mytext='Elements for scatter plot comparison: '+', '.join(plotelems)
    a=tk.Label(root, text=mytext).grid(row=0, column=0)
    a=tk.Label(root, text='Enter threshold for outliers (higher yields more)').grid(row=1, column=0)
    b=tk.Entry(root, textvariable=thresh).grid(row=1, column=1)
    a=tk.Label(root, text='Errors for plots (x, y, xy or none)').grid(row=2, column=0)
    b=tk.Entry(root, textvariable=errbars).grid(row=2, column=1)

    # radiobuttons for how to join separate datasets
    radio1 = tk.Radiobutton(root, text='Join on filenumber & areanumber', value='fa', variable = jointype).grid(row=3, column=0)
    radio1 = tk.Radiobutton(root, text='Join on sample', value='sample', variable = jointype).grid(row=4, column=0)
    # basis or at percent radiobutton
    radio1 = tk.Radiobutton(root, text='Compare at. %', value='atper', variable = plottype).grid(row=3, column=1)
    radio1 = tk.Radiobutton(root, text='Compare basis', value='basis', variable = plottype).grid(row=4, column=1)

    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    d=tk.Button(root, text='Plot')
    d.bind('<Button-1>', plot)
    d.grid(row=5, column=0)

    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=5, column=1)

    root.mainloop()
        
    mychoice=choice.get()
    if mychoice=='plot':
        # Set up kwargs for scatter comp plot 
        kwargs={}
        if errbars.get()!='': # use newly assigned values
            kwargs.update({'errbars':errbars.get()}) 
        try:
            kwargs.update({'thresh':float(thresh.get())})
        except:
            print('Entered threshold must be a float')
            kwargs.update({'thresh':0.1})
        # plot type
        if plottype.get()=='basis':
            kwargs.update({'basis':True}) # use metals basis
        else:
            kwargs.update({'basis':False}) # use at. 5
        # data join type
        if plottype.get()=='fa':
            kwargs.update({'joinlist':['Filenumber','Areanumber']}) # default data join technique (most strict)
        elif plottype.get()=='sample':
            kwargs.update({'joinlist':['Sample']}) # default data join technique (most strict)
        scattercompplot(comp1, comp2, elemlist, **kwargs)
    return kwargs
    
# Possible legacy functions
def cloneparamrows(df):
    ''' Make param log entry for for each areanum - used by calccomposition to correctly process spe files with multiple spatial areas
    passed df is usually list of spe files
    this solves problem that AugerParamLog only has one entry (despite possibly having multiple distinct areas with different spectra'''
    df['Areanumber']=1 # set existing entries as area 1
    mycols=df.dtypes.index
    newrows=pd.DataFrame(columns=mycols)
    for index, row in df.iterrows():
        numareas=int(df.loc[index]['Areas'])
        for i in range(2,numareas+1):
            newrow=df.loc[index] # clone this row as series
            newrow=newrow.set_value('Areanumber',i)
            newrows=newrows.append(newrow)
    df=pd.concat([df,newrows], ignore_index=True) # merge new rows with existing ones
    df=df.sort_values(['Filenumber','Areanumber'])
    mycols=['Filenumber','Areanumber','Filename','Evbreaks']
    df=df[mycols]
    return df
  
def uniquify(mylist):
    ''' Removes exact duplicates from nested lists '''
    tuplelist=[tuple(l) for l in mylist]
    seen = set()
    seen_add = seen.add
    uniquelist=[x for x in tuplelist if not (x in seen or seen_add(x))]
    uniquelist=[list(l) for l in uniquelist]
    return uniquelist