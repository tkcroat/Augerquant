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
import re
from matplotlib.backends.backend_pdf import PdfPages
import scipy
import scipy.stats # load in this sequence to get linregress working
import math
from statsmodels.formula.api import ols # ordinary least squares
import ternary # python-ternary plotting module

def uniquify(mylist):
    tuplelist=[tuple(l) for l in mylist]
    seen = set()
    seen_add = seen.add
    uniquelist=[x for x in tuplelist if not (x in seen or seen_add(x))]
    uniquelist=[list(l) for l in uniquelist]
    return uniquelist

def cloneparamrows(df):
    ''' Make param log entry for for each areanum - used by calccomposition to correctly process spe files with multiple spatial areas
    passed df is usually list of spe files
    this solves problem that AugerParamLog only has one entry (despite possibly having multiple distinct areas with different spectra'''
    df['Areanumber']=1 # set existing entries as area 1
    mycols=df.dtypes.index
    newrows=pd.DataFrame(columns=mycols) # blank df for new entries
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
    
def getplotboundaries(Augerfile, plotelems, AESquantparams, colname='Energy'):
    ''' Gets typical boundary of plots of given line from quantparams from plotelems but remove duplicates
    or remove if no data in range
    do not return range for duplicates (ie. Fe2 is in Fe plot range, Ca in C plot range so return 0s 
    defaults to checking for vals in energy column but can check peaks or other cols
    can also just pass string with ev range '''
    plotranges=[] # returns list of length 2 lists for valid elements
    
    for i, elem in enumerate(plotelems):
        if '-' in elem: # direct ev range specified in length 1 list (not list of elements)
            plotranges.append([int(elem.split('-')[0]),int(elem.split('-')[1])])
            return plotranges
        thiselemdata=AESquantparams[(AESquantparams['element']==elem)]
        if len(thiselemdata)==1:
            thisrange=thiselemdata.iloc[0]['plotrange']
            try:
                match= re.finditer(r'\d+', thisrange)
                if match:
                    thisrange=[m.group(0) for m in match] # parse range into lower upper
                    thisrange=[int(i) for i in thisrange]
                    Augerslice=Augerfile[(Augerfile[colname]>thisrange[0]) & (Augerfile[colname]<thisrange[1])]
                    if not Augerslice.empty:
                        plotranges.append(thisrange)
            except:
                pass # combining multiple lines in single window isn't an error
                # print('Duplicate plot range for ', elem)
        # knock out duplicated plot ranges (i.e. just plot Fe and Fe2 together)
        # Duplicate lines to knock out are set to identical plotranges in AESquantparams        
        else:
            print ('Problem finding plot range for element ', elem)
    # now knock out any duplicate plot ranges, such as C and Ca, Fe and Fe2
    # this is controllable via plotrange specified in AESquantparams
    plotranges=uniquify(plotranges) # knocks out duplicated ranges in my nested lists
    return plotranges
    
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
            # typically tokens are plain old integers
            if int(i)<=numareas: # skips if value inadvertently exceeds # of areas
                myareas.add(int(i))
        except: # if not, then it might be a range
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # we have items seperated by a dash
                    # try to build a valid range
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
    return myareas

def getelemenergy(plotelems, plotrange, AESquantparams, deriv=False):
    ''' Pass plotted data range as tuple or string and list of desired elements for plotting, return paired list of elemental and ideal energies 
    only returns element if in range of the given plot; can return either direct peak ideal position or deriv 
    peak ideal position (using deriv switch)
    negpeak is energy of peak deriv peak whereas direct peak is negpeak- integpeak (expressed as eV shift'''
    # elems=[str(s) for s in plotelems.split(' ')]
    elemlines=[] # energies list
    # Need to return empty list of plotelems only contains an ev range 
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
        tempstring=tempstring.split(']')[0] # remove brackets
        evbreaks=[int(s) for s in tempstring.split(',')] # turn string into list of ints 
    # find energy x values from known evbreak indices (Augerslice.index matches original index #s)
    energyvals=[] # new list for corresponding energy x vals within range
    for i in range(0,len(Augerfile)):
        index=Augerfile.index[i]
        if index in evbreaks:
            energyvals.append(Augerfile.iloc[i]['Energy']) # find energy vals assoc with index # 
    return energyvals
    
def setplotrange(plotrange, Augerfile):
    ''' Set range of plot based on element, numerical range or default to max range of spectrum  
    commonly called by Auger plot functions below'''
    if '-' in plotrange: # determine range for plot
        myplotrange=(int(plotrange.split('-')[0]),int(plotrange.split('-')[1]))
    elif plotrange=='C':
        myplotrange=(236,316)  
    elif plotrange=='Ca':
        myplotrange=(236,336)   
    elif plotrange=='O':
        myplotrange=(470,540)
    elif plotrange=='Fe':
        myplotrange=(560,747)  
    elif plotrange=='Mg':
        myplotrange=(1145,1225)
    elif plotrange=='Al':
        myplotrange=(1350,1430)
    elif plotrange=='Si':
        myplotrange=(1570,1650)
    else: # defaults to full data range
        lower=Augerfile.Energy.min()        
        upper=Augerfile.Energy.max()
        myplotrange=(lower, upper)
    return myplotrange

def findelemlines(plotelems,myplotrange):
    ''' Return energies of lines mentioned in plotelems string from main, return only if within plot rage (range passed as tuple)''' 
    plotelems=plotelems.replace('Fe','Fe Fe1 Fe2') # choose all Fe lines
    elems=[str(s) for s in plotelems.split(' ')]
    elemlines=[] # energies list    
    AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Params\\AESquantparams.csv', encoding='utf-8')
    AESquantparams=AESquantparams[AESquantparams['element'].isin(elems)] # select rows in element list
    AESquantparams=AESquantparams[(AESquantparams['negpeak']>myplotrange[0]) & (AESquantparams['negpeak']<myplotrange[1])]
    elemlines=AESquantparams['negpeak'].tolist()
    return elemlines

def maketitlestring(plotpts):
    ''' extract label for plot containing element peak amplitude, noise amplitudes (at lower and higher eV), and significance (elem amplitude/noise amplitude)
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

def get_plotkwargs(mykwargs, elem):
    ''' gets subset of keyword arguments that are needed as plot arguments
    pass elem to construct x and y err col names which may be needed 
    
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
    #mycols=mycols.append('Element')
    # Prepare dataframes for full returned datasets
    outliers=pd.DataFrame(columns=mycols) # empty dataframe for outlying points
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
        kwarg1=get_plotkwargs(kwargs, elem)
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
    # TODO fix error Nonetype has no attribute is_bbox
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
        
    # create plot
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
    
def AESplot1(Params, Peaks,plotrange, plotelems=''): # single file already selected by number
    '''Plot a single file over desired range (after dataframe slicing via filenumber 
    plotelems kwarg is list of elements who's ideal position will be shown as blue vert line on plot'''
    AugerFileName=Params.Filename # filename if logmatch is series
    numareas=int(Params.Areas)
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    myplotrange=setplotrange(plotrange, Augerfile) # passes ev range or element or defaults to max range (returns tuple)
    
    # Things added to plots (breaks in multiplex scans, points used for S7D7 amplitudes)
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
    # add elemental peaks as vert lines
    if plotelems!='':
        elemlines = findelemlines(plotelems,myplotrange) # find energy vals of peaks within plotted range
           
    cols=divmod(numareas,2)[0]+ divmod(numareas,2)[1]
    if numareas>1:
        rows=2
    else:
        rows=1        
    fig, axes = plt.subplots(nrows=rows, ncols=cols) # axes is array
    
    # slice to desired energy plotrange 
    Augerslice=Augerfile[(Augerfile['Energy']>myplotrange[0]) & (Augerfile['Energy']<myplotrange[1])]
    
    for i in range(0,numareas): # myareas is a set and can be non-continuous
        if (i+1)%2==1:
            rownum=0
        else:
            rownum=1
        colnum=int((i+1)/2.1)
        colname='S7D7'+str(i+1)
        plotpts=Peaks[(Peaks['Areanumber']==(i+1))] # pull peak max/min points from only this spatial area
        if numareas==1:
            Augerslice.plot(x='Energy', y=colname, ax=axes) # single plot axes has no [#,#]
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes, color='r') # plotting of max/min pts for smdif quant
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes, color='r') 
            for i, val in enumerate(energyvals):
                axes.axvline(x=val, color='r') 
            for i, val in enumerate(elemlines):
                axes.axvline(x=val, color='b') 
        else: # case for >1 plots axes array
            Augerslice.plot(x='Energy', y=colname, ax=axes[rownum,colnum]) # single plot axes has no [#,#]            
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[rownum,colnum], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[rownum,colnum], color='r') 
            # add evbreaks as red vertical lines
            for i, val in enumerate(energyvals):
                axes[rownum,colnum].axvline(x=val, color='r') 
            for i, val in enumerate(elemlines):
                axes[rownum,colnum].axvline(x=val, color='b') 
    return

def reportderivcnt(paramlog, plotelems, AESquantparams, Smdifdf=False, backfitdf=False, PDFname='SDplots_report.pdf'):
    ''' Comparison plots for both derivative and counts itself (don't have to worry about axes w/o indexing)
    plots selected filenumber + areanumber
    plots all spe files, associated quant points, and labels selected elemental lines
    both backfitlog and smdiflogs are optional params
    '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            areanum=int(paramlog.loc[index]['Areanumber'])  
            try:
                Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            except:
                print(AugerFileName,' skipped ... not found.')
                continue
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.loc[index] # grab row for this spe file as Series            
            if type(Smdifdf)==pd.core.frame.DataFrame:
                thisfilepeaks=Smdifdf[(Smdifdf['Filename']==AugerFileName)&(Smdifdf['Areanumber']==areanum)] # retrieve assoc. subset of peaks data
                addsmdif=True # flag to include smooth-diff points 
            # SKIP energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # determine size of plot for this filenumber (same for all areas); plotranges combines some plotelems (i.e. C and Ca together)
            plotranges=getplotboundaries(Augerfile, plotelems, AESquantparams) # returns plot ranges for all regions with data from plotelems
            # set plot rows and columns
            numrows=2
            numcols=len(plotranges)
            if type(backfitdf)==pd.core.frame.DataFrame: # for plotting background fit points integral method
                thisfilebackpts=backfitdf[(backfitdf['Filename']==AugerFileName) & (backfitdf['Areanumber']==areanum)]
                plotbackpts=True
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
                indexptslist=[int(i) for i in indexptslist] # make sure all index #s are ints 
                indexptslist.sort()
            else:
                plotbackpts=False
            try:
                fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
                mytitle=AugerFileName.replace('.csv','') +' area #'+str(areanum)                
                plt.suptitle(mytitle)
                for j, bounds in enumerate(plotranges):
                    [lower, upper]=bounds                    
                    thiscol=j
                    Augerslice=Augerfile[(Augerfile['Energy']>lower) & (Augerfile['Energy']<upper)] # already known that this isn't empty
                    Augerslice.plot(x='Energy', y='S7D7'+str(areanum), ax=axes[0,thiscol]) # deriv in upper plot
                    Augerslice.plot(x='Energy', y='Counts'+str(areanum), ax=axes[1,thiscol]) # counts in lower 
                    # Section for labeling plotelements
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
                    # Section for adding smooth-diff quant data
                    if addsmdif:        
                        plotpts=thisfilepeaks[(thisfilepeaks['Peakenergy']>lower) & (thisfilepeaks['Peakenergy']<upper)]
                        if not plotpts.empty:
                            try:
                                plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,thiscol], color='r')
                                plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,thiscol], color='r')                                        
                                titlestring=maketitlestring(plotpts)
                                axes[0,thiscol].set_title(titlestring, fontsize=10)
                            except:
                                print('Problem adding points from smdif quant calcs for ', AugerFileName,'area # ', areanum )
                    # add red vert line at multiplex energy break if present
                    # removed... however evbreaks could be retrieved from AugerParamLog if desired
                    '''
                    for l, val in enumerate(energyvals):
                        if val > lower and val < upper: 
                            axes[0,thiscol].axvline(x=val, color='r') # on deriv plot
                            axes[1,thiscol].axvline(x=val, color='r') # on counts plot 
                    '''
                    # Now plot counts and background fits in bottom row
                    if plotbackpts==True: # flag to show points from which background was determined 
                        # Now add scatter plot points at fit region boundaries
                        backpts=Augerslice[Augerslice.index.isin(indexptslist)] # gets background fitted pts but only from this data slice
                        if not backpts.empty: # show fitted pts from counts
                            backpts.plot.scatter(x='Energy', y='Counts'+str(areanum), ax=axes[1,thiscol])
                        Augerslice.plot(x='Energy', y='Backfit'+str(areanum), ax=axes[1,thiscol]) 
                    # now label elements for counts plot
                    elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=False) # can pass plot range as lower,upper tuple
                    # list of tuples with energy,elemname
                    for k, elemtuple in enumerate(elemlines):
                        # elemtuple[0] is energy and [1] is element symbol
                        # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                        try:
                            axes[1,thiscol].axvline(x=elemtuple[0], color='b') # O line
                            yval=(Augerslice['Counts'+str(areanum)].max()-Augerslice['Counts'+str(areanum)].min())*0.9+Augerslice['Counts'+str(areanum)].min()
                            axes[1,thiscol].text(elemtuple[0],yval, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                        except:
                            print('Problem labeling elements') 
                pdf.savefig(fig)
                plt.close(fig)
            except:
                print('Unknown problem plotting', AugerFileName,' area #', areanum)
    plt.ion()
    return

def reportderivcntall(paramlog,  plotelems, AESquantparams, Smdifdf=False, backfitdf=False, PDFname='this_report.pdf'):
    ''' Comparison plots for both derivative and counts itself (don't have to worry about axes w/o indexing)
    plots selected filenumber for all area (looped)
    plots all spe files, associated quant points, and labels selected elemental lines
    '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])  
            try:
                Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            except:
                print(AugerFileName,' skipped ... not found.')
                continue
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            # Params=paramlog.loc[index] # grab row for this spe file as Series
            if type(Smdifdf)==pd.core.frame.DataFrame:
                addsmdif=True # flag to plot sm-diff quant points 
                thisfilepeaks=Smdifdf[Smdifdf['Filename']==AugerFileName] # retrieve assoc. subset of peaks data
            # SKIP energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # determine size of plot for this filenumber (same for all areas); plotranges combines some plotelems (i.e. C and Ca together)
            plotranges=getplotboundaries(Augerfile, plotelems, AESquantparams) # returns plot ranges for all regions with data from plotelems
            for i in range(0,numareas):
                areanum=i+1
                
                # set plot rows and columns
                numrows=2
                numcols=len(plotranges)
                if type(backfitdf)==pd.core.frame.DataFrame: # for plotting background fit points integral method
                    thisfilebackpts=backfitdf[(backfitdf['Filename']==AugerFileName) & (backfitdf['Areanumber']==areanum)]
                    plotbackpts=True
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
                    indexptslist=[int(i) for i in indexptslist] # make sure all index #s are ints 
                    indexptslist.sort()
                else:
                    plotbackpts=False
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
                        # Section for labeling plotelements
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
                        # Section for adding smooth-diff quant data
                        if addsmdif:
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
                        # removed... however evbreaks could be retrieved from AugerParamLog if desired
                        '''
                        for l, val in enumerate(energyvals):
                            if val > lower and val < upper: 
                                axes[0,thiscol].axvline(x=val, color='r') # on deriv plot
                                axes[1,thiscol].axvline(x=val, color='r') # on counts plot 
                        '''
                        # Now plot counts and background fits in bottom row
                        if plotbackpts==True:
                            # Now add scatter plot points at fit region boundaries
                            backpts=Augerslice[Augerslice.index.isin(indexptslist)] # gets background fitted pts but only from this data slice
                            if not backpts.empty: # show fitted pts from counts
                                backpts.plot.scatter(x='Energy', y='Counts'+str(areanum), ax=axes[1,thiscol])
                            Augerslice.plot(x='Energy', y='Backfit'+str(areanum), ax=axes[1,thiscol]) 
                        # now label elements for counts plot
                        elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=False) # can pass plot range as lower,upper tuple
                        # list of tuples with energy,elemname
                        for k, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                            try:
                                axes[1,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                yval=(Augerslice['Counts'+str(areanum)].max()-Augerslice['Counts'+str(areanum)].min())*0.9+Augerslice['Counts'+str(areanum)].min()
                                axes[1,thiscol].text(elemtuple[0],yval, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                            except:
                                print('Problem labeling elements')
                    # now hide empty subplots
                    for i in range(0,numrows*numcols):
                        if i>len(plotranges)-1:
                            thisrow=i//numcols
                            thiscol=i%numcols
                            axindex=thisrow, thiscol # tuple to index axes 
                            axes[axindex].set_visible(False)         
                    pdf.savefig(fig)
                    plt.close(fig)
                except:
                    print('Unknown problem plotting', AugerFileName,' area #', areanum)
    plt.ion()
    return

def reportSD(paramlog, Smdifpeakslog, plotelems, AESquantparams, PDFname='SDplots_report.pdf'):
    ''' generalized newer version of SDmajor, plots all spe files, associated quant points, and labels selected elemental lines  
    '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])  
            try:
                Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            except:
                print(AugerFileName,' skipped... not found.')
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

def reportcountsback(paramlog, plotelems, AESquantparams, backfitdf=False, PDFname='countsback_report.pdf'):
    ''' Plot of list of passed elements 
    pass list of files and selected background regions from automated fitting
    background fits themselves stored with auger csv files
    optional pass of backfitlog (w/ points defining region boundary for background fitting useful for troubleshooting fits)
    evbreaks from multiplex plotted as red lines (often different dwell times for different elements)
    plotback switch -- whether or not to plot 
    '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])            
            try:
                Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            except:
                print(AugerFileName,' skipped ... not found.')
                continue
            Params=paramlog.loc[index] # grab row for this spe file as Series
            # filenumber=Params.Filenumber # retrieve filenumber            
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # determine size of plot for this filenumber (same for all areas)
            plotranges=getplotboundaries(Augerfile, plotelems, AESquantparams) # returns plot ranges for all regions with data from plotelems
            # boundaries of backfit range from backfitlog are helpful for plotting(lower1 & 2 and upper 1&2 which are index #s) 
            if type(backfitdf)==pd.core.frame.DataFrame:
                thisfilebackpts=backfitdf[backfitdf['Filename']==AugerFileName]
                plotbackpts=True
            else:
                plotbackpts=False
            for i in range(0,numareas): # create separate plot page for each area 
                areanum=i+1
                if plotbackpts==True: # this gets all the lower1, lower2, upper1, upper2 index point boundaries
                    indexptslist=[]
                    thisareabackpts=thisfilebackpts[(thisfilebackpts['Areanumber']==areanum)] # get subset of peaks for this area number 
                    thisarr=thisareabackpts.Lower1.unique()
                    thislist=np.ndarray.tolist(thisarr)
                    indexptslist.extend(thislist)
                    thisarr=thisareabackpts.Lower2.unique()
                    thislist=np.ndarray.tolist(thisarr)
                    indexptslist.extend(thislist)
                    thisarr=thisareabackpts.Upper1.unique()
                    thislist=np.ndarray.tolist(thisarr)
                    indexptslist.extend(thislist)
                    thisarr=thisareabackpts.Upper2.unique()
                    thislist=np.ndarray.tolist(thisarr)
                    indexptslist.extend(thislist)
                    indexptslist=[int(i) for i in indexptslist] # make sure all index #s are ints 
                    indexptslist.sort()
                # set plot row and column for this element range (from plotelems -- plotranges)
                if len(plotranges)==1:
                    numcols=1
                    numrows=1
                else:
                    numcols=2 # 
                    numrows=math.ceil(len(plotranges)/2)
                # new plot for each spatial area
                try:
                    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by 3 axes array
                    colname='Counts'+str(areanum)
                    mytitle=AugerFileName +' area #'+str(areanum)
                    plt.suptitle(mytitle)
                    # now loop over the elemental plot ranges
                    for j, bounds in enumerate(plotranges):
                        [lower, upper]=bounds                    
                        thisrow=j%numrows
                        thiscol=j//numrows
                        axindex=thisrow, thiscol
                        Augerslice=Augerfile[(Augerfile['Energy']>=lower) & (Augerfile['Energy']<=upper)] # already known that this isn't empty
                        Augerslice.plot(x='Energy', y=colname, ax=axes[axindex]) # plot counts
                        if plotbackpts==True:
                            # Now add scatter plot points at fit region boundaries
                            backpts=Augerslice[Augerslice.index.isin(indexptslist)] # gets background fitted pts but only from this data slice
                            if not backpts.empty: # show fitted pts from counts
                                backpts.plot.scatter(x='Energy', y=colname, ax=axes[axindex])
                            backname='Backfit'+str(areanum) # don't need separate empty test since df already checked for empty in getplotboundaries
                            Augerslice.plot(x='Energy', y=backname, ax=axes[thisrow,thiscol]) 
                            
                        # Section for labeling plotelements
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
                        # add red vertical lines at multiplex energy breaks 
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

def reportpeaksall(spelist, plotelems, AESquantparams, PDFname='Peaks_report.pdf'):
    ''' Generalized version of reportpeaksmajor working from spelist (loop over spatial areas 
    plots subtracted data for all spe files
    '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in spelist.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=spelist.loc[index]['Filename']
            numareas=int(spelist.loc[index]['Areas'])  
            try:
                Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            except:
                print(AugerFileName,' skipped ... not found.')
                continue
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=spelist.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
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

def reportpeaks(dfselect, plotelems, AESquantparams, PDFname='Peaks_report.pdf'):
    ''' Peaks report for subset of files (i.e. outliers) and only plots listed filenum and areanum (no all areas loop
    plots subtracted data for all spe files
    '''
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        for index,row in dfselect.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=dfselect.loc[index]['Filename']
            areanum=dfselect.loc[index]['Areanumber']  
            try:
                Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            except:
                print(AugerFileName,' skipped ... not found.')
                continue
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=dfselect.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # determine size of plot for this filenumber (same for all areas)
            plotranges=getplotboundaries(Augerfile, plotelems, AESquantparams, colname='Peaks1') # returns plot ranges for data-containing regions
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

