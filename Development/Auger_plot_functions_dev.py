# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 12:07:56 2016

@author: tkc
"""

#%% Load modules
import matplotlib.pyplot as plt
import re, os, glob, sys, csv # already run with functions 
import math
import pandas as pd
import numpy as np
import scipy
import scipy.stats
from matplotlib.backends.backend_pdf import PdfPages # Only needed for plotting
from collections import defaultdict
from math import factorial # used by Savgol matrix
from itertools import cycle

import ternary # python-ternary module
#%% Auger plotting functions under development 

kwargs.update({'groupby':['Filenumber','Areanumber']})
Smdifpeakslog.columns
grouped=Smdifpeakslog.groupby(['Filenumber','Areanumber'])

grouplists=[]
for i, col in enumerate(groupcols):
    grouplists[i].append()
groupcols=kwargs.get('groupby',[])

flist=np.ndarray.tolist(Smdifpeakslog['Filenumber'].unique())
alist=np.ndarray.tolist(Smdifpeakslog['Areanumber'].unique())


for key, group in grouped:
    
    print(key[0],key[1])

list(dict(list(grouped)))

# FUNCTION TO MODIFY LEGEND/TITLE
# tk Interactive function to change legends on active AES plot
L=plt.legend()
L.get_texts()[0].set_text('10 nA')
L.get_texts()[1].set_text('5 nA')
L.get_texts()[2].set_text('1 nA')

def changeplot_tk():
    ''' Alter legend and/or change title '''
    L=plt.legend()
    L.get_texts()[0].set_text('10 nA')
    L.get_texts()[1].set_text('5 nA')
    L.get_texts()[2].set_text('1 nA')
    L.get_texts()[0].set_text('1 nA')
    L.get_texts().set_text('1 nA')
    L.set_text('1 nA')
    type(L.get_texts())
    plt.suptitle('FIB area 8 (no Ar sputtering)')
    

def reportSD(spelist, plotrange, Smdifpeakslog, AESquantparams, PDFname='SDplot_report', **kwargs):
    ''' generalized newer version of SDmajor, plots all spe files, associated quant points, and 
    labels selected elemental lines
    PDFname -  altered/passed by tk batch PDf reporting interface
    
    kwargs:
        plotelems -  for element/peak labelling (optional and separate from plotrange)
        smdifpeaks - bool... mark points used for amplitude calculation (negpeak and pospeak)
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
            Params=spelist.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            thisfilepeaks=Smdifpeakslog[(Smdifpeakslog['Filename']==AugerFileName)] # retrieve assoc. subset of peaks data
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            # determine size of plot for this filenumber (same for all areas)
            plotranges=getplotboundaries(Augerfile, plotrange, AESquantparams) # returns plot ranges for all regions with data from plotrange

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
                        # Section for labeling plotelements (optional)
                        if 'plotelems' in kwargs:
                            plotelems=kwargs.get('plotelems',[]) # list of element peaks for labeling
                            elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=True) # can pass plot range as lower,upper tuple
                            # list of tuples with energy,elemname
                            for k, elemtuple in enumerate(elemlines):
                                # elemtuple[0] is energy and [1] is element symbol
                                # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                try:
                                    axes[axindex].axvline(x=elemtuple[0], color='b') # O line
                                    yval=(Augerslice[colname].max()-Augerslice[colname].min())*0.9+Augerslice[colname].min()
                                    axes[axindex].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18)
                                except:
                                    print('Problem labeling elements')
                        # Section for adding smooth-diff quant data (optional)
                        if 'smdifpeaks' in kwargs:
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
                                                # Section for labeling plotelements
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


def findgroupnum():
    ''' Figure out color and or marker for this '''

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
        groupby -- columns for plot grouping
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
    # convert passed areas string 
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
    # add elemental peaks as vert lines
    plotelems=kwargs.get('plotelems',[])
    if len(plotelems)>0:
        # find energy vals of peaks within plotted range
        elemlines = findelemlines(plotelems, xrange, AESquantparams) 
           
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
    Handling of color, line or symbol type cycling using groupby kwarg
    first grouped col uses default color cycling
    # TODO does colorlist represent different areas or different filenumbers?
    # Default is to use color for areas (different sample regions) as filenumbers are often 
    temporal reruns (time slices of same underlying spectrum possibly w/ eneryg shifts)

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
    if len(plotelems)>0:
        for i, val in enumerate(elemlines):
            axes.axvline(x=val, color='b') 
    return


def convertspelist(spelist):
    '''Convert spelist file log which has number of areas but not areanumber 
    to file pointer dataframe with filenumber, areanumber and filename, cleaner than
    previous clone param rows method '''
    mycols=['Filenumber','Areanumber','Filename','Evbreaks']
    areasdf=pd.DataFrame(columns=mycols)
    for index, row in spelist.iterrows():
        tempser=pd.Series(index=mycols)
        tempser.Filename=spelist.loc[index]['Filename']
        tempser.Filenumber=int(spelist.loc[index]['Filenumber'])
        tempser.Areanumber=int(spelist.loc[index]['Areas'])
        tempser.Evbreaks=int(spelist.loc[index]['Evbreaks'])
        areasdf=areasdf.append(tempser, ignore_index=True)
    return areasdf


def reportpeaks(spelist, plotelems, numplots, AESquantparams, Integquantlog, PDFname='Peaks_report.pdf'):
    ''' Generalized version of reportpeaksmajor working from spelist (loop over spatial areas 
    plots subtracted data for all spe files;  plot 
    integquantlog often has Gaussian fit attempts
    subplots are usually different elemental regions
    '''
    # single plot multiple areas 
    plt.ioff()
    with PdfPages(PDFname) as pdf:
        areasdf=cloneparamrows(spelist) # changes to iteration over areas not over spe files
        for i in range(0,len(areasdf)+1,numplots):
            thisplot=areasdf[i:i+numplots] # spatial areas to plot on this page
            spefilelist=np.ndarray.tolist(thisplot.Filename.unique()) # these go on single plot
            # should work for deriv or peaks data (direct counts not so much)
            Augerfiles=[] # list of dataframes
            for j, file in enumerate(spefilelist):
                try:
                    Augerfile=pd.read_csv(file) # reads entire spectra into df (all areas)
                except:
                    pass
        # TODO this looks screwed up for some reason
        for index,row in spelist.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=spelist.loc[index]['Filename']
            numareas=int(spelist.loc[index]['Areas'])
            # Grab first n areas (numplots var) and put together on same page with legend (different colors)
            
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

ratiolist='Mg/Si, Fe/Si, S/Si' # Pass as single string

ratiolist=['MgSi', 'FeSi', 'SSi']
elempairs=parseratios(ratiolist)

def parseratios(ratiolist):
    '''Convert list of ratios into list of lists '''

    return elempairs

def duplicatesplot(d, ratiolist,joinlist=['Sample']):
    ''' Find only duplicated measurements and plot these ratios to look for consistency '''
    # compare at % (basis would give same result though for ratios)
    ratiolist=ratiolist.split(',') # split into pairs
    ratiolist=[[x.split('/')[0],x.split('/')[1]] for x in ratiolist]

def plotduplicatesternary(duplicates, ternelems):
    ''' Take duplicated compositional data and plot on ternary diagram in different colors '''
        # Calculate ternary quantities for all in duplicates dataframe  
    hyphensplit = re.compile('(\+[a-zA-Z]+)').split
    ternlist=[part for img in ternelems for part in hyphensplit(img) if part]
    ternlist=[x.replace('+','') for x in ternlist]    
    try:
        duplicates['Tbasis']=duplicates[ternlist].sum(axis=1) # sum a list of columns 
    except:
        print('Failed sum.. missing data for given element?')
    # calculate T1, T2, T3 values for ternary plot (auto normalized to unity)
    for i, val in enumerate(ternelems):
        num=str(i+1)
        if '+' not in val:
            duplicates['T'+num]=duplicates[val]/duplicates['Tbasis']
        else:
            elems=[str(s) for s in val.split('+')]
            duplicates['T'+num]=duplicates[elems].sum(axis=1)/duplicates['Tbasis']
    # colorlist=['black','red','blue','green','cyan','magenta','yellow']
    # 26 distinct colors from Green-Armytage work 
    rgbacolors=[(240,163,255,1),(0,117,220,1),(153,63,0,1),(76,0,92,1),(25,25,25,1),(0,92,49,1),(43,206,72,1),(255,204,153,1),(128,128,128,1),(148,255,181,1),(143,124,0,1),(157,204,0,1),(194,0,136,1),(0,51,128,1),(255,164,5,1),(255,168,187,1),(66,102,0,1),(255,0,16,1),(94,241,242,1),(0,153,143,1),(224,255,102,1),(116,10,255,1),(153,0,0,1),(255,255,128,1),(255,255,0,1),(255,80,5,1)]
    # now parse compositional data into unique duplicates
    samplelist=duplicates.Sample.unique()
    samplelist=np.ndarray.tolist(samplelist)
    # Prepare ternary plot
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

    # data should be list of tuples or list of numpy arrays    
    samplelist=samplelist[:26] # only 26 available distinct colors
    for i, sample in enumerate(samplelist):
        match=duplicates[duplicates['Sample']==sample]
        ternpts=[] # holds list of tuples for ternary plot 
        for index, row in match.iterrows(): # deal with each spectrum in pair (or multiple)
            # need to calculate denom for ternary plot
            ternpts.append((duplicates.loc[index]['T1'],duplicates.loc[index]['T2'],duplicates.loc[index]['T2'])) # append this data 3-tuple
        tax.scatter(ternpts, marker='s', s=40, color=rgbacolors[i], label=sample)  # s is point size
    return duplicates # return with ternary compositional data


test= [[float(x/255) for x in l] for l in rgbacolors] # divide al
test2=[

tax.scatter(ternpts, marker='s', s=40, color=rgbacolors[i], label=sample)
tax.scatter(ternpts, marker='s', s=40, color=(0,0,255,1), label=sample)
        
    # can also specify colors as rgba tuples (r,g,b,a) all between zero and 1; pass colormap=False to heatmap()
def scatterratioplot(compdata,  ratiolist, thresh=0.1, joinlist=['Sample']):
    '''Compare elemental ratios using 2 different versions of composition calculation
        uses inner merge to select only subset with values from each df'''

    # compare at % (basis would give same result though for ratios)
    ratiolist=ratiolist.split(',') # split into pairs
    ratiolist=[[x.split('/')[0],x.split('/')[1]] for x in ratiolist]
    # Create subplot for each pair 
    numcols=min(len(ratiolist),2) # 1 or 2 columns
    numrows=math.ceil(len(ratiolist)/2)
    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False)
    # Merge dfs with comp1 and comp2 using inner join
    df=pd.merge(compdata, compdata, how='inner', on=joinlist, suffixes=('','b'))
    # inner merge but throw
    for j, elems in enumerate(ratiolist):
        [elem1, elem2]=elems 
        thisrow=(j)//2 # select this subplot (corrected for zero based indexing)
        thiscol=(j)%2
        rationame=elem1+elem2+'ratio'
        errname='err'+elem1+elem2
        elem1='%'+elem1 # convert to at % (has most complete propagated errors)
        elem2='%'+elem2
        # Drop those with blank error col for these elements  (4 cols due to 2 elements per comp dataset)
        thisdf=df.replace([np.inf, -np.inf], np.nan).dropna(subset=['err'+elem1,'err'+elem2,'err'+elem1+'b','err'+elem2+'b'])
        # Compute this elemental ratio using comp1 data
        df[rationame]=df[elem1]/comp1[elem2]
        df[rationame+'b']=df[elem1+'b']/comp1[elem2+'b'] # same for second set
        # error propagation err%Mg col is absolute error (not fractional)
        df[errname]=0 # blank error columns
        df[errname+'b']=0 
        for index, row in df.iterrows():
            fracterr1=df.loc[index]['err'+elem1]/df.loc[index][elem1] # elem1 has absolute err so convert to fractional
            fracterr2=df.loc[index]['err'+elem2]/df.loc[index][elem2]
            thiserr=np.sqrt(fracterr1**2+fracterr2**2) # percent error in calculated ratio
            df=df.set_value(index, errname, df[rationame]*thiserr) # set absolute err for this ratio
            # now propagate errors for comp2 data (b as suffix)
            fracterr1=df.loc[index]['err'+elem1+'b']/df.loc[index][elem1+'b'] # elem1 has absolute err so convert to fractional
            fracterr2=df.loc[index]['err'+elem2+'b']/df.loc[index][elem2+'b']
            thiserr=np.sqrt(fracterr1**2+fracterr2**2) # percent error in calculated ratio
            df=df.set_value(index, errname+'b', df[rationame+'b']*thiserr) # set absolute err for this ratio
        df.plot.scatter(x=rationame, y=rationame+'b', xerr=errname, yerr=errname+'b', ax=axes[thisrow,thiscol])
        # Slope,intercept=np.polyfit(data1, data2, 1)  numpy version
        xdata=df[rationame].as_matrix()
        ydata=df[rationame+'b'].as_matrix()
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xdata, ydata) # imported from scipy.stats
        theseoutliers=returnoutliers(xdata, ydata, elem1+'/') # return the outliers from these scatter plots
        # set x range for linear plot 
        text1=str(round(slope,2))+' *x +' + str(round(intercept,2))
        text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,3))
        xmax=max(max(data1),max(data2))*1.1 # set to slightly larger than max of dataset
        x=np.linspace(0,xmax,100) # setting range for 
        axes[thisrow, thiscol].text(0.025,0.9, text1, fontsize=12, transform=axes.transAxes)
        axes[thisrow, thiscol].text(0.025,0.8, text2, fontsize=12, transform=axes.transAxes)
        plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        
    # TODO organize returned dataset
    return df


def fileslice(filelist, AugerParamLog,Smdifpeakslog):
    '''Using shortened input list, return peaks and params for that subset
    string is of format filelist='138, 143, 198a2'''
    # pull tuples from filelist string (filenum and desired area)
    strlist=[str(s) for s in filelist.split(',')] # split into various ints
    fileareas=[] # list of tuples with filenum and area
    filenums=[]  # list of desired integer file numbers and desired area
    for i in range(0,len(strlist)): # construct tuple of filenum + area
        if 'a' not in strlist[i]:
            filenums.append(int(strlist[i])) # add to filtered list
            fileareas.append((int(strlist[i]),1)) # filenumber and area 1 for tuples list
        if 'a' in strlist[i]: # area is specified 
            tempstr=strlist[i][0:strlist[i].index('a')]
            area=int(strlist[i][strlist[i].index('a')+1:])
            filenums.append(int(tempstr))
            fileareas.append((int(tempstr),area))     
    Params=AugerParamLog.loc[AugerParamLog.Filenumber.isin(filenums)]
    Peaks=Smdifpeakslog.loc[Smdifpeakslog.Filenumber.isin(filenums)]
    return Params, Peaks, fileareas
    
def AESplotstack(Params,Peaks, plotrange, plotelems=plotelems):
    ''' Plot the Auger counts and the smooth-diff spectra for a given filenumber (details in Params& Peaks already sliced dataframes)
    plotrange and plotelems (vertical blue lines) are parsed to correct energy values
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    numareas=int(Params.Areas)
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    myplotrange=setplotrange(plotrange, Augerfile) # passes ev range or element or defaults to max range (returns tuple)
    
    # Things added to plots (breaks in multiplex scans, points used for S7D7 amplitudes)
    energyvals=AESplot.findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
    # add elemental peaks as vert lines
    if plotelems!='':
        elemlines = findelemlines(plotelems,myplotrange) # find energy vals of peaks within plotted range
        
    cols=numareas # 2 plots per area for counts & smooth-diff 
    rows=2     
    fig, axes = plt.subplots(nrows=rows, ncols=cols) # axes is array
    
    # slice spectrum and peaks to desired energy plotrange 
    Augerslice=Augerfile[(Augerfile['Energy']>myplotrange[0]) & (Augerfile['Energy']<myplotrange[1])]
    Peaks=Peaks[(Peaks['Peakenergy']>myplotrange[0]) & (Peaks['Peakenergy']<myplotrange[1])]
    
    for i in range(0, numareas): # myareas is a set and can be non-continuous
        colnum=i         
        plotpts=Peaks[(Peaks['Areanum']==(i+1))] # pull points from only this spatial area        
        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,colnum], color='r') # scatter plot from df
        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,colnum], color='r')                 
        colname='Counts'+str(i+1)        
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,colnum]) 
        colname='S7D7'+str(i+1)        
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,colnum])
        titlestring=AESplot.maketitlestring(plotpts)
        axes[1,colnum].set_title(titlestring, fontsize=10)
        # add evbreaks as red vertical lines
        for j, val in enumerate(energyvals):
            axes[0,colnum].axvline(x=val, color='r') # on counts plot
            axes[1,colnum].axvline(x=val, color='r') # on smooth-diff plot
        for j, val in enumerate(elemlines):
            axes[0,colnum].axvline(x=val, color='b') # on counts plot
            axes[1,colnum].axvline(x=val, color='b') # on smooth-diff plot
            
# is this obsolete... use AESplot1
def plotAESreg(Params, Peaks, plotrange='', areas='', markpeak=1):   
    ''' Utility function to quickly look at peak region out of entire Auger spectrum df 
    Params is df with selected file(s), number of areas
    kwargs list: 
    1) range=, string with element name or direct hyphenated range in eV (ie. 100-300); defaults to entire data range
    2) markpeak=1 or 0; if 1 shows points from which amplitude was determined and adds peak ampl to plot
    3) areanum= ; defaults to area 1 or all areas?
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    numareas=int(Params.Areas)    
    myareas=parseareas(areas, numareas) # set of areas for plotting
        
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    myplotrange=setplotrange(plotrange, Augerfile) # passes ev range or element or defaults to max range (returns tuple)
    # Slice of df from full range Augerfile (better solution to y autoscale problem)
    Augerslice=Augerfile[(Augerfile['Energy']>myplotrange[0]) & (Augerfile['Energy']<myplotrange[1])]
    
    energyvals=AESplot.findevbreaks(Params, Augerslice) # get x energy vals for evbreaks in this multiplex (if they exist) as float
    
    # set plot dimension based on # of areas
    cols=divmod(len(myareas),2)[0]+ divmod(len(myareas),2)[1] #num of subplot cols, round odds up to next even
    # cols=divmod(myareas,2)[0]+ divmod(myareas,2)[1]
    fig, axes = plt.subplots(nrows=2, ncols=cols) # axes is array
    
    if markpeak==1: # can easily deal with entire range
        quantpoints=Peaks[(Peaks['Peakenergy']>myplotrange[0]) & (Peaks['Peakenergy']<myplotrange[1])]
        # don't need to slice by area here... do that when plotting
        quantpoints=findpeakdata(logmatch, elem, myareas) # return list of tuples with (areanum, index, value) in quantpoints 
    
    # default here to plot of all areas
    # create data plots + evbreak lines + x,y data from 
    whichplots=[]
    for i,areanum in enumerate(myareas): # myareas is a set and can be non-continuous
        if (i+1)%2==1:
            rownum=0
        else:
            rownum=1
        colnum=int((i+1)/2.1)
        temptuple=(areanum, rownum, colnum)
        whichplots.append(temptuple)
        colname='S7D7'+str(areanum)
        Augerslice.plot(x='Energy', y=colname, ax=axes[rownum,colnum])
        
    # add vert lines for all axes in array
    for i in range(0,axes.shape[0]):
        for j in range(0,axes.shape[1]):
            axes[i,j].axvline(x=val, color='r')

    # add scatter points
    if markpeak==1:
        for i, areanum in enumerate(myareas):
            thisplot=whichplots[i] # tuple associating spatial area number with position in plot
            area=thisplot[0]
            rownum=thisplot[1]
            colnum=thisplot[2]
            print('areanum, Area, row, column are ', areanum, area, rownum, colnum)
            mypoints=quantpoints[(quantpoints['Areanum']==areanum)]     
            mypoints.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[rownum,colnum])
            mypoints.plot.scatter(x='Peakenergy'-'Peakwidth', y='Posintensity', ax=axes[rownum,colnum])    
    
    quantpoints.plot(x='Peakenergy')
    axes[0,0].scatter
    fig, axs = plt.subplots(nrows=2, ncols=cols) # axes is array
    for i in np.nditer(axes):
        print (axes)
    axarr[0,0].plot
    # plot multiplex breaks as vertical line (find evbreaks)   
    
    for i, val in enumerate(energyvals): 
        print(val)
        plt.axvline(x=val, color='r') # vertical red line to active subplot w/o rescale    

   # add a bunch of X-Y scatter points to same graph (structured array or pandas
    plt.scatter(x,y) # 
    vlines()
    Augerfile.plot(x='Energy', y='S7D71')

def plotAESmultiple(logmatch, plotrange='', areanum='', peaks=1): 
    # Pass spectra from multiple spe/csv files and plot in same window
    
    # TODO finish plotting of multiple spectra
    subplots_adjust(hspace=0.000)
    number_of_subplots=len(myareas)
    
    for i in range(0,len(logmatch)):
        # basically same as single plot version  
        # Set up plot for n spatial areas 
    for numarea in range(1, len(myareas)+1):
        colname='S7D7',str(numarea)
        Augerfile.plot(x='Energy', y=colname, xlim=myplotrange)
    for i in range(0, len(quantpoints)): # add peak points to plot