# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:11:09 2016
Contains all functions related to peak finding in smooth-differentiated Auger spectra (which are called by Auger_quant_batch.py)

findindices - takes user's list of desired elements and background regions, checks each spe file and
returns ideal positions of those peaks (usually index # of peak, not peak energy)
Any desired peaks that are out of the data's range are knocked out here (so desired element list can be made more inclusive without worrying about functions returning errors or garbage)

findnegpeak & findpospeak -- after Auger file csv is sliced into region right around peak, these find and 
return minima, maxima and their positions

createpeakdf -- called from various location to create empty dataframes with the correct columns for storing peak data

smdifquant (called from main)- finds and returns elemental peaks in each given spectrum; logmatch argument contains selected row from Augerparams.csv which has everything you need to know about a given spe (esp. path to file which is then loaded)

smdifbackground (called from main) - finds and returns magnitude of background fluctuations (noise) which can be used to determine the significance of a given peak's amplitude

Addbacktopeaks (called from main) - for each elemental peak it finds the closest regions where background noise was measured and adds that info to peak dataframe

checkpeakboundaries - whenever the original spectrum is sliced into energy subregions, this ensures that a boundary (i.e. between different multiplex regions) is not included
avoids garbage peaks in smooth-diff spectra related to discontinuities in multiplex scans


Precondition: Auger_batch_import already run on spe files in batch mode, leading to creation of 
Augerparamlog.csv
@author: tkc
"""
#%%
import re, csv, datetime, os
import pandas as pd
import numpy as np
from collections import defaultdict
import tkinter as tk
#%%

# def checkindices(Elements, Backregs): # maybe faster version of findindices if same energy structure in multiplex?

def smdiffquant_gui(spelist, Elements, AESquantparams, Smdifpeakslog):
    ''' tk Front-end for batch quant... 
    1) filter by filenumber available (all areas) 
    2) reprocess (bool)
    3) apply predetermined large shift (for samples with charging)  '''
    root = tk.Tk()
    root.title("Run quant on derivative spectra ")
    # Filtering by filenumber or filename or combo files only
    filterstr=tk.StringVar() # comma separated or range of filenumbers for plot or string for sample name
    filtercol=tk.StringVar()  # to which column is optional filter applied
    filtercol.set('Filenumber') # set default to filename (not sample)
    combobool=tk.BooleanVar() # Select only combine-averaged files
    # Apply custom shift 
    customshiftbool=tk.BooleanVar()
    shift=tk.IntVar()
    choice=tk.StringVar() # for action button

    # string showing active elements
    elemstr=tk.StringVar()
    mytext=', '.join(Elements) # elements for labelling 
    elemstr.set(mytext)
    # reprocess existing files or not
    reprocess=tk.BooleanVar()  
    # autosave altered file or not
    autosave=tk.BooleanVar()
    autosave.set(True)
    # Optional filtering of chosen spectra by filenumber or string
    rownum=0
    tk.Label(root, text='Filenumber(s) or sample name filter').grid(row=rownum, column=0)
    tk.Radiobutton(root, text='Filenumber(s) filter ', value='Filenumber', variable = filtercol).grid(row=rownum, column=1)
    tk.Checkbutton(root, variable=reprocess, text="Choose combine-averaged files").grid(row=rownum, column=2)
    rownum+=1
    tk.Radiobutton(root, text='Filter on sample column', value='Filename', variable = filtercol).grid(row=rownum, column=1)
    tk.Entry(root, textvariable=filterstr).grid(row=rownum, column=0)
    rownum+=1
    tk.Label(root, text="Elements for quant").grid(row=rownum, column=0)
    tk.Entry(root, textvariable=elemstr).grid(row=rownum, column=1)
    rownum+=1
    tk.Checkbutton(root, variable=reprocess, text="Reprocess all spectra").grid(row=rownum, column=0)
    rownum+=1
    tk.Checkbutton(root, variable=autosave, text="autosave smdifquant log file?").grid(row=rownum, column=0)
    rownum+=1
    tk.Checkbutton(root, variable=customshiftbool, text='Apply custom shift?').grid(row=rownum, column=0)
    tk.Entry(root, textvariable=shift).grid(row=rownum, column=1)
    rownum+=1
    
    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def quantify(event):
        choice.set('quantify')        
        root.destroy()  
        
    def repickelems(event):
        ''' trigger pickelemsgui again (from within this tk interface) '''
        Elements=repickelems(AESquantparams)
        # Elements=AESutils.pickelemsGUI(AESquantparams)
        print('Reselected elements:')
        print(', '.join(Elements))
        # TODO fix this... passing back unaltered 
        mytext=', '.join(Elements)
        elemstr.set(mytext)
        
    a=tk.Button(root, text='Repick elements')
    a.bind('<Button-1>', repickelems)
    a.grid(row=rownum, column=0)

    a=tk.Button(root, text='Abort')
    a.bind('<Button-1>', abort)
    a.grid(row=rownum, column=1)
    
    a=tk.Button(root, text='Quantify')
    a.bind('<Button-1>', quantify)
    a.grid(row=rownum, column=2)
    
    root.mainloop()
    
    if choice.get()=='quantify':
        kwargs={} # reconstruct with new choices
        if filtercol.get()=='Filenumber' and filterstr.get()!='':
            filenums=parsefilenums(filterstr.get())  # list of ints 
            spelist=spelist[spelist['Filenumber'].isin(filenums)]
            # store filenums used in kwargs (for possible next run)
            kwargs.update({'fileset':', '.join([str(i) for i in filenums])})
        elif filtercol.get()=='Filename' and filterstr.get()!='':
            filstr=filterstr.get()
            # Assumes filtering on string in sample column
            spelist=spelist[spelist[filtercol.get()].str.contains(filstr)]
            if len(spelist)==0:
                print('No files with ',filstr, 'in sample name')
                return
        if combobool.get():
            # Choose only combine averaged files (named 100103 or comparable)
            spelist=spelist[spelist['Filenumber']>10000]
        if reprocess.get():
            kwargs.update({'reprocess':True})
        if shift.get()!=0:
            # add custom offset to quant params (but don't change core version)
            Aquant=AESquantparams.copy()
            Aquant['negpeak']+=int(shift.get())
            Aquant['pospeak']+=int(shift.get())
        else:
            Aquant=AESquantparams.copy()
            # alter all negpeak in AESquantparams 
        Smdifpeakslog=smdifbatchquant(spelist, Elements, Aquant, Smdifpeakslog, **kwargs)
        if autosave.get(): # flag to automatically write changes to log
            Smdifpeakslog.to_csv('Smdifpeakslog.csv', index=False)
        return Smdifpeakslog

    return Smdifpeakslog

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
   
def compareavg_subs(Smdifcomp, Smdifcompsubs):
    '''Compare compositions on same samples derived from avg-combined files named 173175
    and the sub-spe files 173,174,175, single line for each'''
    combofiles=Smdifcomp[Smdifcomp['Filenumber']>3000] # this is core of returned df
    tempcols=np.ndarray.tolist(combofiles.columns.unique())
    elemlist=[col for col in tempcols if '%' in col]
    elemlist=[el.replace('%','') for el in elemlist]
    # make merge1, merge2 and merge3 cols for pd.merge
    # Create copies of Smdifcompsubs and rename Filenumber to merge1,2,3
    # This allows pd merge on merge col, areanumber (since sample might not be unique)

    # determine # of underlying sub spe files (assumed same for all in combofiles)
    filestr=str(int((combofiles.iloc[0]['Filenumber'])))
    if len(filestr)==3: # e.g. 221 for combined file is avg-combination of 2 through 21
        firstnum=int(filestr[0:1])
        lastnum=int(filestr[1:])
    else:
        firstnum=int(filestr[0:int((len(filestr)/2))])
        lastnum=int(filestr[int((len(filestr)/2)):])
    for i in range(1,(lastnum-firstnum+2)):
        nummerges=lastnum-firstnum+1
        colname='Merge'+str(i)
        combofiles[colname]='' # add correct number of merge columns
    for index, row in combofiles.iterrows():        
        filestr=str(int((combofiles.loc[index]['Filenumber'])))
        firstnum=int(filestr[0:int((len(filestr)/2))])
        lastnum=int(filestr[int((len(filestr)/2)):])
        # maybe need to associate file append string to source file
        filerange=[]
        filerange.extend(range(firstnum,lastnum+1))
        # deal with possible irregular combo file structure (assumes same number of averaged sub spe files)        
        if len(filerange)!=nummerges:
            print('Different number of average sub spefiles for ', filestr)
            if len(filerange)!=nummerges: # truncate if too long to avoid error
                filerange=filerange[0:nummerges+1]
        # Now assign correct filenumber for each merge
        for i, val in enumerate(filerange):
            colname='Merge'+str(i+1)
            combofiles=combofiles.set_value(index, colname, val) # assign correct filenumber to each merge (cols already created)
    # now read for n merges with compositional data from 
    # prepare compsubs for merge
    dropcols=['Project', 'Filename', 'FilePath', 'Sample', 'Comments']
    dropcollist=[s for s in Smdifcompsubs.dtypes.index if s in dropcols] # ensure the drop col is actually present
    Smdifcompsubs=Smdifcompsubs.drop(dropcollist, axis=1)
    for i in range(1,nummerges+1):
        colname='Merge'+str(i)
        tempdf=Smdifcompsubs
        tempdf=tempdf.rename(columns={'Filenumber':colname})
        combofiles=pd.merge(combofiles,tempdf, how='left', on=[colname,'Areanumber'],  suffixes=('',str(i)))
    # create and return a subset of compositional comparison
    mycols=['Filenumber','Areanumber', 'Sample','Filename','AESbasis']

    # now make average and stdev from composition of sub spes for comparison with average
    # Now compute avg and stdev of elemental bases and at. % compositions
    for i, elem in enumerate(elemlist):
        mycols.extend([elem+'ampl','%'+elem]) # add to truncated output
        numrange=[str(i) for i in range(1,nummerges+1)]
        collist=[elem+'ampl'+val for val in numrange]
        # average amplitudes for each element
        newcol=elem+'amplavg'
        mycols.extend([newcol])
        combofiles[newcol]=combofiles[collist].mean(axis=1)
        newcol=elem+'amplstdev'
        mycols.extend([newcol])
        combofiles[newcol]=combofiles[collist].std(axis=1)
        # now compute average at.% for each element
        collist=['%'+elem+val for val in numrange]
        newcol='%'+elem+'avg'
        mycols.extend([newcol])
        combofiles[newcol]=combofiles[collist].mean(axis=1)
        newcol='%'+elem+'stdev'
        mycols.extend([newcol])
        combofiles[newcol]=combofiles[collist].std(axis=1)
        # also keep S1, S2, S3 (or whatever depending on nummerges)
        for j in range(1,nummerges+1):
            mycols.append(elem+'ampl'+str(j))
    # Output a smaller subset of 
    compslice=combofiles[mycols]
    return compslice,combofiles

def getpeakstats(df, Elements):
    ''' Returns summary statistics for peak amplitude and shifts from smooth-diff quant method '''
    mycols=['Element','Count','Ampmean','Ampstdev','Ampmin','Ampmax','Shiftmean','Shiftstdev','Shiftmin','Shiftmax']
    peakstats=pd.DataFrame(columns=mycols)
    for i, elem in enumerate(Elements):
        thisrow=pd.Series(index=mycols)
        thisel=df[df['PeakID']==elem]
        if len(thisel)>0: # it's possible that element wasn't included in smdiffquant run
            thisrow.Element=elem
            thisrow.Count=len(thisel)
            thisrow.Ampmean=int(thisel.Amplitude.mean())
            thisrow.Ampstdev=int(thisel.Amplitude.std())
            thisrow.Ampmin=int(thisel.Amplitude.min())
            thisrow.Ampmax=int(thisel.Amplitude.max())
            thisrow.Shiftmean=thisel.Shift.mean()
            thisrow.Shiftstdev=thisel.Shift.std()
            thisrow.Shiftmin=thisel.Shift.min()
            thisrow.Shiftmax=thisel.Shift.max()
            peakstats=peakstats.append(thisrow, ignore_index=True)
    return peakstats
 
def organizecolumns(df1,mycols):
    ''' Pass df and template (list of desired columns in desired order) and return reorganized newdf
    '''
    cols1=df1.columns.tolist()
    newdf=df1 # avoids modification of passed df
    uniquelist=[i for i in cols1 if i not in mycols]
    for i,colname in enumerate(uniquelist): # remove cols from df1 that are absent from df2
        # newdf.drop(colname, axis=1, inplace=True) # this modifies both passed and returned dfs
        newdf=newdf.drop(colname, axis=1)
    newdf=newdf[mycols] # reorder columns based on template df
    return newdf
    
def parseelemlist(elemlist):
    '''Find and separate multielement peaks to be averaged (e.g. Fe2 & Fe) from longer string of element peaks
    e.g. splits "Mg Fe Fe2 Si" into "Mg Si" and "{Fe,[Fe,Fe2]} dictionary'''
    # Strip numbers from strings within list  
    newlist=[re.match('\D+',i).group(0) for i in elemlist]
    
    # find duplicated peaks (multiple peaks per element)
    Multielem = defaultdict(list)
    for i, item in enumerate(newlist):
        Multielem[item].append(i)
    Multielem = {k:v for k,v in Multielem.items() if len(v)>1} # dictionary with duplicated item and list with indices
    
    duplist=list(Multielem.values()) # get list  
    duplist=[item for sublist in duplist for item in sublist] # single list with positions of duplicated elements
    
    # now alter multipeak elements list to give dict with element and then list of peak for that element    
    for key,value in Multielem.items():
        templist=value # dictionary value is list of elem peak index positions
        peaklist=[]
        for i, index in enumerate(templist): # create new list with original elem peak from index positions
            peaklist.append(elemlist[index])
        # now replace list of index positions with elempeak names
        Multielem.update({key:peaklist}) # key will be multipeak element string i.e. "Fe"
    # finally construct new single elements list with multipeak ones removed (handle each separately)
    newelemlist=[]
    for i in range(0,len(elemlist)):
        if i not in duplist:
            newelemlist.append(elemlist[i])
    return newelemlist, Multielem
    
def parseelem2(elemlist, Multielem):
    ''' After multielement peaks removed, also move secondary peaks used as primary to dict (handle separately)
    e.g. splits "S Mg Fe2 Si" into "S Mg Si" and "{Fe,[Fe2]} dictionary; same structure and df output 
    for averaging of Fe, Fe2, or straight Fe2 or straight Fe'''
    
    newelemlist=[]
    for i, elem in enumerate(elemlist):
        if re.search(r'\d',elem): # has number
            match=re.search(r'\d',elem)
            newkey=elem[0:match.start()] # This is the element name (as opposed to peak name)
            templist=[] # peakIDs added as list (of length 1)
            templist.append(elem) # list containing single string (keeps identical data structure)
            # split at the number to get associated element
            Multielem.update({newkey:templist}) # add to existing dictionary for separate handling
        else:
            newelemlist.append(elemlist[i]) # just copy over 
    return newelemlist, Multielem # return altered element list and multielem dictionary    
    
def thresholdtest(df, threshold):
    '''Pass df holding smdifpeak info on single peak; return df if ratio above threshold, return empty frame if below
	if below threshold, set adjamp to zero (but leave raw amplitude at same value '''
    # avgbackamplitude is average measure of noise level in smoothed-differentiated spectrum 
    avgbackamp=(df.iloc[0]['Lowbackamplitude']+df.iloc[0]['Highbackamplitude'])/2
    if df.iloc[0]['Amplitude']<threshold*avgbackamp:
        df=df.set_value(df.index[0],'Adjamp', 0) # zero out adjusted amplitude
        return df
    else:
        return df # return empty df if amplitude below threshold

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
    return df

# Testing df=spelist.copy()  elemlist=Peaks
def calccomposition(df, Smdifpeakslog, elemlist, threshold=0.0):
    '''Calculate elemental composition of given subset of files (generally an spelist) based on input element list 
    can also pass a later list of files (i.e. outliers) that are already split by filenumber/areanumber    
    returns raw amplitudes, adjusted amplitudes and at % calculations for each in element list (above threshold)
    should elements be eliminated if amplitude is less than 2x that of noise background?
    threshold - ratio of element peak to noise peak (0 means no threshold applied '''

    elemlist, multipeaklist = parseelemlist(elemlist) # list of single peak elements and dict with multipeaks
    # check if any of the single peaks are secondary (i.e. quant on Fe2 not main Fe)
    elemlist, multipeaklist= parseelem2(elemlist, multipeaklist)
    
    # to process compositions from multiple areas, clone rows from spe log (one for each areanum)
    if 'Areas' in df: # indicates augerparamlog with single row per multi-area spectrum
        df=cloneparamrows(df)
    df=df.reset_index(drop=True)
    df['AESbasis']=0.0 # resets to zero if already present from calcamplitude
    mycols=['Filenumber', 'Project', 'Filename', 'FilePath', 'Sample', 'Comments',
            'Phase','AESbasis','Areanumber']
    for i, col in enumerate(mycols):
        if col not in df.columns:
            df[col]=''
    for i, elem in enumerate(elemlist):  # add columns for raw peak amplitudes and significance (smooth-differentiation method)
        df[elem+'ampl']=0.0 # add col for each element to df
        df['sig'+elem]=0.0 # significance (based on ratio of peak amplitude over average noise amplitude)
        mycols.append(elem+'ampl')
        mycols.append('sig'+elem)
    for i,elem in enumerate(list(multipeaklist.keys())): # do same for multipeak elements
        df[elem+'ampl']=0.0
        df['sig'+elem]=0.0 
        mycols.append(elem+'ampl')
        mycols.append('sig'+elem)
    for i, elem in enumerate(elemlist):  # add columns for basis
        df[elem]=0.0 # add col for each element to df
        mycols.append(elem)
    for i,elem in enumerate(list(multipeaklist.keys())): # get elements (keys) from dict
        df[elem]=0.0
        mycols.append(elem)
    for i, elem in enumerate(elemlist):  # now add at.% columns (e.g. %S, %Mg)
        colname='%'+elem # at % columns named %S, %Mg, etc.
        mycols.append(colname)  # add to column list template
        df[colname]=0.0
    for i,elem in enumerate(list(multipeaklist.keys())): # add multipeak elements
        colname='%'+elem # at % columns named %S, %Mg, etc.
        mycols.append(colname)  # add to column list template
        df[colname]=0.0
    for index, row in df.iterrows(): # loop through each desired spectrum
        filename=df.loc[index]['Filename']
        areanum=df.loc[index]['Areanumber']
        match=Smdifpeakslog[Smdifpeakslog['Filename']==filename] # find smdif data for this filenumber
        # TODO add loop through different areas 
        match=match[match['Areanumber']==areanum] # matches correct spatial area number
        basis=0.0 #
        # Testing elem=elemlist[0]
        for j, elem in enumerate(elemlist): # handle the single peak elements
            temp=match[match['PeakID']==elem] # finds entry for this element 
            if len(temp)==1: # peakid might be skipped (len 0) if out of data range
                avgnoise=(abs(temp.iloc[0]['Lowbackamplitude'])+abs(temp.iloc[0]['Highbackamplitude']))/2 #
                thissig=temp.iloc[0]['Amplitude']/avgnoise # ratio of peak amplitude to a measure of spectral noise amplitude
                df=df.set_value(index, 'sig'+elem, thissig) # store this value (even if removed from basis)
                temp=thresholdtest(temp, threshold) # zeros out adjamp if below noise threshold 
                colname=elem+'ampl'
                df=df.set_value(index, colname, temp.iloc[0]['Amplitude']) # copy raw amplitude of this element
                df=df.set_value(index, elem, temp.iloc[0]['Adjamp']) # copy adjusted amplitude of this element
                basis+=temp.iloc[0]['Adjamp'] # add this element's value to AES basis
        # now handle the multipeak elements (get average value from both peaks)
        for key, value in multipeaklist.items(): # key is element (aka colname in df), value is list of peaks in Smdifpeakslog
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # this is number of lines that are average (i.e. 2 for Fe&Fe2)            
            avgval=0.0 # working value for averaged adjamplitude
            avgvalraw=0.0 # working value for raw amplitude (just average even though this is a bit goofy)
            avgrat=0.0 # for significance ratio if using multiple lines again take the average
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                temp=match[match['PeakID']==peak] # finds entry for this peak
                if len(temp)==1: # peakid might be skipped (len 0) if out of data range
                    avgnoise=(abs(temp.iloc[0]['Lowbackamplitude'])+abs(temp.iloc[0]['Highbackamplitude']))/2 #
                    avgrat+=temp.iloc[0]['Amplitude']/avgnoise # ratio of peak amplitude to a measure of spectral noise amplitude
                    temp=thresholdtest(temp, threshold) # zeros out adjamp if below noise threshold                
                    avgval+=temp.iloc[0]['Adjamp']
                    avgvalraw+=temp.iloc[0]['Amplitude']
                else:
                    numlines=numlines-1 # if peak is zeroed out and not added, this reduces # peaks in average
            if numlines>0: # avoid divbyzero if peak is too small
                avgval=avgval/numlines # this is now average basis for given element
                avgvalraw=avgvalraw/numlines # this is now average basis for given element
                avgrat=avgrat/len(templist) # ratio calc done for all peaks (none knocked out)
            df=df.set_value(index, key, avgval) # copy adjusted amplitude of this element
            df=df.set_value(index, 'sig'+key, avgrat) # copy average ratio of peak relative to noise amplitude of this element
            df=df.set_value(index, key+'ampl', avgvalraw) # copy adjusted amplitude of this element
            # add value from this element to AESbasis
            basis+=avgval
        df=df.set_value(index, 'AESbasis', basis) # write total basis value to df 
        # now compute at.% for each listed element
        for j, elem in enumerate(elemlist):
            colname='%'+elem
            ratio=df.loc[index][elem]/df.loc[index]['AESbasis']
            df=df.set_value(index, colname, ratio)
        # also calculate for elements w/ multiple peaks (if present)
        for key, value in multipeaklist.items(): 
            colname='%'+key
            ratio=df.loc[index][key]/df.loc[index]['AESbasis']
            df=df.set_value(index, colname, ratio)
        # end of loop calculation for each spectrum 
                
    # organize data based on mycols template
    df=organizecolumns(df,mycols)
    return df

def calcamplitude(df, AESquantparams):
    '''For each elemental peak in smdifpeaks log, calculate basis with k-factor and mass
    result stored in adjamp column and used for subsequent compositional determinations
    can change AESquantresults and recalc at any time ''' 
    if 'Adjamp' not in df:
        df['Adjamp']=0.0 # new column for adjusted amplitude (if not already present)
    # what's the fastest way to do this?
    # loop for each element, mask df, get appropriate k-factor & mass
    df=df.reset_index(drop=True) # go ahead and reset index
    elemlist=np.ndarray.tolist(df.PeakID.unique()) # list of unique elements from df
    for i,elem in enumerate(elemlist):
        match=AESquantparams[(AESquantparams['element']==elem)]
        match=match.reset_index(drop=True)
        kfactor=match.iloc[0]['kfactor'] # kfactor and mass for this element/peak
        mass=match.iloc[0]['mass'] 
        elemmask=(df['PeakID']==elem) # mask for this element in loop 
        for j in range(0,len(df)): # loop and set adjamplitude to amp*kfact/mass
            if elemmask[j]==True: # row has this element
                newval=df.iloc[j]['Amplitude']*kfactor/mass
                df=df.set_value(j,'Adjamp',newval)
    return df


def getbestduplicate(df, Smdifpeakslog, elemlist):
    ''' Used to choose a single best spectrum for each sample based on elemlist basis (drop inferior ones based on value in sortparam)
    Pass df and column name for dupl. search and params sort param for finding best value
    '''
    # choose best on spectral area #1 and ignore others for now
    #TODO modify to perform separate AES basis calculation for multiple areas
    df=df.reset_index(drop=True)
    df['AESbasis']=0.0
    # for each spectrum in spelist, calculate underlying AES basis from Smdifpeakslog (filenumber match)
    for i in range(0,len(df)):
        filenum=df.iloc[i]['Filenumber']
        match=Smdifpeakslog[Smdifpeakslog['Filenumber']==filenum] # find smdif data for this filenumber
        match=match[match['Areanumber']==1] # select area 1
        AESbasis=0.0 # init temp variable for basis calculation
        for j, elem in enumerate(elemlist):
            temp=match[match['PeakID']==elem] # finds entry for this element 
            if len(temp)==1: # avoid problems if certain element is not in range of scan
                AESbasis+=temp.iloc[0]['Adjamp'] # grab and sum adjusted amplitude for this single match
        df=df.set_value(i,'AESbasis', AESbasis) # write this summed value to df col before sort
    
    df=df.sort_values(['Sample','AESbasis'], ascending=False)    
    df=df.drop_duplicates(['Sample'], keep='first') # keeps best duplicate of multiples on same sample name
    return df
    
def findindices(Elements, Backregs, logmatch, AESquantparams):
    ''' takes element strings and energies of background regs and returns tuple for each elem symbol containing all params necessary to find each Auger peak from given spe file 
    also returns 2-tuple with energy val and index of chosen background regions
    ''' 
    # load standard Auger values from AESquantparams (csv encoding must be utf-8)
    AugerFileName=logmatch.Filename # full filename w/o path stored in 3rd col
    # get energy values
    Energyvals = [] #empty list for energy values
    values = csv.reader(open(AugerFileName, 'r'))
    try: # this can occasionally fail
        for row in values: # read energy column (first col position)
            Energyvals.append(row[0]) # get first energy column only
    except:
        print('Problem reading energy column for ', AugerFileName)
        Backregdata=[]
        return Energyvals, Backregdata # return empty lists and test
    Energyvals=Energyvals[1:] # remove header value
    Energyvals=[float(i) for i in Energyvals] # convert energy list to floats
    Elemdata=[] # returns list of length5 tuples for all elements
    for elem in Elements:
        # find row in AESquantparams for this element
        thiselemdata=AESquantparams[(AESquantparams['element']==elem)]
        thiselemdata=thiselemdata.squeeze() # series with this elements params
        negpeakev=thiselemdata.negpeak # ideal energy value of negative Auger peak in smooth-diff spectrum
        pospeakev=thiselemdata.pospeak # same for associated positive peak (at lower energy than negative one)
        width=int(thiselemdata.searchwidth) # search width used to find actual peak in real data
        kfact=thiselemdata.kfactor # typical sensitivity k-factor associated with element
        # find index # for ideal neg and pos peaks... use lambda funct.
        #   min(Energyvals, key=lambda x:abs(x-negpeakev)) gives value but not index #
        
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-negpeakev)) # tuple with index of closest and closest value
        negpeak=temptuple[0] # first of tuple is index #
        peakinrange=temptuple[1]-negpeakev # should be ~0 if desired peak is in data range
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-pospeakev)) # tuple with index of closest and closest value
        pospeak=temptuple[0] 
        # return list of length numelements with 5-tuple for each containing 1) element symbol, index # of ideal energy for 2)negpeak and 3) pospeak (not energy value itself!), 4)associated width to search for peak (element dependent) 5)sensitivity kfactor.. 
        # Skip typical width for peak of given element but can look up later if necessary
        if abs(peakinrange)<1: # Must skip entire desired element here if it's out of range of the data in this particular spe        
            elemtuple=(elem, negpeak, pospeak, width, kfact) # add tuple with info for this element     
            Elemdata.append(elemtuple)
        else:
            print('Warning: No quant for ',elem,' for ',AugerFileName, 'data not collected in this energy range.')
    Backregdata=[] # returns list of length2 tuples (only energy val and associated index number), default search width of 10 channels
    for reg in Backregs: # reg will be an ev value
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-reg)) # tuple with closest index # and associated closest value
        regindex=temptuple[0] # index number of this energy value for given spe file (not straightforward for multiplex with unknown eV gaps)
        if temptuple[1]!=reg:
            actualreg=temptuple[1] # overwrite desired energy (in case out of range)
            print('Warning..background region shifted to', actualreg, 'for file ', AugerFileName,': desired energy region out of data range')
            regtuple=(actualreg, regindex) # reassign in case desired background region is out-of-data range (just pick closest region)
        else:
            regtuple=(reg, regindex) 
        Backregdata.append(regtuple)
    #Elemdata is a list (of length number of elements) containing length5 tuples (tuple contain everything needed by smdifquant to search/find that element)
    return Elemdata, Backregdata

def findnegpeak(peakregion):  # pass slice of pandas Series within searchwidth of ideal peak
    ''' finds negative value and peak position in smooth-diff spectra around selected element'''
    negval=peakregion.min()  # min value from series slice around this peak
    negindex=peakregion.idxmin() # gives index of minimum from series
    return negval, negindex

def findpospeak(peakregion): # pass slice of pandas Series within searchwidth of ideal peak
    ''' finds positive value and peak position in smooth-diff spectra around selected element; pospeak to left in energy from negpeak in S7D7 smoothdiff spectrum'''
    posval=peakregion.max()  # min value from series slice around this peak
    posindex=peakregion.idxmax() # gives index of minimum from series (not energy value but index from original dataset)
    return posval, posindex
    
def checkpeakboundaries(lower,upper, index, evbreaks):
    ''' modifies boundaries of background peak search regions if these cross and eV break in multiplex scan
    this typically shouldn't occur unless you make energy regions too small in multiplex scans
    '''
    ends=[i for i in evbreaks if i > lower-3 and i < upper+3]
    # must be min 3 channels from endpoint (since S7D7 set to zero in this range)
    if ends: # list is not empty (aka selected background region crosses an ev break)
        if ends[0]>=index: # modify upper bound
            upper=ends[0]-3
            if upper-lower<10: # also change lower in case resulting region is too small
                lower=upper-10                
        else: # modify lower bound
            lower=ends[0]+3
            if upper-lower<10:
                upper=lower+10
        # print('Background region indices shifted to ', lower, 'and ', upper)
    return lower, upper

## Currently peak searches assume 1 eV/channel ... need to rescale search widths 
def addbacktopeaks(Smdifpeaks, Smdifback):
    '''For each elemental peak add the measure of background noise in the closest regions (from correct spatial area)
    '''
    for i in range(0,len(Smdifpeaks)):
        areanum=Smdifpeaks.Areanumber[i] #multiple spatial areas from same spe are passed together
        thisarea=Smdifback[(Smdifback['Areanumber']==areanum)] # Slice to select only background fits for correct area number
        peakval=Smdifpeaks.Peakindex[i]
        thisarea=thisarea.reset_index() # index must be reset for below argsort to work (avoids nan return)
        matches=thisarea.ix[(thisarea['Peakindex']-peakval).abs().argsort()[:2]] # then find two closest values to correct peakindex
        matches=matches.sort_values(['Peakenergy'], ascending=True) # sort matches in order of peak index
        Smdifpeaks.Lowback[i]=matches.iloc[0]['Peakenergy'] # energy val and not index 
        Smdifpeaks.Lowbackamplitude[i]=matches.iloc[0]['Amplitude']
        Smdifpeaks.Highback[i]=matches.iloc[1]['Peakenergy'] # energy val and not index 
        Smdifpeaks.Highbackamplitude[i]=matches.iloc[1]['Amplitude']
    return Smdifpeaks    
        
def smdifquant(logmatch, Elemdata):
    ''' Quant routine for finding peaks in smooth-differentiated spectrum (opens source spectrum as Augerfile, finds all peaks using Elemdata, 
    returns dataframe with position/amplitude/etc; desired elements out of data range are skipped (in prior findindices function)
    peakenergy- actual negpeak location in eV, shift- deviation from ideal value, (negint/posint/amplitude)- all related to magnitude of smooth-diff peak, peakwidth- width of peak as measured by derivative method (i.e. negative peak at left edge and pospeak at right; 
    5 column with *back* are measures of noise amplitude in differentiated spectrum
    '''
    #create Smdifpeaks dataframe for temp storage of each peak's params
    AugerFileName=logmatch.Filename  # get Auger filename from Series... should be .csv file in current directory
    numareas=int(logmatch.Areas) # must be integer for looping
    dim=numareas*len(Elemdata)# can't write to non-existant df row so set # of rows as numareas*len(Elements)
    mycols=['Project','Filepath','Date','Sample','Filename','Filenumber','Areanumber','Peakenergy','Peakindex','PeakID','Shift',
    'Negintensity','Posintensity','Pospeak','Amplitude','Peakwidth','Lowback','Lowbackamplitude','Highback','Highbackamplitude',
    'Avgbackamplitude','Quantdetails']    
    Smdifpeaks=pd.DataFrame(index=np.arange(0,dim), columns=mycols) # creates a temporary dataframe with all columns necessary for holding smooth-diff peak data
    dfrow=0 # keep track of row # for smdifpeaks dataframe
    Augerfile=pd.read_csv(AugerFileName) # read entire spectra into df and then slice
    
    # Since column structure in CSV is not 100% certain, extract each area's data by name (e.g. counts1, S7D71, counts2)
    # add error handling in case converted csv doesn't exist (unlikely since created with batch run) 
    
    for areanum in range(1,numareas+1): # loop over each separate area in spe
        # energy column not necessary (already have index #s)
        colname='S7D7'+str(areanum)
        thisS7D7=Augerfile[colname] # pulls smooth-diff column from correct area (Series)
        # loop through all smooth-diff peaks for each element in this spatial area
        for i, (elem, negpeak, pospeak, width, kfact) in enumerate(Elemdata):
            peakregion=thisS7D7[negpeak-width:negpeak+width+1] # select region around peak with specified search width from AESquantparams.csv
            negval, negindex=findnegpeak(peakregion) 
            peakregion=thisS7D7[pospeak-width:pospeak+width+1]            
            posval,posindex=findpospeak(peakregion) # return int and position (via shift from ideal)
            # calculate and assign all peak parameters to row in peak df
            Smdifpeaks.iloc[dfrow]['PeakID']=elem            
            Smdifpeaks.iloc[dfrow]['Peakindex']=negindex # index num of peak and not eV (negpeak is ideal index #)
            Smdifpeaks.iloc[dfrow]['Peakenergy']=Augerfile.iloc[negindex][0] # gets associated energy of true peak from df imported from csv
            Smdifpeaks.iloc[dfrow]['Areanumber']=areanum
            Smdifpeaks.iloc[dfrow]['Shift']=negindex-negpeak # ideal position minus found position (in eV assuming 1eV/channel)
            Smdifpeaks.iloc[dfrow]['Negintensity']=negval
            Smdifpeaks.iloc[dfrow]['Posintensity']=posval
            Smdifpeaks.iloc[dfrow]['Pospeak']=Augerfile.iloc[negindex][0]-(negindex-posindex) # energy val of corresponding positive peak
            Smdifpeaks.iloc[dfrow]['Amplitude']=posval-negval
            Smdifpeaks.iloc[dfrow]['Peakwidth']=negindex-posindex
            # construct quant details string to keep track of assorted smooth-diff quant params
            now=datetime.datetime.now()
            thisdate=now.strftime('%m/%d/%Y')
            tempstring=thisdate+' '+ elem+' '+str(negpeak)+' '+str(negpeak)+' '+str(width)+' '+'%.2f' % round(kfact, 2)
            Smdifpeaks.iloc[dfrow]['Quantdetails']=tempstring
            dfrow=dfrow+1 # increment to next data row
            # find and assign associated energy
            
    # assign params that are common to all areas/all peaks into rows of df
    for i in range(0,dfrow):
        Smdifpeaks.iloc[i]['Project']=logmatch.Project
        Smdifpeaks.iloc[i]['Filepath']=logmatch.FilePath
        Smdifpeaks.iloc[i]['Date']=logmatch.Date
        Smdifpeaks.iloc[i]['Sample']=logmatch.Sample
        Smdifpeaks.iloc[i]['Filename']=logmatch.Filename
        Smdifpeaks.iloc[i]['Filenumber']=logmatch.Filenumber   
    return Smdifpeaks # df with smooth-diff peaks for all areas/ all elements

def smdifbackground(logmatch,Backregdata):
    ''' Searches peak-free regions in background of smooth-diff spectrum and returns the largest negative peak with 10 channels of energy regions chosen in Backregs
    returns energy value (in eV not index #) and associated peak-to-peak amplitude of background fluctuations  
    since just returning magnitude of noise fluctuation, pospeak is not required to be lower energy than negpeak
    '''
    searchwidth=6 # default of +/-6 channels for background peak search is similar to most elemental peaks
    pospeakwidth=12 # search for associated positive "peak" with 15 channels/eV to left of negative peak
    AugerFileName=logmatch.Filename  # get Auger filename from Series... should be .csv file in current directory
    numareas=int(logmatch.Areas) # number of spatial areas within spe... must be integer for looping
    dim=numareas*len(Backregdata)# can't write to non-existant df row so set # of rows as numareas*len(Elements)
    mycols=['Project','Filepath','Date','Sample','Filename','Filenumber','Areanumber','Peakenergy','Peakindex','PeakID','Shift',
    'Negintensity','Posintensity','Pospeak','Amplitude','Peakwidth','Lowback','Lowbackamplitude','Highback','Highbackamplitude',
    'Avgbackamplitude','Quantdetails']    
    Smdifback=pd.DataFrame(index=np.arange(0,dim), columns=mycols) # dataframe of same type for background region "peaks" (aka noisy fluctuations) in smooth-diff
    dfrow=0 # keep track of row # for smdifpeaks dataframe
    Augerfile=pd.read_csv(AugerFileName) # read entire spectra into df, slice it later
    evbreaks=logmatch.Evbreaks # needed to avoid endpoints or energy breaks during peak searches
    if type(evbreaks)==str: # create list of ints from loaded string of format [0, 115, 230, 361]
        tempstring=evbreaks.split('[')[1]
        tempstring=tempstring.split(']')[0] # remove brackets
        evbreaks=[int(s) for s in tempstring.split(',')] # turn string into list of ints 
        
    for areanum in range(1,numareas+1): # loop over each separate area in spe
        # energy column not necessary (already have index #s)
        colname='S7D7'+str(areanum)
        thisS7D7=Augerfile[colname] # pulls smooth-diff column from correct area (Series)
        # loop through all smooth-diff peaks for each element in this spatial area
        for i, (energy, index) in enumerate(Backregdata):
            lower=index-searchwidth
            upper=index+searchwidth+1
            # ensure that background regions do not cross a energy break in multiplex scans (which could induce goofy peak)
            lower, upper = checkpeakboundaries(lower,upper, index, evbreaks)
            peakregion=thisS7D7[lower:upper] # select region around peak using index #s not eV
            try:
                negval, negindex=findnegpeak(peakregion) # returns neg peak amplitude and index of position
                Smdifback.iloc[dfrow]['Shift']=negindex-index # shift from specified position (probably just random)
                Smdifback.iloc[dfrow]['Negintensity']=negval
                Smdifback.iloc[dfrow]['Peakindex']=negindex # index num of peak and not eV (negpeak is ideal index #)
                Smdifback.iloc[dfrow]['Peakenergy']=Augerfile.iloc[negindex][0] # associated energy of true peak
            except:
                print ('No negative peak found for indices', lower, 'to', upper)
            # find associated positive "peak" (no ideal position since these are background fluctuations)
            lower=negindex-pospeakwidth
            upper=negindex
            lower, upper = checkpeakboundaries(lower,upper, index, evbreaks)
            peakregion=thisS7D7[lower:upper]
            try:            
                posval,posindex=findpospeak(peakregion) # return int and position (via shift from ideal)
                Smdifback.iloc[dfrow]['Posintensity']=posval
            except:
                print ('No positive peak found for indices', lower, 'to', upper) 
            # calculate and assign all peak parameters to row in peak df
            Smdifback.iloc[dfrow]['Areanumber']=areanum
            try:
                Smdifback.iloc[dfrow]['Pospeak']=Augerfile.iloc[negindex][0]-(negindex-posindex) # energy val of corresponding positive peak
                Smdifback.iloc[dfrow]['Amplitude']=abs(posval-negval) # just noise measure so can be negative
                Smdifback.iloc[dfrow]['Peakwidth']=negindex-posindex
            except:
                print('Problem with noise negpeak or pospeak.')
            # construct quant details string to keep track of assorted smooth-diff quant params
            dfrow=dfrow+1 # increment to next data row
            # find and assign associated energy  
    return Smdifback # df with smooth-diff peaks for all areas/ all elements

def smdifbatchquant(myfiles, Elements, AESquantparams, Smdifpeakslog, Backregs=[121,200,405,800,1475,1850], **kwargs):
    ''' Batch quantification of all peaks in Elements list and noise amplitude at all chosen background regions (Backregs) 
    returns df with peak positions, amplitudes, width, energy shift, etc.
    kwargs:
        reprocess: flag to redo all deriv quant (default is False)
    '''   
    # Backregs was an attempt to monitor signal to noise in differentiated spectra but it's not very good
    # often insufficient energy range to do this well in multiplex spectra
    # should replace this with Ogliore/Gazda 2015 error estimation methods 
    Elemdata=[] # list of tuples with elem symbol and other params needed to find smoothed-differentiated peak in given spectrum
    Backregdata=[] # list of energy value for background peak checks and associated index value (for easier slicing)
    if kwargs.get('reprocess', False):
        # do not keep any prior data runs if reprocess all set to true
        Smdifpeakslog=pd.DataFrame() # empty dataframe to hold all peaks, all spe files
    for index, row in myfiles.iterrows():
        # get ith row from parameters log for subset of selected spe files (i.e. from myfiles)
        logmatch=myfiles.loc[index] #contains row with filename and all other parameters from a given spectra 
        logmatch=logmatch.squeeze() # convert/flatten to Series
        # work with indices or directly with energy vals?  indices probably faster for smooth-deriv method
        filename=logmatch.Filename
        if not os.path.isfile(filename):
            print ('File ', filename, 'not found.')
            continue
        if not kwargs.get('reprocess', False): # check if already present in smdiflog, if so skip reprocessing
            try:
                match=Smdifpeakslog[Smdifpeakslog['Filename']==filename]
                if len(match)!=0:
                    priorelem=np.ndarray.tolist(match.PeakID.unique())
                    missingelems=[elem for elem in Elements if elem not in priorelem]
                    print('Entire file ', filename, 'skipped/already processed... ensure same AESquantparams ')
                    if len(missingelems)>0:
                        missingstr=",".join(missingelems)
                        print('Elements ', missingstr, ' not processed for ', filename)
                    continue # skip entire spe file on the assumption that it's already been processed
            except: # handles empty dataframe problem (no prior run)
                pass
        # Elemdata has index #s of ideal positions for this energy list
        Elemdata, Backregdata=findindices(Elements, Backregs, logmatch, AESquantparams) # passes ideal position
        if len(Elemdata)==0:
            print('Problem finding indices for ', filename)
            continue
        # S7D7 quant routine pass 1) File params pulled from log as Series (incl. spectra filename) 2) element list for quant 
        Smdifpeaks=smdifquant(logmatch, Elemdata) # temporary df for smooth-differentiated quant results from each peak of single spe file
        
        # Smdifback is a temp dataframe for noise amplitude in spe at various energy regions (selected by eV in Backregs list)
        Smdifback=smdifbackground(logmatch,Backregdata) # Add appropriate noise peak magnitudes to each elemental peak (use energy region near each peak's energy)
        # find and append background region above and below peak and append (also test to ensure they exist)
        Smdifpeaks=addbacktopeaks(Smdifpeaks, Smdifback) # for each row/peak in SmdifQuantresult find and insert nearest background regions, overwrite existing df with more complete one
        
        # now that dataset is complete append smoothdiff data from single spe to longer master list
        Smdifpeakslog=Smdifpeakslog.append(Smdifpeaks)
        # Set ignore_index?  reindexes the whole df each time
        print(filename, 'processed.')
    Smdifpeakslog=Smdifpeakslog.sort_values(['Filenumber','Areanumber'])
    # Unlikely to occur, but drop duplicates of exact same peak
    Smdifpeakslog=Smdifpeakslog.drop_duplicates(subset=['Filenumber', 'Areanumber','PeakID'])
    mycols=['Project','Filepath','Date','Sample','Filename','Filenumber','Areanumber','Peakenergy','Peakindex','PeakID','Shift',
    'Negintensity','Posintensity','Pospeak','Amplitude','Peakwidth','Lowback','Lowbackamplitude','Highback','Highbackamplitude',
    'Avgbackamplitude','Quantdetails']
    Smdifpeakslog=Smdifpeakslog[mycols] # reorder and drop extra columns (may drop adjamp but easily recalculated)
    # Go ahead and run calc amplitude 
    Smdifpeakslog=calcamplitude(Smdifpeakslog, AESquantparams)
    return Smdifpeakslog # not autosaved (optional autosave in tk gui)
