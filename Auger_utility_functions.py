# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:11:09 2016
Assorted Auger utility functions 

@author: tkc
"""
#%%
import pandas as pd
import numpy as np
import os, glob,sys, shutil, re # already run with functions 
import datetime
#%%
''' RUNNING BELOW FILE MANAGEMENT UTILITIES
destdir='H:\\Research_data\\Stardust\\C2010W\\Auger\\18Nov15\\sub\\'
movefiles(spelist,destdir)

excludelist=combinelist[combinelist['Areas']==2] # finds filenames with screwed up import
AugerParamLog=removefromlog(excludelist, AugerParamLog) # removes log entries from any dataframe by filename if in excludelist

AugerParamLog=removelistdups(AugerParamLog,'Evbreaks')
'''

def loadsubfiles():
    ''' Load of standard sub spe files (before combination via averaging from working Auger data directory '''
    if os.path.isfile('sub\\Smdifpeakslog_subs.csv'):
        Smdifpeakslogsubs=pd.read_csv('sub\\Smdifpeakslog_subs.csv', encoding='cp437')
    else:
        print('Smdifpeakslogsubs not found.')
        Smdifpeakslogsubs=pd.DataFrame()
    if os.path.isfile('sub\\Backfitlog_subs.csv'):
        Backfitlogsubs=pd.read_csv('sub\\Backfitlog_subs.csv', encoding='cp437')
    else:
        print('Backfitlogsubs not found.')
        Backfitlogsubs=pd.DataFrame()
    if os.path.isfile('sub\\Integquantlog.csv'):
        Integquantlogsubs=pd.read_csv('sub\\Integquantlog_subs.csv', encoding='cp437')
    else:
        print('Integquantlogsubs not found.')
        Integquantlogsubs=pd.DataFrame()
    return Smdifpeakslogsubs, Backfitlogsubs, Integquantlogsubs

def loadmainfiles():
    ''' Load of standard files from main Auger data directory '''
    if os.path.isfile('Augerparamlog.csv'):
        AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
        spelist=AugerParamLog[(AugerParamLog['Areas']>=1)] 
        excludemask=spelist['Comments'].str.contains('exclude', case=False, na=False)
        spelist=spelist.loc[~excludemask]
        spelist=spelist.sort_values(['Filenumber'], ascending=True)
    else:
        print('Augerparamlog not found.')
        AugerParamLog=pd.DataFrame()
        spelist=pd.DataFrame()
    if os.path.isfile('Smdifpeakslog.csv'):
        Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437')
    else:
        print('Smdifpeakslog not found.')
        Smdifpeakslog=pd.DataFrame()
    if os.path.isfile('Integquantlog.csv'):
        Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
    else:
        print('Integquantlog not found.')
        Integquantlog=pd.DataFrame()
    if os.path.isfile('Backfitlog.csv'):
        Backfitlog=pd.read_csv('Backfitlog.csv', encoding='cp437')
    else:
        print('Backfitlog not found.')
        Backfitlog=pd.DataFrame()
    if os.path.isfile('C:\\Users\\tkc\\Documents\\Python_Scripts\\AESquantparams.csv'):
        AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\AESquantparams.csv', encoding='cp437')
    else:
        print('AESquantparams not found.')
        AESquantparams=pd.DataFrame()
    return AugerParamLog, spelist, Smdifpeakslog, Integquantlog, Backfitlog, AESquantparams
    
def dropexcluded(df, spelist):
    '''Only returns subset of df with filenumber in list passed by second argument'''
    filelist=spelist.Filenumber.unique()
    df=df[df.Filenumber.isin(filelist)]
    return df
    
def writecomps(df, AESquantparams, Elems):
    '''Store composition results in auto-named, autosaved xls along with subset of elements used and 
    second tab with k-factors employed, lines of element used'''
    now=datetime.datetime.now()
    elemstr=Elems[0]+'ampl'
    if elemstr in df: # comp derived from deriv data
        savename='smdiffcomp_'+datetime.date.strftime(now, "%d%b%y")+'.xlsx'
        sheetname='smdiffcomp'
    else: # comp from integ method
        savename='integcomp_'+datetime.date.strftime(now, "%d%b%y")+'.xlsx'
        sheetname='integcomp'
    if os.path.exists(savename): # prompt for overwritting existing xls of same name
        print('Enter Y to overwrite file', savename)
        overwrite=input()
        if overwrite!='Y':
            print('Exiting without overwrite of existing file')
            return
    writer=pd.ExcelWriter(savename, engine='openpyxl', datetime_format='mm/dd/yy')
    # Write elements to comment (writer is )
    # elemstr='_'.join(Elems)
    # writer["A1"].comment=elemstr
    df.to_excel(writer,sheetname,index=False) # this overwrites existing file
    # write date and elements list into comment
    AESq=AESquantparams[AESquantparams.element.isin(Elems)]# subset of elements used in this quant (can be Fe or Fe+Fe2, etc.)
    AESq.to_excel(writer,'Kfactors', index=False)    
    writer.save() # saves xls file with all modified data
    return
    
def updatepaths(AugerParamLog):
    '''Checks parameter log to make sure path to Auger csv files in AugerParamLog is correct 
    single log holds those in dir (combined or unique files) and those in subdir (source files
    before average via combination)'''
    path=os.getcwd() #
    basepath='Research_data'+path.split('Research_data')[1]
    for index, row in AugerParamLog.iterrows():
        name=AugerParamLog.loc[index]['Filename']
        if 'csv' not in name: # jpgs derived from sem and map files will not be in //sub
            AugerParamLog=AugerParamLog.set_value(index,'FilePath',basepath)
            continue
        if os.path.exists(name):
            AugerParamLog=AugerParamLog.set_value(index,'FilePath',basepath)
        elif os.path.exists('sub\\'+name): # check for source file in /sub directory
            AugerParamLog=AugerParamLog.set_value(index,'FilePath',basepath+'\\sub')
        else: # file not found in directory or subdirectory (remove from AugerParamLog)            
            print(name,' removed from log... not found in dir or subdir')
            AugerParamLog.drop(AugerParamLog.index[index], inplace=True)
    return AugerParamLog
    

def assembledataset(paramloglist, integloglist, smdifloglist):
    '''Construct master paramlog, integlog, and smdiflog for list of directories 
    used to create master data sets for later processing, plotting, etc.'''
    # Structure of created files sometimes changes over time... ensure consistency
    mycols=['Filenumber', 'Project', 'Filename', 'FilePath', 'Sample', 'Comments',
       'Date', 'FieldofView', 'Type', 'Energy', 'GunVoltage', 'nA', 'Areas',
       'Cycles', 'Timestep', 'Details', 'Evbreaks', 'Acqtime', 'Scanarea', 'X',
       'Y', 'Z', 'Tilt', 'Rotation', 'ImageshiftX', 'ImageshiftY']
    Masterparamlog=pd.DataFrame(columns=mycols)
    # current df structure for integration results
    mycols2=['Filenumber', 'Filename', 'Filepath', 'Sample', 'Comments', 'Areanumber', 'Element', 'Integcounts', 
    'Backcounts', 'Significance', 'Xc', 'Width', 'Peakarea', 'Y0','Rsquared','Numchannels']
    Masterinteglog=pd.DataFrame(columns=mycols2) # empty frame
    # current df structure for smdiff quant results (see auger_smdiff_functions)
    mycols3=['Project','Filepath','Date','Sample','Filename','Filenumber','Areanumber','Peakenergy','Peakindex','PeakID','Shift',
    'Negintensity','Posintensity','Pospeak','Amplitude','Peakwidth','Lowback','Lowbackamplitude','Highback','Highbackamplitude',
    'Avgbackamplitude','Quantdetails']
    Mastersmdiflog=pd.DataFrame(columns=mycols3) # empty smdif frame
    for i, logfile in enumerate(paramloglist):
        thisparam=pd.read_csv(logfile, encoding='cp437')
        # check for and print differences between sets of columns
        cols1=list(thisparam.columns.values)
        uniquecols=[col for col in cols1 if col not in mycols]
        if len(uniquecols)>0:
            print('Unique columns in ', logfile,': ', ','.join(uniquecols))
        missingcols=[col for col in mycols if col not in cols1]
        if len(missingcols)>0:
            print('Missing columns in ', logfile,': ', ','.join(missingcols))
        Masterparamlog=pd.concat([Masterparamlog,thisparam], ignore_index=True)
    for i, logfile in enumerate(integloglist):
        thisinteg=pd.read_csv(logfile, encoding='cp437')
        cols1=list(thisinteg.columns.values)
        uniquecols=[col for col in cols1 if col not in mycols2]
        if len(uniquecols)>0:            
            print('Unique columns in ', logfile,': ', ','.join(uniquecols))
        missingcols=[col for col in mycols2 if col not in cols1]
        if len(missingcols)>0:
            print('Missing columns in ', logfile,': ', ','.join(missingcols))
        Masterinteglog=pd.concat([Masterinteglog,thisinteg], ignore_index=True)
    for i, logfile in enumerate(smdifloglist):
        thissmdif=pd.read_csv(logfile, encoding='cp437')
        cols1=list(thissmdif.columns.values)
        uniquecols=[col for col in cols1 if col not in mycols3]
        if len(uniquecols)>0:
            print('Unique columns in ', logfile,': ', ','.join(uniquecols))
        missingcols=[col for col in mycols3 if col not in cols1]
        if len(missingcols)>0:
            print('Missing columns in ', logfile,': ', ','.join(missingcols))
        Mastersmdiflog=pd.concat([Mastersmdiflog,thissmdif], ignore_index=True)
    return Masterparamlog, Masterinteglog, Mastersmdiflog

def checkparamlog(AugerParamLog, makeentry=False):
    ''' Checks the Auger parameters logbook against actual csv spectral files, correct path if necessary
    prints out filenumbers that have a problem to console
    ''' 
    #TODO fix this so it works with new blanklog containing filenames (due to degenerate filenumber problem)
    # find csv file in current and /sub directory 
    spelist=[] # list of filenumbers of spe files in directory
    sublist=[] # csv files that are in sub directory
    allcsvfiles=glob.glob('**/*.csv', recursive=True)
    for i, name in enumerate(allcsvfiles): # deals with 3 cases (spe, sem or map)
        if "sub\\" in name:
                tempstring=name.split('sub\\')[1]
                match=re.search(r'.\d+.csv',tempstring) # Auger spectra should have format like filename.100.csv
                if match: # check if it's an Auger spectrum                 
                    sublist.append(tempstring)
        else:            
            match=re.search(r'.\d+.csv',name) # Auger spectra should have format like filename.100.csv
            if match:
                spelist.append(name) # add to spe files in dir list
    # Filter AugerParamLog into spe files in main vs those in sub
    spelog=AugerParamLog[(AugerParamLog['Areas']>=1)]
    excludemask=spelog['FilePath'].str.contains('sub', case=False, na=False)
    sublogfiles=spelog.loc[excludemask]
    spelog=spelog.loc[~excludemask]
    logfilelist=spelog.Filename.unique()
    logfilelist=np.ndarray.tolist(logfilelist)
    sublogfilelist=sublogfiles.Filename.unique()
    sublogfilelist=np.ndarray.tolist(sublogfilelist)
    # compare spes in main
    missingdata=[i for i in logfilelist if i not in spelist] # list comprehension for missing data file (but present in excel logbook)
    missingentry=[i for i in spelist if i not in logfilelist]     
    for i, val in enumerate(missingdata):
        print ('Data file number ', val, ' mentioned in AugerParamlog but missing from directory')
    for i, val in enumerate(missingentry):
        print ('Data file number ', val, ' present in directory but missing from AugerParamlog')
    missingdata=[i for i in sublogfilelist if i not in sublist] # list comprehension for missing data file (but present in excel logbook)
    missingentry=[i for i in sublist if i not in sublogfilelist]     
    for i, val in enumerate(missingdata):
        print ('Data file number ', val, ' mentioned in AugerParamlog but missing from directory')
    for i, val in enumerate(missingentry):
        print ('Data file number ', val, ' present in directory but missing from AugerParamlog')
    if makeentry: # not yet implemented and shouldn't happen... fix manually for now?
        pass
        # AugerParamLog=makemissingentries(missingentry, AugerParamLog)
    return AugerParamLog

def smooth7cnts(df, areanum, evbreaks):
    '''create smooth differentiated column from counts using S7D7 PHI algorithm (Multipak tables A-5 and A-1
    version for rerun on combined spectra with internal ev breaks''' 
    countname='Counts'+str(areanum)
    smcountname='Smcounts'+str(areanwum)
    counts=df[countname].tolist() # convert to list
    numpts=len(counts)
    smooth=[0]*numpts # empty list of correct length for smoothed data
    # smoothing of endpoints according to Multipak algorithm appendix table A-5    
    for i in range(0,numpts): # special cases for endpoints (within 3 of an evbreak)
        diff=i-min(evbreaks, key=lambda x:abs(x-i)) # distance from closest evbreak index # in list            
        if diff==0:
            if i==numpts-1: #last point
                smooth[i]=(2*counts[i]+2*counts[i-1]+1)/4 # additional special case for last point
            else: # first point
                smooth[i]=(2*counts[i]+2*counts[i+1]+1)/4 # all others at exact breaks can use value and adj higher value
        elif abs(diff)==1:  # works for +1 or -1 from nearest break
            smooth[i]=(1*counts[i-1]+2*counts[i]+1*counts[i+1]+1)/4
        elif abs(diff)==2:
            smooth[i]=(-3*counts[i-2]+12*counts[i-1]+17*counts[i]+12*counts[i+1]+-3*counts[i+2]+1)/35
        else:
            smooth[i]=(-2*counts[i-3]+3*counts[i-2]+6*counts[i-1]+7*counts[i]+6*counts[i+1]+3*counts[i+2]-2*counts[i+3]+1)/21
    df[smcountname]=smooth # add smoothed data as new dataframe column
    return df

def addsmoothloop(spelist):
    ''' Add smoothed counts data column for each area '''   
    for i in range(0,len(spelist)):
        # get ith row from parameters log for subset of selected spe files (i.e. from spelist)
        logmatch=spelist.iloc[i] #contains row with filename and all other parameters from a given spectra 
        logmatch=logmatch.squeeze() # convert/flatten to Series
        numareas=int(logmatch.Areas) # get # of spatial areas for this spe
        # load Auger spe file of interest here
        AugerFileName=logmatch.Filename  # get Auger filename from Series
        Augerfile=pd.read_csv(AugerFileName) # read entire spectra into df
        savefile=False
        if 'Smcounts1' not in Augerfile: # same evbreaks for all areas
            savefile=True
            evbreaks=logmatch.Evbreaks # needed for possible savgol smooth-diff
            tempstring=evbreaks.split('[')[1] # remove brackets from list
            tempstring=tempstring.split(']')[0]
            evbreaks=[int(s) for s in tempstring.split(',')] # convert string to list of break index values
        # now loop through any areas within this spectrum (typically only 1 area)
        for areanum in range(1,numareas+1): # loop over each separate area in spe
            # Now check to ensure this Augerfile has all necessary columns for this area
            # print('Processing area ', areanum)  TESTING
            colname='Counts'+str(areanum)
            if colname not in Augerfile:# this shouldn't happen
                print(colname, ' not present in file ', AugerFileName)
                continue # skip to next area
            smoothname='Smcounts'+str(areanum) # Sav-gol 2nd deriv column used to guide selection of fitting regions
            if smoothname not in Augerfile: # returns df with this Savgol column added
                Augerfile=smooth7cnts(Augerfile, areanum, evbreaks)  # FUNCT pass full spectrum for given area (saved below)   
        # direct save of modified auger csv with new linear background fits (after all areas processed)
        if savefile: # only save if new smoothed counts have been added
            Augerfile.to_csv(AugerFileName, index=False)
    return
    
def compareduplicatecomps(df, elemlist):
    '''Return subset with same sample name but different file number (avoids multiple spatial areas from same crater which may be heterogeneous)'''
    samplelist=df.Sample.unique()
    samplelist=np.ndarray.tolist(samplelist)
    # make dataframe for output
    mycols=['Sample', 'Filenumber', 'Filename','Areanumber', 'AESbasis']
    commoncols=['Sample', 'Filenumber', 'Filename','Areanumber', 'AESbasis']
    for i, elem in enumerate(elemlist):
        mycols.append('%'+elem)
        mycols.append('err%'+elem)
        mycols.append('minerr'+elem)
        mycols.append('maxerr'+elem)
    duplicatecompdata=pd.DataFrame(columns=mycols)
    for i, sample in enumerate(samplelist):
        # check for multiple spatial areas and process separately
        match=df[df['Sample']==sample]
        arealist=match.Areanumber.unique()
        arealist=np.ndarray.tolist(arealist)
        for j, areanum in enumerate(arealist):
            thissample=pd.DataFrame(index=np.arange(0,1),columns=mycols) # single row for this sample
            match2=match[match['Areanumber']==areanum]
            for j, elem in enumerate(elemlist):
                colname='%'+elem
                vals=match2[colname].unique() # gets atomic percents
                thissample=thissample.set_value(0,'%'+elem,vals.mean()) # set mean value as at. %
                thissample=thissample.set_value(0,'err%'+elem,vals.std()) # set stdev as err in at. %
                colname='err%'+elem # Same for errors 
                vals=match2[colname].unique() # gets error values in at % calculations (in case any are fubar)
                thissample=thissample.set_value(0,'minerr'+elem,vals.min())
                thissample=thissample.set_value(0,'maxerr'+elem,vals.max())
            # Copy value or values for each common parameter 
            for j,colname in enumerate(commoncols):
                templist=match2[colname].unique()
                templist=np.ndarray.tolist(templist)
                templist=set(templist) # eliminate duplicates
                templist=list(templist)
                if colname=='Filenumber' or colname=='Areanumber':
                    templist=[int(i) for i in templist] # convert floats to ints
                templist=[str(s) for s in templist] # convert to string for output into df
                templist=', '.join(templist) # convert to string
                thissample=thissample.set_value(0,colname,templist)
            duplicatecompdata=duplicatecompdata.append(thissample) # add to multirowed df
    return duplicatecompdata	

def findduplicates(df):
    '''Return subset with same sample name but different file number (avoids multiple spatial areas from same crater which may be heterogeneous)'''
    samplelist=df.Sample.unique()
    samplelist=np.ndarray.tolist(samplelist)
    keeplist=[] # index # of 
    for i, sample in enumerate(samplelist):
        match=df[df['Sample']==sample]
        temp=match.Filenumber.unique()
        if len(temp)>1: # keep them
            indlist=match.index.unique()
            indlist=np.ndarray.tolist(indlist)
            keeplist.extend(indlist) # add these to list of values to keep
    dups=df[df.index.isin(keeplist)] # actual duplicates list 
    return dups		
		
def keepbestvals(df, spectralregs):
    '''For duplicate energy values, keep the one with largest # of sweeps '''
    # generate temporary index ranges for spectralregs (can't just use eneryg values)
    start=0
    spectralregs['Start']=0 # for index #s
    spectralregs['End']=0
    for index,row in spectralregs.iterrows():
        lower=spectralregs.loc[index]['Lower']
        upper=spectralregs.loc[index]['Upper']
        thisrange=int(upper-lower)
        spectralregs=spectralregs.set_value(index,'Start',start)
        spectralregs=spectralregs.set_value(index,'End',start+thisrange)
        start=start+thisrange+1 # adjust for next loop
    dupl=df.duplicated(['Energy'], keep=False) # only duplicate vals and keep both
    dupl=df.loc[dupl]
    energyvals=dupl.Energy.unique()
    energyvals=np.ndarray.tolist(energyvals)
    removelist=[]
    for i, val in enumerate(energyvals):
        thismatch=dupl[dupl['Energy']==val]
        if len(thismatch)!=2:
            print('Unknown error in energy duplicate elimination')
            continue
        else: # pick counts value with highest number of sweeps (best value)
            index1=thismatch.index[0]
            index2=thismatch.index[1]
            specmatch1=spectralregs[(spectralregs['Start']<=index1)&(spectralregs['End']>=index1)]
            specmatch2=spectralregs[(spectralregs['Start']<=index2)&(spectralregs['End']>=index2)]
            
            try:
                if specmatch1.iloc[0]['Sweeps']>=specmatch2.iloc[0]['Sweeps']:
                    # first is best value... remove it from dupl (which will be used as knockout df)
                    removelist.append(index1)
                else:
                    removelist.append(index2)
            except:
                print('Problem with duplicate removal of ', index1, index2)
    df=df[-df.index.isin(removelist)]
    print (len(removelist), ' duplicated energy values removed from multiplex')
    return df

def sortmultiplex(multiplex, evbreaks):
    ''' Rearrange multiplex if taken out of order and adjust evbreaks accordingly
    evbreaks holds index #s of multiplex breaks so must be altered along with multiplex sort
    only called if multiplex is not monotonically increasing '''
    energylist=[]
    for i, val in enumerate(evbreaks):
        energylist.append(multiplex.loc[val]['Energy'])
        if i>0 and i<len(evbreaks)-1:
            energylist.append(multiplex.loc[val+1]['Energy']) # starting ev for next spectral region
    multiplex=multiplex.sort_values(['Energy'])
    multiplex=multiplex.drop_duplicates(['Energy'])
    multiplex=multiplex.reset_index(drop=True)
    matches=multiplex[multiplex['Energy'].isin(energylist)]
    # Some can be internal data breaks if multiplex regions overlap
    newevbreaks=matches.index.tolist()
    newevbreaks=[int(i) for i in newevbreaks if i-1 not in newevbreaks] # remove adjacent values
    return multiplex, newevbreaks

def avglogentry(logmatches):
    ''' create a logbook entry for a combine-averaged spectrum '''    
    firststr=str(logmatches.iloc[0][0]) # averaged file name is numbered "firstnum"lastnum"
    laststr=str(logmatches.iloc[len(logmatches)-1][0])
    avgnum=int(firststr+laststr)
    avgentry=logmatches.iloc[[0]] # avg entry has same values as first file for cols 1,3-15, 17-24
    avgentry=avgentry.squeeze() # convert to series (avoids warnings regarding changing a slice of the df)
    # now change the small number of header params that will differ between single spe and averaged one    
    avgentry.iloc[0]=avgnum # assign filenum for averaged as firstnum-lastnum 
    tempval=avgentry.iloc[2]
    tempval=tempval.replace(firststr+'.',firststr+laststr+'.') # avoids error if number is elsewhere in name string
    avgentry.iloc[2]=tempval # reassign to correct filename
    avgentry.iloc[16]*=len(logmatches) # multiply acquisition time by number of files
    return avgentry # pandas series with all params and of correct dimension for append to main df log

def makecomboentries(combinelist, AugerParamLog):
    ''' Make log entry for average combined file using this loop and avglogentry function
    for index, row in combinelist.iterrows(): # large loop through each match
    this is normally done within combinespeloop '''
    for index, row in combinelist.iterrows(): # large loop through each match
        firstfile=int(combinelist.loc[index]['Filenumber']) # iloc grabs correct row from list of matches
        lastfile=int(combinelist.loc[index]['Lastnumber'])
        # look up filenames and # of areas associated with the above file range (all should have same # of area)
        logmatches=AugerParamLog[(AugerParamLog['Filenumber']>=firstfile) & (AugerParamLog['Filenumber']<=lastfile)]
        # Start of function combinespectra(logmatches) ... return new line for AugerFileParams logbook
        # new entry into Augerparamlog  
        avgentry = avglogentry(logmatches) # create new logbook entry (series) for averaged spectrum (mostly from firstfile's info)
        # append new Series entry to end of AugerParamLog
        AugerParamLog=AugerParamLog.append(avgentry)
    return AugerParamLog

def checklog(filelist, AugerParamLog):
    ''' Pass list of csv files in directory and checks the user Auger parameter log matches the actual data file list from directory
    prints out filenumbers that have a problem to console'''   
    spelog=AugerParamLog[(AugerParamLog['Areas']>=1)]
    loglist=spelog.Filename.unique()
    loglist=np.ndarray.tolist(loglist)
    missingdata=[i for i in loglist if i not in filelist] # list comprehension for missing data file (but present in excel logbook)
    missingentry=[i for i in filelist if i not in loglist] # list comprehension for data file with missing log entry (could also convert to sets and compare)
    for i, val in enumerate(missingdata):
        print ('Data file number ', val, ' present in Auger params log but missing from directory')
    for i, val in enumerate(missingentry):
        print ('Data file number ', val, ' present in directory but missing from Auger params log')
    # check for duplicate entries in logbook
    myset=set([x for x in loglist if loglist.count(x) > 1])
    for i in myset:
        print('Duplicate entry for file number', i, ' in Auger params log.')
    return

def removelistdups(df, colname):
    ''' Remove duplicates from any list stored within df column and return '''
    for index,row in df.iterrows():
        thisstr=df.loc[index][colname]
        if str(thisstr)!='nan':
            tempstring=thisstr.split('[')[1]
            tempstring=tempstring.split(']')[0] # remove brackets
            evbreaks=[int(s) for s in tempstring.split(',')] # turn string into list of ints
            evbreaks=set(evbreaks)
            evbreaks=list(evbreaks)
            evbreaks.sort()
            df=df.set_value(index,colname,evbreaks) # write it back
    return df

def deletefiles(df):
    ''' Grab csv of each file in current df (e.g. after slicing) and copy to select subfolder 
    also save current filelog as AugerParamLog after path modification'''        
    for index, row in df.iterrows():
        filename=df.loc[index]['Filename']
        if os.path.exists(filename):
            os.remove(filename)
            print(filename, 'deleted.')
    return

def movefiles(df, destdir):
    ''' Moves files named in dataframe to selected location '''
    for index, row in df.iterrows():
        filename=df.loc[index]['Filename']
        if os.path.exists(filename):
            newname=destdir+filename
            shutil.move(filename, newname)
            print(filename, 'moved to ', destdir)
        else: 
            print(filename, 'not found.')
    return


def copyexclusion(df1,df2):
    '''If AugerFileParams shows exclusion of an entire spectrum, copy 'excluded' into comments of other logs 
    for each peak '''
    #TODO finish this comment copying function
    return


def copyproblemcomments(df1,df2):
    '''If there's a problem noted in comments of integlog, copy over to smdiflog (or vice versa) '''
    problemlist=df1['Comments'].str.contains('problem', case=False, na=False)
    problems=df1.loc[problemlist]
    if 'Comments' not in df2:
        df2['Comments']='' # add comments col (sometimes missing from smdifpeakslog)
    for index, row in problems.iterrows():
        filenum=problems.loc[index]['Filenumber']
        elem=problems.loc[index]['Element']
        comment=problems.loc[index]['Comments']
        match=df2[(df2['Filenumber']==filenum) & (df2['PeakID']==elem)]
        if len(match)==1: # found associated peak
            # append or add comment to other log
            if str(match.iloc[0]['Comments'])=='' or str(match.iloc[0]['Comments'])=='nan':
                # set value directly for df2 using index of match
                df2=df2.set_value(match.index[0],'Comments',comment)
            else:
                currcomm=match.iloc[0]['Comments']
                newcomm=comment+' '+currcomm
                df2=df2.set_value(match.index[0],'Comments',newcomm)
    return df2
                


def removefromlog(df1, df2):
    '''If filename is in df1, remove these entries from df2.. then manually save ''' 
    removal=df1.Filename.unique()
    removal=np.ndarray.tolist(removal)
    df2=df2[-df2.Filename.isin(removal)] # inverse of isin
    return df2

def dropcolumns(df1,df2):
    ''' Pass two dfs with df2 being the template.. extra unnecessary columns dropped from df1
    inplace=True modifies both passed and returned df  '''
    cols1=df1.columns.tolist()
    cols2=df2.columns.tolist()
    newdf=df1 # avoids modification of passed df
    uniquelist=[i for i in cols1 if i not in cols2]
    for i,colname in enumerate(uniquelist): # remove cols from df1 that are absent from df2
        # newdf.drop(colname, axis=1, inplace=True) # this modifies both passed and returned dfs
        newdf=newdf.drop(colname, axis=1)
    return newdf
    
def outputduplicates(df, colname):
    '''Pass df and column name to search for duplicates (as string); outputs duplicates to console'''
    df=df.reset_index(drop=True)
    tempdf=df.duplicated([colname], keep=False) # marks all duplicates as true (even first last)
    # get list of duplicated filenumber
    for i in range(0,len(tempdf)):
        if tempdf[i]==True:
            print('Duplicated spectrum for sample: ', df.iloc[i]['Sample'])
    return    

def copyselectfiles(df,foldername):
    ''' Grab csv of each file in current df (e.g. after slicing) and copy to select subfolder 
    also save current filelog as AugerParamLog after path modification'''    
    
    if not os.path.exists(foldername): # create subdirectory for raw spe/sem/map files & csv sub-files (when combined)
        os.makedirs(foldername)
    # just use current path 
    curpath=df.iloc[0]['FilePath'] # should be same for entire folder
    newpath=curpath + "\\" + foldername
    df=df.replace(curpath, newpath) # changes all instances of current path in df
    for i in range(0,len(df)):
        filename=df.iloc[i]['Filename']
        newname=foldername + "\\" + filename
        shutil.copy2(filename, newname) # copy file
    # Save parameters list for this subset
    paramfile=foldername+ "\\" +'Augerparamlog.csv'
    df.to_csv(paramfile,index=False) # save to sub
    return

def truncatenumber(AugerParamLog,Smdifpeakslog):
    ''' Shortens long combo file name to 3 digit version'''
    for i in range(0,len(AugerParamLog)):
        val=str(AugerParamLog.iloc[i]['Filenumber'])
        if len(val)==6 or len(val)==7:
            newval=int(val[0:2])
        if len(val)==8:
            newval=int(val[0:3])
    AugerParamLog.iloc[i]['Filenumber']=newval
    for i in range(0,len(AugerParamLog)):
        val=str(Smdifpeakslog.iloc[i]['Filenumber'])
        if len(val)==6 or len(val)==7:
            newval=int(val[0:2])
        if len(val)==8:
            newval=int(val[0:3])
    Smdifpeakslog.iloc[i]['Filenumber']=newval
    return AugerParamLog, Smdifpeakslog    
