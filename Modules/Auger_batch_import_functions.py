# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 12:51:54 2016

@author: tkc
"""
import re, struct, os, glob # already run with functions
import pandas as pd
from PIL import Image, ImageDraw, ImageFont # needed for jpg creation
import numpy as np # needed for image arrays
# no need to clean or find/replace from .spe or .sem file headers (mostly : delimited)
from math import factorial
from skimage.feature import register_translation

def findshift(imagenum1, imagenum2, paramlog):
    ''' Pass pre and post- images, determine stage drift in microns and any uncorrected pixel shift and print it
    error is returned from register_translation '''
    match=paramlog[paramlog['Filenumber']==imagenum1]
    
    if len(match)==1:
        fname1=match.iloc[0]['Filename'].replace('.sem','.jpg')
        image1=Image.open(fname1)
        imshiftx1=match.iloc[0]['ImageshiftX']
        imshifty1=match.iloc[0]['ImageshiftY']
        FOV=match.iloc[0]['FieldofView']
    match=paramlog[paramlog['Filenumber']==imagenum2]
    if len(match)==1:
        fname2=match.iloc[0]['Filename'].replace('.sem','.jpg')
        image2=Image.open(fname2)
        imshiftx2=match.iloc[0]['ImageshiftX']
        imshifty2=match.iloc[0]['ImageshiftY']
    # for scikit, image sections need to be numpy arrays not images
    image1= np.array(image1)
    image2= np.array(image2)
    driftx=imshiftx2-imshiftx1 # drift in micron corrected with image registration
    drifty=imshifty2-imshifty1
    
    # Use scikit image register translation to find pixel shift
    pixshift, error, diffphase = register_translation(image1, image2)
    # shift gives pixel shift that was uncorrected
    shift=FOV*pixshift/512 # uncorrected shift in microns
    print('X Drift =', str(driftx), ' microns.  Y drift =', str(drifty), 'microns. Uncorr shift =', shift)
    return shift, error
    
def annotatejpg(jpgimage, areamatch):
    '''Pass Auger sem image and info about spatial areas and create annotated jpg with superimposed ROIs'''
    draw=ImageDraw.Draw(jpgimage) # single draw instance to label above image
    ttfont=ImageFont.truetype('arial.ttf', size=20)
    for index, row in areamatch.iterrows():
        areanum=int(areamatch.loc[index]['Areanumber'])
        x1=areamatch.loc[index]['X1']
        y1=areamatch.loc[index]['Y1']
        x2=areamatch.loc[index]['X2']
        y2=areamatch.loc[index]['Y2']
        fname=areamatch.loc[index]['Filename']
        draw.rectangle((x1,y1,x2,y2), outline='red') # red rectangle at specified position
        message=str(areanum) # label only with area number
        draw.text((x2+2,y2+2),message, font=ttfont, fill='red')
    annotjpgname=fname.replace('.spe','_areas.jpg')
    jpgimage.convert('RGB').save(annotjpgname)
    jpgimage.close()
    return

def annotateone(jpgfname, spename, SpatialAreas):
    '''Annotate a selected jpg image with spatial areas corresponding to selected spe file (pass both by filename)
    '''
    # now get spatial areas info from spatial areas log
    spename=spename.replace('.csv','.spe') # .spe extension used in spatialareas log
    areamatch=SpatialAreas[SpatialAreas['Filename']==spename]
    # Drop areanum duplicates (sometimes duplicate entries from rerun of import scripts)
    areamatch=areamatch.drop_duplicates(subset=['Areanumber'])
    
    try:
        jpgimage=Image.open(jpgfname)
        annotatejpg(jpgimage, areamatch)
    except:
        print('Problem creating jpg with annotated areas for spe file', spename)
    return

def makeannotatedjpg(AugerParamLog, SpatialAreas):
    ''' Searches log for spe and prior sem file to create image with annotated spatial areas
    For each spe, assume prior sem is image; pass/find number of areas and slice of spatialareas for this fnumber'''
    for i in range(0,len(AugerParamLog)):
        if '.csv' in AugerParamLog.iloc[i]['Filename']:
            spenum=AugerParamLog.iloc[i]['Filenumber']
            # only create annotation if prior pre-image is taken (assumes same field as spe)
            jpgmatch=AugerParamLog[AugerParamLog['Filenumber']==spenum-1]
            if len(jpgmatch)==1 and '.sem' in jpgmatch.iloc[0]['Filename']:
                jpgfname=jpgmatch.iloc[0]['Filename'].replace('.sem','.jpg')                
                # now get spatial areas info from spatial areas log
                areamatch=SpatialAreas[SpatialAreas['Filenumber']==spenum]
                # Drop areanum duplicates (sometimes duplicate entries from rerun of import scripts)
                areamatch=areamatch.drop_duplicates(subset=['Areanumber'])                
                try:
                    jpgimage=Image.open(jpgfname)
                    annotatejpg(jpgimage, areamatch)
                except:
                    print('Problem creating jpg with annotated areas for spe file', spenum)
    return 

def makeblanklog(filelist):
    ''' Make blank Excel log matching existing set of files (in cases where one was not created during data collection
    default naming scheme is prename.filenumber.sem or whatever but some changed to string_filenumber.sem '''
    mycols=['Project', 'Filename','Filenumber', 'Lastnumber', 'Sample', 'Comment', 'X', 'Y', 'Mag']    
    Augerlogbook=pd.DataFrame(columns=mycols) # blank df 
    # get project name from directory 
    fullpath=os.getcwd()
    pattern=re.compile(r'(\\)')
    match=re.finditer(pattern, fullpath)    
    indices=[m.start(0) for m in match] # digging out path name 
    projname=fullpath[indices[-1]+1:] # get project name (last folder) from full path
    for i, filename in enumerate(filelist):
        Samplelogrow=pd.DataFrame(index=np.arange(0,1), columns=mycols) # single df row for given file 
        match=re.finditer(r'\.',filename)
        starts=[m.start(0) for m in match]
        if len(starts)==2: # default formatted to *.filenumber.spe  
            try:
                filenumber=int(filename[starts[0]+1:starts[1]])               
                thisext=filename[starts[1]+1:] # file extension
                Samplelogrow=Samplelogrow.set_value(0,'Filenumber', filenumber) 
                Samplelogrow=Samplelogrow.set_value(0,'Filename', filename) # include full filename in case filenumber is not unique 
                if thisext=='spe':
                    Samplelogrow=Samplelogrow.set_value(0,'Mag', '.spe')
            except:
                print('Problem extracting filenumber.')
        if len(starts)==1:
            try:
                thisname=filename.split('.')[0]
                thisext=filename.split('.')[1]
                filenumber=int(re.findall(r'(\d+)', thisname)[-1]) # last string of numbers is filenumber
                Samplelogrow=Samplelogrow.set_value(0,'Filenumber', filenumber)
                Samplelogrow=Samplelogrow.set_value(0,'Filename', filename) 
                if thisext=='spe':
                    Samplelogrow=Samplelogrow.set_value(0,'Mag', '.spe')
            except:
                print('Problem extracting filenumber.')      
        Samplelogrow=Samplelogrow.set_value(0,'Project', projname)        
        # leave lastnumber, sample, comments, x,y mag all blank... just temporary before pass to Augerparamlog
        Augerlogbook=pd.concat([Augerlogbook,Samplelogrow], ignore_index=True)
    csvname=projname+'Auger_logbook.csv'
    Augerlogbook=Augerlogbook.sort_values(['Filenumber'])
    Augerlogbook=Augerlogbook[mycols] # reorder columns to standard
    Augerlogbook.to_csv(csvname, index=False)
    print('Blank logbook created for project ', projname, '; Sample names and comments can be manually entered .')
    return Augerlogbook

def openorcreatelogbook(filelist):
    ''' Looks for existing csv or xls log file ... if not found makes new one by calling makeblanklog ''' 
    logfile=glob.glob('*Auger_log*') # find ... Auger_logbook.
    if len(logfile)==1: # found logbook
        name=logfile[0]
        if '.xls' in name: # open log tab of existing excel file
            Augerlogbook=pd.read_excel(name, sheetname='Log')        
        if '.csv' in name: # open csv
            Augerlogbook=pd.read_csv(name)
    elif len(logfile)==0:
        # detect naming scheme (i.e. look for duplicate filenumber)
        Augerlogbook=makeblanklog(filelist)
    else:
        print('Error: There can be only one Auger logbook in this folder.')
        return 
    return Augerlogbook

def avglogentry(logmatches):
    ''' Create a logbook entry for a combine-averaged spectrum '''
    numfiles=len(logmatches) # number of combined files
    firststr=str(logmatches.iloc[0][0]) # averaged file name is numbered "firstnum"lastnum"
    laststr=str(logmatches.iloc[len(logmatches)-1][0])
    avgnum=int(firststr+laststr)
    avgentry=logmatches.iloc[[0]] # avg entry has same values as first file for cols 1,3-15, 17-24
    avgentry=avgentry.squeeze() # convert to series (avoids warnings regarding changing a slice of the df)
    # now change the small number of header params that will differ between single spe and averaged one
    avgentry=avgentry.set_value('Filenumber',avgnum) # reset filenumber
    avgentry.Cycles*=numfiles # adjust number of cycles due to combine-averaging
    avgentry.Acqtime*=numfiles # adjust effective aquisition time
    tempval=avgentry.Filename
    tempval=tempval.replace(firststr+'.',firststr+laststr+'.') # avoids error if number is elsewhere in name string
    avgentry=avgentry.set_value('Filename',tempval) # reassign with correct filename
    # Add avg 2-20 string to comments
    tempstr=str(logmatches.iloc[0]['Comments'])
    if tempstr!='' and tempstr!='nan':
        avgentry=avgentry.set_value('Comments',tempstr+' avg '+firststr+'-'+laststr)
    else:
        avgentry=avgentry.set_value('Comments','avg '+firststr+'-'+laststr)  
    avgentry=avgentry.set_value('Filename',tempval) # reassign with correct filename
    return avgentry # pandas series with all params and of correct dimension for append to main df log

def checklogfile(filelist, Augerlogbook):
    ''' Checks the user Auger logbook Excel for consistency with the actual data file list from directory
    prints out filenumbers that have a problem to console
    ''' 
    #TODO fix this so it works with new blanklog containing filenames (due to degenerate filenumber problem)
    spe=[] # list of filenumbers of spe files in directory
    semmap=[] # list of filenumbers for sem or map files in directory
    for i, name in enumerate(filelist): # deals with 3 cases (spe, sem or map)
        match=re.search(r'\d+.spe',name)
        if match:
            tempstring=match.group(0)
            num=int(tempstring.split('.')[0])
            spe.append(num) # add to spe files in dir list
        match=re.search(r'\d+.sem',name)
        if match:
            tempstring=match.group(0)
            num=int(tempstring.split('.')[0])
            semmap.append(num) # add to spe files in dir list        
        match=re.search(r'\d+.map',name)
        if match:
            tempstring=match.group(0)
            num=int(tempstring.split('.')[0])
            semmap.append(num) # add to spe files in dir list       
    logspecombine=[] # combineable spe files from excel logbook
    logother=[] # other file from Excel logbook (sem, map or single spe)  
    alllog=[]
    combinelist=Augerlogbook[(Augerlogbook['Lastnumber']>0)] # get file ranges to combine 
    tempdf=Augerlogbook.replace(np.nan, 0)    
    singlelist=tempdf[(tempdf['Lastnumber']==0)]
    for i in range(0,len(singlelist)):
        logother.append(singlelist.iloc[i][1])
    for i in range(0,len(combinelist)):
        first=int(combinelist.iloc[i][1])
        last=int(combinelist.iloc[i][2])
        for i in range(first,last+1):            
            logspecombine.append(i)
    alllog=logspecombine+logother
    alldir=spe+semmap # all filenumbers in directory
    for i,val in enumerate(logspecombine): # ensure all are in spe list
        if val not in spe:
            print('File combination error in Excel logfile.  File # ',val,' is not a combinable spe file.')
    missingdata=[i for i in alllog if i not in alldir] # list comprehension for missing data file (but present in excel logbook)
    missingentry=[i for i in alldir if i not in alllog] # list comprehension for data file with missing log entry (could also convert to sets and compare)
    for i, val in enumerate(missingdata):
        print ('Data file number ', val, ' mentioned in logbook but missing from directory')
    for i, val in enumerate(missingentry):
        print ('Data file number ', val, ' present in directory but missing from logbook')
    # check for duplicate entries in logbook
    myset=set([x for x in alllog if alllog.count(x) > 1])
    for i in myset:
        print('Duplicate entry for file number', i, ' in Excel logbook.')
    return 

def movespes(filelist, AugerParamLog):
    ''' Move average-combined spectra to subfolder named 'sub', split AugerParamlogs to reflect this
    '''
    # split AugerParamLog between combined csv subfiles and everything else
    newparamlist=AugerParamLog[AugerParamLog['Filename'].isin(filelist)]
    AugerParamLog=AugerParamLog[~AugerParamLog['Filename'].isin(filelist)] # main log with csv subs removed (not in filelist)
    newparamlist=newparamlist.reset_index(drop=True)
    AugerParamLog=AugerParamLog.reset_index(drop=True)
    # Append sub to path for each entry
    for i in range(0,len(newparamlist)):
        tempstr=newparamlist.iloc[i]['FilePath']
        newfolder=tempstr+'\sub'
        newparamlist=newparamlist.set_value(i,'FilePath', newfolder) # update path
        filename=newparamlist.iloc[i]['Filename'] # existing file name
        newname='sub/'+filename
        os.rename(filename,newname) # moves csv file in question to folder named sub
    # move csv subfiles to sub folder (append to or create )
    if os.path.exists('sub/Augerparamlog_subs.csv'): # append if file already exists
        with open('sub/Augerparamlog_subs.csv','a') as filedata: # appends so could have duplicated data rows
            newparamlist.to_csv(filedata, header=False, index=False) # should append to existing
    else:
        newparamlist.to_csv('sub/Augerparamlog_subs.csv',index=False) # creates new subs param file
    return AugerParamLog # send back w/ subfiles removed

def autocombinespe(AugerParamLog):
    ''' automatic method of finding combinable measurements (same X,Y, areas, basename) and combine all areas separately
    via averaging 
    naming convention of basename.1.csv, basename.2.csv, basename.3.csv etc to basename.13.csv
    '''
    if 'Basename' not in AugerParamLog: # in case basename not already added
        AugerParamLog['Basename']=''
        for index, row in AugerParamLog.iterrows():
            fname=AugerParamLog.loc[index]['Filename']
            basename=fname.split('.')[0]
            AugerParamLog=AugerParamLog.set_value(index,'Basename',basename)
    AugerParamLog.Comments=AugerParamLog.Comments.fillna(value='')
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
                from operator import itemgetter
                from itertools import groupby                
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
    return AugerParamLog
    
def combinespelist(filelist, AugerParamLog, csvname='', movefiles=False):
    ''' Manual method of combining files via averaging; Make filelist based on naming rules or however'''
    match=AugerParamLog.Filename.isin(filelist) # passing entire AugerParamlog so knock out those not in combine filelist
    logmatches=AugerParamLog.loc[match] # selects above list of files
    # Find first and last in filenumber range (files can be in random order)
    logmatches=logmatches.sort_values(['Filenumber']) # sort by filenumber
    # firstfile=logmatches.Filenumber.min() # 
    lastfile=logmatches.Filenumber.max()
    if csvname=='': # name automatically based on first-last filenumber
        filename=logmatches.iloc[0]['Filename'] # lowest numbered file after sort 
        csvname=filename.replace('.csv',str(lastfile)+'.csv') # same prefix but filenumber is firstnumlastnum concatenated
    if not os.path.isfile(csvname):  # skip combined file creation if it already exists (however log is still regenerated)    
        avgcombinespe(logmatches,csvname) # combines above files and makes new averaged csv 
        print('Average-combined file ', csvname, ' created.')
        # new entry into Augerparamlog  
        avgentry = avglogentry(logmatches) # create new logbook entry (series) for averaged spectrum (mostly from firstfile's info)
        # append new Series entry to end of AugerParamLog
        AugerParamLog=AugerParamLog.append(avgentry)
    else:
        print(csvname,' already exists')
    if movefiles==True: # moves all csv subfiles (those combined) into /sub directory
        # need different method without combine log file
        AugerParamLog=movespes(filelist, AugerParamLog) # shuffle csv sub files to subfolder and split Paramslog accordingly
    return AugerParamLog
    
def movecombinedspe(combinelist, AugerParamLog):
    '''if spectrum is average-combined, move these subfiles to 'sub' directory and create
    new associated AugerParamLog... get file numbers from combinelist '''
    # construct int list of file numbers from combinelist
    matchlist=[] # list of integer filenumbers for df mask (to create separate AugerParamLog for sub csvs)
    for i in range(0,len(combinelist)):
        filenumber=int(combinelist.iloc[i]['Filenumber'])
        lastnumber=int(combinelist.iloc[i]['Lastnumber'])
        for num in range(filenumber,lastnumber+1):
            matchlist.append(num)
    # split AugerParamLog between combined csv subfiles and everything else
    newparamlist=AugerParamLog[AugerParamLog['Filenumber'].isin(matchlist)]
    AugerParamLog=AugerParamLog[~AugerParamLog['Filenumber'].isin(matchlist)] # main log with csv subs removed
    # TODO is this selection method screwed up 
    newparamlist=newparamlist.reset_index(drop=True)
    AugerParamLog=AugerParamLog.reset_index(drop=True)
    # append sub to path for each entry
    for i in range(0,len(newparamlist)):
        tempstr=newparamlist.iloc[i]['FilePath']
        newfolder=tempstr+'\sub'
        newparamlist.set_value(i,'FilePath', newfolder)
        filename=newparamlist.iloc[i]['Filename'] # existing file name
        newname='sub/'+filename
        os.rename(filename,newname) # moves csv file in question to folder named sub
    # move csv subfiles to sub folder
    newparamlist.to_csv('sub/Augerparamlog_subs.csv',index=False) # save to sub
    return AugerParamLog # send back w/ subfiles removed
    
def avgcombinespe(logmatches, newcsvname): # logmatches is slice of Auger params dataframe with filenames, params, etc.
    ''' Takes consecutive spe files (chosen by Augerlogbook) and average-combines them into single new average dataset 
    pass name to be used for new avg-combined file'''    
    firstfile=int(logmatches.iloc[0]['Filenumber'])
    lastfile=int(logmatches.iloc[len(logmatches)-1]['Filenumber'])
    numfiles=lastfile-firstfile+1
    # Test to ensure that x, y,and # areas are same for files to be combined
    xval=logmatches.iloc[0]['X']
    yval=logmatches.iloc[0]['Y']
    numareas=int(logmatches.iloc[0]['Areas'])
    df=logmatches[(logmatches['X']==xval) & (logmatches['Y']==yval) & (logmatches['Areas']==numareas)]
    if len(df)!=len(logmatches):
        print('Check data log for error... X, Y and # of areas in combined files are not the same (atypical for AutoTool loop acquired spectra)' )
        print(logmatches) # just output df slice with offensive data
    # throw all n files into different columns of same frame? 
    filename=logmatches.iloc[0]['Filename']
    if os.path.exists(filename):
        newdf=pd.read_csv(filename) # imports firstfile
    elif os.path.exists('sub\\'+filename): # check for source file in /sub directory
        newdf=pd.read_csv('sub\\'+filename)
    else: # can't find file 
        print(filename,' is missing!')
        #return # if base file not found, end it
    evbreaks=logmatches.iloc[0]['Evbreaks'] # load appropriate evbreaks list from params log (difficult to determine with header info)    
    if type(evbreaks)==str: # create list of ints from string of format [0, 115, 230, 361]
        tempstring=evbreaks.split('[')[1]
        tempstring=tempstring.split(']')[0] # remove brackets
        evbreaks=[int(s) for s in tempstring.split(',')] # turn string into list of ints
        myset=set(evbreaks) # remove any evbreak duplicates (happens in some lists for unknown reason)
        evbreaks=list(myset)
        evbreaks.sort()
    for i in range(1,numfiles): #loop sums values from other spectra into newdf
        filename=logmatches.iloc[i]['Filename']
        if os.path.exists(filename):
            df=pd.read_csv(filename) # load nth spectrum
        elif os.path.exists('sub\\'+filename): # check for source file in /sub directory
            df=pd.read_csv('sub\\'+filename)
        else:
            print(filename,' is missing!')
            continue
        # add counts from each area (aka counts1, counts2, etc.) from nth file to original sheet      
        for j in range(1,numareas+1): # counts1 in [1], counts 2 in [3], counts3 in [5]
            colname='Counts'+str(j) # i is areanumber
            newdf[colname]=newdf[colname]+df[colname] # add values from next file to this area
    # now take average and recompute S7D7 and savgol (for all counts columns/ spatial areas)
    for i in range(1,numareas+1): # now compute average counts from above sum
        colname='Counts'+str(i) # 
        newdf[colname]=newdf[colname]/numfiles # counts/sec avg so divide by # of summed files
        # now compute S7D7 for average-combined counts (again for each area)
        thisser=newdf[colname] # grab these counts as a series 
        counts=thisser.tolist() # list with appropriate counts values
        smoothcol, smdiffcol=smoothdiffS7D7b(counts, evbreaks) # returns new S7D7 smoothed values
        newdf['S7D7'+str(i)]=pd.Series(smdiffcol)# turn sm-diff list to column, assigns to S7D71, S7D72, columns etc. 
        newdf['Smcounts'+str(i)]=pd.Series(smoothcol)
        thissavgol=makesavgolbreak(counts, evbreaks) # returned as list of values
        colname='Savgol'+str(i)
        newdf[colname]=pd.Series(thissavgol) # copy list directly to df column
    # Save the average dataset to new csv named firstfile_lastfileavg.csv
    print(str(firstfile)+'to '+ str(lastfile)+'combined via averaging to '+newcsvname)
    newdf.to_csv(newcsvname, index=False) # output survey data frame to csv file
    return

def combinespeloop(combinelist, AugerParamLog, movefiles=False):
    ''' Combine files  '''    
    for index, row in combinelist.iterrows(): # large loop through each match
        firstfile=int(combinelist.loc[index]['Filenumber']) # iloc grabs correct row from list of matches
        lastfile=int(combinelist.loc[index]['Lastnumber'])
        # look up filenames and # of areas associated with the above file range (all should have same # of area)
        logmatches=AugerParamLog[(AugerParamLog['Filenumber']>=firstfile) & (AugerParamLog['Filenumber']<=lastfile)]
        # Start of function combinespectra(logmatches) ... return new line for AugerFileParams logbook
        # Check if combine-averaged csv file has already been created for this set of spe files
        try:
            filename=logmatches.iloc[0]['Filename'] # get first file name
        except:
            print('Log is screwed up for ', firstfile,' and/or ', lastfile)
            continue
        newcsvname=filename.replace('.csv',str(lastfile)+'.csv') 
        if not os.path.isfile(newcsvname):  # skip combined file creation if it already exists (however log is still regenerated)    
            avgcombinespe(logmatches, newcsvname) # combines above files and makes new averaged csv 
            # new entry into Augerparamlog  
            avgentry = avglogentry(logmatches) # create new logbook entry (series) for averaged spectrum (mostly from firstfile's info)
            # append new Series entry to end of AugerParamLog
            AugerParamLog=AugerParamLog.append(avgentry)
        else:
            print(newcsvname,' already exists')
    if movefiles==True: # moves all csv subfiles (those combined) into /sub directory
        AugerParamLog=movecombinedspe(combinelist, AugerParamLog) # shuffle csv sub files to subfolder and split Paramslog accordingly
    return AugerParamLog
    
    
def getparamfromlog(filename, filenumber, df):
    ''' Search through logbook and return 1) project name, 2) sample name 3) nA 4) comments from df (which is Augerlogbook)
    must keep in mind that for Autotool looped data with multiple filenumbers need to find given # in range
    pass AugerFileName since '''
    # pass Augerlogbook dataframe as df and search it
    # can use .all() .any() or .loc() .. SQL-like splitting of dataframes
    logrowmatch=df[(df['Filenumber']==filenumber)] # matches all except multiple .spe files collected with Autotool
    if logrowmatch.empty: # finds 
        logrowmatch=df[(df['Filenumber']<filenumber) & (df['Lastnumber']>=filenumber)] # throws runtimewarning but works.. maybe install issue
    if len(logrowmatch)==1: # unique match on filenumber
        logrowmatch=logrowmatch.iloc[0] # convert to series (single row)
        return logrowmatch # df but only with single matching row
    elif len(logrowmatch)>1: # must have restarted numbering ... instead look for filename
        if 'Filename' in df: # those created with makeblanklog do contain filename
            logmatch=df[(df['Filename']==filename)] # look for passed filename
            if len(logmatch)==1:
                logmatch=logmatch.iloc[0] # convert to series (single row)
                return logmatch # df but only with single matching row
    elif len(logrowmatch)==0:
        print('No match for filenumber '+str(filenumber) +' in log file!')
        return # return nothing in case of no match    
    return

def makesavgol(counts):
    '''Perform python smooth-diff used to guide selection of background regions
    perform this in chunks between evbreaks, works for survey or multiplex
    returns list with smooth-diff columns
    '''
    savgollist=[] # Empty list to hold savgol data
    thisreg=np.asarray(counts) # convert list to array
    window_size=11
    deriv=2 
    order=2 # order of savgol fit 
    rate=1
    order_range = range(order+1) # range object
    half_window = (window_size -1) // 2 # type int
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    # b is matrix 3 by window size
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv) # series as long as array
    # linalg.pinv gets pseudo-inverse of a matrix (window-sized series)
    # .A of any matrix returns it as ndarray object 
    
    # Pad the signal at the extremes with values taken from the signal itself
    firstvals = thisreg[0] - np.abs( thisreg[1:half_window+1][::-1] - thisreg[0] )
    lastvals = thisreg[-1] + np.abs(thisreg[-half_window-1:-1][::-1] - thisreg[-1])
    thisreg = np.concatenate((firstvals, thisreg, lastvals))
    # Now convolve input signal and sav-gol processing 1D array)
    thisreg=np.convolve( thisreg, m[::-1], mode='valid')
    thisreg=thisreg.tolist() # convert to list 
    savgollist.extend(thisreg) # append to full spectrum list
    while len(savgollist)<len(counts): # can be one element too short 
        savgollist.append(0)
    return savgollist #  returns savitsky-golay smooth diff over same full region

def makesavgolbreak(counts, evbreaks):
    '''Perform python smooth-diff used to guide selection of background regions
    perform this in chunks between evbreaks, works for survey or multiplex
    returns list with smooth-diff columns
    '''
    savgollist=[] # Empty list to hold savgol data
    for i in range(0,len(evbreaks)-1):
        thisreg=counts[evbreaks[i]:evbreaks[i+1]] # slice into separate multiplex regions and process separately
        thisreg=np.asarray(thisreg) # convert list to array
        window_size=11
        deriv=2 
        order=2 # order of savgol fit 
        rate=1
        order_range = range(order+1) # range object
        half_window = (window_size -1) // 2 # type int
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        # b is matrix 3 by window size
        m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv) # series as long as array
        # linalg.pinv gets pseudo-inverse of a matrix (window-sized series)
        # .A of any matrix returns it as ndarray object 
        
        # Pad the signal at the extremes with values taken from the signal itself
        firstvals = thisreg[0] - np.abs( thisreg[1:half_window+1][::-1] - thisreg[0] )
        lastvals = thisreg[-1] + np.abs(thisreg[-half_window-1:-1][::-1] - thisreg[-1])
        thisreg = np.concatenate((firstvals, thisreg, lastvals))
        # Now convolve input signal and sav-gol processing 1D array)
        thisreg=np.convolve( thisreg, m[::-1], mode='valid')
        thisreg=thisreg.tolist() # convert to list 
        savgollist.extend(thisreg) # append to full spectrum list
    while len(savgollist)<len(counts): # can be one element too short 
        savgollist.append(0)
    return savgollist #  returns savitsky-golay smooth diff over same full region 
    
def smoothdiffS7D7(cnts):
    ''' create smooth differentiated column from counts using S7D7 PHI algorithm (Multipak tables A-5 and A-1''' 
    numpts=len(cnts)
    endpts=[0,len(cnts)-1] # legacy way of defining endpoints and internal breaks (although those no longer present)
    smooth=[0]*numpts # empty list of correct length for smoothed data
    smoothdiff=[0]*numpts # 7 pt diff of above smoothed data
    # smoothing of endpoints according to Multipak algorithm appendix table A-5    
    for i in range(0,numpts): # special cases for endpoints (within 3 of an evbreak)
        diff=i-min(endpts, key=lambda x:abs(x-i)) # distance from closest evbreak index # in list            
        if diff==0:
            if i==numpts-1: #last point
                smooth[i]=(2*cnts[i]+2*cnts[i-1]+1)/4 # additional special case for last point
            else: # first point
                smooth[i]=(2*cnts[i]+2*cnts[i+1]+1)/4 # all others at exact breaks can use value and adj higher value
        elif abs(diff)==1:  # works for +1 or -1 from nearest break
            smooth[i]=(1*cnts[i-1]+2*cnts[i]+1*cnts[i+1]+1)/4
        elif abs(diff)==2:
            smooth[i]=(-3*cnts[i-2]+12*cnts[i-1]+17*cnts[i]+12*cnts[i+1]+-3*cnts[i+2]+1)/35
        else:
            smooth[i]=(-2*cnts[i-3]+3*cnts[i-2]+6*cnts[i-1]+7*cnts[i]+6*cnts[i+1]+3*cnts[i+2]-2*cnts[i+3]+1)/21

    # same structure to perform differentiation on smoothed datalist above
    for i in range(0,numpts): # special cases for endpoints (within 3 of an evbreak)
        diff=i-min(endpts, key=lambda x:abs(x-i)) # distance from closest evbreak index # in list
        if abs(diff)<=2:
            smoothdiff[i]=0  # just zero out endpoints (old code moved to python software dev xls cell)
        else:
            smoothdiff[i]=(-3*smooth[i-3]-2*smooth[i-2]-1*smooth[i-1]+1*smooth[i+1]+2*smooth[i+2]+3*smooth[i+3]+1)/28
    return smooth, smoothdiff
    
def smoothdiffS7D7b(counts, evbreaks):
    '''create smooth differentiated column from counts using S7D7 PHI algorithm (Multipak tables A-5 and A-1
    version for rerun on combined spectra with internal ev breaks''' 
    
    numpts=len(counts)
    smooth=[0]*numpts # empty list of correct length for smoothed data
    smoothdiff=[0]*numpts # 7 pt diff of above smoothed data
    diffcol=[0]*numpts # testing
    # smoothing of endpoints according to Multipak algorithm appendix table A-5    
    for i in range(0,numpts): # special cases for endpoints (within 3 of an evbreak)
        diff=i-min(evbreaks, key=lambda x:abs(x-i)) # distance from closest evbreak index # in list            
        if diff==0:
            if i==numpts-1: #last point
                smooth[i]=(2*counts[i]+2*counts[i-1]+1)/4 # additional special case for last point
                diffcol[i]=diff
            else: # first point
                smooth[i]=(2*counts[i]+2*counts[i+1]+1)/4 # all others at exact breaks can use value and adj higher value
                diffcol[i]=diff
        elif abs(diff)==1:  # works for +1 or -1 from nearest break
            smooth[i]=(1*counts[i-1]+2*counts[i]+1*counts[i+1]+1)/4
            diffcol[i]=diff
        elif abs(diff)==2:
            smooth[i]=(-3*counts[i-2]+12*counts[i-1]+17*counts[i]+12*counts[i+1]+-3*counts[i+2]+1)/35
            diffcol[i]=diff
        else:
            smooth[i]=(-2*counts[i-3]+3*counts[i-2]+6*counts[i-1]+7*counts[i]+6*counts[i+1]+3*counts[i+2]-2*counts[i+3]+1)/21
            diffcol[i]=diff

    # same structure to perform differentiation on smoothed datalist above
    for i in range(0,numpts): # special cases for endpoints (within 3 of an evbreak)
        diff=i-min(evbreaks, key=lambda x:abs(x-i)) # distance from closest evbreak index # in list
        if abs(diff)<=2:
            smoothdiff[i]=0  # just zero out endpoints (old code moved to python software dev xls cell)
        else:
            smoothdiff[i]=(-3*smooth[i-3]-2*smooth[i-2]-1*smooth[i-1]+1*smooth[i+1]+2*smooth[i+2]+3*smooth[i+3]+1)/28
    return smooth, smoothdiff # return both smoothed and smooth-differentiated columns
        
def makesurvey(bindata, csvdataname, extractargs):
    ''' function to extract binary spe survey data to dataframe for csv export
    called from processAuger if file is survey spe and csv doesn't exist
    needed variables for survey: starting eV (aka lower); evstep, # of points, # of areas'''
    
    # can read numbers of areas and # of pts from binary header (after text header) but pass them anyway
    survey=pd.DataFrame() # data frame for energy, and counts col for each area
     # create energy values list and add to survey dataframe
    startev=extractargs.get('startev')
    evstep=extractargs.get('evstep')
    numpts =extractargs.get('numpts')
    numareas=extractargs.get('numareas')
    energy=[]
    # evbreaks=[0, numpts-1] # indicates to smoothdiff that there are no energy breaks/discontinuities
    
    for i in range(0,numpts): # generate energy values for survey
        energy.append(startev+evstep*i)
    survey['Energy']=energy
    
    # read in the binary data from bindata and convert
    # structure is 4 byte floats of length numpts starting at byte 112
    # binary header contains numareas, numspectralregions, numpts but just pass these from text header
    # if multiple spatial areas another list of length numpts just added to end
    startbyte=112 # doesn't seem to vary for surveys   
    mystruct=struct.Struct('f')
    # TEST for correct byte location... should have leading zeros for bytes 100:111 
    for i in range(startbyte-12,startbyte,4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval) # little endian encoding
        if unpackbyte[0]!=0:
            print('Possible binary read error... leading zeros are missing for file', csvdataname)

    # read in and convert all numeric values
    counts=[] # single list for all Auger counts values in file
    for i in range(startbyte,len(bindata),4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval)
        counts.append(unpackbyte[0])
    # if numareas>1, split up counts list 
    if len(counts)!=numpts*numareas:
        print('Possible binary data reading error: Data string length does not match expected value for ',csvdataname)
    # if numareas> 1, split counts list into n separate areas
    if numareas>1:
        tempvals=counts[0:numpts]
        survey['Counts1']=tempvals # add counts from 1st area to survey dataframe
        smoothdat, S7D71=smoothdiffS7D7(tempvals) # no evbreaks for surveys
        survey['S7D71']=S7D71
        survey['Smcounts1']=smoothdat
        Surveysavgol=makesavgol(tempvals) # call function to make normal Savgol 
        survey['Savgol']=Surveysavgol # add column with 2nd deriv, savgol poly 2, 11 pts 
        for i in range(2,numareas+1):
            tempvals=counts[numpts*(i-1):numpts*i]
            tempcolname='Counts'+str(i) # name for nth column is counts2, counts3, etc.
            tempsmname='Smcounts'+str(i)
            tempSDname='S7D7'+str(i)
            tempsavgolname='Savgol'+str(i)
            tempsmooth, tempS7D7=smoothdiffS7D7(tempvals) # perform and return 7pt smooth and smooth-diff S7D7 for each
            survey[tempcolname]=tempvals # add nth area data to dataframe
            survey[tempSDname]=tempS7D7 # add col with Multipak S7D7
            survey[tempsmname]=tempsmooth # add smoothed col
            tempsavgol=makesavgol(tempvals) # call Savgol 2nd deriv w/ smoothing
            survey[tempsavgolname]=tempsavgol # add column with Savgol 
            
    else: # single area case
        survey['Counts1']=counts # for single area  name survey as counts1 anyway to avoid problems
        smooth, S7D7=smoothdiffS7D7(counts) # perform smooth diff using S7D7 Multipak algorithm
        survey['S7D71']=S7D7 # name as S7D71 even for single region
        survey['Smcounts1']=smooth # name as S7D71 even for single region
        thissavgol=makesavgol(counts) # add sav-gol smooth diff data
        survey['Savgol1']=thissavgol
        
    survey.to_csv(csvdataname, index=False) # output survey data frame to csv file      
    return

def makemultiplex(bindata,  csvdataname, energy, evbreaks, numareas, spectralregs):
    ''' function to create a csv file containing the extracted multiplex spectral values
    list of energy values are passed from SpectralRegions function (energy x values are not continuous for multiplex)
    bindata is rest of binary after text header extraction
	evbreaks are index# of boundaries between spectral regions '''
    mystruct=struct.Struct('f')
    # binary header of variable length.. just count backward using known data length
    startbyte=len(bindata)-4*numareas*len(energy)
    
    # TEST for correct byte location... data regions begins right after last occurence of triple zero 4 byte float values.. for survey usually found in bytes 100:111 
    for i in range(startbyte-12,startbyte,4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval) # little endian encoding
        if unpackbyte[0]!=0:
            print('Possible binary read error... leading zeros are missing for ',csvdataname)

    # create data frame for multiplex and assign energy values   
    multiplex=pd.DataFrame() # data frame for energy, and counts col for each area   
    multiplex['Energy']=energy # energy values found in SpectralRegions and passed as argument
    
    # read in and convert all numeric values
    
    alldata=[] # single list for all Auger counts values in file
    for i in range(startbyte,len(bindata),4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval)
        alldata.append(unpackbyte[0])
    
    if len(alldata)!=len(energy)*numareas: # number of data values expected for each area
        print('Possible binary data reading error: Data string length does not match expected value for ',csvdataname)
    # multiplex file structure has same energy region of all areas bundled together (split counts into spectral regions based on evbreaks)
    # so organization is multiplex spectral region 1 (all areas), spectral region 2 (all areas), etc.  
    datachunks=[]
    for i in range(0,len(evbreaks)-1):
        datachunks.append(evbreaks[i+1]-evbreaks[i])
    datachunks[0]=datachunks[0]+1 # adjust for length of first regions
    datachunks=[i*numareas for i in datachunks] # total lengths of each spectral region
    databoundaries=[]
    for i,val in enumerate(datachunks):
        temp=datachunks[0:i]
        databoundaries.append(sum(temp))
    databoundaries.append(sum(datachunks))
    specregs=[] # list of lists containing counts values for each spectral region
    # now split data into single spectral region with all spatial areas
    for i in range(0,len(databoundaries)-1):
        specregs.append(alldata[databoundaries[i]:databoundaries[i+1]])
    # Now construct counts for each area
    counts=[[] for x in range(0,numareas)] # list of empty lists one for each spatial area
    smoothdata=[[] for x in range(0,numareas)] # for 7 pt smoothing of counts for each area
    s7d7data=[[] for x in range(0,numareas)] # for s7d7 of each area
    savgoldata=[[] for x in range(0,numareas)] # for savgol deriv of each area
    
    for i,vals in enumerate(specregs):
        # vals is list of count values for this spectral region (for each spatial area back to back)
        numvals=int(datachunks[i]/numareas) # number of values in single area/single spectral region    
        for j in range(0, numareas):
            counts[j].extend(vals[j*numvals:(j+1)*numvals]) #
            thischunk=vals[j*numvals:(j+1)*numvals]
            thissmooth, thiss7d7=smoothdiffS7D7(thischunk)        
            s7d7data[j].extend(thiss7d7) # sent as single region w/o internal data breaks 
            smoothdata[j].extend(thissmooth) # sent as single region w/o internal data breaks 
            thissavgol= makesavgol(thischunk)           
            savgoldata[j].extend(thissavgol) # compute savgol and add to list
    # now assign these values to multiplex columns 
    for i in range(1,numareas+1):
        cntsname='Counts'+str(i) # name for nth column is counts2, counts3, etc.
        SDname='S7D7'+str(i)
        savgolname='Savgol'+str(i)
        multiplex[cntsname]=counts[i-1] # assign counts to frame (and switch from 0 based indexing)
        multiplex[SDname]=s7d7data[i-1]
        multiplex[savgolname]=savgoldata[i-1] # add savgol data to new column    
    # solve multiplex out of order problem
    if not pd.algos.is_lexsorted([multiplex.Energy.values]): # data and evbreaks out of order problem
        multiplex, evbreaks=sortmultiplex(multiplex, evbreaks, spectralregs)
    multiplex.to_csv(csvdataname, index=False) # output data frame to csv file            
    return evbreaks 

def keepbestvals(df, spectralregs):
    '''For duplicate energy values in multiplex scans, keep the one with largest # of sweeps '''
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
    
def sortmultiplex(multiplex, evbreaks, spectralregs):
    ''' Rearrange multiplex if taken out of order and adjust evbreaks accordingly
    evbreaks holds index #s of multiplex breaks so must be altered along with multiplex sort
    only called if multiplex is not monotonically increasing '''
    energylist=[]
    for i, val in enumerate(evbreaks):
        energylist.append(multiplex.loc[val]['Energy'])
        if i>0 and i<len(evbreaks)-1:
            energylist.append(multiplex.loc[val+1]['Energy']) # starting ev for next spectral region
    # now if duplicates exist, keep the best values (larger 3 of sweeps)
    multiplex= keepbestvals(multiplex, spectralregs) # knocks out duplicates based on # of sweeps
    multiplex=multiplex.sort_values(['Energy'])
    multiplex=multiplex.reset_index(drop=True)
    matches=multiplex[multiplex['Energy'].isin(energylist)]
    # Some can be internal data breaks if multiplex regions overlap
    newevbreaks=matches.index.tolist()
    newevbreaks=[int(i) for i in newevbreaks if i-1 not in newevbreaks] # remove adjacent values
    return multiplex, newevbreaks
   
def makejpg(bindata, jpgname, resolution, fieldofview):
    ''' Extract raw binary SEM image (512x512 or whatever), convert to jpg and saves with scaled dpi 
    such that 10 micron field of view image measures 10x10cm in photoshop; only called by processSEM if jpg doesn't already exist '''
              
    # Read in the binary data from bindata and convert
    # structure is 4 byte floats of length resolution^2 starting at byte 112
    startbyte=112 # doesn't seem to vary for surveys   
    mystruct=struct.Struct('f')
    # TEST for correct byte location... should have leading zeros for bytes 100:111 for raw sem image
    for i in range(startbyte-12,startbyte,4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval) # little endian encoding
        if unpackbyte[0]!=0:
            print('Possible binary read error... leading zeros are missing for ',jpgname)

    # read in and convert all numeric values
    intensity=[] # single list for all pixel intensity values in file
    for i in range(startbyte,len(bindata),4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval)
        intensity.append(unpackbyte[0])
    if len(intensity)!=resolution**2:
        print('Possible binary data reading error: Data string length for intensities does not match expected length for ',jpgname)
    imagearray=np.array(intensity).reshape(resolution,resolution) # create numpy image array (typically 52x512)
    imagearray=np.rot90(imagearray, k=3) # rotate by 90CW
    imagearray=np.fliplr(imagearray) # flip horizontal
    # testable with imshow(imagearray) after pylab import
    img=Image.fromarray(imagearray) # convert numpy array to PIL format
    thisdpi=(int((512*2.54)/(fieldofview)),int((512*2.54)/(fieldofview))) # set dpi tuple based on field of view (so that 10micron image measures 10cm)
    img.convert('RGB').save(jpgname,dpi=thisdpi) # convert to RGB format and save at appropriate resolution
    return 

def makemaps(bindata, jpgbasename, resolution, fieldofview, mapelements):
    '''Extracts raw image data and creates jpg for each mapped element:  pass binary data string, base name for jpgs, resolution, field of view and 
    list of mapped elements '''   
    # find startbyte by just subtracting known data length from the end    
    startbyte=len(bindata)-4*len(mapelements)*resolution**2
    maplength=resolution**2 # size of each elements map in bytes
    mystruct=struct.Struct('f')
    # TEST for correct byte location... should have leading zeros for bytes 100:111 for raw sem image
    for i in range(startbyte-12,startbyte,4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval) # little endian encoding
        if unpackbyte[0]!=0:
            print('Possible binary read error... leading zeros are missing for ',jpgbasename)
        # read in and convert all numeric values
    intensity=[] # single list for all pixel intensity values in file
    for i in range(startbyte,len(bindata),4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval)
        intensity.append(unpackbyte[0])
    thisdpi=(int((512*2.54)/(fieldofview)),int((512*2.54)/(fieldofview))) # set dpi tuple based on field of view (so that 10micron image measures 10cm)    
    # extract each map of length maplength and save as 'jpgbasename"_elem.jpg
    for i, elem in enumerate(mapelements):
        tempint=intensity[maplength*i:maplength*(i+1)]
        imagearray=np.array(tempint).reshape(resolution,resolution) # create numpy image array (typically 52x512)
        imagearray=np.rot90(imagearray, k=3) # rotate by 90CW
        imagearray=np.fliplr(imagearray) # flip horizontal
        # testable with imshow(imagearray) after pylab import
        img=Image.fromarray(imagearray) # convert numpy array to PIL format    
        img.convert('RGB').save(jpgbasename+'_'+elem+'.jpg',dpi=thisdpi) # convert to RGB format and save at appropriate resolution
        
def SpectralRegions(filenumber, AugerFileName, numcycles, timestep, header): # gets spectral region info from header (called from processAuger below)
    ''' Subroutine to pull multiplex spe params (e.g. elements and spectral regions), calculate acquisition time for multiplex .spe file '''
    tempstring=header.split('NoSpectralReg: ')[1] # unlike NoSpectralRegFull inactives are already removed
    match=re.search(r'\d+',tempstring)
    numdefregions=int(match.group(0)) # number of defined regions (can be larger than # active regions)
    SpectralRegs=pd.DataFrame(columns=['Filenumber','Filename','Numcycles','Timestep','Element','Sweeps','Evstep','Lower','Upper','Time'])
    numregions=0 # active regions (as opposed to nominal defined regions)
    timeperarea=0 #initialize as zero
    details='' # string for element list
    energy=[] # list for energy x values
    evbreaks=[0] # region boundaries needed for smoothdiff include first
    for i in range(numdefregions):
        tempstr=tempstring.split('SpectralRegDef: ')[i+1] # should work to split 
        element=(tempstr.split(' ')[2]) # name of elemental line
        numpts=int(tempstr.split(' ')[4])        
        evstep=float(tempstr.split(' ')[5]) #eV/step        
        startev=float(tempstr.split(' ')[6]) # starting eV
        endev=float(tempstr.split(' ')[7]) # ending eV
        for j in range(0,numpts): # generate energy values for survey
            energy.append(startev+evstep*j)
        evbreaks.append(len(energy)-1) # gives indices of discontinuities in energy scan (needed for smoothdiffS7D7)   
        tempstr=tempstr.split('SpectralRegDef2: ')[1]
        sweeps=int(tempstr.split(' ')[2]) # number of sweeps through element        
        time=numcycles*timestep*sweeps*(endev-startev)/1000 # acquisition time in seconds for this multiplex region        
        timeperarea+=time # add on acquisition time for this element eV range only
        # add ith new data row to spatialareas dataframe
        # make a list and then append as row? 
        SpectralRegs.loc[i]=[filenumber, AugerFileName, numcycles, timestep, element, sweeps, evstep, startev, endev, time]
        numregions+=1
        details=details+str(sweeps)+'x'+element+' ' # add current region to detailed string
        
    with open('multiplexspectralregionslog.csv','a') as filedata:
        SpectralRegs.to_csv(filedata, header=False, index=False) # directly write-appended to existing log 
    return details, evbreaks, timeperarea, energy, SpectralRegs
    # returns 1) details string summarizing multiplex conditions, 2) index #s of breaks/discontiuities in energy scan
    # 3) acq time per area and 4) energy is list of x vals, 5) df chunk for handling duplicate values
    
def SpectralMapRegions(filenumber, AugerFileName, timestep, header): 
    ''' Subroutine to pull mapping params and calculate acquisition time for .map file; all saved in param log csv'''
    # determine if 2pt map or 3 pt map from header
    tempstring=header.split('ScanMode: ')[1]
    tempstring=tempstring.split('point')[0]
    ptsmode=int(tempstring) # should be 2 or 3
    # determine number of pixels
    tempstring=header.split('NoMapPixelsXY: ')[1]
    tempstring=tempstring.split(' ')[0]
    pixels=int(tempstring)
    numcycles=0 # cycles not relevant for map
    tempstring=header.split('NoSpectralReg: ')[1] # unlike NoSpectralRegFull inactives are already removed
    match=re.search(r'\d+',tempstring)
    numdefregions=int(match.group(0)) # number of defined regions (can be larger than # active regions)
    SpectralRegs=pd.DataFrame(columns=['Filenumber','Filename','Numcycles','Timestep','Element','Sweeps','Evstep','Lower','Upper','Time'])
    numregions=0 # active regions (as opposed to nominal defined regions)
    timeperarea=0 #initialize as zero
    details='' # string for element list
    for i in range(numdefregions):
        tempstr=tempstring.split('SpectralRegDef: ')[i+1] # should work to split 
        element=(tempstr.split(' ')[2]) # name of elemental line
        # evstep=float(tempstr.split(' ')[5]) 
        evstep='' # skip eV per step for map 
        startev=float(tempstr.split(' ')[6]) # for map 3 pt use low background energy
        endev=float(tempstr.split(' ')[7]) # for map 3 pt use high background energy
        tempstr=tempstr.split('SpectralRegDef2: ')[1]
        sweeps=int(tempstr.split(' ')[2]) # number of sweeps through element (same format in map or Auger multiplex)       
        time=ptsmode*timestep*sweeps*pixels**2/1000 # acquisition time in seconds for this multiplex region        
        # acquisition time for maps sum(sweeps for all EL)*2pt or 3pt * resolution**2 /(1000*60)    
        timeperarea+=time # add on acquisition time for this element eV range only
        # add ith new data row to spatialareas dataframe
        # make a list and then append as row? 
        SpectralRegs.loc[i]=[filenumber, AugerFileName, numcycles, timestep, element, sweeps, evstep, startev, endev, time]
        numregions+=1
        details=details+str(sweeps)+'x'+element+' ' # add current region to detailed string
            
    with open('multiplexspectralregionslog.csv','a') as filedata:
        SpectralRegs.to_csv(filedata, header=False, index=False)    
    return details, timeperarea
    
def SpatialAreas(header,fieldofview, filenumber, AugerFileName):
    ''' SpatialAreas outputs the areas defined in Auger and calculates/returns number of areas and total area''' 
    tempstring=header.split('NumSpatialAreas: ')[1]
    match=re.search(r'\d+',tempstring)
    numdefareas=int(match.group(0)) # number of defined areas (can be larger than # active areas)
    SpatialAreas=pd.DataFrame(columns=['Filenumber','Filename','Areanumber','X1','Y1','X2','Y2','Area']) # easier to use frame than dict
    numareas=0
    totalarea=0 #initialize as zero
    for i in range(numdefareas):
        tempstr=tempstring.split('SpatialArea: ')[i+1] # should work to split 
        # Check if area is active        
        active=int(tempstr.split(' ')[1]) # 1 if active and 0 if inactive
        if active==1:
            X1=float(tempstr.split(' ')[3])
            Y1=float(tempstr.split(' ')[4])
            X2=float(tempstr.split(' ')[5])
            tempstr=tempstr.split(' ')[6]
            Y2=float(tempstr.split('\n')[0])
            area=fieldofview**2/(512**2)*(X2-X1)*(Y2-Y1) # add area of current box to total area     
            totalarea+=area
            # add ith new data row to spatialareas dataframe
            # make a list and then append as row? 
            SpatialAreas.loc[i]=[filenumber, AugerFileName, i+1, X1, Y1, X2, Y2, area]
            numareas+=1
    ## open CSV in append mode and add new X1,Y1, etc. incrementally to spatialareaslog
    with open('spatialareaslog.csv','a') as filedata: # appends so could have duplicated data rows
        SpatialAreas.to_csv(filedata, header=False, index=False) # creates file if it doesn't exist
    return numareas, totalarea
    
def processAuger(AugerFileName):   
    '''processAuger pulls 20 important parameters from spe file header and saves in csv param log'''
    with open(AugerFileName, 'rb') as file:
        filedata = file.read()
    end=filedata.find(b'EOFH')
    headerdata=filedata[0:end+6] # works to cut off binary part    
    header=headerdata.decode(encoding='cp437') # more generic encoding than utf-8
    # temporary dict for holding important file parameters
    Augerparams={'Filenumber':'','Project':'','Filename':'','FilePath':'','Sample':'','Comments':'','Date':'','FieldofView':'','Type':'','Energy':'','GunVoltage':'','nA':'','Areas':'','Cycles':'','Timestep':'','Details':'','Evbreaks':'','Acqtime':'','Scanarea':'','X':'','Y':'','Z':'','Tilt':'','Rotation':'','ImageshiftX':'','ImageshiftY':''}
    tempstring=AugerFileName.replace('spe','csv') # save csv name not original spe
    Augerparams.update({'Filename':tempstring})
   # pull acquisition date 
    tempstring=header.split('FileDate: ')[1]
    match=re.search(r'(\d+\s\d+\s\d+)',tempstring)
    tempstring=match.group(0)
    if match:
        Augerparams.update({'Date':tempstring.split(' ')[1]+'/'+tempstring.split(' ')[2]+'/'+tempstring.split(' ')[0]}) # group0 returns matching string

    tempstring=header.split('AcqFilename: ')[1] #find file number
    tempstring=tempstring.split('.')[1] # convert to int to avoid later comparison problems
    if 'Tmp' in tempstring: # avoid error with Tmp1.sem problem caused by PHI duplicate numbering
        tempstring=tempstring.split('Tmp')[0] # remove Tmp string       
    tempint=int(tempstring) # file number according to header
    templist=re.findall(r'\d+', AugerFileName)
    tempint2=int(templist[-1])  # last number in AugerFileName is filenumber
    Augerparams.update({'Filenumber':tempint})    
    if tempint!=tempint2:
        Augerparams.update({'Filenumber':tempint2}) # value in filename probably more reliable   
        print('Suspected acqfilename error in header of file ', str(tempint),'. Manually edit before combining.')
    tempstring=header.split('FieldOfView: ')[1] # find field of view
    match=re.search(r'\d+\.\d+',tempstring) 
    if match:
        tempstr=match.group(0)
        fieldofview=float(tempstr)
        Augerparams.update({'FieldofView':fieldofview}) # group0 returns matching string

    tempstring=header.split('ImageShift: ')[1] # find image shift x and y
    match=tempstring.split(' ')[0]
    tempstring=tempstring.split(' ')[1]
    Augerparams.update({'ImageshiftX':float(match)}) # group0 returns matching string
    Augerparams.update({'ImageshiftY':float(tempstring)})

    tempstring=header.split('EBeamEnergy: ')[1] # find kV
    tempstring=tempstring.split(' ')[0]
    Augerparams.update({'Energy':tempstring}) # group0 returns matching string

    tempstring=header.split('GunLensVoltage: ')[1] # find gunvoltage (related to current)
    tempstring=float(tempstring.split(' ')[0])
    Augerparams.update({'GunVoltage':tempstring})

    tempstring=header.split('NumCycles: ')[1] # find # cycles (works for either survey or multiplex)
    match=re.search(r'\d+',tempstring)
    if match:
        Augerparams.update({'Cycles':int(match.group(0))})
    # find spectral type (Survey or Multiplex)
    match=re.search('MultNumCycles:',header)
    if match:
        Augerparams.update({'Type':'Multiplex'}) 
        spectraltype='multiplex'
    else:
        Augerparams.update({'Type':'Survey'})
        spectraltype='survey' 
    tempstring=header.split('TimePerStep:')[1] # find time per step
    match=re.search(r'\d+\.\d+',tempstring)  
    if match:
        Augerparams.update({'Timestep':float(match.group(0))})        
    # find X,Y,Z, Rotation, Tilt (can be negative)
    tempstring=header.split('DefectPosX:')[1]
    xval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'X':xval})
    
    tempstring=header.split('DefectPosY:')[1] # 
    yval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Y':yval})  
 
    tempstring=header.split('DefectPosZ:')[1] # 
    zval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Z':zval})  

    tempstring=header.split('DefectPosTilt:')[1] # 
    tilt=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Tilt':tilt})  

    tempstring=header.split('DefectPosRotation:')[1] # 
    rot=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Rotation':rot})      

    tempstring=header.split('DefectPosTilt:')[1] # 
    match=re.search(r'\d+\.\d+',tempstring)  
    if match:
        Augerparams.update({'Tilt':float(match.group(0))})  

    tempstring=header.split('DefectPosRotation:')[1] # 
    match=re.search(r'\d+\.\d+',tempstring)  
    if match:
        Augerparams.update({'Rotation':float(match.group(0))})  

    # pull info about spectral regions (pass file number and header)
    numcycles=(Augerparams.get('Cycles'))
    filenumber=(Augerparams.get('Filenumber'))
    timestep=float((Augerparams.get('Timestep')))
    if spectraltype=='multiplex': # extract and return multiplex params
        details, evbreaks, timeperarea, energy, spectralregs =SpectralRegions(filenumber, AugerFileName, numcycles, timestep, header)
        # returns details string on multiplex conditions for inclusion in parameters log file
        # return single list of energy value for multiplex of separate regions (energy x values not continuous)
        # this list (energy) is then used for all areas by makemultiplex() 
        # returns calculated time per area used along with numareas for total acquisition time
    if spectraltype=='survey': # determine survey range
        tempstring=header.split('Sur1 ')[1]
        numpts=int(tempstring.split(' ')[1])        
        evstep=float(tempstring.split(' ')[2])       
        startev=float(tempstring.split(' ')[3]) # starting eV
        endev=float(tempstring.split(' ')[4]) # ending eV
        timeperarea=numcycles*timestep/1000*(endev-startev)/evstep #acquisition time for single area in survey mode
        details=' '        
        details=str(numcycles)+'x '+str(startev)+'-'+str(endev)
        # starting eV, eV/step and # of points needed for cnts vs energy
        extractargs={'startev':startev,'evstep':evstep,'numpts':numpts}
        evbreaks=[0, numpts-1] # list of evbreak index numbers saved to main param log 
    # pull # areas and compute actual scanned area (SpatialAreas append these to separate areas log)
    numareas, totalarea = SpatialAreas(header,fieldofview, filenumber, AugerFileName)     
    # strip header and extra columns from HMR scan data and save
    Augerparams.update({'Scanarea':totalarea,'Areas':numareas,'Details':details})
    Augerparams.update({'Areas':numareas})
    acqtime=timeperarea*numareas  # gives overall time for multiplex acquisition
    Augerparams.update({'Acqtime':acqtime})
    # if not done already,extract binary spe data 
    csvdataname=AugerFileName.split('.spe')[0]+'.csv'
    
    # Call to binary data parser to create survey or multiplex csv data files
    # this is done automatically unless csv file already exists
    # only creates CSV file... nothing returned
    if not os.path.isfile(csvdataname): # will not recreate csv if it exists
        bindata=filedata[end+6:] # gets bytes holding actual data
        if spectraltype=='survey':         
            extractargs.update({'numareas':numareas}) # numareas, numpoints, ev/step and starting ev        
            makesurvey(bindata, csvdataname, extractargs)
        if spectraltype=='multiplex':
            # for multiplex only need to pass 1) binary data string 2) constructed x values list and 3) number of areas            
            evbreaks=makemultiplex(bindata, csvdataname, energy, evbreaks, numareas, spectralregs)
    Augerparams.update({'Evbreaks':evbreaks}) # for multiplex calculated in SpectralRegions (or modified after sort), as list not string
    newname='sub/'+ AugerFileName
    os.rename(AugerFileName,newname) # move raw sem file to subfolder sub
    return Augerparams

def processSEM(AugerFileName): 
    ''' pulls header info parameters from .sem image files and saves in master param log'''
    with open(AugerFileName, 'rb') as file:
        filedata = file.read()
    end=filedata.find(b'EOFH')
    headerdata=filedata[0:end+6] # works to cut off binary part    
    header=headerdata.decode(encoding='cp437') # more generic encoding than utf-8
    # temporary dict for holding important file parameters
    Augerparams={'Filenumber':'','Project':'','Filename':'','FilePath':'','Sample':'','Comments':'','Date':'','FieldofView':'','Type':'','Energy':'','GunVoltage':'','nA':'','Areas':'','Cycles':'','Timestep':'','Details':'','Acqtime':'','Scanarea':'','X':'','Y':'','Z':'','Tilt':'','Rotation':'','ImageshiftX':'','ImageshiftY':''}
    
    Augerparams.update({'Filename':AugerFileName})
    path=os.getcwd()
    Augerparams.update({'FilePath':path})
   # pull acquisition date 
    tempstring=header.split('FileDate: ')[1]
    match=re.search(r'(\d+\s\d+\s\d+)',tempstring)
    tempstring=match.group(0)
    if match:
        Augerparams.update({'Date':tempstring.split(' ')[1]+'/'+tempstring.split(' ')[2]+'/'+tempstring.split(' ')[0]}) # group0 returns matching string
    # get filenumber from AugerFileName not from AcqFilename in header... avoids rename file (e.g.) Tmp problems
    tempstring=header.split('AcqFilename: ')[1] #find file number
    tempstring=tempstring.split('.')[1] # convert to int to avoid later comparison problems
    if 'Tmp' in tempstring: # avoid error with Tmp1.sem problem caused by PHI duplicate numbering
        tempstring=tempstring.split('Tmp')[0] # remove Tmp string   
    tempint=int(tempstring)
    Augerparams.update({'Filenumber':tempint})
    
    tempstring=header.split('FieldOfView: ')[1] # find field of view
    match=re.search(r'\d+\.\d+',tempstring) 
    if match:
        tempstr=match.group(0)
        fieldofview=float(tempstr)
        Augerparams.update({'FieldofView':fieldofview}) # group0 returns matching string

    tempstring=header.split('ImageShift: ')[1] # find image shift x and y
    match=tempstring.split(' ')[0]
    tempstring=tempstring.split(' ')[1]
    Augerparams.update({'ImageshiftX':float(match)}) # group0 returns matching string
    Augerparams.update({'ImageshiftY':float(tempstring)})

    tempstring=header.split('EBeamEnergy: ')[1] # find kV
    tempstring=tempstring.split(' ')[0]
    Augerparams.update({'Energy':tempstring}) # group0 returns matching string

    tempstring=header.split('GunLensVoltage: ')[1] # find gunvoltage (related to current)
    tempstring=float(tempstring.split(' ')[0])

    Augerparams.update({'GunVoltage':tempstring})
    
    Augerparams.update({'Areas':np.nan}) # set to nan so sorting works correctly
    # find X,Y,Z, Rotation, Tilt
    tempstring=header.split('DefectPosX:')[1]
    xval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'X':xval})

    tempstring=header.split('DefectPosY:')[1] # 
    yval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Y':yval})  
 
    tempstring=header.split('DefectPosZ:')[1] # 
    zval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Z':zval})  

    tempstring=header.split('DefectPosTilt:')[1] # 
    tilt=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Tilt':tilt})  

    tempstring=header.split('DefectPosRotation:')[1] # 
    rot=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Rotation':rot})  

    tempstring=header.split('NoMapPixelsXY: ')[1]
    tempstring=tempstring.split(' ')[0]
    resolution=int(tempstring) # need resolution for jpg file creation (must be int)
    tempstr=tempstring+'x'+tempstring+' SE image'
    Augerparams.update({'Details':tempstr})  
    # call makejpg  pass: ) name for jpg file 2) field of view for jpg scaling 3) bindata
    jpgname=AugerFileName.replace('.sem','.jpg')
    if not os.path.isfile(jpgname): # creates jpg if not done already
        bindata=filedata[end+6:] # gets bytes holding actual data
        makejpg(bindata, jpgname, resolution, fieldofview) # extracts and saves RGB jpg file
    del(filedata) #binary byte object    
    del(header) # above truncated and converted to header string
    #  should spatial areas log and spectral details logs also be passed to main
    newname='sub/'+ AugerFileName
    os.rename(AugerFileName,newname) # move raw sem file to subfolder sub
    return Augerparams
    
def processMap(AugerFileName):   
    ''' pulls header info parameters from Auger .map files and saves in master param log'''
    with open(AugerFileName, 'rb') as file:
        filedata = file.read()
    end=filedata.find(b'EOFH')
    headerdata=filedata[0:end+6] # works to cut off binary part    
    header=headerdata.decode(encoding='cp437') # more generic encoding than utf-8
    # temporary dict for holding important file parameters
    Augerparams={'Filenumber':'','Project':'','Filename':'','FilePath':'','Sample':'','Comments':'','Date':'','FieldofView':'','Type':'','Energy':'','GunVoltage':'','nA':'','Areas':'','Cycles':'','Timestep':'','Details':'','Acqtime':'','Scanarea':'','X':'','Y':'','Z':'','Tilt':'','Rotation':'','ImageshiftX':'','ImageshiftY':''}
    
    Augerparams.update({'Filename':AugerFileName})
   # pull acquisition date 
    tempstring=header.split('FileDate: ')[1]
    match=re.search(r'(\d+\s\d+\s\d+)',tempstring)
    tempstring=match.group(0)
    if match:
        Augerparams.update({'Date':tempstring.split(' ')[1]+'/'+tempstring.split(' ')[2]+'/'+tempstring.split(' ')[0]}) # group0 returns matching string

    tempstring=header.split('AcqFilename: ')[1] #find file number
    tempstring=tempstring.split('.')[1] # convert to int to avoid later comparison problems
    if 'Tmp' in tempstring: # avoid error with Tmp1.sem problem caused by PHI duplicate numbering
        tempstring=tempstring.split('Tmp')[0] # remove Tmp string   
    tempint=int(tempstring)
    Augerparams.update({'Filenumber':tempint})
    
    tempstring=header.split('FieldOfView: ')[1] # find field of view
    match=re.search(r'\d+\.\d+',tempstring) 
    if match:
        tempstr=match.group(0)
        fieldofview=float(tempstr)
        Augerparams.update({'FieldofView':fieldofview}) # group0 returns matching string
        Augerparams.update({'Scanarea':fieldofview**2}) # simple determination of scanned area for maps 
    tempstring=header.split('ImageShift: ')[1] # find image shift x and y
    match=tempstring.split(' ')[0]
    tempstring=tempstring.split(' ')[1]
    Augerparams.update({'ImageshiftX':float(match)}) # group0 returns matching string
    Augerparams.update({'ImageshiftY':float(tempstring)})

    tempstring=header.split('EBeamEnergy: ')[1] # find kV
    tempstring=tempstring.split(' ')[0]
    Augerparams.update({'Energy':tempstring}) # group0 returns matching string

    tempstring=header.split('GunLensVoltage: ')[1] # find gunvoltage (related to current)
    tempstring=float(tempstring.split(' ')[0])
    Augerparams.update({'GunVoltage':tempstring})

    tempstring=header.split('TimePerStep: ')[1] # find time per step
    match=re.search(r'\d+\.\d+',tempstring)  
    if match:
        Augerparams.update({'Timestep':float(match.group(0))})    
        # find X,Y,Z, Rotation, Tilt
    tempstring=header.split('DefectPosX:')[1]
    xval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'X':xval})

    tempstring=header.split('DefectPosY:')[1] # 
    yval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Y':yval})  
 
    tempstring=header.split('DefectPosZ:')[1] # 
    zval=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Z':zval})  

    tempstring=header.split('DefectPosTilt:')[1] # 
    tilt=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Tilt':tilt})  

    tempstring=header.split('DefectPosRotation:')[1] # 
    rot=float(tempstring.split('\r\n')[0])
    Augerparams.update({'Rotation':rot})  

    # pull info about spectral regions (pass file number and header)
    # numcycles=0 # irrelevant so set to zero for mapsww
    filenumber=(Augerparams.get('Filenumber'))
    timestep=float((Augerparams.get('Timestep')))
    details, timeperarea =SpectralMapRegions(filenumber, AugerFileName, timestep, header)
    # no need for spatial areas determination for maps (entire field of view is used)
    # acquisition time for maps sum(sweeps for all EL)*2pt or 3pt * resolution**2 /(1000*60)
    Augerparams.update({'Details':details})

    # get necessary params for jpg element intensity map
    tempstring=header.split('NoMapPixelsXY: ')[1]
    tempstring=tempstring.split(' ')[0]
    resolution=int(tempstring) # need resolution for jpg file creation (must be int)

    # construct list of mapped elements from details (remove 'x's and numbers, leave spaces)
    tempstring=details.replace('x','') 
    mapelements= ''.join([i for i in tempstring if not i.isdigit()]) # list comprehension to remove numbers
    mapelements=mapelements.strip() # remove rogue leading/trailing spaces
    mapelements=[str(s) for s in mapelements.split(' ')] # list of mapped elements from above string
    jpgbasename=AugerFileName.replace('.map','')
    testjpgname=jpgbasename+'_'+mapelements[0]+'.jpg'
    if not os.path.isfile(testjpgname): # test if jpgs from maps were previously 
        bindata=filedata[end+6:] # gets bytes holding actual data    
        makemaps(bindata, jpgbasename, resolution, fieldofview, mapelements) # call function to make jpg intensity map for each elem
    acqtime=timeperarea  # gives overall time for map (only 1 area)
    Augerparams.update({'Acqtime':acqtime})    
    del(filedata) #binary byte object    
    del(header) # above truncated and converted to header string
    #  should spatial areas log and spectral details logs also be passed to main
    newname='sub/'+ AugerFileName
    os.rename(AugerFileName,newname) # move raw sem file to subfolder sub
    return Augerparams

def Augerbatchimport(filelist, Augerlogbook):
    ''' Main import processing loop for spe, map and sem images '''
    if not os.path.exists('sub'): # create subdirectory for raw spe/sem/map files & csv sub-files (when combined)
        os.makedirs('sub') # the process functions will move raw files here
    # Set up log file for multiplex spectral regions
    if not os.path.isfile('multiplexspectralregionslog.csv'):
        SpectralRegs=pd.DataFrame(columns=['Filenumber','Filename','Numcycles','Timestep','Element','Sweeps','Evstep','Lower','Upper','Time'])
        with open('multiplexspectralregionslog.csv','w') as filedata:
            SpectralRegs.to_csv(filedata, header=True, index=False)
    # Set up log for spatial areas info
    if not os.path.isfile('spatialareaslog.csv'):
        SpatialAreas=pd.DataFrame(columns=['Filenumber','Filename','Areanumber','X1','Y1','X2','Y2','Area']) # easier to use frame than dict
        with open('spatialareaslog.csv','w') as filedata:
            SpatialAreas.to_csv(filedata, header=True, index=False)
    # load AugerParamLog if prior run was made
    mycols=['Filenumber','Project','Filename','FilePath','Sample','Comments','Date','FieldofView','Type','Energy','GunVoltage','nA','Areas','Cycles','Timestep','Details','Evbreaks','Acqtime','Scanarea','X','Y','Z','Tilt','Rotation','ImageshiftX','ImageshiftY']
    if os.path.isfile('Augerparamlog.csv'):
        AugerParamLog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
    else: # create blank one
        AugerParamLog=pd.DataFrame(columns=mycols) # Create blank df log for Auger params
    # check for processed spe files in sub and add these to processed list
    if os.path.isfile('Augerparamlog_subs.csv'):
        AugerParamLogsub=pd.read_csv('Augerparamlog_subs.csv', encoding='cp437')
        AugerParamLog=pd.concat(AugerParamLog,AugerParamLogsub,ignore_index=True)
    # Process only files not previously processed and loaded into AugerParamLog
    unprocessed=[] # list for unprocessed files
    for i,AugerFileName in enumerate(filelist):
        AugerFileName=AugerFileName.replace('.spe','.csv')
        processed=AugerParamLog.Filename.unique()
        if AugerFileName not in processed:
            AugerFileName=AugerFileName.replace('.csv','.spe') # convert back to spe 
            unprocessed.append(AugerFileName)                    
    for i,AugerFileName in enumerate(unprocessed): # Now deal with only unprocessed files
        if AugerFileName.endswith('avg.spe'): # skip averaged spectra
            continue
        AugerFileParams={} # temp dictionary for this filenumber's Auger file params
        if AugerFileName.endswith('.spe'):
            try:            
                AugerFileParams = processAuger(AugerFileName) # retrieve IM file params as dictionary
                print(AugerFileName, ' processed.')
            except:
                print('Problem processing spe file ', AugerFileName)
                continue # move to next
        if AugerFileName.endswith('.sem'):
            try:            
                AugerFileParams = processSEM(AugerFileName) # retrieve params from sem file
                print(AugerFileName, ' processed.')
            except:
                print('Problem processing SEM file ', AugerFileName)
                continue # move to next
        if AugerFileName.endswith('.map'):
            try:
                AugerFileParams = processMap(AugerFileName) # retrieve params from sem file
                print(AugerFileName, ' processed.')
            except:
                print('Problem processing map file ', AugerFileName)
                continue # move to next            
        filepath=os.getcwd()
        AugerFileParams.update({'FilePath':filepath}) # same filepath for all 
        # Get associated project name, sample name, current (nA) and comments from auger filelog (based on filenumber)
        # names in imported xls log file must match the other fields in AugerFileParams
        filenumber=int(AugerFileParams.get('Filenumber'))  
        logrowmatch=getparamfromlog(AugerFileName, filenumber, Augerlogbook) # Use filename to avoid filenumber reset errors
        try:
            AugerFileParams.update({'Project':logrowmatch.get('Project'),'Sample':logrowmatch.get('Sample'),'Comments':logrowmatch.get('Comments'),'nA':logrowmatch.get('nA') })    
        except: 
            print("Couldn't find sample, comment, etc. for ", AugerFileName)
        
        # pass temp dictionary to dataframe row through use of Series (problems with passing scalars directly)
        AugerParams=pd.Series(AugerFileParams) # make series from dictionary 
        AugerParamLog=AugerParamLog.append(AugerParams, ignore_index=True) # append series as dataframe 
        # end of loop through unprocessed files
    AugerParamLog.Filenumber=AugerParamLog.Filenumber.apply(int) # converts all filenumbers to int (sometimes end up as str for unknown reason)
    AugerParamLog['Comments']=AugerParamLog['Comments'].replace(np.nan,'', regex=True) # avoids string search filtering errors
    AugerParamLog['Sample']=AugerParamLog['Sample'].replace(np.nan,'', regex=True) # avoids string search filtering errors
    AugerParamLog=AugerParamLog.sort_values(['Filenumber'], ascending=True) # sort by filenumber
    AugerParamLog=AugerParamLog.reset_index(drop=True)
    AugerParamLog=AugerParamLog[mycols] # put in original order
    return AugerParamLog

