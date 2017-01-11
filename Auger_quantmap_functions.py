    # -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 17:10:19 2016

@author: tkc
"""
import pandas as pd
import numpy as np
import shutil, sys, fileinput, os, math
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont # needed for jpg creation

#%%

# Quant map data processing done after normal batch import process (working on csv files)

def plotmaps(elementmaps, Elements, savename=''):
    ''' Plot arbitrary number of element maps (passed as list of numpy arrays) into single figure '''
    # determine which of the elementmaps will be plotted 
    nummaps=0    
    for i, [elem, thismap] in enumerate(elementmaps):
        if elem in Elements:
            nummaps+=1
    # Determine shape of figure
    if nummaps<=3:
        numrows=1
    else:
        numrows=2
    numcols=math.ceil(nummaps/numrows)
    
    # fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False)
    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, squeeze=False)

    for i, [elem, thismap] in enumerate(elementmaps): # thismap is element string followed by 512 x 512 array
        thisrow=i//numcols
        thiscol=i%numcols
        axindex=thisrow, thiscol # tuple to index axes 
        axes[axindex].set_aspect('equal')
        axes[axindex].set_title(elem)
        axes[axindex].imshow(thismap, cmap='hot') # plots element map to correct subplot
    fig.tight_layout()
    if savename!='':
        fig.savefig(savename) # optional saving of figure
    return

def overlayimage(elementmaps, elem, imagefile, savename='', alphaval=0.3):
    '''Find given element's quantmap (among list of created elementmaps) and overlay on pre-QM SE image, variable alpha value and 
    optional savename'''
    try:
        image=Image.open(imagefile)
    except:
        print ("Can't find ", imagefile)
        return
    plt.imshow(image)
    for i, (el, elemmap) in enumerate(elementmaps):
        if el==elem:
            thismap=elemmap
        else:
            print("Can't find ", elem, "map")
            return
    plt.imshow(thismap, alpha=alphaval)
    if savename!='':
        plt.savefig(savename) 
    return    
    
def createampmaps(filenum, Elements, Smdifpeakslog, SpatialAreas):
    '''Create elemental amplitude maps from QM elemental amplitudes (smooth-diff amplitude for each line in Smdifpeakslog and area information in 
    spatialareas '''
    elementmaps=[] # list of 512x512 arrays (oversampling actual 10x10 or 20x20 but more convenient)
    thisfiledat=Smdifpeakslog[Smdifpeakslog['Filenumber']==filenum]
    thisfileareas=SpatialAreas[SpatialAreas['Filenumber']==filenum]
    for i, elem in enumerate(Elements):
        thiselemdat=thisfiledat[thisfiledat['PeakID']==elem]
        if len(thiselemdat)==0: # check to ensure this element is present in quant map (multi-area multiplex)
            print(elem, ' not mapped in filenumber', filenum)
            continue
        elemarray = np.zeros((512,512)) # new np array for this element's map
        
        for index, row in thiselemdat.iterrows(): # loops through each areanumber
            areanum=thiselemdat.loc[index]['Areanumber']
            elemamp=thiselemdat.loc[index]['Amplitude'] # smooth- diff amplitude for this pixel range
            thispixarea=thisfileareas[thisfileareas['Areanumber']==areanum] # finds boundaries for this areanumber
            if len(thispixarea)!=1: # 
                print('Multiple area matches for ', str(filenum),' area ', str(areanum))
                continue
            else: # gets boundaries in 512 x 512 elemarray for this elemamp
                x1=int(thispixarea.iloc[0]['X1'])
                x2=int(thispixarea.iloc[0]['X2'])
                y1=int(thispixarea.iloc[0]['Y1'])
                y2=int(thispixarea.iloc[0]['Y2'])
                # now assign value to specified range within array
                for k in range(x1,x2):
                    for l in range(y1,y2):
                        elemarray[l][k]=elemamp # indexing is [col][row] not [row][col]
        elementmaps.append([elem, elemarray]) # add to list of element arrays
    return elementmaps # list of elem, array for all passed elements 

def showsubareas(df, SpatialAreas, image='tr184.147.jpg', savestr='Mgrichregion', label=True):
    '''After combination of QM data, find and show selected set of areanumbers (i.e Mgrich or whatever) 
    pass df subset based on some peak criteria'''
    # Put in a check to ensure this is a single file number/ single image
    filenum=int(df.iloc[0]['Filenumber'])
    filematch=SpatialAreas[SpatialAreas['Filenumber']==filenum]
    areas=np.ndarray.tolist(df.Areanumber.unique())
    areas=[int(i) for i in areas]
    areamatches=filematch[filematch['Areanumber'].isin(areas)]
    savename=savestr+str(filenum)+'.jpg'
    # open image
    try:
        jpgimage=Image.open(image)
    except:
        print("Couldn't open image ", image)
        return
    annotatejpg(jpgimage, areamatches, savename, label) # automatically saves image with correct ROIs
    return

def annotatejpg(jpgimage, areamatch, savename, label=True):
    '''Pass Auger sem image and info about spatial areas and create annotated jpg with superimposed ROIs'''
    draw=ImageDraw.Draw(jpgimage) # single draw instance to label above image
    ttfont=ImageFont.truetype('arial.ttf', size=20)
    for index, row in areamatch.iterrows():
        areanum=int(areamatch.loc[index]['Areanumber'])
        x1=areamatch.loc[index]['X1']
        y1=areamatch.loc[index]['Y1']
        x2=areamatch.loc[index]['X2']
        y2=areamatch.loc[index]['Y2']
        draw.rectangle((x1,y1,x2,y2), outline='red') # red rectangle at specified position
        if label:
            message=str(areanum) # label only with area number
            draw.text((x2+2,y2+2),message, font=ttfont, fill='red')
    jpgimage.convert('RGB').save(savename)
    jpgimage.close()
    return

def renumberQMareas(df, fnums,QMname=''):
    '''Duplicate spatial area entries in spatial areas to match those in renumbered QM file (i.e. area 1 of arrayfile2 is area 21'''
    startfile=int(fnums.split('-')[0])
    endfile=int(fnums.split('-')[1])
    filenumber=int(str(startfile)+str(endfile))
    filenums=list(range(startfile,endfile+1))
    areas=df[df['Filenumber'].isin(filenums)]
    areas=areas.sort_values(['Filenumber','Areanumber'])
    if QMname=='': # use default naming scheme if no new filename is passed
        # file format should be name.filenumber.csv...
        filename=areas.iloc[0]['Filename']
        if 'spe' in filename: # usually spe extensiion in spatialareaslog
            QMname=filename.replace('.spe',str(endfile)+'.spe') # name for quant map combined areas csv file
        else:
            QMname=filename.replace('.csv',str(endfile)+'.csv') # name for quant map combined areas csv file
    # just do a straight renumbering then concat and return
    for i in range(0, len(areas)):
        areas=areas.set_value(areas.index[i],'Filenumber', filenumber) # sets all to combined filenum
        areas=areas.set_value(areas.index[i],'Filename', QMname) # sets all to combined filenum
        areas=areas.set_value(areas.index[i],'Areanumber', i+1) # sets all to combined filenum
    df=pd.concat([df,areas],ignore_index=True)
    df.to_csv('spatialareaslog.csv', index=False)
    return df
    
def combineQMdata(df,fnums,QMname=''):
    ''' Pass sequences of spe files each containing 20 spatial areas, combine into single csv file; using
    autonumbering convention w/ unique filenumbers'''
    startfile=int(fnums.split('-')[0])
    endfile=int(fnums.split('-')[1])
    filenums=list(range(startfile,endfile+1))
    spefiles=df[df['Filenumber'].isin(filenums)]
    spefiles=spefiles.sort_values(['Filenumber'])
    if QMname=='': # use default naming scheme if no new filename is passed
        # file format should be name.filenumber.csv...
        filename=spefiles.iloc[0]['Filename']
        QMname=filename.replace('.csv',str(endfile)+'.csv') # name for quant map combined areas csv file
    # Check X, Y, multiplex details to ensure no error in pass file number list
    tempdf=spefiles[(spefiles['X']==spefiles.iloc[0]['X']) & (spefiles['Y']==spefiles.iloc[0]['Y']) & (spefiles['Details']==spefiles.iloc[0]['Details'])]
    if len(tempdf)!=len(spefiles):
        print('X, Y, multiplex details are not consistent... error in passed filelist?' )
        return
    # Working on converted csv files 
    filelist=np.ndarray.tolist(spefiles.Filename.unique()) # list of spe files for this quant map    
    firstfname=filelist[0]
    if os.path.exists(firstfname):
        newdf=pd.read_csv(firstfname)
    else:
        print("Error: Can't find ", firstfname)
    # keep only counts and S7D7 columns ... drop smcounts, savgol
    mycols = [c for c in newdf.columns if c.lower()[:6] != 'savgol'] # case insensitive
    newdf=newdf[mycols]
    mycols = [c for c in newdf.columns if c.lower()[:8] != 'smcounts'] # case insensitive
    newdf=newdf[mycols]
    currentarea=int(spefiles.iloc[0]['Areas'])+1 # for area renumbering (start is 1 more than #areas in file 1
    try:
        for i, spename in enumerate(filelist[1:]): # loop through 2nd through nth files
            if os.path.exists(spename):
                tempdf=pd.read_csv(spename) # loads current dataset (length match effectively checked already)
                mycols = [c for c in tempdf.columns if c.lower()[:6] != 'savgol'] # case insensitive
                tempdf=tempdf[mycols]
                mycols = [c for c in tempdf.columns if c.lower()[:8] != 'smcounts'] # case insensitive
                tempdf=tempdf[mycols]           
                thisfile=spefiles[spefiles['Filename']==spename]
                numnewareas=int(thisfile.iloc[0]['Areas'])
                for j in range(1,numnewareas+1):
                    oldcntname='Counts'+str(j)
                    oldS7D7name='S7D7'+str(j)
                    newcntname='Counts'+str(currentarea)
                    newS7D7name='S7D7'+str(currentarea)
                    tempdf.rename(columns={oldcntname:newcntname}, inplace=True)
                    tempdf.rename(columns={oldS7D7name:newS7D7name}, inplace=True)
                    currentarea=currentarea+1 # counter for areas
                # add new spatial areas to master file (merge on matching energy column)
                newdf=pd.merge(newdf,tempdf,how='left', on=['Energy']) 
    except:
        print('Problem combining spatial areas from all files.. missing files?')
        return 
    # Automatically save quantmap csv file with all areas from all sub-spe files    
    newentry=makeQMlogentry(spefiles, QMname) # new AugerParamLog entry for combined quant map 
    df=pd.concat([df,newentry], ignore_index=True) # append new entry to AugerParamLog
    df.to_csv('AugerParamLog.csv', index=False) # autosave of AugerParamLog with QM appended
    newdf.to_csv(QMname, index=False) # direct save of spe with all areas
    return df # returns modified AugerParamLog  (also added quant map file entry is appended and autosaved

def makeQMlogentry(spefiles, QMname):
    '''Make entry for Augerparam log for quant map dataset '''    
    newentry=spefiles.iloc[[0]] # first row as template
    newentry=newentry.reset_index(drop=True) # index carries no important info here
    # change  filename, areas, acqtime, scanarea to reflect entire quant map conditions
    newentry=newentry.set_value(0,'Areas',spefiles.Areas.sum()) # sum areas from each sub-spe file
    newentry=newentry.set_value(0,'Acqtime',spefiles.Acqtime.sum())
    newentry=newentry.set_value(0,'Scanarea',spefiles.Scanarea.sum())
    newentry=newentry.set_value(0,'Filename',QMname) # use new pass filename
    # Filenumber is firstnumlastnum combination
    newfilenum=int(str(spefiles.Filenumber.min())+str(spefiles.Filenumber.max())) # combine end to end
    newentry=newentry.set_value(0,'Filenumber',newfilenum)    
    newentry=newentry.set_value(0,'Filename',QMname)
    # Append quantmap string to current comments
    newentry.Comments=newentry.Comments.apply(str) # convert to string (in case all nan)
    tempstr=str(newentry.iloc[0]['Comments']) 
    if tempstr!='' and tempstr!='nan':
        newentry=newentry.set_value(0,'Comments',tempstr+' quantmap')
    else:
        newentry=newentry.set_value(0,'Comments','quantmap')      
    return newentry # returned as single rowed df 
    
def writeautotool(df, autotoolname):
    ''' Write of standard autotool loop for quantmap 
    weird encoding so just read and modify existing pff file '''
    # TODO pass as arg if others use this
    datastr=''
    datastr+='[AutoTool]\nTask Count='
    datastr+=str(len(df)) # 
    datastr+='\n'
    for index, row in df.iterrows():
        datastr+='Task '
        datastr+='%d' % index
        command=df.loc[index]['Command']
        datastr+='='
        datastr+=command
        datastr+='\n'
    datastr+='Data Count='
    datastr+=str(len(df)) # 
    datastr+='\n'    
    for index, row in df.iterrows():
        datastr+='Data '
        datastr+='%d' % index
        datastr+='='
        val=df.loc[index]['Data']
        if str(val)!='nan':
            datastr+=str(val) # could be int in some cases
        datastr+='\n'
        # Write this chunk of files to .phi spatial areas file (done w/ file replace method since encoding is weird unknown type)
    shutil.copyfile('C:\\Users\\tkc\\Documents\\Python_Scripts\\Utilities\\spatial_areas_sample_min.phi',autotoolname)
    for line in fileinput.input(autotoolname, inplace=1):
        sys.stdout.write(datastr)
    return

def makeautotool(filelist, multacq='C:\\Temp\\QM_multiplex.phi'):
    '''Generates df with Autotool commands and data values (for generating Autotool phi file)'''
    df=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Utilities\\QM_autotool.csv')
    atframe=df.loc[0:2]
    # filelist.sort(reverse=True) # sort just screws things up
    mycolumns=['Command','Data']
    for i, file in enumerate(filelist):        
        newrow=pd.DataFrame(index=np.arange(0,2), columns=mycolumns)
        newrow=newrow.set_value(0,'Command','AES:Load Area Define Setting...')
        newrow=newrow.set_value(0,'Data','C:\\Temp\\'+file)
        newrow=newrow.set_value(1,'Command','AES:Multiplex Acquire')
        atframe=pd.concat([atframe,newrow], ignore_index=True)
    atframe=pd.concat([atframe,df.loc[[3]]], ignore_index=True)
    return atframe

def makesquarearray(margin, arraysize, basename):
    ''' Divide up 512x512 pixels in map into n areas and format correctly for spatial areas phi files
    (which are loaded using Autotool loops into PHI Smartsoft); margin is % of field to skip mapping if desired '''
    pix=512
    width=(pix*(1-margin)/arraysize) # exact width (height is same)
    startxy=int(pix*margin/2) # split margin between top/bottom, left/right
    mycols=['Xindex','Yindex','Areanumber','PHIname','Subnumber','X1','Y1','X2','Y2', 'Width', 'Height']
    dim=arraysize**2
    square=pd.DataFrame(index=np.arange(0,dim), columns=mycols)
    for index,row in square.iterrows():
        xindex=index//arraysize # remainder is row
        yindex=index%arraysize # mod is correct column
        square.loc[index]['Xindex']=xindex # remainder is row
        square.loc[index]['Yindex']=yindex # mod is correct column
        left=int(width*yindex+startxy) # left-right position depending on column 
        square.loc[index]['X1']=left
        right=int(width*yindex+startxy+width)
        square.loc[index]['X2']=right
        top=int(width*xindex+startxy)
        square.loc[index]['Y1']=top # top-bottom position depending on row
        bottom=int(width*xindex+startxy+width)
        square.loc[index]['Y2']=bottom
        square.loc[index]['Width']=right-left # variations due to rounding error
        square.loc[index]['Height']=bottom-top
        # max of 20 areas allowed per spatial area .phi file
        square.loc[index]['Areanumber'] # true area number describing pix position after file combination
        square.loc[index]['Subnumber']=index%20 # Before combination w/ 20 areas per file
        filenum=index//20+1
        square.loc[index]['PHIname']=basename+str(filenum)
    filelist=square.PHIname.unique()
    filelist=np.ndarray.tolist(filelist)
    for i, fname in enumerate(filelist):
        thisfile=square[square['PHIname']==fname]
        writeAESareas(thisfile, fname) # writes each list of 25 areas to separate .phi text file
    atframe=makeautotool(filelist, multacq='QM_multiplex')
    # instead of C:\Temp copy multiplex and spatial areas files to Smartsoft settings folders
    return square, atframe

def writeAESareas(df, PHIname):
    ''' Write of stage positions to .phi spatial areas file in chunks of 25 positions
    Some weird encoding so just read and modify existing pff file '''    
    datastr=''
    datastr+='[SpatialArea]\nArea Count='
    datastr+=str(len(df)) # number of areas
    datastr+='\n'
    for i in range(0,len(df)):
        datastr+='Area Active '
        datastr+='%d' % i
        datastr+='=True\n'
    for i in range(0,len(df)):
        datastr+='Area Mode '
        datastr+='%d' % i
        datastr+='=Area\n'
    for i in range(0,len(df)):
        datastr+='Area Left '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['X1']
        datastr+='%d' % val
        datastr+='\n'
    for i in range(0,len(df)):
        datastr+='Area Top '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Y1']
        datastr+='%d' % val
        datastr+='\n' 
    for i in range(0,len(df)):
        datastr+='Area Right '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['X2']
        datastr+='%d' % val
        datastr+='\n'         
    for i in range(0,len(df)):
        datastr+='Area Bottom '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Y2']
        datastr+='%d' % val
        datastr+='\n'      
    for i in range(0,len(df)):
        datastr+='Area Width '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Width']
        datastr+='%d' % val
        datastr+='\n'      
    for i in range(0,len(df)):
        datastr+='Area Height '
        datastr+='%d' % i
        val=df.iloc[i]['Height']
        datastr+='='
        datastr+='%d' % val
        datastr+='\n' 
    # Write this chunk of files to .phi spatial areas file (done w/ file replace method since encoding is weird unknown type)
    filename=PHIname+'.phi'
    shutil.copyfile('C:\\Users\\tkc\\Documents\\Python_Scripts\\Utilities\\spatial_areas_sample_min.phi',filename)
    for line in fileinput.input(filename, inplace=1):
        sys.stdout.write(datastr)
    return

def writemultiplex(PHIname, dwelltime=20, numcycles=1, reginterval=2, regmode='Areas'):
    ''' Write of multiplex settings file (max 20 regions)
    Some weird encoding so just read and modify existing pff file 
    image registration choices are cycles or areas, '''
    df=pd.read_csv('C:\\Users\\tkc\Documents\\Python_Scripts\\Utilities\\multiplex_acq.csv')
    # Eliminate any duplication in scan range (check/adjust lower upper ranges)
    for i in range(0, len(df)-1):
        priorend=df.iloc[i]['Upper']
        nextstart=df.iloc[i+1]['Lower']
        if nextstart<priorend:
            df=df.set_value(df.index[i+1],'Lower',priorend+1)
    # TODO pass as arg if others use this 
    datastr=''
    datastr+='[MultiplexAcq]\nRegister Image=True\nImage Registration Interval='
    datastr+=str(reginterval)
    datastr+='\n'
    datastr+='Image Registration Mode='
    datastr+=regmode
    datastr+='\nTime Per Step (ms)='
    datastr+='%0.1f' % dwelltime
    datastr+='\nNegative Values=Allow\nNumber of Regions='
    datastr+=str(len(df)) # number of multiplex regions
    datastr+='\nAtomic Number List Count='
    datastr+=str(len(df))
    datastr+='\n'
    for i in range(0,len(df)):
        datastr+='Atomic Number List '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['AtomNum']
        datastr+='%d' % val
        datastr+='\n'
    datastr+='Element List Count='
    datastr+=str(len(df))
    datastr+='\n'
    for i in range(0,len(df)):
        datastr+='Element List '
        datastr+='%d' % i
        datastr+='='
        strval=df.iloc[i]['Elem']
        datastr+=strval
        datastr+='\n'
    datastr+='Active Count='
    datastr+=str(len(df))
    datastr+='\n'
    for i in range(0,len(df)):
        datastr+='Active '
        datastr+='%d' % i
        datastr+='=True\n'  # won't be present in Df if not active
    datastr+='Sweeps Count='
    datastr+=str(len(df))
    datastr+='\n' 
    for i in range(0,len(df)):
        datastr+='Sweeps '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Sweeps']
        datastr+='%d' % val
        datastr+='\n'
    datastr+='Lower Acq Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Lower Acq '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Lower']
        datastr+='%0.1f' % val # float with tenths place precision
        datastr+='\n'
    datastr+='Upper Acq Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Upper Acq '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Upper']
        datastr+='%0.1f' % val
        datastr+='\n'        
    datastr+='Acq Range Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Acq Range '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Range']
        datastr+='%0.1f' % val
        datastr+='\n'    
    datastr+='Lower Analysis Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Lower Analysis '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Lowpeak']
        datastr+='%0.1f' % val
        datastr+='\n'     
    datastr+='Upper Analysis Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Upper Analysis '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Hipeak']
        datastr+='%0.1f' % val
        datastr+='\n'    
    datastr+='eV Per Step Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='eV Per Step '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['EVstep']
        datastr+='%0.1f' % val
        datastr+='\n'   
    datastr+='Peak Energy Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Peak Energy '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Peak']
        datastr+='%0.1f' % val
        datastr+='\n'  
    datastr+='Background 1 Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Background 1 '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Back1']
        datastr+='%0.1f' % val
        datastr+='\n'  
    datastr+='Background 2 Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Background 2 '
        datastr+='%d' % i
        datastr+='='
        val=df.iloc[i]['Back2']
        datastr+='%0.1f' % val
        datastr+='\n'  
    datastr+='Number Of Channels Count='
    datastr+=str(len(df))
    datastr+='\n'     
    for i in range(0,len(df)):
        datastr+='Number Of Channels '
        datastr+='%d' % i
        datastr+='=1 thru 8\n'
    datastr+='Number of Cycles='
    datastr+=str(numcycles)
    datastr+='\n'
    # Write this chunk of files to .phi spatial areas file (done w/ file replace method since encoding is weird unknown type)
    shutil.copyfile('C:\\Users\\tkc\\Documents\\Python_Scripts\\Utilities\\spatial_areas_sample_min.phi',PHIname)
    for line in fileinput.input(PHIname, inplace=1):
        sys.stdout.write(datastr)
    return