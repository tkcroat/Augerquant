    # -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 17:10:19 2016

@author: tkc
"""
import pandas as pd
import numpy as np
import shutil, sys, fileinput, os, math, re
import matplotlib.pyplot as plt
import struct
from PIL import Image, ImageDraw, ImageFont # needed for jpg creation
from matplotlib.backends.backend_pdf import PdfPages
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules')
from Auger_batch_import_functions import SpectralRegions
# functions used by plotpixel
from Auger_plot_functions import combineranges, findelemlines
import tkinter as tk
import glob
from operator import itemgetter
from itertools import groupby
from scipy.ndimage import filters  # uniform filter and others
from decimal import Decimal
import datetime
import random
from scipy.signal import argrelmin
from tkinter import filedialog

#%%

plt.rcParams.update({'font.size':12})
plt.rcParams.update({'legend.fontsize':12})

# Quant map data processing done after normal batch import process (working on csv files)

def multi_shift_setup(Elements, AESquantparams, QMpixarray):
    ''' Interactive set up parameters for quantmap 
    differs from QMmultiplexsetup -- allows different energy shifts for different X,Y pixels within map
    max of 20 spectral regions for multiplex means limit of only 6 elems 
    Saves custom QMmultiplex.phi file (normally needed)
    (incl low and high energy regions for each)
    Params:  array size (typically max 100x100 or Autotool will crash)
        dwelltime - in ms
        multiplex filename -- contains multiplex scan ranges (can directly alter these files) 
        margin - percent of image field at margins (not scanned)
    kwargs:
        image registration
        
        '''
    # Get range of peak shifts directly from QMpixarray (shifts determined from O map)
    if len(Elements)>6:
        print('Limit of 6 elements for conventional multiplex setup') 
        print('manually combine or revise multiplex regions')
    root = tk.Tk()
    root.title("QM quant map setup of Autotool & multiplex files")
    '''
    mypath=tk.StringVar()
    mypath= filedialog.askdirectory()
    tk.Label(root, 'Save directory')
    '''
    dwell=tk.IntVar()
    dwell.set(20)
    numcycles=tk.IntVar()
    numcycles.set(3)    
    peakshiftmin=tk.IntVar()
    # use mean of first 20 areas as min value
    peakshiftmin.set(int(QMpixarray.iloc[0:20]['Shift'].mean()))
    peakshiftmax=tk.IntVar()
    # Use mean of # areas in last spatial area file as max (usually 20 but could be less)
    if len(QMpixarray)%20==0:
        numfiles=20
    else:
        numfiles=len(QMpixarray)%20
    peakshiftmax.set(int(QMpixarray.iloc[-numfiles:-1]['Shift'].mean()))
    peakshiftincr=tk.IntVar() # if adjusting 
    peakshiftincr.set(4)
    peakwidth=tk.IntVar()
    peakwidth.set(7)
    regintmult=tk.IntVar() # image registration interval (if done w/in multiplex)
    regintmult.set(0)
    regintAT=tk.IntVar()  # frequency of image registration embedded in autotool
    regintAT.set(2)
    phiname=tk.StringVar()
    phiname.set('QMmultiplex.phi')
    timeest=tk.DoubleVar()
    choice=tk.StringVar() # abort or generate
    
    # Initialize w/ shift of zero
    # elem name, lower region center, peak center, upper reg center, # sweeps
    elemparams=getQMparams(Elements, AESquantparams,0)
    rownum=0
    # First col are energy region/element labels (not alterable)
    tk.Label(root, text='Region').grid(row=rownum, column=0)
    tk.Label(root, text='Peak').grid(row=rownum, column=1)
    tk.Label(root, text='Sweeps').grid(row=rownum, column=2)
    tk.Label(root, text='Width').grid(row=rownum, column=3)
    rownum+=1
    # Display of array size and other params
    tk.Label(root, text='Dwell time (ms)').grid(row=0, column=4)
    tk.Entry(root, textvariable=dwell).grid(row=0, column=5)
    tk.Label(root, text='# of cycles').grid(row=1, column=4)
    tk.Entry(root, textvariable=numcycles).grid(row=1, column=5)
    # Display peak shifts determined from QMpixarray
    tk.Label(root, text='Peak shift min').grid(row=2, column=4)
    tk.Entry(root, textvariable=peakshiftmin).grid(row=2, column=5)
    tk.Label(root, text='Peak shift max').grid(row=3, column=4)
    tk.Entry(root, textvariable=peakshiftmax).grid(row=3, column=5)
    tk.Label(root, text='Peak shift increment').grid(row=4, column=4)
    tk.Entry(root, textvariable=peakshiftincr).grid(row=4, column=5)
    tk.Label(root, text='Half width').grid(row=5, column=4)
    tk.Entry(root, textvariable=peakwidth).grid(row=5, column=5)
    tk.Label(root, text='Image reg multiplex interval').grid(row=6, column=4)
    tk.Entry(root, textvariable=regintmult).grid(row=6, column=5)
    tk.Label(root, text='Image reg Autotool interval').grid(row=7, column=4)
    tk.Entry(root, textvariable=regintAT).grid(row=7, column=5)
    tk.Label(root, text='Multiplex name').grid(row=8, column=4)
    tk.Entry(root, textvariable=phiname).grid(row=8, column=5)
    tk.Label(root, text='Time (hrs)').grid(row=9, column=4)
    tk.Entry(root, textvariable=timeest).grid(row=9, column=5)

    # tk var not needed for peak/ element regions ... just display
    peaks=[]
    sweeps=[]
    widths=[]

    # Create lists of tk string variables
    for i in range(0,3*len(Elements)):
        peaks.append(tk.DoubleVar())
        sweeps.append(tk.IntVar())
        widths.append(tk.IntVar())

    def initvals():
        ''' update all string variable values based on shift, global width'''
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
            peaks[3*i].set(elemparams[i][1])
            widths[3*i].set(7)
            sweeps[3*i].set(1)
            peaks[3*i+1].set(elemparams[i][2])
            widths[3*i+1].set(7)
            sweeps[3*i+1].set(elemparams[i][4])
            peaks[3*i+2].set(elemparams[i][3])
            widths[3*i+2].set(7)
            sweeps[3*i+2].set(1)
    
    def updatevals(event):
        ''' update peak positions and widths based on peak shift and global width
        use min shift in cases where range of shifts will be used '''

        for i, elem in enumerate(elemparams):
            peaks[3*i].set(elemparams[i][1]+peakshiftmin.get())
            widths[3*i].set(peakwidth.get())
            peaks[3*i+1].set(elemparams[i][2]+peakshiftmin.get())
            widths[3*i+1].set(peakwidth.get())
            peaks[3*i+2].set(elemparams[i][3]+peakshiftmin.get())
            widths[3*i+2].set(peakwidth.get())
        # acquisition time est. in hours
        numchannels=0
        for i in range(0,len(widths)):
            numchannels+=int(sweeps[i].get())*(int(widths[i].get())*2+1)
        timeest.set(int(numcycles.get())*numchannels*int(dwell.get())/1000)

    initvals() # should initialize values with zero shift
    
    # Display peak params
    for i, elem in enumerate(Elements):
        # first list element regions (i.e. SiL, Si, SiH)
        tempstr=elem+'L'
        tk.Label(root, text=tempstr).grid(row=rownum, column=0)
        tk.Entry(root, textvariable=peaks[3*i], width=7).grid(row=rownum, column=1)
        tk.Entry(root, textvariable=sweeps[3*i], width=7).grid(row=rownum, column=2)
        tk.Entry(root, textvariable=widths[3*i], width=7).grid(row=rownum, column=3)
        rownum+=1
        # now handle peak line itself
        tk.Label(root, text=elem).grid(row=rownum, column=0)
        tk.Entry(root, textvariable=peaks[3*i+1], width=7).grid(row=rownum, column=1)
        tk.Entry(root, textvariable=sweeps[3*i+1], width=7).grid(row=rownum, column=2)
        tk.Entry(root, textvariable=widths[3*i+1], width=7).grid(row=rownum, column=3)
        rownum+=1
        tempstr=elem+'H'
        tk.Label(root, text=tempstr).grid(row=rownum, column=0)
        tk.Entry(root, textvariable=peaks[3*i+2], width=7).grid(row=rownum, column=1)
        tk.Entry(root, textvariable=sweeps[3*i+2], width=7).grid(row=rownum, column=2)
        tk.Entry(root, textvariable=widths[3*i+2], width=7).grid(row=rownum, column=3)
        rownum+=1

    def abort(event):
        choice.set('abort')        
        root.destroy()  

    def plot(event):
        choice.set('generate')        
        root.destroy()  
        
    a=tk.Button(root, text='Update')
    a.bind('<Button-1>', updatevals)
    a.grid(row=rownum, column=0)

    a=tk.Button(root, text='Abort')
    a.bind('<Button-1>', abort)
    a.grid(row=rownum, column=1)
    
    a=tk.Button(root, text='Generate')
    a.bind('<Button-1>', plot)
    a.grid(row=rownum, column=2)
    
    root.mainloop()
    
    if choice.get()=='generate':
        mkwargs={}
        # Pass params to write multiplex
        phibasename=phiname.get()
        dwell=dwell.get()
        numcycles=numcycles.get()
        # handle peak shifts built in to multiplex
        # optional embedding of image registration within multiplex
        if regintmult.get()>0: 
            # only for embed into multiplex itself (normally use zero and embed imagereg into autotool)
            mkwargs.update({'regmode':'Areas'})
            mkwargs.update({'reginterval':regintmult.get()})
        # Note ... usually faster to embed image registrations in Autotool loop
        peaks=convert_tklist(peaks)
        sweeps=convert_tklist(sweeps)
        widths=convert_tklist(widths)
        # add multiplexname col to QMpixarray
        if 'Multiplexname' not in QMpixarray.columns:
            QMpixarray['Multiplexname']=''
        # Figure out binning of peak shifts
        peakmin=peakshiftmin.get()
        peakmax=peakshiftmax.get()
        peakincr=peakshiftincr.get()
        shiftvals=[i for i in range(peakmin, peakmax+1, peakincr)]
        # Multiplex breaks need to line up with existing spatial area file breaks (each 20 regions)
        multibreaks=[0] # area number of break
        multinames=[] # names of multiplex files for loading
        for i, shift in enumerate(shiftvals):
            peakshifted=[val+shift-shiftvals[0] for val in peaks]
            phiname=phibasename.split('.')[0]+'_'+str(shift)+'.phi'
            # Start at min shift value and loop
            multdf=makemultdf(elemparams, peakshifted, sweeps, widths)
            # Save modified multiplex scan for QM
            writemultiplex(multdf, phiname, dwell, numcycles, **mkwargs)
            # Assign multiplex to QMpixarray (then use for autotool creation)
            match=QMpixarray[(QMpixarray['Shift']>=shift) & (QMpixarray['Shift']<=shift+peakincr)]
            multibreaks.append(match.index.max()+1) # area numbers (i.e index # -1)
            multinames.append(phiname)
        # align multiplex breaks w/ 20 spatial area breaks (set row to nearest 20)
        multibreaks=[(i//20) for i in multibreaks]
        fbreaks=[20*i for i in multibreaks] # need for write of multiplex file to QMpixarray
        for i in range(0, len(fbreaks)-1):
            match=QMpixarray[(QMpixarray['Areanumber']>fbreaks[i]) & (QMpixarray['Areanumber']<=fbreaks[i+1])]
            for index, row in match.iterrows():
                QMpixarray=QMpixarray.set_value(index,'Multiplexname', multinames[i])
        # Make updated associated autotool file
        filelist=np.ndarray.tolist(QMpixarray.PHIname.unique()) # list of spatial area files
        # two types of regint are possible (in multiplex or in autotool)    
        atkwargs={}
        atkwargs.update({'multibreaks': multibreaks}) # spatial areas at which to insert multiplex breaks
        atkwargs.update({'multinames': multinames})
        if regintAT.get()!=0:
            atkwargs.update({'regint': int(regintAT.get())})
        # for changing multiplex, pass breakpoints between different multiplex and associated names
        atframe=makeautotool(filelist, multacq='QM_multiplex', **atkwargs)
        ATname='QM_autotool_shifted.phi'
        writeautotool(atframe, ATname) 
        # Now write non-sequential spatial areas (in sets of 20)
        for i, fname in enumerate(filelist):
            thisfile=QMpixarray[QMpixarray['PHIname']==fname]
            writeAESareas(thisfile, fname) 
    return QMpixarray

def QMmultiplex_setup(Elements, AESquantparams):
    ''' Interactive set up parameters for quantmap 
    max of 20 spectral regions for multiplex means limit of only 6 elems 
    Saves custom QMmultiplex.phi file (normally needed)
    (incl low and high energy regions for each)
    Params:  array size (typically max 100x100 or Autotool will crash)
        dwelltime - in ms
        multiplex filename -- contains multiplex scan ranges (can directly alter these files) 
        margin - percent of image field at margins (not scanned)
    kwargs: 
        image registration
        
        '''
    if len(Elements)>6:
        print('Limit of 6 elements for conventional multiplex setup') 
        print('manually combine or revise multiplex regions')
    root = tk.Tk()
    root.title("QM quant map setup of Autotool & multiplex files")
    '''
    mypath=tk.StringVar()
    mypath= filedialog.askdirectory()
    tk.Label(root, 'Save directory')
    '''
    dwell=tk.IntVar()
    dwell.set(20)
    numcycles=tk.IntVar()
    numcycles.set(3)
    peakshift=tk.IntVar() # globally appliable peak shift (and no local shift option)
    peakshift.set(0)
    peakwidth=tk.IntVar() # globally appliable peak shift
    peakwidth.set(7)
    regint=tk.IntVar() # image registration interval (if done w/in multiplex)
    regint.set(0)
    phiname=tk.StringVar()
    phiname.set('QMmultiplex.phi')
    timeest=tk.DoubleVar()
    choice=tk.StringVar() # abort or generate
    
    # Initialize w/ shift of zero
    # elem name, lower region center, peak center, upper reg center, # sweeps
    elemparams=getQMparams(Elements, AESquantparams,0)
    rownum=0
    # First col are energy region/element labels (not alterable)
    tk.Label(root, text='Region').grid(row=rownum, column=0)
    tk.Label(root, text='Peak').grid(row=rownum, column=1)
    tk.Label(root, text='Lower').grid(row=rownum, column=2)
    tk.Label(root, text='Upper').grid(row=rownum, column=3)
    tk.Label(root, text='Sweeps').grid(row=rownum, column=4)
    tk.Label(root, text='Half-width').grid(row=rownum, column=5)
    rownum+=1
    # Display of array size and other params
    tk.Label(root, text='Dwell time (ms)').grid(row=0, column=6)
    tk.Entry(root, textvariable=dwell).grid(row=0, column=7)
    tk.Label(root, text='# of cycles').grid(row=1, column=6)
    tk.Entry(root, textvariable=numcycles).grid(row=1, column=7)
    tk.Label(root, text='Peak shift').grid(row=2, column=6)
    tk.Entry(root, textvariable=peakshift).grid(row=2, column=7)
    tk.Label(root, text='Half width').grid(row=3, column=6)
    tk.Entry(root, textvariable=peakwidth).grid(row=3, column=7)
    tk.Label(root, text='Image reg interval').grid(row=4, column=6)
    tk.Entry(root, textvariable=regint).grid(row=4, column=7)
    tk.Label(root, text='Multiplex name').grid(row=5, column=6)
    tk.Entry(root, textvariable=phiname).grid(row=5, column=7)
    tk.Label(root, text='Time (hrs)').grid(row=6, column=6)
    tk.Entry(root, textvariable=timeest).grid(row=6, column=7)

    # tk var not needed for peak/ element regions ... just display
    peaks=[]
    sweeps=[]
    widths=[]
    lowers=[] # lower ev of scan range
    uppers=[]

    # Create lists of tk string variables
    for i in range(0,3*len(Elements)):
        peaks.append(tk.DoubleVar()) # center of scan window
        lowers.append(tk.DoubleVar()) # lower limit of scan window
        uppers.append(tk.DoubleVar())
        sweeps.append(tk.IntVar())
        widths.append(tk.IntVar())

    def recalc():
        ''' Update time estimates after any value changes (not tied to button) but called from 
        within every method
        '''
        # acquisition time est. in hours
        numchannels=0
        for i in range(0,len(widths)):
            numchannels+=int(sweeps[i].get())*(int(widths[i].get())*2+1)
        timeest.set(int(numcycles.get())*numchannels*int(dwell.get())/1000)

    def updatewidth():
        ''' Use current lowers/uppers values to calculate each local width '''
        # first reset centers/ranges to initial values 
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
            if int(uppers[3*i].get() - lowers[3*i].get())/2 > 0:
                widths[3*i].set(int(uppers[3*i].get() - lowers[3*i].get())/2 )
            else:
                widths[3*i].set(0)
            if int(uppers[3*i+1].get() - lowers[3*i+1].get())/2 >0:   
                widths[3*i+1].set(int(uppers[3*i+1].get() - lowers[3*i+1].get())/2)
            else:
                widths[3*i+1].set(0)
            if int(uppers[3*i+2].get() - lowers[3*i+2].get())/2 >0:
                widths[3*i+2].set(int(uppers[3*i+2].get() - lowers[3*i+2].get())/2)
            else:
                widths[3*i+2].set(0)
    
    def adjustpeakcenter():
        ''' Squeezing of background regions changes effective center  '''
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
            peaks[3*i].set((uppers[3*i].get() + lowers[3*i].get())/2) # lower window for element
            peaks[3*i+2].set((uppers[3*i+2].get() + lowers[3*i+2].get())/2)
 
    def checkoverlaps():
        ''' After any changes adjust lower and upper backgrounds if overlapping w/ main region '''
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
            if uppers[3*i].get() >= lowers[3*i+1].get():
                uppers[3*i].set( lowers[3*i+1].get() - 1)
            if lowers[3*i+2].get() <= uppers[3*i+1].get():
                lowers[3*i+2].set(uppers[3*i+1].get() + 1)
        # TODO width of OL or OH cannot be negative (constraint in updatewidth?)
        # Changing lowers/ uppers change resulting width
        updatewidth()
        # changing lowers uppers for lowback and highback regions can change peakcenter
        adjustpeakcenter()
        
    def initvals():
        ''' update all string variable values based on shift, global width'''
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
            peaks[3*i].set(int(elemparams[i][1])) # lower window for element
            lowers[3*i].set(int(elemparams[i][1]-7))
            uppers[3*i].set(int(elemparams[i][1]+7))
            widths[3*i].set(7)
            sweeps[3*i].set(1)
            peaks[3*i+1].set(int(elemparams[i][2])) # peak itself
            lowers[3*i+1].set(int(elemparams[i][2]-7))
            uppers[3*i+1].set(int(elemparams[i][2]+7))
            widths[3*i+1].set(7)
            sweeps[3*i+1].set(int(elemparams[i][4]))
            peaks[3*i+2].set(int(elemparams[i][3])) # upper window for element
            lowers[3*i+2].set(int(elemparams[i][3]-7))
            uppers[3*i+2].set(int(elemparams[i][3]+7))
            widths[3*i+2].set(7)
            sweeps[3*i+2].set(1)
        checkoverlaps()
        recalc()

    initvals() # initialize values with zero shift from AESquantparams

    # Display peak params (centers/ lower ev/ upper ev/ width/
    for i, elem in enumerate(Elements):
        # first list element regions (i.e. SiL, Si, SiH)
        tempstr=elem+'L'
        tk.Label(root, text=tempstr).grid(row=rownum, column=0)
        tk.Entry(root, textvariable=peaks[3*i], width=7).grid(row=rownum, column=1)
        tk.Entry(root, textvariable=lowers[3*i], width=7).grid(row=rownum, column=2)
        tk.Entry(root, textvariable=uppers[3*i], width=7).grid(row=rownum, column=3)
        tk.Entry(root, textvariable=sweeps[3*i], width=7).grid(row=rownum, column=4)
        tk.Entry(root, textvariable=widths[3*i], width=7).grid(row=rownum, column=5)
        rownum+=1
        # now handle peak line itself
        tk.Label(root, text=elem).grid(row=rownum, column=0)
        tk.Entry(root, textvariable=peaks[3*i+1], width=7).grid(row=rownum, column=1)
        tk.Entry(root, textvariable=lowers[3*i+1], width=7).grid(row=rownum, column=2)
        tk.Entry(root, textvariable=uppers[3*i+1], width=7).grid(row=rownum, column=3)
        tk.Entry(root, textvariable=sweeps[3*i+1], width=7).grid(row=rownum, column=4)
        tk.Entry(root, textvariable=widths[3*i+1], width=7).grid(row=rownum, column=5)
        rownum+=1
        tempstr=elem+'H'
        tk.Label(root, text=tempstr).grid(row=rownum, column=0)
        tk.Entry(root, textvariable=peaks[3*i+2], width=7).grid(row=rownum, column=1)
        tk.Entry(root, textvariable=lowers[3*i+2], width=7).grid(row=rownum, column=2)
        tk.Entry(root, textvariable=uppers[3*i+2], width=7).grid(row=rownum, column=3)
        tk.Entry(root, textvariable=sweeps[3*i+2], width=7).grid(row=rownum, column=4)
        tk.Entry(root, textvariable=widths[3*i+2], width=7).grid(row=rownum, column=5)
        rownum+=1

    def applyglobals(event):
        ''' Update each peak positions and widths based on global peak shift and global width '''

        for i, elem in enumerate(elemparams):
            peaks[3*i].set(int(elemparams[i][1]+peakshift.get()))
            lowers[3*i].set(int(elemparams[i][1]-peakwidth.get() +peakshift.get()))
            uppers[3*i].set(int(elemparams[i][1] + peakwidth.get() +peakshift.get()))
            widths[3*i].set(peakwidth.get())
            peaks[3*i+1].set(int(elemparams[i][2]+peakshift.get()))
            lowers[3*i+1].set(int(elemparams[i][2]-peakwidth.get() +peakshift.get()))
            uppers[3*i+1].set(int(elemparams[i][2] + peakwidth.get() +peakshift.get()))
            widths[3*i+1].set(peakwidth.get())
            peaks[3*i+2].set(int(elemparams[i][3]+peakshift.get()))
            lowers[3*i+2].set(int(elemparams[i][3]-peakwidth.get() +peakshift.get()))
            uppers[3*i+2].set(int(elemparams[i][3] + peakwidth.get() +peakshift.get()))
            widths[3*i+2].set(peakwidth.get())
        checkoverlaps()
        recalc() # update time estimate
        
    def updatepeaks(event):
        ''' Use local widths (and global shift) to update peak centers and ranges 
        no vals for local energy shift but could change peaks/lowers/uppers '''
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
            peaks[3*i].set(int(elemparams[i][1] + peakshift.get() )) # lower window for element
            lowers[3*i].set(int(elemparams[i][1]+ peakshift.get() - widths[3*i].get()))
            uppers[3*i].set(int(elemparams[i][1]+ peakshift.get() + widths[3*i].get()))
            peaks[3*i+1].set(int(elemparams[i][2] + peakshift.get() )) # peak itself
            lowers[3*i+1].set(int(elemparams[i][2] + peakshift.get() - widths[3*i+1].get()))
            uppers[3*i+1].set(int(elemparams[i][2] + peakshift.get() + widths[3*i+1].get()))
            peaks[3*i+2].set(int(elemparams[i][3]  + peakshift.get())) # upper window for element
            lowers[3*i+2].set(int(elemparams[i][3]+ peakshift.get() - widths[3*i+2].get()))
            uppers[3*i+2].set(int(elemparams[i][3]+ peakshift.get() + widths[3*i+2].get()))
        
        checkoverlaps()
        # acquisition time est. in hours
        recalc()

    def updatewidths(event):
        ''' Use current lowers/uppers values to calculate each local width '''
        # first reset centers/ranges to initial values 
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
            widths[3*i].set(int(uppers[3*i].get() - lowers[3*i].get())/2 )
            widths[3*i+1].set(int(uppers[3*i+1].get() - lowers[3*i+1].get())/2)
            widths[3*i+2].set(int(uppers[3*i+2].get() - lowers[3*i+2].get())/2)
        checkoverlaps()
        # acquisition time est. in hours
        recalc()

    def recalcupdate(event):
        ''' Update time estimates, widths, adjust overlap boundaries
        '''
        checkoverlaps()
        recalc()
        
    def abort(event):
        choice.set('abort')        
        root.destroy()  

    def plot(event):
        choice.set('generate')        
        root.destroy()  
    
    tk.Label(root, text='Known values + global adjustments').grid(row=rownum, column=0, columnspan=2)
    a=tk.Button(root, text='Apply global shift/width')
    a.bind('<Button-1>', applyglobals)
    a.grid(row=rownum, column=2)
    rownum+=1

    tk.Label(root, text='Make local peak changes').grid(row=rownum, column=0, columnspan=2)
    a=tk.Button(root, text='Change ranges using widths')
    a.bind('<Button-1>', updatepeaks)
    a.grid(row=rownum, column=2)

    a=tk.Button(root, text='Update widths from ranges')
    a.bind('<Button-1>', updatewidths)
    a.grid(row=rownum, column=3)

    a=tk.Button(root, text='Recalc/update')
    a.bind('<Button-1>', recalcupdate)
    a.grid(row=rownum, column=4)
    rownum+=1

    tk.Label(root, text='Create custom multiplex scan').grid(row=rownum, column=0, columnspan=2)    
    a=tk.Button(root, text='Abort')
    a.bind('<Button-1>', abort)
    a.grid(row=rownum, column=2)
    
    a=tk.Button(root, text='Generate')
    a.bind('<Button-1>', plot)
    a.grid(row=rownum, column=3)
    
    root.mainloop()
    
    if choice.get()=='generate':
        mkwargs={}
        # Pass params to write multiplex
        phiname=phiname.get()
        dwell=dwell.get()
        numcycles=numcycles.get()
        # optional embedding of image registration within multiplex
        if regint.get()>0: #
            mkwargs.update({'regmode':'Areas'})
            mkwargs.update({'reginterval':regint.get()})
        # Note ... usually faster to embed image registrations in Autotool loop
        peaks=convert_tklist(peaks)
        sweeps=convert_tklist(sweeps)
        widths=convert_tklist(widths)
        lowers=convert_tklist(lowers)
        uppers=convert_tklist(uppers)
        multdf=makemultdf(elemparams, peaks, lowers, uppers, sweeps, widths)
        # Save modified multiplex scan for QM
        writemultiplex(multdf, phiname, dwell, numcycles, **mkwargs)
    return multdf

def makerectarray(rectsize, startpix, arraywidth, basename, **kwargs):
    ''' Divide up 512x512 pixels in map into n areas and format correctly for spatial areas phi files
    (which are loaded using Autotool loops into PHI Smartsoft); 
    Mapping proceed horizontally (x horizontal, y vertical) 
    ARGS:
    rectsize -
    arraywidth -- # of pixels wide
    margin - % of image field that is unmapped (20% means 10% of 512 field at both edges (aka 51 pixels))
        is unmapped
    basename - basename for area definition files 
        e.g. "50x50array20m" basemname makes files 50x50array20m1, 50x50array20m2, etc.
        
    KWARGS: 'regint' - interval at which to build in image registration into autotool loop; val of 1 means
           register every 20 pixels (since each area file holds 20 defined spatial areas);  passed to makeautotool
           * this is best way to incorporate image registration in quantmap process... more flexible
           interval allowed; if instead one builds in image reg into multiplex itself, one has to 
           run image registration much more frequently which unnecessarily increases acquisition time
        'writeareas' (bool) - write new spatial area definition files (old ones can be reused if same arraysize 
             and margin)
        'writeAutotool' --  write of new Autotool sequence (can be reused if same arraysize/regint)
        scrambled -- randomize scan to minimize charging 
    '''
    arrayheight=math.ceil((rectsize[1]/rectsize[0])*arraywidth)
    dim=arraywidth*arrayheight # length of rect 
    mycols=['Xindex','Yindex','Areanumber','PHIname','Subnumber','X1','Y1','X2','Y2', 'Width', 'Height']
    # Rect is the pixarray file correlating 101.spe with associated pixel in quantmap 
    rect=pd.DataFrame(index=np.arange(0,dim), columns=mycols)
    ''' width/height of single scan pixel in terms of 512x512 pixel field
    can be non-integer (in which case pixels in 512 scan can be slightly variable width 
    depending on rounding '''
    scanpixsize=(512*(rectsize[0])/arraywidth) 
    
    # x is horizontal axis and mapping proceeds by going across top row
    for index,row in rect.iterrows():
        xindex=index//arraywidth # remainder is row (0th is top row)
        yindex=index%arraywidth # mod is correct column (0th is left column)
        rect.loc[index]['Xindex']=xindex # remainder is row
        rect.loc[index]['Yindex']=yindex # mod is correct column
        left=int(scanpixsize*yindex+startpix[0]) # left-right position depends on column 
        rect.loc[index]['X1']=left
        right=int(scanpixsize*(yindex+1)+startpix[0])
        rect.loc[index]['X2']=right
        top=int(scanpixsize*xindex+startpix[1]) 
        rect.loc[index]['Y1']=top # top-bottom position depends on row
        bottom=int(scanpixsize*(xindex+1)+startpix[1])
        rect.loc[index]['Y2']=bottom
        rect.loc[index]['Width']=right-left # variations due to rounding error
        rect.loc[index]['Height']=bottom-top
        # true area number describing pix position after file combination
        rect.loc[index]['Areanumber']=index+1
        # max of 20 areas allowed per spatial area .phi file
        rect.loc[index]['Subnumber']=index%20 # Before combination w/ 20 areas per file
        filenum=index//20+1 # filenumber of multiplex 
        rect.loc[index]['PHIname']=basename+str(filenum)
    filelist=rect.PHIname.unique()
    filelist=np.ndarray.tolist(filelist)
    # map/move beam non-sequentially to minimize local charging
    if 'scrambled' in kwargs:
        areanums=np.ndarray.tolist(rect.Areanumber.unique())
        areanums=np.random.permutation(areanums)
        rect['Areanumber']=pd.Series(areanums)
        rect=rect.sort_values(['Areanumber'])
        # need to reassign subnumber and PHIname
        rect=rect.reset_index(drop=True)
        for index,row in rect.iterrows():
            rect=rect.set_value(index,'Subnumber', index%20)
            rect=rect.set_value(index,'PHIname', basename+str(index//20+1))
    if 'writeareas' in kwargs:
        for i, fname in enumerate(filelist):
            thisfile=rect[rect['PHIname']==fname]
            writeAESareas(thisfile, fname) # Writes each list of 20 areas to separate .phi text file
    if 'writeAutotool' in kwargs:
        atkwargs={}
        if 'regint' in kwargs: # pass image registration interval 
            atkwargs.update({'regint':kwargs.get('regint',0)})
        # TODO preload multiplex file and skip adding multacq step?
        atframe=makeautotool(filelist, multacq='QM_multiplex', **atkwargs)
        # Use autonaming convention for autotool file
        ATname='QM_autotool_rect'+str(arraywidth)+'x'+str(arraywidth)+'.phi'
        writeautotool(atframe, ATname) 
    # instead of C:\Temp copy multiplex and spatial areas files to Smartsoft settings folders
    return rect

def QMrectarray_setup():
    ''' Interactive set up of pixarray and autotool file for quantmap if using new params

    Params: array size (typically max 100x100 or Autotool will crash)
            margin - percent of image field at margins (not scanned)
    kwargs: 
        image registration
        scrambled -- jump around w/ beam to avoid excessive charging 
    '''
    root = tk.Tk()
    root.title("QM rectangular array and custom Autotool setup")
    arrwidth=tk.IntVar()
    arrwidth.set(100)
    rectwidth=tk.IntVar() # width of scan area in 512 pix coords
    rectwidth.set(0)
    rectheight=tk.IntVar()
    rectheight.set(0)
    xstartpix=tk.IntVar() # left edge of scan in 512 pix coords
    xstartpix.set(0)
    ystartpix=tk.IntVar() # top edge of scan in 512 pix coords
    ystartpix.set(0)
    regint=tk.IntVar()
    regint.set(2)
    atbool=tk.BooleanVar() # make new Autotool file
    atbool.set(False)
    areabool=tk.BooleanVar() # make new spatial area files
    areabool.set(False)
    scramblebool=tk.BooleanVar() 
    scramblebool.set(False)
    basename=tk.StringVar()
    basename.set('rect_array')
    choice=tk.StringVar() # abort or generate
    
    # TODO set up label and entry stuff
    rownum=0
    tk.Label(root, text='Array width (~100 max)').grid(row=rownum, column=2)
    tk.Entry(root, textvariable=arrwidth, width=14).grid(row=rownum, column=0)
    rownum+=1
    tk.Label(root, text='Rect width and height (in scan pix)').grid(row=rownum, column=2)
    tk.Entry(root, textvariable=rectwidth, width=7).grid(row=rownum, column=0)
    tk.Entry(root, textvariable=rectheight, width=7).grid(row=rownum, column=1)
    rownum+=1
    tk.Label(root, text='Starting pixel x,y coord (in scan pix)').grid(row=rownum, column=2)
    tk.Entry(root, textvariable=xstartpix, width=7).grid(row=rownum, column=0)
    tk.Entry(root, textvariable=ystartpix, width=7).grid(row=rownum, column=1)
    rownum+=1

    tk.Label(root, text='Image reg interval (Autotool)').grid(row=rownum, column=2)
    tk.Entry(root, textvariable=regint, width=7).grid(row=rownum, column=0)
    rownum+=1
    tk.Label(root, text='Basename for files').grid(row=rownum, column=2)
    tk.Entry(root, textvariable=basename, width=10).grid(row=rownum, column=0)
    rownum+=1
    tk.Checkbutton(root, variable=atbool, text='Make new Autotool file?').grid(row=rownum, column=0)
    rownum+=1
    tk.Checkbutton(root, variable=scramblebool, text='randomly move beam to minimize charging?').grid(row=rownum, column=0)
    rownum+=1
    tk.Checkbutton(root, variable=areabool, text='Create new spatial area files?').grid(row=rownum, column=0)
    rownum+=1
    tk.Label(root, text='If same parameters, existing Autotool and spatial area files can be reused.').grid(row=rownum, column=0, columnspan=3)
    rownum+=1    
    def abort(event):
        choice.set('abort')
        root.destroy()

    def create(event):
        choice.set('create')        
        root.destroy()  

    a=tk.Button(root, text='Abort')
    a.bind('<Button-1>', abort)
    a.grid(row=rownum, column=0)
    
    a=tk.Button(root, text='Create')
    a.bind('<Button-1>', create)
    a.grid(row=rownum, column=1)
    
    root.mainloop()
    
    if choice.get()=='create':
        kwargs={}
        startpix=(xstartpix.get(), ystartpix.get()) # coords of upper left corner of scan
        rectsize=(rectwidth.get()/512, rectheight.get()/512) # convert to fraction of area 
        basename=basename.get() # convert to string
        # Always create square pixarray file
        if atbool.get(): # optional create new AutoTool file
            kwargs.update({'writeAutotool':True})
        if areabool.get(): # optional create new spatial areas file
            kwargs.update({'writeareas':True})        
        if regint.get()>0:
            kwargs.update({'regint':regint.get()})
        if scramblebool.get(): # scramble spatial areas/ hop beam to avoid local charging
            kwargs.update({'scrambled':True})    
        # create pixarray (and optional create of Autotool and spatial areas files) 
        pixarray=makerectarray(rectsize, startpix, arrwidth.get(), basename, **kwargs)
        return pixarray
    else:
        return pd.DataFrame()
    return

def QMarray_setup():
    ''' Interactive set up of pixarray and autotool file for quantmap if using new params

    Params: array size (typically max 100x100 or Autotool will crash)
            margin - percent of image field at margins (not scanned)
    kwargs: 
        image registration
        scrambled -- jump around w/ beam to avoid excessive charging 
    '''
    root = tk.Tk()
    root.title("QM pixel array and custom Autotool setup")
    arrsize=tk.IntVar() # Number of pixels for QM scan (usually differs from 512 pix of imaging)
    arrsize.set(100)
    margin=tk.DoubleVar()
    margin.set(0.2)
    regint=tk.IntVar()
    regint.set(2)
    atbool=tk.BooleanVar() # make new Autotool file
    atbool.set(False)
    areabool=tk.BooleanVar() # make new spatial area files
    areabool.set(False)
    scramblebool=tk.BooleanVar() 
    scramblebool.set(False)
    basename=tk.StringVar()
    basename.set('100x100array20marg')
    choice=tk.StringVar() # abort or generate
    
    # TODO set up label and entry stuff
    rownum=0
    tk.Label(root, text='Array Size (~100 x 100 max)').grid(row=rownum, column=1)
    tk.Entry(root, textvariable=arrsize).grid(row=rownum, column=0)
    rownum+=1
    tk.Label(root, text='Unmapped margin (fraction)').grid(row=rownum, column=1)
    tk.Entry(root, textvariable=margin).grid(row=rownum, column=0)
    rownum+=1
    tk.Label(root, text='Image reg interval (Autotool)').grid(row=rownum, column=1)
    tk.Entry(root, textvariable=regint).grid(row=rownum, column=0)
    rownum+=1
    tk.Label(root, text='Basename for files').grid(row=rownum, column=1)
    tk.Entry(root, textvariable=basename).grid(row=rownum, column=0)
    rownum+=1
    tk.Checkbutton(root, variable=atbool, text='Make new Autotool file?').grid(row=rownum, column=0)
    rownum+=1
    tk.Checkbutton(root, variable=scramblebool, text='hop beam to minimize charging?').grid(row=rownum, column=0)
    rownum+=1
    tk.Checkbutton(root, variable=areabool, text='Create new spatial area files?').grid(row=rownum, column=0)
    rownum+=1
    tk.Label(root, text='If same parameters, existing Autotool and spatial area files can be reused.').grid(row=rownum, column=0, columnspan=3)
    rownum+=1    
    def abort(event):
        choice.set('abort')
        root.destroy()

    def create(event):
        choice.set('create')        
        root.destroy()  

    a=tk.Button(root, text='Abort')
    a.bind('<Button-1>', abort)
    a.grid(row=rownum, column=0)
    
    a=tk.Button(root, text='Create')
    a.bind('<Button-1>', create)
    a.grid(row=rownum, column=1)
    
    root.mainloop()
    
    if choice.get()=='create':
        kwargs={}
        arraysize=arrsize.get()
        margin=margin.get()
        basename=basename.get() # convert to string
        # Always create square pixarray file
        if atbool.get(): # optional create new AutoTool file
            kwargs.update({'writeAutotool':True})
        if areabool.get(): # optional create new spatial areas file
            kwargs.update({'writeareas':True})        
        if regint.get()>0:
            kwargs.update({'regint':regint.get()})
        if scramblebool.get(): # scramble spatial areas/ hop beam to avoid local charging
            kwargs.update({'scrambled':True})    
        # create pixarray (and optional create of Autotool and spatial areas files) 
        pixarray=makesquarearray(margin, arraysize, basename, **kwargs)
        return pixarray
    else:
        return pd.DataFrame()
    return



def convert_tklist(tklist):
    ''' convert list of tk inter variables to normal python list '''
    newlist=[]
    for i, val in enumerate(tklist):
        newlist.append(val.get())
    return newlist

def getQMparams(Elements, AESquantparams, shift):
    ''' retrieve direct peak location, low & high energy regions (background below and 
    above peak), and default # of sweeps for passed elements '''
    elemparams=[]
    for i, elem in enumerate(Elements):
        match=AESquantparams[AESquantparams['element']==elem]
        if len(match)==1:
            lowback=round(match.iloc[0]['QMlow']+shift,1)
            hiback=round(match.iloc[0]['QMhigh']+shift,1)
            # Direct peak location val is relative to negpeak
            peak=round(match.iloc[0]['negpeak']+match.iloc[0]['integpeak']+shift,1)
            # Default # of sweeps for actual peak
            sweeps=int(match.iloc[0]['QMsweeps'])
            atno=int(match.iloc[0]['atno']) # atomic number
            elemparams.append([elem, lowback, peak, hiback, sweeps, atno])
    return elemparams

def findelems(Elements, xrange, AESquantparams):
    ''' Find list of elements in this plotrange (for lookup of shifts, sums  
    returns element and its ideal direct peak counts position '''
    peaks=AESquantparams[AESquantparams['element'].isin(Elements)] # select rows in element list
    peaks=peaks[(peaks['negpeak']>xrange[0]) & (peaks['negpeak']<xrange[1])]
    peaks['integpeak']=peaks['integpeak']+peaks['negpeak']
    elems=[]
    for index, row in peaks.iterrows():
        elems.append([peaks.loc[index]['element'],peaks.loc[index]['integpeak']])    
    return elems

def showpeakampls(amplmapdict, Elements):
    ''' Plot/show all summed images from spectral image in same window 
    '''
    numcols=min(len(Elements),2) # 1 or 2 columns
    numrows=math.ceil(len(Elements)/2)    
    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
    plotnum=0
    fig.tight_layout()
    for i, elem in enumerate(Elements):
        image=amplmapdict.get(elem) # element name is dictionary key
        # Last element is the subtracted data 
        thisrow=(plotnum)//2 # correct zero based indexing
        thiscol=(plotnum)%2
        axindex=thisrow,thiscol        
        axes[axindex].imshow(image, cmap='hot')
        # Label subplot
        axes[axindex].set_title(elem)
        plotnum+=1
    return 

def showintegimages(sumdict, Elemdata):
    ''' Plot/show all summed images from spectral image in same window 
    '''
    numcols=min(len(Elemdata),2) # 1 or 2 columns
    numrows=math.ceil(len(Elemdata)/2)    
    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
    plotnum=0
    fig.tight_layout()
    for key, val in sumdict.items():
        image=sumdict.get(key)
        # Last element is the subtracted data 
        image=image[:,:,2]
        thisrow=(plotnum)//2 # correct zero based indexing
        thiscol=(plotnum)%2
        axindex=thisrow,thiscol        
        axes[axindex].imshow(image, cmap='hot')
        # Label subplot
        axes[axindex].set_title(key)
        plotnum+=1
    return 

                
def showshiftimages(shiftdict, peakstats):
    ''' Plot/show all shift images from spectral image in same window 
    peakstats has:  elem, avg shift, stdev of shift, min and max shift '''
    numcols=min(len(shiftdict),2) # 1 or 2 columns
    numrows=math.ceil(len(shiftdict)/2)    
    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
    plotnum=0
    fig.tight_layout()
    for key, val in shiftdict.items():
        image=shiftdict.get(key)
        thisrow=(plotnum)//2 # correct zero based indexing
        thiscol=(plotnum)%2
        axindex=thisrow,thiscol        
        axes[axindex].imshow(image, cmap='hot')
        thispeak=peakstats.get(key)
        thislabel= key + ': '+ str(thispeak[0])+' +/- '+str(thispeak[1])+ ' eV'
        # Label subplot
        axes[axindex].set_title(thislabel)
        plotnum+=1
    return 

def getregboundaries(spectralregs, Elements):
    ''' Get energy boundaries for each elemental subregion (incl. adjacent low/high background) '''
    plotranges=[]
    for i, elem in enumerate(Elements):
        matchL=spectralregs[spectralregs['Element']==elem+'L']
        matchH=spectralregs[spectralregs['Element']==elem+'H']
        if len(matchL)==1 and len(matchH)==1:
            lower=matchL.iloc[0]['Lower']
            upper=matchH.iloc[0]['Upper']
            plotranges.append([lower,upper])
        else: # handle single wide spectral range scan
            match=spectralregs[spectralregs['Element']==elem]
            if len(match)==1:
                lower=match.iloc[0]['Lower']
                upper=match.iloc[0]['Upper']
                plotranges.append([lower,upper])            
            else:
                print("Couldn't find spectral range for", elem)
    return plotranges

def pixelreport_tk(specimage, pixlist, energy, Elemdata, spectralregs, amplmaps, integmaps, shiftmaps, AESquantparams, **kwargs):
    ''' tk interface for args/kwargs of quantmap pixel subset plots 
	 filtering by filenumber, filename incorporated in tk interface 
    all args/dataframes must be passed through to plot functions 
    kwargs: fileset, areas, xrangestr -- usually just holds values entered during prior run 
    report will use all areas '''
    # first print out existing info in various lines
    root = tk.Tk()
    root.title('Generate PDF plot report of quantmap pixels')
    xrangestr=tk.StringVar()  # energy range in eV 
    xrangestr.set(kwargs.get('xrangestr',''))
    plottype=tk.StringVar() # column for plot report/ radio1 choice
    plottype.set('Both')
    smdifbool=tk.BooleanVar() # Bool for plotting of peak locations in smdifpeakslog
    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    backptbool=tk.BooleanVar() 
    elemstr=tk.StringVar()
    maxareas=tk.StringVar()
    maxareas.set('1')
    pixliststr=tk.StringVar()
    # Unpack list of lists from pixlist
    pixlist=[str(sublist[0])+','+str(sublist[1]) for sublist in pixlist for item in sublist]
    pixlist=set(pixlist)
    pixlist=list(pixlist)
    pixliststr.set(';'.join(pixlist))    
    Elements=[i[0] for i in Elemdata]
    mytext=', '.join(Elements) # elements for labelling 
    elemstr.set(mytext)
    plotneighborsbool=tk.BooleanVar()
    plotelemsbool=tk.BooleanVar()  # optional labelling of elements
    plotelemsbool.set(True) # default to true
    now=datetime.datetime.now()
    PDFname='Countsback_report'+'_'+now.strftime('%d%b%y')+'.pdf'
    PDFnamestr=tk.StringVar()  # energy range in eV 
    PDFnamestr.set(PDFname)
    chargeshift=tk.IntVar() # optional shift of element peaks from charging
    chargeshift.set(0)
    choice=tk.StringVar()  # plot or abortw
    # Functions to enable/disable relevant checkboxes depending on radiobutton choice
    def Countopts():
        ''' Disable irrelevant checkboxes '''
        smdiff.config(state=tk.DISABLED)
        backfit.config(state=tk.NORMAL)
        backpt.config(state=tk.NORMAL)
        PDFnamestr.set('Counts_report'+'_'+now.strftime('%d%b%y')+'.pdf')
    def Derivopts():
        ''' Disable irrelevant checkboxes '''
        smdiff.config(state=tk.NORMAL)
        backfit.config(state=tk.DISABLED)
        backpt.config(state=tk.DISABLED)
        PDFnamestr.set('Deriv_report'+'_'+now.strftime('%d%b%y')+'.pdf')
    def Bothopts():
        ''' Disable irrelevant checkboxes '''
        smdiff.config(state=tk.NORMAL)
        backfit.config(state=tk.NORMAL)
        backpt.config(state=tk.NORMAL)
        PDFnamestr.set('Countsderiv_report'+'_'+now.strftime('%d%b%y')+'.pdf')
    # Choose counts, deriv, both or peaks plot
    rownum=0
    tk.Label(root, text='Report name').grid(row=rownum, column=0)
    tk.Entry(root, textvariable=PDFnamestr, width=25).grid(row=rownum, column=1)    
    rownum+=1
    tk.Radiobutton(root, text='Counts', value='Counts', variable = plottype, 
                   command=Countopts).grid(row=rownum, column=0)
    tk.Radiobutton(root, text='Derivative', value='S7D7', variable = plottype, 
                   command=Derivopts).grid(row=rownum, column=1)
    tk.Radiobutton(root, text='Both', value='Both', variable = plottype, 
                   command=Bothopts).grid(row=rownum, column=2)
    rownum+=1
    backfit=tk.Checkbutton(root, variable=backfitbool, text='Plot background fits?')
    backfit.grid(row=rownum, column=0) # can't do immediate grid or nonetype is returned
    backpt=tk.Checkbutton(root, variable=backptbool, text='Plot background points?')
    backpt.grid(row=rownum, column=1) 
    smdiff=tk.Checkbutton(root, variable=smdifbool, text='Plot smdiff peak locations?')
    smdiff.grid(row=rownum, column=2)
    rownum+=1

    # Mark element positions (option)
    tk.Checkbutton(root, variable=plotelemsbool, 
                   text='Label element peaks?').grid(row=rownum, column=0)
    tk.Label(root, text='Elements to plot:').grid(row=rownum, column=1)
    elementry=tk.Entry(root, textvariable=elemstr)
    elementry.grid(row=rownum, column=2)
    rownum+=1
    tk.Label(root,text='List of plotted pixels:').grid(row=rownum, column=0)
    tk.Entry(root, textvariable=pixliststr).grid(row=rownum, column=1)
    tk.Checkbutton(root,variable=plotneighborsbool, 
                   text='Plot single pixel neighbors?').grid(row=rownum, column=2)
    rownum+=1
    # range choice for quant map should just be elements (data from Elemdata)
    tk.Label(root, text='Charging shift').grid(row=rownum, column=0)
    tk.Entry(root, textvariable=chargeshift).grid(row=rownum, column=1)
    tempstr='Number of pixels: '+str(int(len(pixlist)))
    tk.Label(root, text=tempstr).grid(row=rownum, column=2)
    
    rownum+=1

    tk.Label(root, text='Max # areas per page').grid(row=rownum, column=0)
    tk.Entry(root, textvariable=maxareas).grid(row=rownum, column=1)
    rownum+=1
    # Combine averaged files only (checkbox)?
   
    # option to reselect labeled elemental peaks 

    def abort(event):
        choice.set('abort')
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    a=tk.Button(root, text='Plot')
    a.bind('<Button-1>', plot)
    a.grid(row=rownum, column=1)
    
    a=tk.Button(root, text='Abort')
    a.bind('<Button-1>', abort)
    a.grid(row=rownum, column=2)

    root.mainloop()
    mychoice=choice.get()
    
    if mychoice=='plot':
        if chargeshift.get()!=0:
            AESqp=shiftAESparams(AESquantparams, chargeshift.get())
        else:
            AESqp=AESquantparams.copy()
        myPDFname=PDFnamestr.get()
        kwargs={} # reconstruct with new choices
        # set up parameters common to all reporting functions
        if maxareas.get()!='':
            maxareas=int(maxareas.get()) # max # areas on single plot page
            kwargs.update({'maxareaplots': maxareas})
        # Plot range will just be entered elements
        try:
            Elements=elemstr.get().split(',')
            Elements=[i.strip() for i in Elements]
            myelemdata=[l for l in Elemdata if l[0] in Elements]
        except:
            print('Problem extracting elements from string.')
            pass
        if pixliststr.get()!='':
            try: # Multiple values
                pixlist=pixliststr.get().split(';')
            except:
                pixlist=[pixliststr.get()]
            # converts comma sep sublist 1,20 into list [1,20]
            pixlist=[[int(i) for i in i.split(',')] for i in pixlist]
        else: # no pixels specified so report all 
            pixlist=[]
        if plotneighborsbool.get():
            [xcent,ycent]=pixlist[0]
            pixlist=[]
            for X in range(xcent-1,xcent+2):
                for Y in range(ycent-1,ycent+2):
                    if X in range(0, specimage.shape[0]) and Y in range(0, specimage.shape[1]):
                        pixlist.append([X,Y])
        plottype=plottype.get() # counts, deriv or both
        if backfitbool.get():
            kwargs.update({'addbackfit':True}) # plotting of background
        if backptbool.get():
            kwargs.update({'addbackpts':True}) # plot points used to create background
        if smdifbool.get():
            kwargs.update({'smdifpts':True}) # plot pts used to fit background 
            # these are retrieved from backfitlog
        if plotelemsbool.get(): # optional labeling of elements
            kwargs.update({'plotelems':Elements}) # pass elements list 
        # TODO modify so that low/high background regions are also plotted 
        reportpixels(plottype, myelemdata, specimage, pixlist, energy, amplmaps, 
                     integmaps, shiftmaps, AESqp, myPDFname, **kwargs)
    return kwargs

'''TESTING
[elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, 
                idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth, 
                kfact, kfact2, mass]=Elemdata[0]
page=0
maxareas=5
[X,Y]=pixlist[0]
'''
 
def reportpixels(plottype, myelemdata, specimage, pixlist, energy, amplmaps, 
                     integmaps, shiftmaps, AESquantparams, myPDFname, **kwargs):
    ''' Generate PDF plot report of underlying data from a subset of pixels
    Elemdata already has index and ev ranges for each element
    shiftmaps -- 0) deriv-based peak 1) integ-based peak
    amplmaps - 0) amplitude, 1) negpeak value, 2) negpeak ev 3) width (to pospeak)
    integmaps - 0 peak integration 1) peak ev 2) countsmax (at peak) 3) slope
        4)  intercept
    TODO figure out how to use backarray, shiftdict, sumdict 
    kwargs:
        addbackfit
        smdifpts -- pass smdiflog in kwargs and optional plot of negpeak/pospeak '''
    
    plt.ioff()
    with PdfPages(myPDFname) as pdf:
        numareas=len(pixlist)
        maxareas=kwargs.get('maxareaplots',1)
        # determine plot structure (rows and columns)
        if plottype=='Both':
            numrows=2 # counts on top and S7D7 on bottom
            numcols=len(myelemdata) # number of elemental subregions for plotting
        else: # Count or derivative plot only
            numcols=min(len(myelemdata),2) # 1 or 2 columns
            numrows=math.ceil(len(myelemdata)/2)   
        for page in range(0, numareas, maxareas): # single plot for every n areas/pixel positions
            fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), 
                                squeeze=False) # 2 by ? axes array
            firstrun=True # add element label on first pass for single page
            if len(myelemdata)>4:
                plt.tight_layout() # shrinks to fit axis labels
            mylegend=[]
            for i, [X, Y] in enumerate(pixlist[page:page+maxareas]):
                Augerfile=pd.DataFrame()
                Augerfile['Energy']=energy # full energy list over all regions
                rawdata=specimage[X,Y,:]
                Augerfile['Counts']=rawdata
                # Deriv/s7d7 must be separately calculated for each elem region due to evbreaks
                if Augerfile.empty:
                    continue
                '''
                if kwargs.get('backfitpts', False):
                    indexptslist=getbackfitpts(Backfitlog, AugerFileName, areanum)
                '''
                # TODO removing single element from passed list screws up retrieval of 
                # smdifpts and slope/intercept
                mylegend.append('Pixel '+str(X)+', '+ str(Y))
                for j, [elem, order, peakind, lowind, highind, peakrange, lowrange, 
                        hirange, idealev, idealind, idealnegpeak, integwidth, 
                        chargeshift, peakwidth, searchwidth, kfact, kfact2, 
                        mass] in enumerate(myelemdata):                   
                    # if plotting background in counts window, get slope/intercept for this pixel/this element
                    if kwargs.get('addbackfit', False): 
                        [integcnts,integev, countmax,slope, intercept] = integmaps[j][X, Y]
                        xvals=np.linspace(int(lowrange[0]),int(hirange[1]),100)
                    if kwargs.get('smdifpts', False): # unpack values from amplmap
                        [ampl, negval, negev, width]= amplmaps[j][X, Y]
                        posev=negev-width
                        posval=ampl+negval
                        derivx=[posev, negev]
                        derivy=[posval, negval]
                    if plottype=='Counts': # single counts or s7d7/deriv report
                        # work out axes which depend on # of separate spectral lines/regions
                        thisrow=j%numrows
                        thiscol=j//numrows
                        # plot of each element should extend from lowrange[0] to hirange[1].. not just narrower peakrange
                        Augerslice=Augerfile[(Augerfile['Energy']>lowrange[0]) & (Augerfile['Energy']<hirange[1])]
                        Augerslice.plot(x='Energy', y=plottype, ax=axes[thisrow,thiscol])
                        if kwargs.get('addbackfit', False): # optional plot of background fits
                            axes[thisrow,thiscol].plot(xvals, xvals*slope+intercept)
                    elif plottype=='S7D7': # single counts or s7d7/deriv report
                        # work out axes which depend on # of separate spectral lines/regions
                        thisrow=j%numrows
                        thiscol=j//numrows
                        # Need to calculate this derivative in piece-wise manner
                        thiscnts=rawdata[lowind[0]:highind[1]+1]
                        thiss7d7=calcderiv(thiscnts, peakind, lowind, highind, peakrange, lowrange, hirange)
                        xvals=range(int(lowrange[0]),int(hirange[1])+1)
                        axes[thisrow,thiscol].plot(xvals, thiss7d7)
                        # optional plotting of markers w/ negpeak/ pospeak
                        if kwargs.get('smdifpts', False):
                            axes[thisrow,thiscol].scatter(derivx, derivy)
                    elif plottype=='Both': # counts and s7d7/deriv
                        thiscol=j # row 1 is counts, row 2 is deriv
                        thiscnts=rawdata[lowind[0]:highind[1]+1] # parse counts to only this peak
                        thiss7d7=calcderiv(thiscnts, peakind, lowind, highind, peakrange, lowrange, hirange)
                        xvals=range(int(lowrange[0]),int(hirange[1])+1)
                        axes[1,thiscol].plot(xvals, thiss7d7)
                        Augerslice=Augerfile[(Augerfile['Energy']>lowrange[0]) & (Augerfile['Energy']<hirange[1])] # already known that this isn't empty
                        Augerslice.plot(x='Energy', y='Counts', ax=axes[0,thiscol]) # counts in lower 
                        if kwargs.get('addbackfit', False)==True and plottype=='Counts': # optional plot of background fits
                            axes[0,thiscol].plot(xvals, xvals*slope+intercept, color='r')
                        if kwargs.get('smdifpts', False):
                            axes[1,thiscol].scatter(derivx, derivy)
                    if firstrun: # only for first pixel on page
                        # Section for labeling plot elements
                        print('labelling elements for ', elem)
                        if 'plotelems' in kwargs:
                            plotelems=kwargs.get('plotelems',[]) # list of element peaks for labeling
                            if plottype!='Both':
                                if plottype=='Counts': # get ideal peak for direct
                                    print('Elements for labeling  are', plotelems)
                                    elemlines=getelemenergy(plotelems, peakrange, AESquantparams, deriv=False)
                                    print('Element lines are', elemlines)
                                else: # get ideal peak for deriv
                                    print('Elements for labeling  are', plotelems)
                                    elemlines=getelemenergy(plotelems, peakrange, AESquantparams, deriv=True)
                                    print('Elements are', elemlines)
                                # list of tuples with energy,elemname
                                for k, elemtuple in enumerate(elemlines):
                                    # elemtuple[0] is energy and [1] is element symbol
                                    # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                    try:
                                        axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                        axes[thisrow,thiscol].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                                    except:
                                        print('Problem labeling elements')
                            else: # handle both in same plot
                                elemlines=getelemenergy(plotelems, peakrange, AESquantparams, deriv=True) # can pass plot range as lower,upper tuple
                                print('Elements are', elemlines)
                                # list of tuples with energy,elemname
                                for k, elemtuple in enumerate(elemlines):
                                    # elemtuple[0] is energy and [1] is element symbol
                                    # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                    try:
                                        axes[1,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                        axes[1,thiscol].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                                    except:
                                        print('Problem labeling elements')
                                # same labeling of elements in upper plot (w/ peaks slightly shifted)
                                elemlines=getelemenergy(plotelems, peakrange, AESquantparams, deriv=False) # can pass plot range as lower,upper tuple
                                # list of tuples with energy,elemname
                                for k, elemtuple in enumerate(elemlines):
                                    # elemtuple[0] is energy and [1] is element symbol
                                    # axes[thisrow,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                    try:
                                        axes[0,thiscol].axvline(x=elemtuple[0], color='b') # O line
                                        axes[0,thiscol].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                                    except:
                                        print('Problem labeling elements')
                firstrun=False  # label elems only for first pixel on page          
                print('Pixel', str(X), str(Y), 'plotted.')
                axes[0,0].legend(mylegend,loc='best',fontsize=8)
                # except:
                # print('Unknown problem plotting pixel', str(X),', ', str(Y))
            pdf.savefig(fig)
            plt.close(fig)
    plt.ion()
    return

def calcderiv(rawdata, peakind, lowind, highind, peakrange, lowrange, hirange):
    ''' Piece-wise calculation of derivative over element region including possible breaks/
    different # sweeps, etc. 
    regions shouldn't have overlap but may have discontinuities '''
    # Find index #s of any ev breaks (checking for continuity)
    # 4 cases for possible evbreaks
    # rawdata contains local peak from lowind[0] to highind[1]
    # all indices below are relative to lowind[0] (subtract that)
    if lowrange[1]+1!=peakrange[0]:
        if peakrange[1]+1!=hirange[0]:
            evranges=[[0,lowind[1]-lowind[0]], [peakind[0]-lowind[0],peakind[1]-lowind[0]], [highind[0]-lowind[0],highind[1]-lowind[0]]]
        else:
            evranges=[[0,lowind[1]-lowind[0]], [peakind[0]-lowind[0],highind[1]-lowind[0]]]
    elif peakrange[1]+1!=hirange[0]:
        evranges=[[0,peakind[1]-lowind[0]], [highind[0]-lowind[0],highind[1]-lowind[0]]]
    else:
        evranges=[[0,highind[1]-lowind[0]]] # full range w/o data breaks
    thiss7d7=np.zeros(0)
    for i, [start,stop] in enumerate(evranges):
        thissect=smoothdiffS7D7(rawdata[start:stop+1])
        thiss7d7=np.concatenate((thiss7d7, thissect), axis=0)
    if len(rawdata)!=len(thiss7d7):
        print('Problem... data and deriv of data are different lengths', len(rawdata), len(thiss7d7))
    return thiss7d7

def plotpix_tk(specimage, energy, Elemdata, spectralregs, amplmaps, integmaps, shiftmaps, AESquantparams, **kwargs):
    ''' tk interface for plotting multiplex data from selected pixels 
    within spectral image; 
    kwargs: fileset, areas, xrangestr -- usually just holds values entered during prior run  
    shiftdict - contains peak shift 100x100 arrays for all peaks
    sumdict - contains rawcnts, backcnts and subcnts for 100x100 arrays/all peaks
    TODO: Add optional plot of shift for this pixel from shiftdict
    '''
    # first print out existing info in various lines
    Elements=[i[0] for i in Elemdata]
    root = tk.Tk()
    root.title("Multiplex plot of spectral image pixel")
    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    backfitbool.set(True)
    pixstr=tk.StringVar()
    pixstr.set(kwargs.get('pixstr',"(0,0);(1,1)"))
    labelbool=tk.BooleanVar() # label elements
    labelbool.set(True)
    shiftbool=tk.BooleanVar() # display peak shift (TODO remove this?)
    shiftbool.set(True)
    countsbool=tk.BooleanVar() # 
    countsbool.set(True)
    smdifbool=tk.BooleanVar()
    smdifbool.set(True)
    plottype=tk.StringVar()  # counts, deriv or both
    plotelemstr=tk.StringVar()
    mytext='Peaks to be labeled:'+', '.join(Elements)
    plotelemstr.set(mytext)
    choice=tk.StringVar()  # plot or abort
    
    tk.Label(root, text='Enter list of X, Y pixels for plot').grid(row=0, column=0)
    tk.Entry(root, textvariable=pixstr).grid(row=0, column=1)
    # Selection of plot type
    tk.Label(root, text='Plot type').grid(row=1, column=0)
    tk.Radiobutton(root, text='Counts', value='Counts', 
                   variable=plottype).grid(row=1, column=1)
    tk.Radiobutton(root, text='Deriv', value='Deriv', 
                   variable=plottype).grid(row=1, column=2)
    tk.Radiobutton(root, text='Both', value='Both', 
                   variable=plottype).grid(row=1, column=3)

    tk.Checkbutton(root, variable=backfitbool, 
                   text='Plot fitted backgrounds?').grid(row=2, column=0)
    tk.Checkbutton(root, variable=labelbool, 
                   text='Label peak centers?').grid(row=2, column=1)
    tk.Checkbutton(root, variable=shiftbool, 
                   text='Display peak shift?').grid(row=3, column=0)
    tk.Checkbutton(root, variable=countsbool, 
                   text='Display integrated counts info?').grid(row=3, column=1)
    tk.Checkbutton(root, variable=smdifbool, 
                   text='Plot derivative peak?').grid(row=4, column=0)
    # option to reselect labeled elemental peaks 
        
    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=5, column=0)

    d=tk.Button(root, text='Plot')
    d.bind('<Button-1>', plot)
    d.grid(row=5, column=1)

    root.mainloop()
        
    mychoice=choice.get()
    
    if mychoice=='plot':
        # Set up kwargs for plot 
        kwargs={}
        plottype=plottype.get()
        if pixstr.get()!='': # extract list of pixels
            pixlist=pixstr.get().split(';')
            # Convert pixels string to X, Y tuple list
            try:
                Xvals=[int(s.split("(")[1].split(',')[0]) for s in pixlist]
                Yvals=[int(s.split(",")[1].split(')')[0]) for s in pixlist]
                pixlist=[list(i) for i in zip(Xvals,Yvals)] # zip together into single list of x,y pixels
                kwargs.update({'pixstr' : pixstr.get()}) # pass-back string
            except:
                print('Failed conversion of pixel values from', pixstr.get())
            #return kwargs
        if backfitbool.get(): # include linear background fit
            kwargs.update({'addbackfit':True})
        if smdifbool.get(): # label peak centers 
            kwargs.update({'smdifpts':True})
        if shiftbool.get(): # include linear background fit
            kwargs.update({'shiftlabel':True})
        if countsbool.get(): # include text of counts info
            kwargs.update({'countslabel':True})            
        if labelbool.get(): # label peak centers 
            kwargs.update({'plotelems':Elements})

        plotpixels(specimage, plottype, energy, pixlist, Elemdata, spectralregs, 
                   amplmaps, integmaps, shiftmaps, AESquantparams, **kwargs)
        return kwargs
    return
''' TESTING
i=0   (X,Y)=pixlist[i]
[elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, 
                idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth, 
                kfact, kfact2, mass]=Elemdata[1]
'''

def plotpixels(specimage, plottype, energy, pixlist, Elemdata, spectralregs, amplmaps, 
               integmaps, shiftmaps, AESquantparams, **kwargs): # single file already selected by number
    '''Extract and plot multiplex spectrum from list of pixel coords (direct counts 
    not s7d7 deriv)
    specimage - entire 3d ndarray spectral data
    energy - energy x values in eV
    pixels - list of x,y pixel values to extract
    backarray -- 2D array with background slope ([0]) and intercept[1]
    AESquantparams- has ideal peak positions
    
    kwargs:  
        addbackfit-- bool to select backfit plotting
        plotelems - optional plotting of ideal energy of peak centers
        shiftlabel (bool) - plots shift value for pixel from smoothed spectrum
        countslabel - adds text summary of integrated counts, backcnts, subcnts
        TODO charge -- offset energy for mapping due to charging (built-in to multiplex)
        '''
    if plottype=='Both':
        numrows=2 # counts on top and S7D7 on bottom
        numcols=len(Elemdata) # number of elemental subregions for plotting
    else: # Count or derivative plot only
        numcols=min(len(Elemdata),2) # 1 or 2 columns
        numrows=math.ceil(len(Elemdata)/2)   

    fig, axes = plt.subplots(nrows=numrows, ncols=numcols, squeeze=False) # axes is array

    # Use color and linestyles to differentiate filenumber(s) and areanumber(s)
    # mylegend=[]
    colorlist=['b','r','g','c','m','y','k', 'olive','pink','purple']
    filecount=0 # used for colorlist

    for i, (X,Y) in enumerate(pixlist):
        # load appropriated data from this pixel from specimage
        Augerfile=pd.DataFrame()
        Augerfile['Energy']=energy
        rawdata=specimage[X,Y,:]
        Augerfile['Counts']=rawdata
        # if counts and deriv plot all in 2row/1 col window
        for j, [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth, kfact, kfact2, mass] in enumerate(Elemdata):
            Augerslice=Augerfile[(Augerfile['Energy']>=lowrange[0]) & 
                                 (Augerfile['Energy']<=hirange[1])]
            plkwargs={'color':colorlist[filecount]} # custom colorlist
            # loop over subplots (each is a multiplex elemental energy region from AESquantparam)                
            if 'addbackfit' in kwargs: # do for counts and both selections                 
                [integcnts,integev, countmax,slope, intercept] = integmaps[j][X, Y]
                # Create range for slope/intercept plot
                xvals=np.linspace(int(lowrange[0]),int(hirange[1]),100)
            if plottype=='Counts':
                axindex=j%numrows, j//numrows # thisrow, this column
                labelname='Pixel X='+ str(X) +' Y=' +str(Y)
                Augerslice.plot(x='Energy', y='Counts', linestyle='', marker='o', 
                        label=labelname, ax=axes[axindex], **plkwargs) # plot counts
                if kwargs.get('addbackfit', False):
                    axes[axindex].plot(xvals, xvals*slope+intercept, **plkwargs)
            elif plottype=='Deriv':
                axindex=j%numrows, j//numrows # thisrow, this column tuple
                # Need to calculate this derivative in piece-wise manner
                thiscnts=rawdata[lowind[0]:highind[1]+1]
                thiss7d7=calcderiv(thiscnts, peakind, lowind, highind, peakrange, 
                                   lowrange, hirange)
                xvals=range(int(lowrange[0]),int(hirange[1])+1)
                axes[axindex].plot(xvals, thiss7d7)
            elif plottype=='Both': # counts and s7d7/deriv
                thiscol=j # row 1 is counts, row 2 is deriv
                axindex=0, j//numrows
                thiscnts=rawdata[lowind[0]:highind[1]+1] # parse counts to only this peak
                thiss7d7=calcderiv(thiscnts, peakind, lowind, highind, peakrange, 
                                   lowrange, hirange)
                xvals=range(int(lowrange[0]),int(hirange[1])+1)
                axes[1,thiscol].plot(xvals, thiss7d7) # deriv in lower (1)
                Augerslice.plot(x='Energy', y='Counts', ax=axes[0,thiscol]) # counts in lower 
                if kwargs.get('addbackfit', False):
                    axes[0, thiscol].plot(xvals, xvals*slope+intercept, **plkwargs)

            # Optional additional plotting of backgrounds (only for direct counts plotting)
            if 'plotelems' in kwargs:
                Elements=kwargs.get('plotelems',[])
                elemlines = findelemlines(Elements, [lowrange[0], hirange[1]], AESquantparams, peaktype='Counts')
                for i, val in enumerate(elemlines):
                    axes[axindex].axvline(x=val, color='k')
            # Adding optional text labels
            if 'shiftlabel' in kwargs or 'countslabel' in kwargs:
                mytitle=''
                if kwargs.get('shiftlabel', False):
                    # print/label integ and deriv shifts for this pixel
                    [shiftder, shiftint]=shiftmaps[j][X, Y]
                    mytitle+='Shift (deriv, integ)'+ str(int(shiftder))+', '+ str(int(shiftint))
                if kwargs.get('countslabel', False):
                    [integcnts,integev, countmax,slope, intercept] = integmaps[j][X, Y]
                    mytitle+=str(X)+','+str(Y)+'cnts raw:'+ str(int(integcnts))+'\n'
                axes[axindex].set_title(mytitle, fontsize=12)
            # axes[axindex].legend(mylegend, loc='best')
        filecount+=1 # used for custom color list (multiple pixels in same window)

    '''
    
    for axis in ['top', 'bottom','left','right']:
        axes[0,0].spines[axis].set_linewidth(2.5) # change axes thickness
        try: # sometimes plotting both counts and deriv
            axes[1,0].spines[axis].set_linewidth(2.5) 
        except:
            pass
    '''
    return

def getelemdata(SpectralRegs, AESquantparams, **kwargs):
    ''' Return index ranges for each element's peak, low background, and high background
    Elemdata contains: 
        0) elem/peak name 
        1) order in original scan
        2) indices of peak range (list) 
                    
        3) indices of low back range (list) 4) indices of high back range (list)
        5) energy vals of peak range (list) 6) energy vals of low back range (list) 
        7) energy vals of high back range (list)
        8) energy of ideal direct peak (float; set to nan if not found)
        9) corresponding index of ideal peak (nan if not within data range).. this is possible 
        if significant peak shift is occurring and predefined shift can be built in to 
        QM multiplex
        10) ideal negpeak ev (position of negative peak in smooth-diff deriv s7d7)
        11) typical integration peak width
        12) chargeshift -- # of eV of charge-compensation shift applied to scan  
		     figure out by comparison of AESquantparams w/ scan center
        13) peakwidth - element specific negpeak-pospeak
        14) searchwidth - region 
    kwarg: Elements - explicit element list (normally autogenerated)...
        only needed if multiplex setup is unusual
        '''
    Elemdata=[] # list of lists with important element data
    if 'elems' in kwargs:
        Elements=kwargs.get('elems','')
    else:
        # Autogenerate elements within quantmap
        Elements=np.ndarray.tolist(SpectralRegs.Element.unique())
        # works only with OL, O, OH type setup (not continuous single O region)
        Elements=Elements[1::3]
        validelem=np.ndarray.tolist(AESquantparams['element'].unique())
        missing=[elem for elem in Elements if elem not in validelem]
        if len(missing)!=0:
            print(','.join(missing),' are not valid elements/peaks in AESquantparams')      
    # Ensure all are standard peak names
    for i, elem in enumerate(Elements):
        thispeak=SpectralRegs[SpectralRegs['Element']==elem]
        lowback=SpectralRegs[SpectralRegs['Element']==elem+'L']
        hiback=SpectralRegs[SpectralRegs['Element']==elem+'H']
        match=AESquantparams[AESquantparams['element']==elem]
        if len(match)!=1 | len(thispeak)!=1 :
            print('Element', elem,' not found in AESquantparams or associated peak missing')
            continue
        idealev=int(match.iloc[0]['negpeak']+match.iloc[0]['integpeak'])
        idealnegpeak=int(match.iloc[0]['negpeak'])
        width=int(match.iloc[0]['integwidth'])
        # typical distance between neg and pospeaks in smooth-diff data
        peakwidth=int(match.iloc[0]['peakwidth'])
        searchwidth=int(match.iloc[0]['searchwidth'])
        kfact=match.iloc[0]['kfactor']
        kfact2=match.iloc[0]['kfactor2']
        mass=match.iloc[0]['mass']
        
        # Find associated index value in this dataset
        peakrange=[j for j in range(int(thispeak.iloc[0]['Lower']),int(thispeak.iloc[0]['Upper']))]
        minrange=int(thispeak.iloc[0]['Start'])
        maxrange=int(thispeak.iloc[0]['End'])
        minevrange=float(thispeak.iloc[0]['Lower'])
        maxevrange=float(thispeak.iloc[0]['Upper'])
        # Some QM scans have built in charging compensation
        chargeshift=int((minevrange + maxevrange)/2 - idealev)
        try:
            val=peakrange.index(idealev)
            indrange=[j for j in range(int(thispeak.iloc[0]['Start']),int(thispeak.iloc[0]['End']))]
            idealind=indrange[val]
        except:
            idealind=np.nan  # ideal peak not in range (due to chargeshift)
        if len(lowback)==1 & len(hiback)==1: 
            minlow=int(lowback.iloc[0]['Start'])
            maxlow=int(lowback.iloc[0]['End'])     
            minhigh=int(hiback.iloc[0]['Start'])
            maxhigh=int(hiback.iloc[0]['End'])
            minevlow=float(lowback.iloc[0]['Lower'])
            maxevlow=float(lowback.iloc[0]['Upper'])     
            minevhigh=float(hiback.iloc[0]['Lower'])
            maxevhigh=float(hiback.iloc[0]['Upper'])
        else: # handle continuous scan case 
            # often used in cases of significant charging so below vals aren't so valid anyway
            minlow=np.nan
            maxlow=np.nan
            minhigh=np.nan
            maxhigh=np.nan
            minevlow=np.nan
            maxevlow=np.nan    
            minevhigh=np.nan
            maxevhigh=np.nan

        Elemdata.append([elem, i, [minrange, maxrange], [minlow, maxlow], [minhigh, maxhigh], 
            [minevrange, maxevrange], [minevlow, maxevlow], [minevhigh, maxevhigh], idealev, idealind, 
            idealnegpeak, width, chargeshift, peakwidth, searchwidth, kfact, kfact2, mass])
    return Elemdata

def uniformfilter(specimage, size=1):
    ''' Uniform filter does adjacent averaging in spatial domain; returns entire
    spectral image but with adjacent pixel spatial averaging 
    other scipy options include: gaussian_filter; convolve2d '''
    specsmooth=np.empty([specimage.shape[0],specimage.shape[1],specimage.shape[2]])
    for z in range(0, specimage.shape[2]):
        slice2d=specimage[:,:, z]
        specsmooth[:,:,z]=filters.uniform_filter(slice2d, size)
    return specsmooth

def smoothdiffS7D7(cnts):
    ''' Create smooth differentiated column from counts using S7D7 PHI algorithm (Multipak tables 
    A-5 and A-1; passing np array (not pandas series)''' 
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

    # Same structure to perform differentiation on smoothed datalist above
    for i in range(0,numpts): # special cases for endpoints (within 3 of an evbreak)
        diff=i-min(endpts, key=lambda x:abs(x-i)) # distance from closest evbreak index # in list
        if abs(diff)<=2:
            smoothdiff[i]=0  # just zero out endpoints (old code moved to python software dev xls cell)
        else:
            smoothdiff[i]=(-3*smooth[i-3]-2*smooth[i-2]-1*smooth[i-1]+1*smooth[i+1]+2*smooth[i+2]+3*smooth[i+3]+1)/28
    smoothdiff=np.asarray(smoothdiff)
    return smoothdiff


''' TESTING
i=0
[elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, idealnegpeak, integwidth, 
     chargeshift, peakwidth, searchwidth, kfact, kfact2, mass]=Elemdata[i]
bigpeak='O'
X=0   Y=0
'''

def findallpeaks(specimage, Elemdata):
    ''' Find charging shift associated w/ some large peak (e.g. O) for all pixels 
    return charge shift and peak amplitude (measure of significance) 
    works if scan regions are sufficiently large to allow smooth-diff peak 
    calcshifts performs similar function on direct peaks (but needs mod for big shifts)
    
    peakind- index # of this peak in np Z direction; peakrange - associated energy range
    lowind, highind, lowrange, hirange - above and below background regions
    Params saved:
        amplmap deriv-related (mostly for spectral plotting)
            [0] ampl and [1] associated negpeak value [2] energy val (not index)
            and smdiff peak width (#ev to popspeak)
        integmap:
            [0] best integcounts value and [1] assoc energy val
            [2] slope and [3] intercept of background fit
        shift map 
        counts max at peak position 
    '''
    numpeaks=3 # max number of negpeaks in deriv spectrum to evaluate

    # Two shift values (deriv and direct peaks) for each element
    shiftmaps=[]
    ''' Parameters for s7d7 derivative peaks:  [0] ampl [1] negval [2] negpeakind
    [3] width (can get pospeak index and val indirectly for use in spectral plots)
    '''
    amplmaps=[]
    ''' Parameters for direct integrations: [0] integcnts [1] peak index (eV value
    available indirectly through shiftmaps) [2] countsmax (unsubtracted int at peak)
    [3] slope [4] intercept
    '''
    integmaps=[]
    for i, [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, 
            idealev, idealind, idealnegpeak, integwidth, chargeshift, peakwidth, 
            searchwidth, kfact, kfact2, mass] in enumerate(Elemdata):
        print('Extracting maps for element', elem)
        if str(lowind[0])=='nan' and str(highind[0])=='nan':
            lowind, lowrange, highind, hirange, peakind, peakrange= fixoddsetup(peakind, peakrange)
        shiftmap=np.empty([specimage.shape[0],specimage.shape[1], 2]) 
        amplmap=np.empty([specimage.shape[0],specimage.shape[1], 4])
        integmap=np.empty([specimage.shape[0],specimage.shape[1], 5])
        for X in range(0,specimage.shape[0]):
            for Y in range(0,specimage.shape[1]):
                # Get raw data associated w/ bigpeak (whole region not just peak subset)
                rawdata=specimage[X,Y,lowind[0]:highind[1]+1]
                if rawdata.max()<0.0001: # data missing.. only zeros for counts
                    shiftmap[X,Y]=np.nan
                    amplmap[X,Y]=np.nan
                    continue
                s7d7=calcderiv(rawdata, peakind, lowind, highind, peakrange, lowrange, hirange)
                # Find all relative minima
                foundind=argrelmin(s7d7) # scipy.stat returned as array (index is reset!)
                mypeaks=pd.DataFrame()
                thisenergy=np.arange(lowrange[0],hirange[1]+1,1.0)
                mypeaks['Energy']=thisenergy[foundind]
                mypeaks['negpeakval']=s7d7[foundind]
                mypeaks['negpeakind']=foundind[0] # relative index (original position lost)
                # Find associated pospeaks
                # row=mypeaks.loc[0]
                for index, row in mypeaks.iterrows():
                    # searching within s7d7 (index zeroed)
                    lowlim=int(row.negpeakind-peakwidth-searchwidth) # again indices relative to lowind[0]
                    uplim=int(row.negpeakind - peakwidth + searchwidth + 1)
                    if lowlim<0:
                        lowlim=0
                    if uplim>len(s7d7)-1:
                        uplim=len(s7d7)
                    try: # This can fail on wide C peaks
                        pospeakval=s7d7[lowlim:uplim].max()
                       # index in terms of original s7d7 dataset
                        pospeakind=np.unravel_index(s7d7[lowlim:uplim].argmax(), s7d7.shape)[0]+lowlim
                        # now calculate amplitude and s7d7 peak width
                        mypeaks=mypeaks.set_value(index,'Ampl',pospeakval-row.negpeakval)
                        # peakwidth = negpeak- pospeak
                        mypeaks=mypeaks.set_value(index,'Peakwidth',row.negpeakind-pospeakind)
                    except:
                        print('No s7d7 amplitude for', elem,'peak', str(index), 'of pixel', str(X),str(Y))
                 # Can consider amplitude, peakwidth and position when looking for real peaks
                mypeaks=mypeaks.sort_values(['Ampl'], ascending=False).head(numpeaks)
                mypeaks=mypeaks.reset_index(drop=True) # necessary?
                # DIRECT peak check..  for top n peaks, also check for direct peak maxima (no background yet)
                peakdata=pd.DataFrame()
                peakdata['Energy']=thisenergy
                peakdata['Counts']=rawdata
                # adjust index to match vals in lowind, highind
                peakdata.index=peakdata.index+lowind[0]
                # Do background fit of peak region
                try:
                    slope, intercept = calcbackground(peakdata, lowind, highind)
                except:
                    print('Slope fit problem for pixel', str(X), str(Y))
                    # print('Lowind/ highind are', print(lowind), print(highind))
                # Add subtracted data
                peakdata['Subdata']=np.nan
                for index, row in peakdata.iterrows():
                    yval=row.Counts-(slope*row.Energy+intercept)
                    peakdata=peakdata.set_value(index,'Subdata', yval)
                # Do direct integration over all suspected peaks
                mypeaks=integpeaks(peakdata, mypeaks, slope, intercept, 
                                   idealnegpeak, idealev, integwidth)
                # Get best peak data from available selections 
                try:
                    mypeak=pickbestderivpeak(mypeaks, idealev, idealnegpeak)
                except:
                    print('Problem with X, Y', str(X), str(Y),'for ', elem)
                # Store needed params in various maps 
                shiftmap[X,Y,0]=mypeak.Energy-idealnegpeak # from deriv
                shiftmap[X,Y,1]=mypeak.Energy2-idealev # from direct
                ''' Parameters for s7d7 derivative peaks:  [0] ampl [1] negval [2] negpeakind
                [3] width (can get pospeak index and val indirectly for use in spectral plots)
                '''
                amplmap[X,Y,0]= mypeak.Ampl
                amplmap[X,Y,1]= mypeak.negpeakval
                amplmap[X,Y,2]= int(mypeak.Energy) # use eV not peak index 
                amplmap[X,Y,3]= int(mypeak.Peakwidth)
                ''' Parameters for direct integrations: [0] integcnts [1] peak index (eV value
                available indirectly through shiftmaps) [2] countsmax (unsubtracted int at peak)
                [3] slope [4] intercept
                '''
                integmap[X,Y,0]= int(mypeak.Integcounts)
                integmap[X,Y,1]= int(mypeak.Energy2) # eV value of direct peak
                integmap[X,Y,2]= int(mypeak.Countsmax)
                integmap[X,Y,3]= slope
                integmap[X,Y,4]= intercept
        amplmaps.append(amplmap)
        shiftmaps.append(shiftmap)
        integmaps.append(integmap)
    return amplmaps, shiftmaps, integmaps

def findnegpeaks(specimage, Elemdata, bigpeak):
    ''' Find charging shift associated w/ some large peak (e.g. O) for all pixels 
    return charge shift and peak amplitude (measure of significance) 
    works if scan regions are sufficiently large to allow smooth-diff peak 
    calcshifts performs similar function on direct peaks (but needs mod for big shifts)
    
    peakind- index # of this peak in np Z direction; peakrange - associated energy range
    lowind, highind, lowrange, hirange - above and below background regions
    
    '''
    #TODO fix handling of empty data (zeros)
    thiselem=[e for e in Elemdata if e[0] == bigpeak]
    [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, 
     idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth]=thiselem[0]
    # Need array with total shift (chargeshift adjusted by actual position relative to center)
    
    if str(lowind[0])=='nan' and str(highind[0])=='nan':
        # modify lowind/highind in case of unconventional setup (remove section from peakind/peakrange)
        numchan=int(0.1*(peakind[1]-peakind[0]))
        lowind[0]=peakind[0]
        lowind[1]=peakind[0]+numchan
        lowrange[0]=peakrange[0]
        # fix me !!!
        lowrange[1]= 600+ 10/100 * 600
        lowrange[1]= peakrange[0]+ (peakind[1]-peakind[0])*(lowind[1]/(peakind[1]-peakind[0]))
        highind[0]=peakind[1]-numchan
        highind[1]=peakind[1]
        hirange[0]= peakrange[0] + (peakind[1]-peakind[0])*(highind[0]/(peakind[1]-peakind[0]))       
        hirange[1]=peakrange[1]
        peakind =[lowind[1]+1,highind[0]-1] # reset peak range to compensate
        peakrange=[lowrange[1]+1, hirange[0]-1 ]
    # Also include peak amplitude
    chargeshiftmap=np.empty([specimage.shape[0],specimage.shape[1]])
    peakamplmap=np.empty([specimage.shape[0],specimage.shape[1]])
    badcount=0
    for X in range(0,specimage.shape[0]):
        for Y in range(0,specimage.shape[1]):
            # Get raw data associated w/ bigpeak (whole region not just peak subset)
            rawdata=specimage[X,Y,lowind[0]:highind[1]+1]
            if rawdata.max()<0.0001: # data missing.. only zeros for counts
                chargeshiftmap[X,Y]=np.nan
                peakamplmap[X,Y]=np.nan
                continue
            s7d7=calcderiv(rawdata, peakind, lowind, highind, peakrange, lowrange, hirange)
            negpeakamp=s7d7.min()
            # Get index of min (shift by lowind[0] to match correct position in unsliced data)
            negpeak=np.unravel_index(s7d7.argmin(), s7d7.shape)[0] + lowind[0] 
            # find energy value and then shift from ideal position 
            if negpeak in range(lowind[0],lowind[1]+1):
                foundpeakev=lowrange[0]+(negpeak-lowind[0])/(lowind[1]-lowind[0]) * (lowrange[1]-lowrange[0])
                thisshift=int(foundpeakev-idealnegpeak)
            elif negpeak in range(peakind[0],peakind[1]+1):
                foundpeakev=peakrange[0]+(negpeak-peakind[0])/(peakind[1]-peakind[0]) * (peakrange[1]-peakrange[0])
                thisshift=int(foundpeakev-idealnegpeak)
            elif negpeak in range(highind[0],highind[1]+1):
                foundpeakev=hirange[0]+(negpeak-highind[0])/(highind[1]-highind[0]) * (hirange[1]-hirange[0])
                thisshift=int(foundpeakev-idealnegpeak)
            else:
                print('negpeak location not found. Ampl:', str(negpeakamp),' neg index:', str(negpeak),'X,Y',str(X), str(Y))
                thisshift=np.nan
            chargeshiftmap[X,Y]=thisshift
            # Pospeak is typically $peakwidth lower than found negpeak above
            lowlim=negpeak-peakwidth-searchwidth-lowind[0] # again indices relative to lowind[0]
            uplim=negpeak-peakwidth+searchwidth+1-lowind[0]
            # Ensure that limits are within data range
            if lowlim<0:
                lowlim=0
            if uplim>len(rawdata)-1:
                uplim=len(rawdata)-1
            posreg=s7d7[lowlim:uplim]
            try:                
                pospeakamp=posreg.max()
                ''' width of peak in smooth-diff data (maybe diagnostic)
                pospeak=np.unravel_index(posreg.argmax(), posreg.shape)[0]+negpeak-peakwidth-searchwidth
                smdifwidth=negpeak-pospeak
                '''
                ampl=pospeakamp-negpeakamp                
                peakamplmap[X,Y]=ampl
            except:
                print('Problem with pixel',str(X), str(Y))
                print('Size of pospeak search region is',len(posreg))
                print('Limits are ', str(lowlim), str(uplim))
                badcount+=1
    return chargeshiftmap, peakamplmap


def get_spectral_regs(QMpixarray):
    ''' Retrieve and save multiplex spectral regions from first spe file 
    saves multiplex details from first spe file to disk; returns spectralregs''' 
    AugerFileName=QMpixarray.iloc[0]['Filename']
    filenumber=AugerFileName.split('.')[1].split('.')[0]
    # check for file in cwd or sub
    if not os.path.exists(AugerFileName):
        AugerFileName='sub\\'+AugerFileName # check in sub directory
        if not os.path.exists(AugerFileName):
            print('Datafile ', AugerFileName,' not found in cwd or sub directories.')
            return pd.DataFrame(), []
    with open(AugerFileName, 'rb') as file:
        filedata = file.read()
    end=filedata.find(b'EOFH')
    headerdata=filedata[0:end+6] # works to cut off binary part    
    header=headerdata.decode(encoding='cp437') # more generic encoding than utf-8
    # get number of cycles and time per step from header
    tempstring=header.split('NumCycles: ')[1] # find # cycles (works for either survey or multiplex)
    match=re.search(r'\d+',tempstring)
    if match:
        numcycles=int(match.group(0))
    tempstring=header.split('TimePerStep:')[1] # find time per step
    match=re.search(r'\d+\.\d+',tempstring) 
    if match:
        timestep=float(match.group(0))
    details, evbreaks, timeperarea, energy, SpectralRegs=SpectralRegions(filenumber, AugerFileName, numcycles, timestep, header)
    # Spectral image data cube signal always sorted by energy (in get_multiplex_data)
    # need to auto-sort here as well
    SpectralRegs=SpectralRegs.sort_values(['Lower'])
    SpectralRegs=SpectralRegs.reset_index(drop=True)
    # make indexrange showing Ev range of elements as index #s with specimage
    # Need to eliminate indices associated with duplicated energy vals
    SpectralRegs['Start']=np.nan
    SpectralRegs['End']=np.nan
    # Need to remove overlapping vals from spectralregs (keep one w/ large # of sweeps)
    for i in range(0, len(SpectralRegs)-1):
        if SpectralRegs.iloc[i+1]['Lower']<=SpectralRegs.iloc[i]['Upper']:
            redval=SpectralRegs.iloc[i]['Upper']-SpectralRegs.iloc[i+1]['Lower']+1
            print('Removed ', str(redval),' overlaping energy values.')
            # remove overlapping vals from inferior one
            if SpectralRegs.iloc[i+1]['Sweeps']>SpectralRegs.iloc[i]['Sweeps']:
                SpectralRegs=SpectralRegs.set_value(i+1,'Lower',SpectralRegs.iloc[i+1]['Lower']+redval)
            else:
                SpectralRegs=SpectralRegs.set_value(i,'Upper',SpectralRegs.iloc[i]['Upper']-redval)
    count=0
    for index, row in SpectralRegs.iterrows():
        thisrange=SpectralRegs.loc[index]['Upper']-SpectralRegs.loc[index]['Lower']+1
        # Zero based indexing
        SpectralRegs=SpectralRegs.set_value(index, 'Start', count) 
        SpectralRegs=SpectralRegs.set_value(index, 'End', count+thisrange-1)
        count+=thisrange
    # Rename certain peaks/elem to align peak naming conventions (generally main peak w/o appended number)
    peakdict={'Mg2':'Mg','Si2':'Si','S1':'S','Fe3':'Fe'}
    for key, val in peakdict.items():
        SpectralRegs['Element']=SpectralRegs['Element'].str.replace(key, val)
    energy.sort() # extracted multiplex spectra always sorted by increasing eV
    return SpectralRegs, energy

''' TESTING
file=spefilelist[0]  i=0
'''

def makespecimage(QMpixarray, spectralregs):
    ''' Make 3D numpy array from stack of multiplex spe files (20 areas per file); 
    pixel position obtained from QMpixarray file
    in final 3D nparray, X, Y are spatial and Z is signal ... electron counts/sec 
    from single multiplex scan '''
    # Get energy range for multiplex from first data file (same for all pixels)
    AugerFileName=QMpixarray.iloc[0]['Filename']
    # get energy values and multiplex breaks directly from spe file
    energy, evbreaks =get_energy_range(AugerFileName) # Values in original order
    xmax=QMpixarray.Xindex.max()
    ymax=QMpixarray.Yindex.max()
    # Make 3D blank numpy array of dimension xmax+1 * ymax+1 * len(energy)- duplicates
    specimage=np.empty([xmax+1,ymax+1,len(set(energy))])

    # Now extract counts from all areas in each multiplex (20 max per spe)
    spefilelist=np.ndarray.tolist(QMpixarray.Filename.unique())
    for i, file in enumerate(spefilelist):
        thisfile=QMpixarray[QMpixarray['Filename']==file]
        numareas=len(thisfile)
        # Dataframe with all spatial areas
        thismult=get_multiplex_data(file, energy, evbreaks, numareas)
        # keep best value (eliminate accidental ev value duplicates in specregs)
        thismult=keepbestvals(thismult, spectralregs)
        if thismult.empty: 
            print('File missing ... no data for ', file)
            continue
        for index, row in thisfile.iterrows():
            xind=thisfile.loc[index]['Xindex']
            yind=thisfile.loc[index]['Yindex']
            thiscolnum=thisfile.loc[index]['Subnumber']+1
            if len(thismult)!=specimage.shape[2]:
                print('spe spectrum has different length than specimage array!')
                continue
            else:
                specimage[xind, yind,:]=thismult['Counts'+str(thiscolnum)]
    energy=list(set(energy)) # now remove eV duplicates
    energy.sort()
    return specimage, energy

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
    df=df[~df.index.isin(removelist)]
    print (len(removelist), ' duplicated energy values removed from multiplex')
    return df

def get_multiplex_data(AugerFileName, energy, evbreaks, numareas):
    ''' Function to extract and return all areas withing multiplex spe file 
    data is sorted by energy
    list of energy values are passed from SpectralRegions function (energy x values are not continuous for multiplex)
    bindata is rest of binary after text header extraction
	 evbreaks are index# of boundaries between spectral regions 
    energy is list of ev values and corresponding index within specimage (element in list) 
    '''
    if not os.path.exists(AugerFileName):
        if os.path.exists('sub//'+AugerFileName):
            AugerFileName='sub//'+AugerFileName # open from sub directory
    try:
        with open(AugerFileName, 'rb') as file:
            filedata = file.read()
    except:
        print(AugerFileName," missing from data directory")
        return pd.DataFrame() # return empty frame?
    end=filedata.find(b'EOFH')
    bindata=filedata[end+6:] # binary data section of file (header removed)
    mystruct=struct.Struct('f')
    # binary header of variable length.. just count backward using known data length
    startbyte=len(bindata)-4*numareas*len(energy)
    
    # TEST for correct byte location... data regions begins right after last occurence of triple zero 4 byte float values.. for survey usually found in bytes 100:111 
    for i in range(startbyte-12,startbyte,4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval) # little endian encoding
        if unpackbyte[0]!=0:
            print('Possible binary read error... leading zeros are missing for ',AugerFileName)

    # create data frame for multiplex and assign energy values   
    multiplex=pd.DataFrame() # data frame for energy, and counts col for each area   
    multiplex['Energy']=energy # energy values found in SpectralRegions and passed as argument
    
    # Read in and convert all numeric values
    alldata=[] # single list for all Auger counts values in file
    for i in range(startbyte,len(bindata),4):
        byteval=bindata[i:i+4]
        unpackbyte=mystruct.unpack(byteval)
        alldata.append(unpackbyte[0])
    
    if len(alldata)!=len(energy)*numareas: # number of data values expected for each area
        print('Possible binary data reading error: Data string length does not match expected value for ',AugerFileName)
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
    
    for i,vals in enumerate(specregs):
        # vals is list of count values for this spectral region (for each spatial area back to back)
        numvals=int(datachunks[i]/numareas) # number of values in single area/single spectral region    
        for j in range(0, numareas):
            counts[j].extend(vals[j*numvals:(j+1)*numvals]) # gets all counts columns
    for i in range(1,numareas+1):
        cntsname='Counts'+str(i) # name for nth column is counts2, counts3, etc.
        multiplex[cntsname]=counts[i-1] # assign counts to frame (and switch from 0 based indexing)
    # Solve multiplex out of order problem
    if not multiplex.Energy.is_monotonic: # energy values out of order problem
        multiplex=multiplex.sort_values(['Energy']) # sort before returning
    return multiplex

def loadQMpixarray():
    ''' Load of standard files from main Auger data directory '''
    pixfiles=glob.glob('*pixarray*')
    if len(pixfiles)==1:
        QMpixarray=pd.read_csv(pixfiles[0])
        if 'Filename' not in QMpixarray:
            print('spe data file names not yet linked with quantmap pixel array.')
    else:
        print("Couldn't locate single pixarray definition file")
    return QMpixarray
''' TESTING 
basename='3Oct17'   startnum=101
'''
def loadspecimage(dir):
    ''' Attempt load of single *specimage* file from directory '''
    npyfiles=glob.glob('*specimage*.npy')
    if len(npyfiles)==1:
        specimage=np.load(npyfiles[0])
    else:
        print("Couldn't locate single spectral image numpy 3D array in current directory")
        return np.array((0,0)) # return empty np array
    return specimage

def linkfilename(QMpixarray, basename, startnum=101):
    ''' Link multiplex spe data file names to associated pixel in quant map,
    startnum is filenumber of first spe; also checks data directory for these filenumbers '''
    if 'Filename' not in QMpixarray:
        QMpixarray['Filename']=''
    else:
        QMpixarray['Filename']=QMpixarray['Filename'].astype(str)
    for index,row in QMpixarray.iterrows():
        QMpixarray=QMpixarray.set_value(index, 'Filename',basename+'.'+str(startnum+index//20)+'.spe')
    # Check for existence of data against list of generated names
    spefiles=np.ndarray.tolist(QMpixarray.Filename.unique())
    datafiles=glob.glob('*.spe')+glob.glob('sub\\*.spe')
    datafiles=[s.replace('sub\\','') for s in datafiles]
    # only considered missing if not in cwd or sub
    missing=[f for f in spefiles if f not in datafiles]
    if len(missing)!=0:
        try:
            fnums=[int(i.split('.')[1].split('.')[0]) for i in missing]
            franges=[] # ranges of files for missing file output
            for key, group in groupby(enumerate(fnums), lambda x: x[0]-x[1]):
                thisgroup=list(map(itemgetter(1), group))
                if len(thisgroup)>1:
                    # more than one consecutive so group as min-max in frange
                    franges.append(str(min(thisgroup))+'-'+ str(max(thisgroup)))
                else:
                    franges.append(str(thisgroup[0])) # single non-consecutive filenumber
            print('Multiplex spe files in Qmpixarray missing from data directory.')
            print('Filenumbers ',', '.join(franges),' are missing.')
        except:
            print('Filenumbers', ', '.join(missing),' are missing.')
        return QMpixarray

def get_energy_range(AugerFileName):
    ''' Get xrange for QM multiplex (should be basically identical for all )
    if multiplex energy range out of order, then rearranged later (evbreaks needed) 
    pulled from processAuger
    energy is list of eV values (and index of list is same index in specimage
    elements
    '''
    if not os.path.exists(AugerFileName):
        if os.path.exists('sub//'+AugerFileName):
            AugerFileName='sub//'+AugerFileName # open from sub directory  
    with open(AugerFileName, 'rb') as file:
        filedata = file.read()
    end=filedata.find(b'EOFH')
    headerdata=filedata[0:end+6] # works to cut off binary part    
    header=headerdata.decode(encoding='cp437') # more generic encoding than utf-8
    # Stripped down version of SpectralRegions function from Auger import
    tempstring=header.split('NoSpectralReg: ')[1] # unlike NoSpectralRegFull inactives are already removed
    match=re.search(r'\d+',tempstring)
    numdefregions=int(match.group(0)) # number of defined regions (can be larger than # active regions)
    energy=[] # list for energy x values
    evbreaks=[0] # region boundaries needed for smoothdiff include first
    for i in range(numdefregions):
        tempstr=tempstring.split('SpectralRegDef: ')[i+1] # should work to split 
        numpts=int(tempstr.split(' ')[4])   
        evstep=float(tempstr.split(' ')[5]) #eV/step        
        startev=float(tempstr.split(' ')[6]) # starting eV
        for j in range(0,numpts): # generate energy values for survey
            energy.append(startev+evstep*j)
            # Even if accidentally duplicated, cannot remove any energy vals here
        evbreaks.append(len(energy)-1) # gives indices of discontinuities in energy scan (needed for smoothdiffS7D7)   
    return energy, evbreaks
  
def showQMregion(jpgfname, margin, AugerParamLog):
    ''' Indicate mapped region for quant map on jpg overlay
    get mag from paramlog and then scale for 1cm = 1 micron '''
    try:
        jpgimage=Image.open(jpgfname)
        draw=ImageDraw.Draw(jpgimage) # single draw instance to label above image
        draw.rectangle((512*margin/2,512*margin/2,512*(1-margin/2),512*(1-margin/2)), outline='red') # red rectangle at specified position
        annotjpgname=jpgfname.replace('.jpg','_map.jpg')
        # Find field of view from AugerParamLog
        match=AugerParamLog[AugerParamLog['Filename']==jpgfname.replace('.jpg','.sem')]
        if len(match)==1:
            fieldofview=match.iloc[0]['FieldofView']
        thisdpi=(int((512*2.54)/(fieldofview)),int((512*2.54)/(fieldofview))) 
        jpgimage.convert('RGB').save(annotjpgname, dpi=thisdpi)
        jpgimage.close()
    except:
        print('Problem creating jpg with overlaid map region for', jpgfname)
    return 

def superimposearr(QMpixarray, allregs=True, crop=True):
    ''' Superimpose pix arrays (boundaries or all boxes) on SE image 
    '''
    fullpath= filedialog.askopenfilename(title='Select SE image', 
                                             filetypes=[("JPG","*.jpg")])
    (directory, filename)=os.path.split(fullpath)
    jpgimage=Image.open(filename)
    draw=ImageDraw.Draw(jpgimage)
    if allregs:
        for index, row in QMpixarray.iterrows():
            x1=QMpixarray.loc[index]['X1']
            y1=QMpixarray.loc[index]['Y1']
            x2=QMpixarray.loc[index]['X2']
            y2=QMpixarray.loc[index]['Y2']
            draw.rectangle((x1,y1,x2,y2), outline='red') # red rectangle at specified position
    else: # Only draw outer boundary
        x1=QMpixarray.X1.min()
        x2=QMpixarray.X2.max()
        y1=QMpixarray.Y1.min()
        y2=QMpixarray.Y2.max()
        draw.rectangle((x1,y1,x2,y2), outline='red')
    jpgname=filename.replace('.jpg','_overlay.jpg')
    jpgimage.convert('RGB').save(jpgname)

    if crop: # separately save cropped QM area
        x1=QMpixarray.X1.min()
        x2=QMpixarray.X2.max()
        y1=QMpixarray.Y1.min()
        y2=QMpixarray.Y2.max()
        jpgimage=jpgimage.crop((x1,y1,x2,y2))
        jpgname=filename.replace('.jpg','_crop.jpg')
        jpgimage.convert('RGB').save(jpgname)
    jpgimage.close()
    return

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
    '''Create elemental amplitude maps from QM elemental amplitudes (smooth-diff amplitude for
    each line in Smdifpeakslog and area information in spatialareas '''
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
    shutil.copyfile('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\spatial_areas_sample_min.phi',autotoolname)
    for line in fileinput.input(autotoolname, inplace=1):
        sys.stdout.write(datastr)
    return

def makeautotool(filelist, multacq='QM_multiplex.phi', **kwargs):
    '''Generates df with Autotool commands and data values (for generating Autotool phi file)
    7/10 Modified with image reg insertion 
    9/5/17 using C:\Temp doesn't work... load in current settings directory (no pre-path)
    kwarg: regint - interval for insertion of image registrations 
        multiplex breaks -- used in shifted situation (series of files at different shifts '''
    df=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\QM_autotool.csv')
    atframe=df.loc[0:1] # register image, take photo, load QMmultiplex.phi
    # filelist.sort(reverse=True) # sort just screws things up
    mycolumns=['Command','Data']
    if 'regint' in kwargs:
        regint=kwargs.get('regint',1)
    if 'multibreaks' in kwargs:
        multibreaks=kwargs.get('multibreaks',[])
        multinames=kwargs.get('multinames',[])
        multacq=multinames[0] # set to first shifted multiplex to load
        if 0 in multibreaks: # first multiplex file loaded w/ multacq
            multibreaks.remove(0)
    # Add first multiplex load (rest are below)
    newrow=pd.DataFrame(index=np.arange(0,1), columns=mycolumns)
    newrow=newrow.set_value(0,'Command','AES:Load Multiplex Setting...')
    newrow=newrow.set_value(0,'Data', multacq)
    atframe=pd.concat([atframe,newrow], ignore_index=True)

    for i, file in enumerate(filelist):
        if 'regint' in kwargs:
            if i%regint==0: # no image reg this cycle 
                # create double rowed frame
                newrow=pd.DataFrame(index=np.arange(0,3), columns=mycolumns)
                newrow=newrow.set_value(0,'Command','AES:Load Area Define Setting...')
                newrow=newrow.set_value(0,'Data', file)
                newrow=newrow.set_value(1,'Command','AES:Multiplex Acquire')
                newrow=newrow.set_value(2,'Command','AES:Register Image')
            else: # has regint but not this cycle
                # don't add path... areas must be in settings/area definition folder 
                newrow=pd.DataFrame(index=np.arange(0,2), columns=mycolumns)
                newrow=newrow.set_value(0,'Command','AES:Load Area Define Setting...')
                newrow=newrow.set_value(0,'Data',file)
                newrow=newrow.set_value(1,'Command','AES:Multiplex Acquire')
        else:
            newrow=pd.DataFrame(index=np.arange(0,2), columns=mycolumns)
            newrow=newrow.set_value(0,'Command','AES:Load Area Define Setting...')
            newrow=newrow.set_value(0,'Data', file)
            newrow=newrow.set_value(1,'Command','AES:Multiplex Acquire')
        atframe=pd.concat([atframe,newrow], ignore_index=True)
        # Now add load of next shifted multiplex file
        if 'multibreaks' in kwargs:
            if i in multibreaks:
                lindex=multibreaks.index(i)
                newrow=pd.DataFrame(index=np.arange(0,1), columns=mycolumns)
                newrow=newrow.set_value(0,'Command','AES:Load Multiplex Setting...')
                # multfile must be in settings/multiplex acquire (no file extensions in Autotool data cols)
                newrow=newrow.set_value(0,'Data', multinames[lindex+1].replace('.phi',''))
                atframe=pd.concat([atframe,newrow], ignore_index=True)
    atframe=pd.concat([atframe,df.loc[[3]]], ignore_index=True) # ending SEM photo
    atframe=pd.concat([atframe,df.loc[[4]]], ignore_index=True) # add beam deflection
    return atframe

def makesquarearray(margin, arraysize, basename, **kwargs):
    ''' Divide up 512x512 pixels in map into n areas and format correctly for spatial areas phi files
    (which are loaded using Autotool loops into PHI Smartsoft); 
    Mapping proceed horizontally (x horizontal, y vertical) 
    ARGS:
    arraysize -- # of pixels in one direction
    margin - % of image field that is unmapped (20% means 10% of 512 field at both edges (aka 51 pixels))
        is unmapped
    basename - basename for area definition files 
        e.g. "50x50array20m" basemname makes files 50x50array20m1, 50x50array20m2, etc.
        
    KWARGS: 'regint' - interval at which to build in image registration into autotool loop; val of 1 means
           register every 20 pixels (since each area file holds 20 defined spatial areas);  passed to makeautotool
           * this is best way to incorporate image registration in quantmap process... more flexible
           interval allowed; if instead one builds in image reg into multiplex itself, one has to 
           run image registration much more frequently which unnecessarily increases acquisition time
        'writeareas' (bool) - write new spatial area definition files (old ones can be reused if same arraysize 
             and margin)
        'writeAutotool' --  write of new Autotool sequence (can be reused if same arraysize/regint)
        scrambled
    '''
    pix=512
    width=(pix*(1-margin)/arraysize) # width/height of scan pixel in terms of 512x512 field
    startxy=int(pix*margin/2) # split margin between top/bottom, left/right
    mycols=['Xindex','Yindex','Areanumber','PHIname','Subnumber','X1','Y1','X2','Y2', 'Width', 'Height']
    dim=arraysize**2
    # square is the pixarray file correlating 101.spe with associated pixel in quantmap 
    square=pd.DataFrame(index=np.arange(0,dim), columns=mycols)
    # x is horizontal axis and mapping proceeds by going across top row
    for index,row in square.iterrows():
        xindex=index//arraysize # remainder is row (0th is top row)
        yindex=index%arraysize # mod is correct column (0th is left column)
        square.loc[index]['Xindex']=xindex # remainder is row
        square.loc[index]['Yindex']=yindex # mod is correct column
        left=int(width*yindex+startxy) # left-right position depends on column 
        square.loc[index]['X1']=left
        right=int(width*yindex+startxy+width)
        square.loc[index]['X2']=right
        top=int(width*xindex+startxy) 
        square.loc[index]['Y1']=top # top-bottom position depends on row
        bottom=int(width*xindex+startxy+width)
        square.loc[index]['Y2']=bottom
        square.loc[index]['Width']=right-left # variations due to rounding error
        square.loc[index]['Height']=bottom-top
        # true area number describing pix position after file combination
        square.loc[index]['Areanumber']=index+1
        # max of 20 areas allowed per spatial area .phi file
        square.loc[index]['Subnumber']=index%20 # Before combination w/ 20 areas per file
        filenum=index//20+1 # filenumber of multiplex 
        square.loc[index]['PHIname']=basename+str(filenum)
    filelist=square.PHIname.unique()
    filelist=np.ndarray.tolist(filelist)
    # map/move beam non-sequentially to minimize local charging
    if 'scrambled' in kwargs:
        areanums=np.ndarray.tolist(square.Areanumber.unique())
        areanums=np.random.permutation(areanums)
        square['Areanumber']=pd.Series(areanums)
        square=square.sort_values(['Areanumber'])
        # need to reassign subnumber and PHIname
        square=square.reset_index(drop=True)
        for index,row in square.iterrows():
            square=square.set_value(index,'Subnumber', index%20)
            square=square.set_value(index,'PHIname', basename+str(index//20+1))
    if 'writeareas' in kwargs:
        for i, fname in enumerate(filelist):
            thisfile=square[square['PHIname']==fname]
            writeAESareas(thisfile, fname) # writes each list of 20 areas to separate .phi text file
    if 'writeAutotool' in kwargs:
        atkwargs={}
        if 'regint' in kwargs: # pass image registration interval 
            atkwargs.update({'regint':kwargs.get('regint',0)})
        # TODO preload multiplex file and skip adding multacq step?
        atframe=makeautotool(filelist, multacq='QM_multiplex', **atkwargs)
        # Use autonaming convention for autotool file
        ATname='QM_autotool_'+str(arraysize)+'x'+str(arraysize)+'m'+str(int(margin*100))+'.phi'
        writeautotool(atframe, ATname) 
    # instead of C:\Temp copy multiplex and spatial areas files to Smartsoft settings folders
    return square

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
    shutil.copyfile('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\spatial_areas_sample_min.phi',filename)
    for line in fileinput.input(filename, inplace=1):
        sys.stdout.write(datastr)
    return

def makemultdf(elemparams, peaks, lowers, uppers, sweeps, widths):
    ''' Reconstruct dataframe holding altered multiplex scan parameters 
    then feed to writemultiplex'''
    mycols=['AtomNum', 'Elem', 'Active', 'Sweeps', 'EVstep', 'Lower', 'Upper',
       'Range', 'Lowpeak', 'Peak', 'Hipeak', 'Back1', 'Back2']
    multdf=pd.DataFrame(columns=mycols)
    atnos=[]
    regs=[]
    for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
        atnos.append(atno)
        atnos.append(atno)
        atnos.append(atno)
        regs.append(elem+'L')
        regs.append(elem)
        regs.append(elem+'H')
    multdf['AtomNum']=atnos
    multdf['Elem']=regs
    multdf['Active']='Y'
    multdf['Sweeps']=sweeps
    multdf['Sweeps']=sweeps
    multdf['EVstep']=1.0
    multdf['Peak']=peaks
    multdf['Lower']=lowers
    multdf['Upper']=uppers
    
    #multdf['Lower']=multdf['Peak']-multdf['Range']
    #multdf['Upper']=multdf['Peak']+multdf['Range']
    # unclear if these params are even used
    multdf['Lowpeak']=multdf['Lower']+2
    multdf['Hipeak']=multdf['Upper']-2
    multdf['Back1']=multdf['Lower']
    multdf['Back2']=multdf['Upper']
    # convert half-widths to full widths
    widths=[2*i+1 for i in widths]
    multdf['Range']=widths
    # Eliminate any overlaps in scan range (remove from L or H not from main peak)
    # overlaps should already be gone b/c of tk gui methods 
    for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(elemparams):
        lowend=multdf.iloc[i]['Upper']
        mainstart=multdf.iloc[i+1]['Lower']
        mainend=multdf.iloc[i+1]['Upper']
        highstart=multdf.iloc[i+2]['Lower']
        # print(str(lowend), str(mainstart), str(mainend),str(highstart))
        if mainstart<lowend:
            multdf=multdf.set_value(multdf.index[i],'Upper', mainstart-1)
        if highstart<mainend:
            multdf=multdf.set_value(multdf.index[i+2],'Lower', mainend+1)        
    return multdf

def writemultiplex(multdf, phiname, dwell, numcycles, **mkwargs):
    ''' Write of multiplex settings file (max 20 regions) after interactive param setting
    Some weird encoding so just read and modify existing pff file 
    image registration choices are cycles or areas
    kwargs:
        regmode - Areas or  (if not present, image registration done using autotool not in multiplex)
        reginterval - 2 (or whatever)'''


    # TODO pass as arg if others use this 
    datastr=''
    datastr+='[MultiplexAcq]\n'
    # image registration option within multiplex (also can be done in autotool at lower frequencies)
    if 'regmode' in mkwargs:
        regmode=mkwargs.get('regmode','Areas')
        reginterval=mkwargs.get('reginterval', 2)
        'Register Image=True\nImage Registration Interval='
        datastr+=str(reginterval)
        datastr+='\n'
        datastr+='Image Registration Mode='
        datastr+=regmode
    else:
        datastr+='Register Image=False'
        # likely can skip other params if this is false
    datastr+='\nTime Per Step (ms)='
    datastr+='%0.1f' % dwell
    datastr+='\nNegative Values=Allow\nNumber of Regions='
    datastr+=str(len(multdf)) # number of multiplex regions
    datastr+='\nAtomic Number List Count='
    datastr+=str(len(multdf))
    datastr+='\n'
    for i in range(0,len(multdf)):
        datastr+='Atomic Number List '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['AtomNum']
        datastr+='%d' % val
        datastr+='\n'
    datastr+='Element List Count='
    datastr+=str(len(multdf))
    datastr+='\n'
    for i in range(0,len(multdf)):
        datastr+='Element List '
        datastr+='%d' % i
        datastr+='='
        strval=multdf.iloc[i]['Elem']
        datastr+=strval
        datastr+='\n'
    datastr+='Active Count='
    datastr+=str(len(multdf))
    datastr+='\n'
    for i in range(0,len(multdf)):
        datastr+='Active '
        datastr+='%d' % i
        datastr+='=True\n'  # won't be present in Df if not active
    datastr+='Sweeps Count='
    datastr+=str(len(multdf))
    datastr+='\n' 
    for i in range(0,len(multdf)):
        datastr+='Sweeps '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Sweeps']
        datastr+='%d' % val
        datastr+='\n'
    datastr+='Lower Acq Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Lower Acq '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Lower']
        datastr+='%0.1f' % val # float with tenths place precision
        datastr+='\n'
    datastr+='Upper Acq Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Upper Acq '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Upper']
        datastr+='%0.1f' % val
        datastr+='\n'        
    datastr+='Acq Range Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Acq Range '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Range']
        datastr+='%0.1f' % val
        datastr+='\n'    
    datastr+='Lower Analysis Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Lower Analysis '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Lowpeak']
        datastr+='%0.1f' % val
        datastr+='\n'     
    datastr+='Upper Analysis Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Upper Analysis '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Hipeak']
        datastr+='%0.1f' % val
        datastr+='\n'    
    datastr+='eV Per Step Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='eV Per Step '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['EVstep']
        datastr+='%0.1f' % val
        datastr+='\n'   
    datastr+='Peak Energy Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Peak Energy '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Peak']
        datastr+='%0.1f' % val
        datastr+='\n'  
    datastr+='Background 1 Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Background 1 '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Back1']
        datastr+='%0.1f' % val
        datastr+='\n'  
    datastr+='Background 2 Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Background 2 '
        datastr+='%d' % i
        datastr+='='
        val=multdf.iloc[i]['Back2']
        datastr+='%0.1f' % val
        datastr+='\n'  
    datastr+='Number Of Channels Count='
    datastr+=str(len(multdf))
    datastr+='\n'     
    for i in range(0,len(multdf)):
        datastr+='Number Of Channels '
        datastr+='%d' % i
        datastr+='=1 thru 8\n'
    datastr+='Number of Cycles='
    datastr+=str(numcycles)
    datastr+='\n'
    # Write this chunk of files to .phi spatial areas file (done w/ file replace method since encoding is weird unknown type)
    shutil.copyfile('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\spatial_areas_sample_min.phi',phiname)
    for line in fileinput.input(phiname, inplace=1):
        sys.stdout.write(datastr)
    print(phiname,' saved to current working directory')
    return

def plot_2_maps(shiftmap, amplmap):
    ''' Side-by-side plot of charging and associated strength of sm-diff amplitude (from which peak
    O peak and subsequent charging were determined; call after findnegpeaks '''
    fig, axes = plt.subplots(nrows=1, ncols=2, squeeze=False)
    axes[0,0].set_title("Charging magnitude (eV)", fontsize=12)
    axes[0,1].set_title("Peak amplitude (derivative)", fontsize=12)
    axes[0,0].imshow(shiftmap) # deriv [0] or integ [1] based 
    axes[0,1].imshow(amplmap) # 0th layer is sm-diff amplitude
    return

def makepixlist(myarr):
    ''' Get list of X,Y pixel coords for masked array pixels '''
    pixlist=[]
    for X in range(0,myarr.shape[0]):
        for Y in range(0,myarr.shape[1]):
            if np.ma.is_masked(myarr[X,Y]):
                pixlist.append([X,Y])
    return pixlist

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

def shiftAESparams(AESquantparams, chargeshift):
    ''' Return copy of AES quant params shifted by uniform amount
    used by tk on-the-fly plotting  '''
    AESqp=AESquantparams.copy()
    # need to alter negpeak, pospeak, QMlow, QMhigh
    AESqp['negpeak']+=chargeshift
    AESqp['pospeak']+=chargeshift
    AESqp['QMlow']+=chargeshift
    AESqp['QMhigh']+=chargeshift
    AESqp=AESqp[pd.notnull(AESquantparams['plotrange'])]
    for index, row in AESqp.iterrows():
        val=row.plotrange
        val=val.replace('[','')
        val=val.replace(']','')
        newlow=int(val.split(',')[0].strip())+chargeshift
        newhigh=int(val.split(',')[1].strip())+chargeshift
        newstr='['+str(newlow)+','+str(newhigh)+']'
        AESqp=AESqp.set_value(index,'plotrange', newstr)
    return AESqp

def replacemaskpix(arr1, arr2):
    ''' Replace masked regions of arr1 with values from arr2 (usually median filtered) '''
    newarr=np.empty([arr1.shape[0],arr1.shape[1]])
    for X in range(0,arr1.shape[0]):
        for Y in range(0,arr1.shape[1]):
            if np.ma.is_masked(arr1[X,Y]):
                newarr[X,Y]=arr2[X,Y]
            else:
                newarr[X,Y]=arr1[X,Y]
    return newarr

def replacemaskpix2(arr1, arr2, arr3):
    ''' Use arr1 unless masked in arr3 (then use arr2 value; normally median, uniform or other filter 
    '''
    newarr=np.empty([arr1.shape[0],arr1.shape[1]])
    for X in range(0,arr1.shape[0]):
        for Y in range(0,arr1.shape[1]):
            if np.ma.is_masked(arr3[X,Y]):
                newarr[X,Y]=arr2[X,Y]
            else:
                newarr[X,Y]=arr1[X,Y]
    return newarr

def plothisto(myarr, nbins=15):
    ''' Make histogram from 2D np values '''
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,9))
    hist, bins = np.histogram(myarr, nbins)
    center=(bins[:-1]+bins[1:])/2
    width = 0.7 * (bins[1] - bins[0])
    axes.bar(center, hist, align='center', width=width)
    return

def reorder_spatial_array(QMpixarray, chargemap, basename):
    ''' using spatial map of charging (modified chargemap 2D np array) assign a charging value to 
    each pixel in mapped array, reorder by charging magnitude, and create corresponding spatial 
    area files, then save to pixarray 
    '''
    for index, row in QMpixarray.iterrows():
        QMpixarray=QMpixarray.set_value(index,'Shift',chargemap[row.Xindex, row.Yindex])
    QMpixarray=QMpixarray.sort_values(['Shift','Yindex','Xindex'])
    # reset index
    QMpixarray=QMpixarray.reset_index(drop=True)
    QMpixarray['Filename']=np.nan
    # reset phiname
    for index, row in QMpixarray.iterrows():
        QMpixarray=QMpixarray.set_value(index,'Subnumber',index%20)
        QMpixarray=QMpixarray.set_value(index,'Areanumber',index+1)
        filenum=index//20+1
        QMpixarray=QMpixarray.set_value(index,'PHIname',basename+str(filenum))
    QMpixarray.to_csv('QMpixarray_shifted.csv', index=False)
    # Write new scrambled spatial area files 
    filelist=np.ndarray.tolist(QMpixarray.PHIname.unique())
    for i, fname in enumerate(filelist):
        thisfile=QMpixarray[QMpixarray['PHIname']==fname]
        writeAESareas(thisfile, fname) # writes each list of 20 areas to separate .phi text file
    # shifted multiplex setup and corresponding autotool creation done w/ tk shifted interface
    return QMpixarray

def pickbestderivpeak(mypeaks, idealev, idealnegpeak):
    ''' Decision algorithms for which amplitude and direct integ result to save
    need to flag marginal results (esp weak s7d7)
    Either deriv or integ (due to backfit errors) can be erroneous
    First 3 rows of mypeaks are best peaks found w/ deriv (and associated integrations)
    last row has integration around subdata maximum

    Choose best version and return series
    '''
    diffthresh=5
    # Allow choice of secondary deriv peaks as small as 50% if consistent w/ direct peak
    amplthresh=0.5 
    
    # Comparison of multiple deriv peaks w/ direct peak integcnts at max
    mypeaks['Normampl']=mypeaks['Ampl']/mypeaks['Ampl'].max()
    # Lastrow integcounts value very likely to be largest (integ around max val)
    peakdiff1=mypeaks.iloc[0]['Energy']+(idealev-idealnegpeak)-mypeaks.iloc[-1]['Energy2']
    peakdiff2=mypeaks.iloc[1]['Energy']+(idealev-idealnegpeak)-mypeaks.iloc[-1]['Energy2']
    peakdiff3=mypeaks.iloc[1]['Energy']+(idealev-idealnegpeak)-mypeaks.iloc[-1]['Energy2']
    # choose superior 2nd peak 
    if len(mypeaks)>1:
        if abs(peakdiff1)> diffthresh and abs(peakdiff2)<diffthresh and mypeaks.iloc[1]['Normampl']>amplthresh:
            mypeaks=mypeaks.set_value(1,'Energy2',mypeaks.iloc[-1]['Energy2'])
            mypeaks=mypeaks.set_value(1,'Countsmax',mypeaks.iloc[-1]['Countsmax'])
            # use slightly superior integ around maximum
            mypeaks=mypeaks.set_value(1,'Integcounts',mypeaks.iloc[-1]['Integcounts'])
            mypeaks=mypeaks.set_value(1,'Peakind2',mypeaks.iloc[-1]['Peakind2'])
            mypeaks=mypeaks[1:2].squeeze()
            return mypeaks
    # choose superior 3rd peak
    if len(mypeaks)>2:
        if abs(peakdiff1)> diffthresh and abs(peakdiff3)<diffthresh and mypeaks.iloc[2]['Normampl']>amplthresh:
            mypeaks=mypeaks.set_value(2,'Energy2',mypeaks.iloc[-1]['Energy2'])
            mypeaks=mypeaks.set_value(2,'Countsmax',mypeaks.iloc[-1]['Countsmax'])
            # use slightly superior integ around maximum
            mypeaks=mypeaks.set_value(2,'Integcounts',mypeaks.iloc[-1]['Integcounts'])
            mypeaks=mypeaks.set_value(2,'Peakind2',mypeaks.iloc[-1]['Peakind2'])
            mypeaks=mypeaks[2:3].squeeze()
            return mypeaks
    # default case (choose first peak)
    mypeaks=mypeaks.set_value(0,'Energy2',mypeaks.iloc[-1]['Energy2'])
    mypeaks=mypeaks.set_value(0,'Countsmax',mypeaks.iloc[-1]['Countsmax'])
    # use slightly superior integ around maximum
    mypeaks=mypeaks.set_value(0,'Integcounts',mypeaks.iloc[-1]['Integcounts'])
    mypeaks=mypeaks.set_value(0,'Peakind2',mypeaks.iloc[-1]['Peakind2'])
    mypeaks=mypeaks[0:1].squeeze()
    return mypeaks

def pickrandompixels(specimage, numpix):
    ''' Grab a set of random pixels from spectral image and return as list
    avoids np.nan pixels from specimage'''
    pixlist=[]
    for i in range(0,numpix):
        randx=random.randrange(0,specimage.shape[0]-1)
        randy=random.randrange(0,specimage.shape[1]-1)
        if specimage[randx, randy, 0]!=0:
            pixlist.append([random.randrange(0,specimage.shape[0]-1), random.randrange(0,specimage.shape[1]-1)])
    return pixlist

def ratioimages(amplmapdict, elem1, elem2):
    ''' Make ratio images from elemental maps '''
    if elem1 and elem2 not in amplmapdict.keys():
        print('Elements', elem1, elem2,' missing from dictionary.')
        return 
    else:
        arr1=amplmapdict.get(elem1)
        arr2=amplmapdict.get(elem2)
        rat=np.divide(arr1, arr2)
        return rat

def calcbackground(peakdata, lowind, highind):
    ''' Linear fit of e- background under peak using low and high background regions
    on single X, Y pixel of data for chosen element...  returns slope, intercept
    single pixel version of calcbackgrounds WORKS... 
    '''
    backrange=[]
    [lowmin, lowmax]=lowind # unpack index range of elem/peak
    [highmin, highmax]=highind
    backrange=[int(i) for i in range(lowmin, lowmax+1)]
    backrange.extend([int(i) for i in range(highmin, highmax+1)])    
    backrange=peakdata[peakdata.index.isin(backrange)]
    data1=backrange['Energy']
    data2=backrange['Counts']
    slope,intercept=np.polyfit(data1, data2, 1)
    return slope,intercept

def fixoddsetup(lowind, lowrange, highind, hirange, peakind, peakrange):
    ''' Standard multiplex setup is 3 regions per peak (low, peak and high)
    if done in single scan some parameters need to be altered '''
    numchan=int(0.1*(peakind[1]-peakind[0]))
    lowind[0]=peakind[0]
    lowind[1]=peakind[0]+numchan
    lowrange[0]=peakrange[0]
    lowrange[1]= peakrange[0]+ (peakind[1]-peakind[0])*(lowind[1]/(peakind[1]-peakind[0]))
    highind[0]=peakind[1]-numchan
    highind[1]=peakind[1]
    hirange[0]= peakrange[0] + (peakind[1]-peakind[0])*(highind[0]/(peakind[1]-peakind[0]))       
    hirange[1]=peakrange[1]
    peakind =[lowind[1]+1,highind[0]-1] # reset peak range to compensate
    peakrange=[lowrange[1]+1, hirange[0]-1 ]
    return lowind, lowrange, highind, hirange, peakind, peakrange

''' TESTING
row=mypeaks.iloc[0]
'''

def integpeaks(peakdata, mypeaks, slope, intercept, idealnegpeak, idealev, integwidth):
    ''' Compute direct integral corresponding to each potential peak (found in s7d7)
    adds integcounts (integration over integwidth at positions suggested by s7d7 peak)
    includes integchan (# of integrated channels since one can hit bounds of data)
     '''

    # For each found negpeak, do integration over integwidth centered on idealnegpeak-idealev
    mypeaks['Integcounts']=np.nan
    mypeaks['Integchan']=np.nan # number of integrated channels (can hit data boundaries)
    mypeaks['Peakind2']=np.nan
    # have to reset index on peakdata since argrelmin index in mypeaks has been reset
    peakdata=peakdata.reset_index(drop=True)
    for index, row in mypeaks.iterrows():
        # Dealing with hitting edges, etc.
        # negpeakind has reset index (0 is first row) whereas peakdata has original index
        centerind=int(row.negpeakind-int(idealnegpeak-idealev))
        lowlim=centerind-int((integwidth-1)/2) # integwidth is full peak integration width
        uplim=centerind+int((integwidth-1)/2)
        # Ensure that limits are within data range
        if lowlim<0:
            lowlim=0
        if uplim>len(peakdata)-1:
            uplim=len(peakdata)-1
        # TODO maybe add back symmetric part on other side
        # Integcounts is sum of subdata over lowlim to uplim
        integcnts=peakdata.loc[lowlim:uplim]['Subdata'].sum()
        mypeaks=mypeaks.set_value(index,'Integcounts', integcnts)
        mypeaks=mypeaks.set_value(index,'Integchan', uplim-lowlim+1)
        mypeaks=mypeaks.set_value(index,'Peakind2', centerind)
        # TODO maybe also include and save backcounts?
    # Also include integration around subdata maximum
    newrow=pd.Series()
    centerind=peakdata['Subdata'].idxmax()
    lowlim=centerind-int((integwidth-1)/2)
    uplim=centerind+int((integwidth-1)/2)
    if lowlim<0:
        lowlim=0
    if uplim>len(peakdata)-1:
        uplim=len(peakdata)-1
    newrow['Integcounts']=peakdata.loc[lowlim:uplim]['Subdata'].sum()
    # Store index in terms of where negpeak would be (shifted from direct peak)
    newrow['Peakind2']=centerind
    newrow['Integchan']=uplim-lowlim+1
    newrow['Energy2']=peakdata.loc[centerind]['Energy']
    # also include raw signal level at max of subdata (peak position) 
    newrow['Countsmax']=peakdata.loc[centerind]['Counts']
    mypeaks=mypeaks.append(newrow, ignore_index=True)
    return mypeaks

def savemaps(amplmaps, shiftmaps, integmaps, basename):
    ''' Save stacked numpy arrays with maps extracted from specimage using 
    find_all_peaks function  '''
    dim1=amplmaps[0].shape[0]
    dim2=amplmaps[0].shape[1]
    temp=np.zeros((dim1,dim2, 4*len(amplmaps)))
    for i, thismap in enumerate(amplmaps):
        temp[:,:,4*i:4*i+4]=thismap
    np.save(basename+'_amplmaps.npy', temp)
    
    # Manual save of stacked numpy arrays
    temp=np.zeros((dim1,dim2, 2*len(shiftmaps)))
    for i, thismap in enumerate(shiftmaps):
        temp[:,:,2*i:2*i+2]=thismap
    np.save(basename+'_shiftmaps.npy', temp)
    
    temp=np.zeros((dim1,dim2, 5*len(integmaps)))
    for i, thismap in enumerate(integmaps):
        temp[:,:,5*i:5*i+5]=thismap
    np.save(basename+'_integmaps.npy', temp)
    return

def loadmaps(dir):
    ''' Load of assorted maps files from current directory 
    these are stack np arrays (different vals for each element)
    same structure as in quantmap class'''
    npyfiles=glob.glob('*amplmaps.npy')
    amplmaps=[]
    if len(npyfiles)==1:
        temp=np.load(npyfiles[0])
        # Blank container for one elems shift files (deriv [0] and direct[1])
        blank=np.zeros((temp.shape[0], temp.shape[1], 4))
        for i in range(0, int(temp.shape[2]/4)):
            blank=temp[:,:, 4*i:4*i+4]
            amplmaps.append(blank)
    else:
        print("Couldn't locate amplitude map files in current directory")
    npyfiles=glob.glob('*shiftmaps.npy')
    shiftmaps=[]
    if len(npyfiles)==1:
        temp=np.load(npyfiles[0])
        # Blank container for one elems shift files (deriv [0] and direct[1])
        blank=np.zeros((temp.shape[0], temp.shape[1], 2))
        for i in range(0, int(temp.shape[2]/2)):
            blank=temp[:,:, 2*i:2*i+2]
            shiftmaps.append(blank)
    else:
        print("Couldn't locate shift map files in current directory")
    npyfiles=glob.glob('*integmaps.npy')
    integmaps=[]
    if len(npyfiles)==1:
        temp=np.load(npyfiles[0])
        # Blank container for one elems shift files (deriv [0] and direct[1])
        blank=np.zeros((temp.shape[0], temp.shape[1], 5))
        for i in range(0, int(temp.shape[2]/5)):
            blank=temp[:,:, 5*i:5*i+5]
            integmaps.append(blank)
    else:
        print("Couldn't locate integcounts map files in current directory")
    npyfiles=glob.glob('*elemmaps.npy')
    # TODO figure out the element map data scheme
    elemmaps=[]
    if len(npyfiles)==1:
        temp=np.load(npyfiles[0])
        # Blank container for one elems shift files (deriv [0] and direct[1])
        blank=np.zeros((temp.shape[0], temp.shape[1], 5))
        for i in range(0, int(temp.shape[2]/5)):
            blank=temp[:,:, 1*i:1*i+1]
            elemmaps.append(blank)
    else:
        print("Couldn't locate elem map files in current directory")
    return amplmaps, shiftmaps, integmaps, elemmaps
