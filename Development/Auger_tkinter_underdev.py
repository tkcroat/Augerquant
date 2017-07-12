# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 10:32:56 2017

@author: tkc
"""
import tkinter as tk
import datetime

def AESreport_tk(spelist, Elements, Smdifflog, Backfitlog, AESquantparams, **kwargs):
    ''' tk interface for args/kwargs of AES report functions
    all args/dataframes must be passed through to plot functions 
    kwargs: fileset, areas, xrangestr -- usually just holds values entered during prior run 
    report will use all areas '''
    # first print out existing info in various lines
    root = tk.Tk()
    root.title('Generate PDF plot report')
    filestr=tk.StringVar() # comma separated or range of filenumbers for plot
    filestr.set(kwargs.get('fileset','')) # get entry from prior run
    xrangestr=tk.StringVar()  # energy range in eV 
    xrangestr.set(kwargs.get('xrangestr',''))
    plottype=tk.StringVar() # column for plot report/ radio1 choice
    rangechoice=tk.StringVar() # string for x data range choice in radio2
    smdifbool=tk.BooleanVar() # Bool for plotting of peak locations in smdifpeakslog
    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    backptbool=tk.BooleanVar() 
    elemstr=tk.StringVar()
    mytext=', '.join(Elements) # elements for labelling 
    elemstr.set(mytext)
    plotelemsbool=tk.BooleanVar()  # optional labelling of elements
    plotelemsbool.set(True) # default to true
    now=datetime.datetime.now()
    PDFname='Countsback_report'+'_'+now.strftime('%d%b%y')+'.pdf'
    PDFnamestr=tk.StringVar()  # energy range in eV 
    PDFnamestr.set(PDFname)
    choice=tk.StringVar()  # plot or abort
    # Functions to enable/disable relevant checkboxes depending on radiobutton choice
    def Countopts():
        ''' Disable irrelevant checkboxes '''
        xrangentry.config(state=tk.NORMAL) # why isn't this working?
        smdiff.config(state=tk.DISABLED)
        backfit.config(state=tk.NORMAL)
        backpt.config(state=tk.NORMAL)
        PDFnamestr.set('Counts_report'+'_'+now.strftime('%d%b%y')+'.pdf')
    def Derivopts():
        ''' Disable irrelevant checkboxes '''
        xrangentry.config(state=tk.NORMAL) 
        smdiff.config(state=tk.NORMAL)
        backfit.config(state=tk.DISABLED)
        backpt.config(state=tk.DISABLED)
        PDFnamestr.set('Deriv_report'+'_'+now.strftime('%d%b%y')+'.pdf')
    def Bothopts():
        ''' Disable irrelevant checkboxes '''
        xrangentry.config(state=tk.NORMAL) 
        smdiff.config(state=tk.NORMAL)
        backfit.config(state=tk.NORMAL)
        backpt.config(state=tk.NORMAL)
        PDFnamestr.set('Countsderiv_report'+'_'+now.strftime('%d%b%y')+'.pdf')
    def Peakopts():
        ''' Disable irrelevant checkboxes '''
        xrangentry.config(state=tk.DISABLED)  # automatically chooses narrow region around peak
        smdiff.config(state=tk.DISABLED)
        backfit.config(state=tk.DISABLED)
        backpt.config(state=tk.DISABLED)
        PDFnamestr.set('Peak_report'+'_'+now.strftime('%d%b%y')+'.pdf')
    # Choose counts, deriv, both or peaks plot
    a=tk.Entry(root, textvariable=PDFnamestr).grid(row=0, column=0)    
    a=tk.Label(root, text='Plot report name and type').grid(row=0, column=1)
    radio1 = tk.Radiobutton(root, text='Counts', value='Counts', variable = plottype, command=Countopts).grid(row=1, column=0)
    radio1 = tk.Radiobutton(root, text='Derivative', value='Deriv', variable = plottype, command=Derivopts).grid(row=1, column=1)
    radio1 = tk.Radiobutton(root, text='Both', value='Both', variable = plottype, command=Bothopts).grid(row=1, column=2)
    radio1 = tk.Radiobutton(root, text='Peaks', value='Peaks', variable = plottype, command=Peakopts).grid(row=1, column=3)
    backfit=tk.Checkbutton(root, variable=backfitbool, text='Plot background fits?')
    backfit.grid(row=2, column=0) # can't do immediate grid or nonetype is returned
    backpt=tk.Checkbutton(root, variable=backptbool, text='Plot background points?')
    backpt.grid(row=2, column=1) 

    smdiff=tk.Checkbutton(root, variable=smdifbool, text='Plot smdiff peak locations?')
    smdiff.grid(row=3, column=0)
    # Mark element positions (option)
    d=tk.Checkbutton(root, variable=plotelemsbool, text='Label element peaks?')
    d.grid(row=3, column=1)
    
    def Elemopts():
        ''' Disable irrelevant checkboxes '''
        # xrangentry.config(state=tk.DISABLED)  # automatically chooses narrow region around peak
        elementry.config(state=tk.NORMAL)
        xrangentry.config(state=tk.DISABLED)
    def Evopts():
        ''' Disable irrelevant checkboxes '''
        # xrangentry.config(state=tk.DISABLED)  # automatically chooses narrow region around peak
        elementry.config(state=tk.DISABLED)
        xrangentry.config(state=tk.NORMAL)
    def Joinopts():
        ''' Disable irrelevant checkboxes '''
        # xrangentry.config(state=tk.DISABLED)  # automatically chooses narrow region around peak
        elementry.config(state=tk.NORMAL)
        xrangentry.config(state=tk.NORMAL)
    def Fullopts():
        ''' Disable irrelevant checkboxes '''
        # xrangentry.config(state=tk.DISABLED)  # automatically chooses narrow region around peak
        elementry.config(state=tk.DISABLED)
        xrangentry.config(state=tk.DISABLED)
    # Radio button for range (use element list, use ev strings or both)
    a=tk.Label(root, text='Choose x data range').grid(row=4, column=0)
    radio2 = tk.Radiobutton(root, text='Plot element ranges', value='elems', variable = rangechoice, command=Elemopts).grid(row=5, column=0)
    radio2 = tk.Radiobutton(root, text='Plot eV range', value='evrange', variable = rangechoice, command=Evopts).grid(row=5, column=1)
    radio2 = tk.Radiobutton(root, text='Join elements and eV range', value='both', variable = rangechoice, command=Joinopts).grid(row=5, column=2)
    radio2 = tk.Radiobutton(root, text='Use full data range', value='full', variable = rangechoice, command=Fullopts).grid(row=5, column=3)

    a=tk.Label(root, text='Elements:').grid(row=6, column=0)
    elementry=tk.Entry(root, textvariable=elemstr)
    elementry.grid(row=6, column=1)

    a=tk.Label(root, text='eV range:').grid(row=6, column=2)
    xrangentry=tk.Entry(root, textvariable=xrangestr)
    xrangentry.grid(row=6, column=3)
    
    # Optional filtering of chosen spectra by filenumber or string
    a=tk.Label(root, text='Optional filtering by filenumber or filename').grid(row=7, column=0)
    a=tk.Label(root, text='Enter filenumbers for plot report: default all').grid(row=8, column=0)
    b=tk.Entry(root, textvariable=filestr).grid(row=8, column=1)
    # Combine averaged files only (checkbox)?
   
    # option to reselect labeled elemental peaks 

    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=9, column=2)

    d=tk.Button(root, text='Plot')
    d.bind('<Button-1>', plot)
    d.grid(row=9, column=1)

    root.mainloop()
        
    mychoice=choice.get()
    
    if mychoice=='plot':
        # TODO optional filtering of spelist according to tk input
        
        # call either reportcountsback, reportSD, reportderivcnt, or reportpeaksall
        myPDFname=PDFnamestr.get()
        kwargs={} # reconstruct with new choices
        # handle plot x range choice (used for all but peaks report)
        plotrange=[]
        if rangechoice.get()=='elems':
            # Unpack elements in case of change
            try:
                Elements=elemstr.get().split(',')
            except:
                pass
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
        # handle specifics depending on plot type choice
        if plotcol.get()=='Counts': # call reportcountsback
            if backfitbool.get():
                kwargs.update({'addbackfit':True}) # plotting of background
            if backptbool.get():
                kwargs.update({'backfitpts':Backfitlog}) # plot pts used to fit background 
            if plotelemsbool.get(): # optional labeling of elements
                kwargs.update({'plotelems':Elements}) # pass elements list 
            reportcountsback(spelist, plotrange, AESquantparams, Backfitlog, PDFname=myPDFname, **kwargs)
        elif plottype.get()=='Deriv':
            if smdifbool.get():
                kwargs.update({'smdifpeaks':True})
            if plotelemsbool.get(): # optional labeling of elements
                kwargs.update({'plotelems':Elements}) # pass elements list 
            reportSD(spelist, plotrange, Smdifpeakslog, AESquantparams, PDFname='SDplot_report', **kwargs)
        elif plottype.get()=='Both':
            if backfitbool.get():
                kwargs.update({'addbackfit':True}) # plotting of background
            if backptbool.get():
                kwargs.update({'backfitpts':Backfitlog}) # plot pts used to fit background 
            if smdifbool.get():
                kwargs.update({'smdifpeaks':True})
            if plotelemsbool.get(): # optional labeling of elements
                kwargs.update({'plotelems':Elements}) # pass elements list 
            reportSD(spelist, plotrange, Smdifpeakslog, AESquantparams, PDFname='SDplot_report', **kwargs)


        elif plottype.get()=='Both':
        if filestr.get()!='': # also put into kwargs to initialize next run
            kwargs.update({'fileset':filestr.get()})
        # get and pass PDF name string as arg (default name is prefilled)
        AESplot1(filestr.get(), spelist, AESquantparams, **kwargs)

    return kwargs


