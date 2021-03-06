# Auger legacy functions


'''TESTING 
i=0   X=0   Y=0  elnum=0
[elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, 
                idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth]=Elemdata[i]
        
'''
def integmaps(specimage, backarray, Elemdata, shiftdict):
    ''' Return integrated counts over each specified peak 
    Charging and any charge compensation 
    specimage - raw spectral image data 
    backarray is slope-intercept of background fits for each pixel 
    Elemdata contains all data regions, integ widths, etc.
    shiftdict hold pixel and peak specific peak shift 
    (calculated using 3x3 smoothed spectral image)
    peak out of range problem handled by double counting same region on other side of peak
    returns: dictionary w/ elem as key and np array containing [0] - total raw counts
        [1] background fit counts [2] subtracted counts integrated over element-specific
        integration width (usually 9eV wide)
    REPLACED BY FINDALLPEAKS METHOD
    '''
    sumdict={}
    for elnum, [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, 
                idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth] in enumerate(Elemdata):
        shiftimage=shiftdict.get(elem,[]) # gets peak centers for this
        xrange=np.arange(peakrange[0], peakrange[1]+1, 1.0) # same xrange for all
        # np array holding integrated counts
        sumarray=np.empty([specimage.shape[0],specimage.shape[1],3])
        for X in range(0,specimage.shape[0]): # iterate over rows
            for Y in range(0,specimage.shape[1]): # iterate over cols
                # Get linear background over elem peak for given pixel
                # Backarray contains slope, intercepts stacked for all elements in elem order
                slope=backarray[X,Y,2*elnum]
                intercept=backarray[X,Y,2*elnum+1]
                backvals=np.zeros(xrange.shape)
                for j in range(0, xrange.shape[0]):
                    backvals[j]=slope*xrange[j]+intercept
                # select data subset using indices of this peak
                rawdata=specimage[X,Y,peakind[0]:peakind[1]+1]
                # Calculated background subtracted data in this region 
                subdata=specimage[X,Y,peakind[0]:peakind[1]+1]-backvals

                thiscent=int((subdata.shape[0]-1)/2)+int(shiftimage[X,Y])
                # Looking for peak integration limits
                low=thiscent-int((integwidth-1)/2)
                high=thiscent+int((integwidth-1)/2)
                # Mirror integration regions if out-of-range (assumes roughly symmetric peak)
                if low<0 | high>subdata.shape[0]:
                    if low<0:
                        # number of missing data channels
                        numaddchan=-low
                        low=0 # reset lower integ limit
                        # compensate w/ double count of high-side counts
                        addintegrange=[high-numaddchan+1,high+1]
                    else:
                        # number of missing data channels
                        numaddchan=high-subdata.shape[0]+1
                        high=subdata.shape[0] # reset upper integ limit
                        # Double count points on opposite side of peak
                        addintegrange=[low,low+numaddchan]
                else:
                    # full peak available in range
                    addintegrange=[]
                # Perform integrations (upper limit not included)
                rawcounts=np.sum(rawdata[low:high])
                backcounts=np.sum(backvals[low:high])
                subcounts=np.sum(subdata[low:high])
                if addintegrange:
                    rawcounts+=np.sum(rawdata[addintegrange[0]:addintegrange[1]])
                    backcounts+=np.sum(backvals[addintegrange[0]:addintegrange[1]])
                    subcounts+=np.sum(subdata[addintegrange[0]:addintegrange[1]])
                    if (rawdata[low:high].shape[0] + rawdata[addintegrange[0]:addintegrange[1]].shape[0]) != integwidth:
                        print('Number of integrated elements not equal to peak integwidth')
                sumarray[X,Y,0]=rawcounts
                sumarray[X,Y,1]=backcounts
                sumarray[X,Y,2]=subcounts
        sumdict.update({elem:sumarray})
    return sumdict

def calcshifts(specimage, Elemdata, AESquantparams):
    ''' Diagnostic to look at peak shifts across spatial array; return stats on all peaks 
    over all pixels in spectral image; looking at direct peak not deriv
    shifts for given element and shift image visualization of peak shifting magnitude
    return:
        numpy arrays for each element in dictionary w/ elem as key (can show position dependent
        charging
        peakstats - dict with each elements avreage, stdev, min max shift
         '''       
    #TODO modify this to handle charging shift determinations?? 
    peakstats={}
    shiftdict={}
    print('Peak shift averages by element')
    print('Element Mean Stdev Min Max')
    for i, [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, 
            idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth] in enumerate(Elemdata):
        [xmin, xmax]=peakind # unpack index range of elem/peak
        # Calculate 2D shift array holding
        shiftarr=np.empty([specimage.shape[0],specimage.shape[1]])
        for x in range(0,specimage.shape[0]):
            for y in range(0,specimage.shape[1]):
                # Index # of true peak
                if str(idealind)=='nan': # for charging samples ideal index can be out of scan range
                    thisshift=peakrange[0]+(peakrange[1]-peakrange[0])*np.argmax(specimage[x,y,xmin:xmax])/(xmax-xmin)-idealev
                else:
                    thisshift=np.argmax(specimage[x,y,xmin:xmax])+xmin-idealind
                shiftarr[x,y]=thisshift
        shiftdict.update({elem:shiftarr})
        # Grab and return mean, stdev, min, max
        peakstats.update({elem: ["%.3f" % shiftarr.mean(), "%.3f" % shiftarr.std(), int(np.min(shiftarr)), int(np.max(shiftarr))]})
        print(elem,'   ', "%.3f" % shiftarr.mean(), "%.3f" % shiftarr.std(),' ', int(np.min(shiftarr)),' ', int(np.max(shiftarr)))
    return shiftdict, peakstats

'''
Testing:
Augerslice.plot.scatter(x='Energy', y='Counts')
'''
def calcbackgrounds(specimage, energy, Elemdata):
    ''' Linear fit of e- background under peak using low and high background regions
    return as numpy slope[0], intercept[1],
    can run on raw spectral image or perform adjacent averaging/filtering first

    '''
    # make 3D blank numpy array of correct X, Y pixel dimension 
    # Z contains slope & intercept of fit for each element from Elemdata
    backarray=np.empty([specimage.shape[0],specimage.shape[1], 2*len(Elemdata)])
    Augerfile=pd.DataFrame()
    # index through the entier 2D image (X and Y pixel positions)
    for X in range(0,specimage.shape[0]):
        for Y in range(0,specimage.shape[1]):
            Augerfile['Energy']=energy
            Augerfile['Counts']=specimage[X,Y]
            # now loop through and fit each element
            for i, [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, 
                    idealind, integwidth, idealnegpeak, chargeshift, peakwidth, searchwidth] in enumerate(Elemdata):
                backrange=[]
                [lowmin, lowmax]=lowind # unpack index range of elem/peak
                [highmin, highmax]=highind
                backrange=[int(i) for i in range(lowmin, lowmax+1)]
                backrange.extend([int(i) for i in range(highmin, highmax+1)])
                Augerslice=Augerfile[Augerfile.index.isin(backrange)]
                data1=Augerslice['Energy']
                data2=Augerslice['Counts']
                slope,intercept=np.polyfit(data1, data2, 1)
                backarray[X,Y,2*i]=slope
                backarray[X,Y,2*i+1]=intercept
    return backarray

	
# old tk interfaces (since combined to single multipurpose reporting function)
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

# old version of scattercompplot before outlier determination/plotting

def plotcntsmajor(Params, areanum):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    calls findevbreaks
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    #numareas=int(Params.Areas)
    #myareas=parseareas(areas, numareas) # set of areas for plotting
    
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max())
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
    colname='Counts'+str(areanum)
    backname='Backfit'+str(areanum)
    
    # find multiplex evbreaks 
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float

    # S region
    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
    if Augerslice.empty==False: # only plottable if data exists in this range          
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region        
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,0])     
        axes[0,0].axvline(x=150.3, color='b') # at ideal S position
               
        for j, val in enumerate(energyvals): # plot multiplex ev breaks in red
            if val > 115 and val < 200: 
                axes[0,0].axvline(x=val, color='r') # on counts plot
                
    # C/Ca region
    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
    if Augerslice.empty==False: # only plottable if data exists in this range       
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,0])         
        axes[1,0].axvline(x=271, color='b') # at ideal C position
        axes[1,0].axvline(x=291, color='b') # at ideal Ca position
        
        # add red vert line at multiplex energy break if present
        for j, val in enumerate(energyvals):
            if val > 225 and val < 320: 
                axes[1,0].axvline(x=val, color='r') # on counts plot
            
    # O regions
    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]     
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region        
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,1])
        axes[0,1].axvline(x=508, color='b') # at ideal O position
        
        for j, val in enumerate(energyvals):
            if val > 475 and val < 535: 
                axes[0,1].axvline(x=val, color='r') # on counts plot
    # Fe region
    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region         
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,1])       
        axes[1,1].axvline(x=595, color='b') # at ideal Fe1 position
        axes[1,1].axvline(x=647.2, color='b') # at ideal Fe2 position
        axes[1,1].axvline(x=702.2, color='b') # at ideal Fe3 position
       
        for j, val in enumerate(energyvals):
            if val > 560 and val < 750: 
                axes[1,1].axvline(x=val, color='r') # on counts plot
    # Mg region
    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region 
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,2])
        axes[0,2].axvline(x=1181.7, color='b') # at ideal Mg position
        
        for j, val in enumerate(energyvals):
            if val > 1140 and val < 1220: 
                axes[0,2].axvline(x=val, color='r') # on counts plot
                 
    # Si region
    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
    if Augerslice.empty==False: # only plottable if data exists in this range 
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,2]) # Si2 region                          
        axes[1,2].axvline(x=1609, color='b') # at ideal Si position  
        
        for j, val in enumerate(energyvals):
            if val > 1550 and val < 1650: 
                axes[1,2].axvline(x=val, color='r') # on counts plot
    return

	# Legacy versions replaced by above generalized functions
def reportpeaksmajor(paramlog, addgauss=True):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass list of files and selected background regions from automated fitting
    '''
    with PdfPages('Peak_report.pdf') as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.iloc[index]['Filename']
            numareas=int(paramlog.iloc[index]['Areas'])            
            Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.iloc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 # 
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
                peakname='Peaks'+str(areanum)
                gaussname='Gauss'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                # S region
                if myplotrange[0] < 135 and  myplotrange[1] > 165:        
                    Augerslice=Augerfile[(Augerfile['Energy']>135) & (Augerfile['Energy']<165)]
                    axes[0,0].axvline(x=150.3, color='b') # at ideal direct peak S position 
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[0,0]) # S region
                        if addgauss==True and gaussname in Augerslice.dtypes.index: # ensure Gaussian fit col exists
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[0,0]) # S region
                # C/Ca region
                if myplotrange[0] < 275 and  myplotrange[1] > 305:         
                    Augerslice=Augerfile[(Augerfile['Energy']>275) & (Augerfile['Energy']<305)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[1,0]) # C/Ca region               
                        if addgauss==True and gaussname in Augerslice.dtypes.index: 
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[1,0]) # C/Ca region                    
                    axes[1,0].axvline(x=291, color='b') # at ideal direct peak Ca position
                # Fe region
                if myplotrange[0] < 625 and  myplotrange[1] > 725:
                    Augerslice=Augerfile[(Augerfile['Energy']>625) & (Augerfile['Energy']<725)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[1,1]) # Fe region
                        if addgauss==True and gaussname in Augerslice.dtypes.index: 
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[1,1])
                    axes[1,1].axvline(x=648, color='b') # at ideal direct peak Fe2 position
                    axes[1,1].axvline(x=702, color='b') # at ideal direct peak Fe3 position
                # Mg region
                if myplotrange[0] < 1165 and  myplotrange[1] > 1200:
                    Augerslice=Augerfile[(Augerfile['Energy']>1165) & (Augerfile['Energy']<1200)]    
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[0,2]) # Mg region
                        if addgauss==True and gaussname in Augerslice.dtypes.index: 
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[0,2])
                    axes[0,2].axvline(x=1185, color='b') # at ideal Mg position
                # Si region
                if myplotrange[0] < 1590 and  myplotrange[1] > 1625:   
                    Augerslice=Augerfile[(Augerfile['Energy']>1590) & (Augerfile['Energy']<1625)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[1,2]) # Si2 region  
                        if addgauss==True and gaussname in Augerslice.dtypes.index: 
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[1,2])
                    axes[1,2].axvline(x=1614, color='b') # at ideal Si position  
                pdf.savefig(fig)
                plt.close('all') # close all open figures
    return
    
def reportSDmajor(paramlog, Smdifpeakslog, PDFname='SDplots_report.pdf'):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    '''
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])            
            Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            mypeaks=Smdifpeakslog[(Smdifpeakslog['Filename']==AugerFileName)] # retrieve assoc. subset of peaks data
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 # 
                Peaks=mypeaks[(mypeaks['Areanumber']==areanum)] # then only this areanum
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
                colname='S7D7'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                # S region
                Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]   
                if not Augerslice.empty: # skip entire plot if no data
                    Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>115) & (Peaks['Peakenergy']<200)]
                    if not plotpts.empty: 
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,0], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,0], color='r')                                        
                        titlestring=maketitlestring(plotpts)
                        axes[0,0].set_title(titlestring, fontsize=10)
                        axes[0,0].axvline(x=154, color='b') # at ideal S position
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[0,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[0,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
                    
                    # add red vert line at multiplex energy break if present                    
                    for j, val in enumerate(energyvals):
                        if val > 115 and val < 200: 
                            axes[0,0].axvline(x=val, color='r') # on counts plot
                # C/Ca region
                Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]   
                if not Augerslice.empty: 
                    Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
                    axes[1,0].axvline(x=276, color='b') # at ideal C position
                    axes[1,0].axvline(x=296, color='b') # at ideal Ca position
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>225) & (Peaks['Peakenergy']<320)]
                    if not plotpts.empty:
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,0], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,0], color='r')                                                                                       
                        titlestring=maketitlestring(plotpts)
                        axes[1,0].set_title(titlestring, fontsize=10)        
                       
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[1,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[1,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

                    # add red vert line at multiplex energy break if present
                    for j, val in enumerate(energyvals): 
                        if val > 225 and val < 320: 
                            axes[1,0].axvline(x=val, color='r') # on counts plot
                        
                # O regions
                Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]
                if not Augerslice.empty:
                    Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region
                    axes[0,1].axvline(x=513, color='b') # at ideal O position
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>475) & (Peaks['Peakenergy']<535)]
                    if not plotpts.empty:
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,1], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,1], color='r')                                                              
                        titlestring=maketitlestring(plotpts)
                        axes[0,1].set_title(titlestring, fontsize=10)                
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[0,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[0,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
                    
                    for j, val in enumerate(energyvals):
                        if val > 475 and val < 535: 
                            axes[0,1].axvline(x=val, color='r') # on counts plot
                # Fe region
                Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
                if not Augerslice.empty:
                    Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region   
                    axes[1,1].axvline(x=600, color='b') # at ideal Fe1 position
                    axes[1,1].axvline(x=654, color='b') # at ideal Fe2 position
                    axes[1,1].axvline(x=707, color='b') # at ideal Fe3 position
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>560) & (Peaks['Peakenergy']<750)]
                    if not plotpts.empty:
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,1], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,1], color='r') 
                        titlestring=maketitlestring(plotpts)
                        axes[1,1].set_title(titlestring, fontsize=10)                
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[1,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[1,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

                    for j, val in enumerate(energyvals):
                        if val > 560 and val < 750: 
                            axes[1,1].axvline(x=val, color='r') # on counts plot
                # Mg region
                Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]    
                if not Augerslice.empty:                    
                    Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region
                    axes[0,2].axvline(x=1185, color='b') # at ideal Mg position
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>1140) & (Peaks['Peakenergy']<1220)]
                    if not plotpts.empty: 
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,2], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,2], color='r')                     
                        titlestring=maketitlestring(plotpts)
                        axes[0,2].set_title(titlestring, fontsize=10)        
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[0,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[0,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

                    for j, val in enumerate(energyvals):
                        if val > 1140 and val < 1220: 
                            axes[0,2].axvline(x=val, color='r') # on counts plot
                # Si region
                Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
                if not Augerslice.empty:                    
                    Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region         
                    axes[1,2].axvline(x=1614, color='b') # at ideal Si position 
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>1550) & (Peaks['Peakenergy']<1650)]
                    if not plotpts.empty: 
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,2], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,2], color='r') 
                        titlestring=maketitlestring(plotpts)
                        axes[1,2].set_title(titlestring, fontsize=10)
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[1,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[1,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
                    for j, val in enumerate(energyvals):
                        if val > 1550 and val < 1650: 
                            axes[1,2].axvline(x=val, color='r') # on counts plot
                    pdf.savefig(fig)
                    plt.close(fig) # closes recently displayed figure (since it's saved in PDF report)
    return

def plotcntsmajor(Params, areanum):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    calls findevbreaks
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    #numareas=int(Params.Areas)
    #myareas=parseareas(areas, numareas) # set of areas for plotting
    
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max())
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
    colname='Counts'+str(areanum)
    backname='Backfit'+str(areanum)
    
    # find multiplex evbreaks 
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float

    # S region
    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
    if Augerslice.empty==False: # only plottable if data exists in this range          
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region        
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,0])     
        axes[0,0].axvline(x=150.3, color='b') # at ideal S position
               
        for j, val in enumerate(energyvals): # plot multiplex ev breaks in red
            if val > 115 and val < 200: 
                axes[0,0].axvline(x=val, color='r') # on counts plot
                
    # C/Ca region
    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
    if Augerslice.empty==False: # only plottable if data exists in this range       
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,0])         
        axes[1,0].axvline(x=271, color='b') # at ideal C position
        axes[1,0].axvline(x=291, color='b') # at ideal Ca position
        
        # add red vert line at multiplex energy break if present
        for j, val in enumerate(energyvals):
            if val > 225 and val < 320: 
                axes[1,0].axvline(x=val, color='r') # on counts plot
            
    # O regions
    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]     
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region        
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,1])
        axes[0,1].axvline(x=508, color='b') # at ideal O position
        
        for j, val in enumerate(energyvals):
            if val > 475 and val < 535: 
                axes[0,1].axvline(x=val, color='r') # on counts plot
    # Fe region
    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region         
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,1])       
        axes[1,1].axvline(x=595, color='b') # at ideal Fe1 position
        axes[1,1].axvline(x=647.2, color='b') # at ideal Fe2 position
        axes[1,1].axvline(x=702.2, color='b') # at ideal Fe3 position
       
        for j, val in enumerate(energyvals):
            if val > 560 and val < 750: 
                axes[1,1].axvline(x=val, color='r') # on counts plot
    # Mg region
    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region 
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,2])
        axes[0,2].axvline(x=1181.7, color='b') # at ideal Mg position
        
        for j, val in enumerate(energyvals):
            if val > 1140 and val < 1220: 
                axes[0,2].axvline(x=val, color='r') # on counts plot
                 
    # Si region
    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
    if Augerslice.empty==False: # only plottable if data exists in this range 
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,2]) # Si2 region                          
        axes[1,2].axvline(x=1609, color='b') # at ideal Si position  
        
        for j, val in enumerate(energyvals):
            if val > 1550 and val < 1650: 
                axes[1,2].axvline(x=val, color='r') # on counts plot
    return
  
def plotSDmajor(Params, Peaks, areanum):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    #numareas=int(Params.Areas)
    #myareas=parseareas(areas, numareas) # set of areas for plotting
    
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max())
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
    colname='S7D7'+str(areanum)
    
    # find multiplex evbreaks 
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float

    # S region
    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
    if Augerslice.empty==False: # only plottable if data exists in this range          
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region        
        axes[0,0].axvline(x=154, color='b') # at ideal S position
        
        plotpts=Peaks[(Peaks['Peakenergy']>115) & (Peaks['Peakenergy']<200)] 
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,0], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,0], color='r')       
            titlestring=maketitlestring(plotpts)
            axes[0,0].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals): # plot multiplex ev breaks in red
            if val > 115 and val < 200: 
                axes[0,0].axvline(x=val, color='r') # on counts plot
                
    # C/Ca region
    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
    if Augerslice.empty==False: # only plottable if data exists in this range       
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
        axes[1,0].axvline(x=276, color='b') # at ideal C position
        axes[1,0].axvline(x=296, color='b') # at ideal Ca position
        
        plotpts=Peaks[(Peaks['Peakenergy']>225) & (Peaks['Peakenergy']<320)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,0], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,0], color='r')
            titlestring=maketitlestring(plotpts)
            axes[1,0].set_title(titlestring, fontsize=10)        
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
        # add red vert line at multiplex energy break if present
        for j, val in enumerate(energyvals):
            if val > 225 and val < 320: 
                axes[1,0].axvline(x=val, color='r') # on counts plot
            
    # O regions
    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]     
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region        
        axes[0,1].axvline(x=513, color='b') # at ideal O position
        
        plotpts=Peaks[(Peaks['Peakenergy']>475) & (Peaks['Peakenergy']<535)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,1], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,1], color='r')                    
            titlestring=maketitlestring(plotpts)
            axes[0,1].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals):
            if val > 475 and val < 535: 
                axes[0,1].axvline(x=val, color='r') # on counts plot
    # Fe region
    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region         
        axes[1,1].axvline(x=600, color='b') # at ideal Fe1 position
        axes[1,1].axvline(x=654, color='b') # at ideal Fe2 position
        axes[1,1].axvline(x=707, color='b') # at ideal Fe3 position

        plotpts=Peaks[(Peaks['Peakenergy']>560) & (Peaks['Peakenergy']<750)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,1], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,1], color='r')                   
            titlestring=maketitlestring(plotpts)
            axes[1,1].set_title(titlestring, fontsize=10)                
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
        
        for j, val in enumerate(energyvals):
            if val > 560 and val < 750: 
                axes[1,1].axvline(x=val, color='r') # on counts plot
    # Mg region
    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region        
        axes[0,2].axvline(x=1185, color='b') # at ideal Mg position
        
        plotpts=Peaks[(Peaks['Peakenergy']>1140) & (Peaks['Peakenergy']<1220)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,2], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,2], color='r') 
            titlestring=maketitlestring(plotpts)
            axes[0,2].set_title(titlestring, fontsize=10)               
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
 
        for j, val in enumerate(energyvals):
            if val > 1140 and val < 1220: 
                axes[0,2].axvline(x=val, color='r') # on counts plot
                 
    # Si region
    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
    if Augerslice.empty==False: # only plottable if data exists in this range 
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region                          
        axes[1,2].axvline(x=1614, color='b') # at ideal Si position  
        
        plotpts=Peaks[(Peaks['Peakenergy']>1550) & (Peaks['Peakenergy']<1650)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,2], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,2], color='r')         
            titlestring=maketitlestring(plotpts)
            axes[1,2].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals):
            if val > 1550 and val < 1650: 
                axes[1,2].axvline(x=val, color='r') # on counts plot
    return

def reportcountsmajor(paramlog, Smdifpeakslog, PDFname='countsplot_report.pdf'):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    REDUNDANT... just use reportcountsback with background switch
    '''
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])            
            Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 # 
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
                colname='Counts'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                # S region
                if myplotrange[0] < 115 and  myplotrange[1] > 200:        
                    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region
                    # TODO look up these constants from AESquantparams
                    axes[0,0].axvline(x=154, color='r') # at ideal S position
                    for j, val in enumerate(energyvals):
                        if val > 115 and val < 200: 
                            axes[0,0].axvline(x=val, color='r') # on counts plot
                # C/Ca region
                if myplotrange[0] < 225 and  myplotrange[1] > 320:         
                    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region               
                    axes[1,0].axvline(x=276, color='r') # at ideal C position
                    axes[1,0].axvline(x=296, color='r') # at ideal Ca position
                    # add red vert line at multiplex energy break if present
                    for j, val in enumerate(energyvals):
                        if val > 225 and val < 320: 
                            axes[1,0].axvline(x=val, color='r') # on counts plot
                        
                # O regions
                if myplotrange[0] < 475 and  myplotrange[1] > 535:
                    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]    
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region
                    axes[0,1].axvline(x=513, color='r') # at ideal O position
                    for j, val in enumerate(energyvals):
                        if val > 475 and val < 535: 
                            axes[0,1].axvline(x=val, color='r') # on counts plot
                # Fe region
                if myplotrange[0] < 560 and  myplotrange[1] > 750:
                    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region                              
                    axes[1,1].axvline(x=600, color='r') # at ideal Fe1 position
                    axes[1,1].axvline(x=654, color='r') # at ideal Fe2 position
                    axes[1,1].axvline(x=707, color='r') # at ideal Fe3 position
                    for j, val in enumerate(energyvals):
                        if val > 560 and val < 750: 
                            axes[1,1].axvline(x=val, color='r') # on counts plot
                # Mg region
                if myplotrange[0] < 1140 and  myplotrange[1] > 1220:
                    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]    
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region
                    axes[0,2].axvline(x=1185, color='r') # at ideal Mg position
                    for j, val in enumerate(energyvals):
                        if val > 1140 and val < 1220: 
                            axes[0,2].axvline(x=val, color='r') # on counts plot
                # Si region
                if myplotrange[0] < 1550 and  myplotrange[1] > 1650:   
                    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region                  
                    axes[1,2].axvline(x=1614, color='b') # at ideal Si position  
                    for j, val in enumerate(energyvals):
                        if val > 1550 and val < 1650: 
                            axes[1,2].axvline(x=val, color='r') # on counts plot
                pdf.savefig(fig)
                plt.close('all') # close all open figures
       
def plotSDmajor(Params, Peaks, areanum):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    #numareas=int(Params.Areas)
    #myareas=parseareas(areas, numareas) # set of areas for plotting
    
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max())
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
    colname='S7D7'+str(areanum)
    
    # find multiplex evbreaks 
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float

    # S region
    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
    if Augerslice.empty==False: # only plottable if data exists in this range          
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region        
        axes[0,0].axvline(x=154, color='b') # at ideal S position
        
        plotpts=Peaks[(Peaks['Peakenergy']>115) & (Peaks['Peakenergy']<200)] 
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,0], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,0], color='r')       
            titlestring=maketitlestring(plotpts)
            axes[0,0].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals): # plot multiplex ev breaks in red
            if val > 115 and val < 200: 
                axes[0,0].axvline(x=val, color='r') # on counts plot
                
    # C/Ca region
    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
    if Augerslice.empty==False: # only plottable if data exists in this range       
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
        axes[1,0].axvline(x=276, color='b') # at ideal C position
        axes[1,0].axvline(x=296, color='b') # at ideal Ca position
        
        plotpts=Peaks[(Peaks['Peakenergy']>225) & (Peaks['Peakenergy']<320)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,0], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,0], color='r')
            titlestring=maketitlestring(plotpts)
            axes[1,0].set_title(titlestring, fontsize=10)        
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
        # add red vert line at multiplex energy break if present
        for j, val in enumerate(energyvals):
            if val > 225 and val < 320: 
                axes[1,0].axvline(x=val, color='r') # on counts plot
            
    # O regions
    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]     
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region        
        axes[0,1].axvline(x=513, color='b') # at ideal O position
        
        plotpts=Peaks[(Peaks['Peakenergy']>475) & (Peaks['Peakenergy']<535)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,1], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,1], color='r')                    
            titlestring=maketitlestring(plotpts)
            axes[0,1].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals):
            if val > 475 and val < 535: 
                axes[0,1].axvline(x=val, color='r') # on counts plot
    # Fe region
    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region         
        axes[1,1].axvline(x=600, color='b') # at ideal Fe1 position
        axes[1,1].axvline(x=654, color='b') # at ideal Fe2 position
        axes[1,1].axvline(x=707, color='b') # at ideal Fe3 position

        plotpts=Peaks[(Peaks['Peakenergy']>560) & (Peaks['Peakenergy']<750)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,1], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,1], color='r')                   
            titlestring=maketitlestring(plotpts)
            axes[1,1].set_title(titlestring, fontsize=10)                
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
        
        for j, val in enumerate(energyvals):
            if val > 560 and val < 750: 
                axes[1,1].axvline(x=val, color='r') # on counts plot
    # Mg region
    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region        
        axes[0,2].axvline(x=1185, color='b') # at ideal Mg position
        
        plotpts=Peaks[(Peaks['Peakenergy']>1140) & (Peaks['Peakenergy']<1220)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,2], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,2], color='r') 
            titlestring=maketitlestring(plotpts)
            axes[0,2].set_title(titlestring, fontsize=10)               
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
 
        for j, val in enumerate(energyvals):
            if val > 1140 and val < 1220: 
                axes[0,2].axvline(x=val, color='r') # on counts plot
                 
    # Si region
    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
    if Augerslice.empty==False: # only plottable if data exists in this range 
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region                          
        axes[1,2].axvline(x=1614, color='b') # at ideal Si position  
        
        plotpts=Peaks[(Peaks['Peakenergy']>1550) & (Peaks['Peakenergy']<1650)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,2], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,2], color='r')         
            titlestring=maketitlestring(plotpts)
            axes[1,2].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals):
            if val > 1550 and val < 1650: 
                axes[1,2].axvline(x=val, color='r') # on counts plot
    return

def fitcubic(df, areanum, elem, AugerFileName):
    '''Pass appropriate chunk from Auger spectral dataframe, perform cubic fit
    return chunk with backfit column added '''
    colname='Counts'+str(areanum)
    backfitname='Backfit'+str(areanum)
    xcol=df['Energy']
    ycol=df[colname] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
		cubicfunc=lambda x, a, b, c, d: a*x**3 + b*x**2 + c*x + d # lambda definition of cubic poly
        A,B,C, D=np.polyfit(xcol, ycol, 3)
    except: # deal with common problems with linregress
        print('Fitting error for', elem, ' in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a','n/a') # return all n/a
        return df, fitparams
    fitparams=(A, B, C, D) # tuple to return coeffs of 2nd order poly fit
    for index,row in df.iterrows():
        xval=df.loc[index]['Energy']
        yval= A * xval**3+ B * xval**2 + C * xval + D
        df=df.set_value(index, backfitname, yval)
    return df, fitparams

def fitCapeak(df, areanum, elem, AugerFileName):
    '''Pass appropriate chunk from Auger spectral dataframe, perform linear fit
    return chunk with backfit column added 
	replaced with curve_fit version'''
    colname='Counts'+str(areanum)
    backfitname='Backfit'+str(areanum)
    xcol=df['Energy']
    ycol=df[colname] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
        A,B,C=np.polyfit(xcol, ycol, 2)
        # diagonal of covariance matrix contains variances for fit params
    except: # deal with common problems with linregress
        print('Fitting error for', elem, ' in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a') # return all n/a
        return df, fitparams
    fitparams=(A, B, C) # tuple to return coeffs of 2nd order poly fit
    for index,row in df.iterrows():
        xval=df.loc[index]['Energy']
        yval= A * xval**2+ B * xval + C
        df=df.set_value(index, backfitname, yval)
    return df, fitparams
	
def scattercompplot(comp1, comp2, elemlist, joinlist=['Sample','Areanum'], basis=False):
    '''Pass two versions of composition calculation (using different lines or whatever) and compare 
    major elements using scatter graphs .. single point for each sample
    uses inner merge to select only subset with values from each df
	use either sample or filenumber'''
    
    elemlist=[re.match('\D+',i).group(0) for i in elemlist] 
    # strip number from peaks like Fe2 if present; columns will be element names (Fe) not peak names (Fe2)
    if basis==False: # use atomic % (which is the default), not basis for each element
        elemlist=['%'+s for s in elemlist]
    numareas=len(elemlist)
    
    # set nrows and ncols for figure of proper size
    cols=divmod(numareas,2)[0]+ divmod(numareas,2)[1]
    if numareas>1:
        rows=2
    else:
        rows=1        
    fig, axes = plt.subplots(nrows=rows, ncols=cols) # axes is array
    # merge dfs with comp1 and comp2 using inner join    
    df=pd.merge(comp1, comp2, how='inner', on=joinlist, suffixes=('','b'))
    mycols=df.dtypes.index # same columns
    outliers=pd.DataFrame(columns=mycols) # empty dataframe for outlying points
    outliers['Element'] # Find outliers from linear fit for each listed element
    for i,elem in enumerate(elemlist):
         # determine which subplot to use
        if (i+1)%2==1:
            rownum=0
        else:
            rownum=1
        colnum=int((i+1)/2.1)       
        xcol=elem
        ycol=elem+'b' # same element from second dataset
        if numareas==1: # deal with single subplot separately
            df.plot.scatter(x=xcol, y=ycol, ax=axes) # single plot axes has no [#,#]
        else:
            df.plot.scatter(x=xcol, y=ycol, ax=axes[rownum,colnum])
            # linear regression: fitting, plot and add labels
        data1=df[elem] # just this data column
        colname=elem+'b'
        data2=df[colname]
        # slope,intercept=np.polyfit(data1, data2, 1)  numpy version
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(data1, data2) # imported from scipy.stats
        # set x range for linear plot 
        text1=str(round(slope,2))+' *x +' + str(round(intercept,2))
        text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,3))
        xmax=max(max(data1),max(data2))*1.1 # set to slightly larger than max of dataset
        x=np.linspace(0,xmax,100) # setting range for 
        if numareas==1: # deal with single subplot separately
            axes.text(0.025,0.9, text1, fontsize=12, transform=axes.transAxes)
            axes.text(0.025,0.8, text2, fontsize=12, transform=axes.transAxes)
            plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        else: # typical multiarea plot
            axes[rownum,colnum].text(0.025,0.9, text1, fontsize=12, transform=axes[rownum,colnum].transAxes)
            axes[rownum,colnum].text(0.025,0.8, text2, fontsize=12, transform=axes[rownum,colnum].transAxes)
            plt.axes(axes[rownum,colnum]) # set correct axes as active
            plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        outliers=returnoutliers(df,slope, intercept, elem, basis=False)
    return df


# old version when entire datachunk including evbreaks was passed
def savgol(counts, evbreaks):
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

# old (not generalized) version of counts report iwth backgrounds 
def reportcountsback(paramlog, plotelems, AESquantparams, plotback=True, PDFname='countsback_report.pdf'):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass list of files and selected background regions from automated fitting
    background fits themselves stored with auger csv files
    optional pass of backfitlog (w/ points defining region boundary for background fitting useful for troubleshooting fits)
    evbreaks from multiplex plotted as red lines (often different dwell times for different elements)
    plotback switch -- whether or not to plot 
    '''
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])            
            Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber          
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 # 
                Peaks=mypeaks[(mypeaks['Area']==areanum)] # get subset of peaks for this area number 
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
                colname='Counts'+str(areanum)
                backname='Backfit'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                # S region
                if myplotrange[0] < 115 and  myplotrange[1] > 200:        
                    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
                    if not Augerslice.empty:      
                        # find elemental lines in this keV range
                        thisrange='115-200'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[0,0].axvline(x=elemtuple[0], color='b') # O line
                            axes[0,0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region                        
                        if bflog!=False: # plotting of Auger background fits from integral quant method
                            Augerslice.plot(x='Energy', y=backname, ax=axes[0,0])   
                            plotpts=Peaks[(Peaks['Element']=='S')]
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # fit boundaries for S
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[0,0])
                        for j, val in enumerate(energyvals): # any possible evbreaks in multiplex
                            if val > 115 and val < 200: 
                                axes[0,0].axvline(x=val, color='r') # on counts plot
                # C/Ca region
                if myplotrange[0] < 225 and  myplotrange[1] > 330:         
                    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<330)]
                    if not Augerslice.empty: 
                        # find elemental lines in this keV range
                        thisrange='225-330'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[1,0].axvline(x=elemtuple[0], color='b') # O line
                            axes[1,0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region
                        if bflog!=False:
                            Augerslice.plot(x='Energy', y=backname, ax=axes[1,0])
                            plotpts=Peaks[(Peaks['Element']=='Ca')]                        
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # fit boundaries for Ca
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[1,0])
                    # add red vert line at multiplex energy break if present
                    for j, val in enumerate(energyvals):
                        if val > 225 and val < 330: 
                            axes[1,0].axvline(x=val, color='r') # on counts plot
                        
                # O regions (Currently skip background fitting for O )
                if myplotrange[0] < 475 and  myplotrange[1] > 535:
                    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]    
                    if not Augerslice.empty: 
                        # find elemental lines in this keV range
                        thisrange='475-535'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[1,0].axvline(x=elemtuple[0], color='b') # O line
                            axes[1,0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region
                        for j, val in enumerate(energyvals):
                            if val > 475 and val < 535: 
                                axes[0,1].axvline(x=val, color='r') # on counts plot
                # Fe region
                if myplotrange[0] < 560 and  myplotrange[1] > 750:
                    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
                    if not Augerslice.empty:
                        # find elemental lines in this keV range
                        thisrange='560-750'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[1,1].axvline(x=elemtuple[0], color='b') # O line
                            axes[1,1].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region
                        if bflog!=False:
                            Augerslice.plot(x='Energy', y=backname, ax=axes[1,1])
                            plotpts=Peaks[Peaks['Element'].str.contains('Fe', na=False, case=False)] # match all Fe peaks
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # should get values from Fe2 and Fe3                
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[1,1])
                        for j, val in enumerate(energyvals):
                            if val > 560 and val < 750: 
                                axes[1,1].axvline(x=val, color='r') # on counts plot
                # Mg region
                if myplotrange[0] < 1140 and  myplotrange[1] > 1220:
                    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]    
                    if not Augerslice.empty:
                        # find elemental lines in this keV range
                        thisrange='1140-1220'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[0,2].axvline(x=elemtuple[0], color='b') # O line
                            axes[0,2].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region
                        if bflog!=False:
                            Augerslice.plot(x='Energy', y=backname, ax=axes[0,2])
                            plotpts=Peaks[(Peaks['Element']=='Mg')]
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # should get values from Fe2 and Fe3                
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[0,2])
                        for j, val in enumerate(energyvals):
                            if val > 1140 and val < 1220: 
                                axes[0,2].axvline(x=val, color='r') # on counts plot
                # Si region
                if myplotrange[0] < 1550 and  myplotrange[1] > 1650:   
                    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
                    if not Augerslice.empty:                        # find elemental lines in this keV range
                        thisrange='1550-1650'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[1,2].axvline(x=elemtuple[0], color='b') # O line
                            axes[1,2].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region
                        if bflog!=False:
                            Augerslice.plot(x='Energy', y=backname, ax=axes[1,2])
                            plotpts=Peaks[(Peaks['Element']=='Si')]
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # should get values from Fe2 and Fe3                
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[1,2])
                        for j, val in enumerate(energyvals):
                            if val > 1550 and val < 1650: 
                                axes[1,2].axvline(x=val, color='r') # on counts plot
                pdf.savefig(fig)
                plt.close('all') # close all open figures
    return

def setplotrange(plotrange, Augerfile):
    ''' Set range of plot based on element, numerical range or default to max range of spectrum  
    commonly called by Auger plot functions below
    LEGACY .. .replaced by getplotboundaries'''
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

# old version of countsbackreport before using tk interface for args/kwargs
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

# Old way with cloned filenumber areanumber rows in paramlog
def reportderivcnt(paramlog, plotelems, AESquantparams, **kwargs):
    ''' Comparison plots for both derivative and counts itself (don't have to worry about axes w/o indexing)
    plots selected filenumber + areanumber
    plots all spe files, associated quant points, and labels selected elemental lines
    Kwargs:   both backfitlog (key bfl) and smdiflogs (key smd) are optional kwarg params 
    pdf rename is also custom (defaults to derivcnts_report_3Mar17)
    
    '''
    plt.ioff()
    now=datetime.datetime.today()
    # Set up default PDF name (custom name is passable via kwargs)
    fname='Derivcnt_report_'+now.strftime('%d%b%y')+'.pdf'
    PDFname=kwargs.get('PDFname',fname)
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows():
            AugerFileName=paramlog.loc[index]['Filename']
            areanum=int(paramlog.loc[index]['Areanumber'])  
            Augerfile=openspefile(AugerFileName)
            if Augerfile.empty: # file not found problem
                continue
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            if 'smd' in kwargs:  
                Smdifdf=kwargs.get('smd','')
                # retrieve assoc. subset of peaks data
                thisfilepeaks=Smdifdf[(Smdifdf['Filename']==AugerFileName)&(Smdifdf['Areanumber']==areanum)] 
            # SKIP energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in multiplex
            # determine size of plot for this filenumber (same for all areas)
            # plotranges fnct sometimes combines certain plotelems (i.e. C and Ca together)
            plotranges=getplotboundaries(Augerfile, plotelems, AESquantparams) # returns plot ranges for all regions with data from plotelems
            # set plot rows and columns
            numrows=2
            numcols=len(plotranges)
            # for plotting background fit points used in integral method
            if 'bfl' in kwargs:  
                backfitdf=kwargs.get('bfl','')
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
            try:
                fig, axes = plt.subplots(nrows=numrows, ncols=numcols, figsize=(16,9), squeeze=False) # 2 by ? axes array
                mytitle=AugerFileName.replace('.csv','') +' area #'+str(areanum)                
                plt.suptitle(mytitle)
                for j, bounds in enumerate(plotranges):
                    [lower, upper]=bounds                    
                    thiscol=j
                    # Augerslice already checked for empty set 
                    Augerslice=Augerfile[(Augerfile['Energy']>lower) & (Augerfile['Energy']<upper)]
                    Augerslice.plot(x='Energy', y='S7D7'+str(areanum), ax=axes[0,thiscol]) # deriv in upper plot
                    Augerslice.plot(x='Energy', y='Counts'+str(areanum), ax=axes[1,thiscol]) # counts in lower 
                    # Section for labeling plotelements (plot range is passable as lower, upper tuple)
                    elemlines=getelemenergy(plotelems, bounds, AESquantparams, deriv=True) 
                    # list of tuples with energy,elemname
                    for k, elemtuple in enumerate(elemlines):
                        # elemtuple[0] is energy and [1] is element symbol
                        try:
                            axes[0,thiscol].axvline(x=elemtuple[0], color='b') # O line
                            axes[0,thiscol].text(elemtuple[0],-250, elemtuple[1],rotation=90, fontsize=18) # use standard -250 y val 
                        except:
                            print('Problem labeling elements')
                    # Section for adding smooth-diff quant data (optional via kwarg)
                    if kwargs.get('addsmdif',False):        
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
                    if 'bfl' in kwargs: # flag to show points from which background was determined 
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
        for index,row in paramlog.iterrows():
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])
            Augerfile=openspefile(AugerFileName)
            if Augerfile.empty:
                continue
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            

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
                    indexptslist=[int(i) for i in indexptslist]
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
