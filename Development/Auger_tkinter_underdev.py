# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 10:32:56 2017

@author: tkc
"""
import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.widgets import LassoSelector
from matplotlib import path

import datetime

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib as mpl 
import matplotlib.pyplot as plt

# still not sure how to 
root=tk.Tk()

figure=plt.figure()

figure = mpl.figure.Figure(figsize=(8,5), dpi=100)
ax = figure.add_subplot(221)
canvas = FigureCanvasTkAgg(figure, root)
canvas.get_tk_widget().grid(column=2, row=1, rowspan=5, sticky="nesw")
canvas.show()
root.mainloop()

# testing QM data structure
myQM=AESquantmap(os.getcwd())


# testing AESdataset file 
myAES=AESdataset()
# Make Augerfile from myAES Auger dataset
Augerfile 
'''TESTING
AESplot_gui(spelist, Elements, AESquantparams)
'''

# Testing QM data structures 
AEStest=AESQMdataset('C:\\Temp\\AugerQM\\14Oct17')
AEStest.aesquantparams
AEStest.spectralregs.Element.unique()
AEStest.spectralregs.Element[0]
AEStest.elemdata
AEStest.specimage.shape
len(AEStest.energy)


def AESplot_gui(spelist, Elements, AESquantparams):
    ''' tk interface for args/kwargs of single AES interactive plot
	filtering by filenumber, filename incorporated in tk interface 
    all args/dataframes must be passed through to plot functions 
    kwargs: fileset, areas, xrangestr -- usually just holds values entered during prior run 
    report will use all areas '''
    # first print out existing info in various lines
    root = tk.Tk()
    root.title('Interactive Auger plot interface')
    # Set up all the tk variables 
    filterstr=tk.StringVar() # comma separated or range of filenumbers for plot or string for sample name
    # filterstr.set(kwargs.get('fileset','')) # get entry from prior run
    areanumstr=tk.StringVar()  # optional choice of subset of area numbers
    # areanumstr.set(kwargs.get('areas',''))

    elemstr=tk.StringVar()
    mytext=', '.join(Elements) # elements for labelling 
    elemstr.set(mytext)
    chargeshift=tk.IntVar()
    chargeshift.set(0)
    plotelemsbool=tk.BooleanVar()  # optional labelling of elements
    plotelemsbool.set(True) # default to true
    choice=tk.StringVar()  # plot or abortw

    # Optional filtering of chosen spectra by filenumber or string
    fr1=tk.Frame(root)
    tk.Label(fr1, text='Filenumber/ filename filter').pack(side=tk.LEFT)
    tk.Entry(fr1, textvariable=filterstr).pack(side=tk.RIGHT)
    fr1.pack(side=tk.TOP)
    fr2=tk.Frame(root)
    tk.Label(fr2, text='Optional area numbers filter').pack(side=tk.LEFT)
    tk.Entry(fr2, textvariable=areanumstr).pack(side=tk.RIGHT)
    fr2.pack(side=tk.TOP)
    
    fr5=tk.Frame(root)
    # Mark element positions (option)
    b=tk.Label(fr5, text='Peak shift from charging')
    c=tk.Entry(fr5, textvariable=chargeshift)
    b.pack(side=tk.LEFT)
    c.pack(side=tk.LEFT)
    fr5.pack(side=tk.TOP)

    fr7=tk.Frame(root)
    a=tk.Checkbutton(fr7, variable=plotelemsbool, text='Label element peaks?')
    a.pack(side=tk.LEFT)
    tk.Label(fr7, text='Elements:').pack(side=tk.LEFT)
    elementry=tk.Entry(fr7, textvariable=elemstr)
    elementry.pack(side=tk.LEFT)
    fr7.pack(side=tk.TOP)
    
    # TODO does this need to be a module... why not called external function? 
    def updateplot(event):
        ''' Call standard AESplotting function AESplot1 but update to active tk window '''
        if filterstr.get()!='':
            # Simultaneously filter by number on filenumber and by string (on filename)
            filenums, filestrs=parsefilefilters(filterstr.get())
            spelist=spelist[spelist['Filenumber'].isin(filenums)]
            if len(filestrs)>0:
                spelist=spelist[spelist['Filename'].str.contains(filestrs[0])]
                if len(filestrs)>1:
                    print("only using first filter string of two entered.")
            # store filenums used in kwargs (for possible next run)
            if len(spelist)==0:
                print('No files remain after', filterstr.get(),' filter applied.')
                return         
            elif len(spelist)>1:
                # this needs to return a single unique spectrum
                filenames=np.ndarray.tolist(spelist.Filename.unique())
                print('Using first of multiple files selected:', ','.join(filenames))
        # Handle shift from charging
        if chargeshift.get()!=0:
            AESqp=shiftAESparams(AESquantparams, int(chargeshift.get()))
        else:
            AESqp=AESquantparams.copy()
        AugerFileName=spelist.iloc[0]['Filename']
        Augerfile=openspefile(AugerFileName)
        # set areanum for plot  
        if areanumstr.get()!='':
            areanum=int(areanumstr.get())
        else:
            areanum=1
        Augerfile.plot(x='Energy', y='Counts'+str(areanum), ax=ax1)
        Augerfile.plot(x='Energy', y='S7D7'+str(areanum), ax=ax2)
        xmin=Augerfile.Energy.min()
        xmax=Augerfile.Energy.max()
        # Optional element labelling w/ shift (using shifted lines)
        elemlines = findelemlines(Elements, [xmin, xmax], AESqp, peaktype='Counts')
        for i, val in enumerate(elemlines):
            ax1.axvline(x=val, color='b')
        elemlines = findelemlines(Elements, [xmin, xmax], AESqp, peaktype='Deriv')
        for i, val in enumerate(elemlines):
            ax2.axvline(x=val, color='b')
        canvas.show()

    def abort(event):
        choice.set('abort')        
        root.destroy()  

    def plot(event):
        choice.set('plot')        
        root.destroy()
    
    # can't use plt.subplots as that generates separate windows
    fig = Figure(figsize=(5, 4), dpi=100)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    # plt auto-creates a separate window
    canvas=FigureCanvasTkAgg(fig, master=root)
    canvas.show()
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=tk.TOP)
    canvas.get_tk_widget().pack(side=tk.TOP)
    canvas.show()
    
    a=tk.Button(root, text='Plot')
    a.bind('<Button-1>', plot)
    a.pack(side=tk.LEFT)

    a=tk.Button(root, text='Update plot', command= lambda: updateplot)
    # a.bind('<Button-1>', updateplot)
    # adding lambda doesn't work
    # a.bind('<Button-1>', lambda: updateplot(spelist, AESquantparams))
    a.pack(side=tk.LEFT)
    
    a=tk.Button(root, text='Abort')
    a.bind('<Button-1>', abort)
    a.pack(side=tk.LEFT)

    root.mainloop()
    print('Charging shift is', str(int(chargeshift.get())))
    return


from tkinter import filedialog

Smdifpeakslog=smdiffquant_gui(spelist, Elements, AESquantparams, Smdifpeakslog)

def printtest():
    ''' '''
    print('Super dupler')
    return

def shiftedplot_tk():
    ''' Pass expected shift, and list of major elements -- plot and adjust 
    use '''
    

def peakdetector_gui():
    
    rownum+=1
    tk.Label(root, text='Legend size').grid(row=rownum, column=0)
    tk.Entry(root, textvariable=legendsize).grid(row=rownum, column=1)
    
        
def select_ROI(nparr):
    '''
    Select a series of ROIs using matplotlib lasso, name them
    
    '''
    root = tk.Tk()
    root.title("Select ROIs")
    ROIname=tk.StringVar() # Temp string for ROI name
    choice=tk.StringVar()  # plot or abort
    
    tk.Label(root, text='Enter ROI name').grid(row=0, column=0)
    tk.Entry(root, textvariable=ROIname).grid(row=0, column=1)
    # list of ROI names and associated binary arrays
    ROIs=[]
    ROIlist=tk.Listbox()
    ROIlist.grid(row=1, column=0)
    
    # Make mpl figure 
    fig = plt.figure()
    ax1 = fig.add_subplot(121) # 1 x 2 grid -first subplot
    ax1.set_title('Lasso selection:')
    ax1.plot()
    ax1.set_xlim([0,nparr.shape[0]])
    ax1.set_ylim([0, nparr.shape[1]])
    ax1.imshow(nparr) # show numpy array
    
    # Make empty np array to hold an ROI    
    ROIarr=np.zeros((nparr.shape[0],nparr.shape[1]))
    # ROIarr=np.zeros((nparr.shape[0],nparr.shape[1]))
    
    ax2 = fig.add_subplot(122) # 1 x 2 grid -second subplot
    ax2.set_title('ROI selection:')
    # mask showing selected array
    msk = ax2.imshow(ROIarr, origin='lower',vmax=1, interpolation='nearest')
    ax2.set_xlim([-1, nparr.shape[0]+1])
    ax2.set_ylim([-1,  nparr.shape[1]+1])
    
    # show mpl figure in tk drawing area 
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.show()
    canvas.get_tk_widget().grid(row=1, column=1)
    # 
    
    def abort(event):
        ''' Exit interact script and return nothing '''
        choice.set('abort')        
        root.destroy()  

    def close(event):
        ''' close and return list of saved ROI names and associated arrays '''
        choice.set('close')
        root.destroy()  
        
    def saveroi(event):
        ''' Save the current np array for this ROI (and select another)'''
        global ROIarr
        thisname=ROIname.get()
        ROIs.append([thisname, ROIarr]) # Add to data save list
        ROIlist.insert(tk.END, thisname) # insert saved ROI into listbox
        ROIname.set('') # reset name to blank

    def resetROI(event):
        ''' Resets current ROI array to zeros'''
        global ROIarr
        ROIarr=np.zeros((nparr.shape[0],nparr.shape[1]))
        msk.set_data(ROIarr) # reset plot at right
        fig.canvas.draw_idle() # redraw
        return ROIarr

    def onselect(verts): # verts are vertices of selected path upon lasso release
        '''Lasso pixel selector for generating ROIs''' 
        global ROIarr, pix
        p = path.Path(verts)
        ind = p.contains_points(pix, radius=1)
        ROIarr= updateArray(ROIarr, ind)
        msk.set_data(ROIarr)
        fig.canvas.draw_idle() # similar to draw but renders figure only when idle

    def updateArray(array, indices):
        '''ROI display tool for quick region of interest selections (used by select ROI)'''
        lin = np.arange(array.size)
        newArray = array.flatten()
        newArray[lin[indices]] = 1
        return newArray.reshape(array.shape) 

    # 
    pix = np.arange(nparr.shape[0])
    xv, yv = np.meshgrid(pix,pix) # just makes rectangular array
    pix = np.vstack( (xv.flatten(), yv.flatten()) ).T # .T returns transpose of array or matrix

    lasso = LassoSelector(ax1, onselect) # associated with subplot 1
    
    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=2, column=0)

    d=tk.Button(root, text='Reset ROI')
    d.bind('<Button-1>', resetROI)
    d.grid(row=2, column=1)
    
    d=tk.Button(root, text='Save ROI')
    d.bind('<Button-1>', saveroi)
    d.grid(row=2, column=2)

    d=tk.Button(root, text='Close')
    d.bind('<Button-1>', abort)
    d.grid(row=2, column=3)

    root.mainloop()
    
    mychoice=choice.get()
    # This basically works OK
    if mychoice=='

thisroi=ROIs[0][1]
plt.imshow(thisroi)

def repick_elems(AESquantparams):
    '''  '''

def plotpix_tk(specimage, energy, Elemdata, backarray, AESquantparams, **kwargs):
    ''' tk interface for plotting multiplex data from selected pixels 
    within spectral image; 
    kwargs: fileset, areas, xrangestr -- usually just holds values entered during prior run  
    TODO: Add optional plot of shift for this pixel from shiftdict
    '''
    # first print out existing info in various lines
    Elements=[i[0] for i in Elemdata]
    root = tk.Tk()
    root.title("Multiplex plot of spectral image pixel(s)")
    backfitbool=tk.BooleanVar() # Bool for plotting background (if counts plot)
    pixstr=tk.StringVar()
    pixstr.set(kwargs.get('pixstr',"(0,0);(1,1)"))
    labelbool=tk.BooleanVar() 
    shiftbool=tk.BooleanVar() 
    plotelemstr=tk.StringVar()
    mytext='Peaks to be labeled:'+', '.join(Elements)
    plotelemstr.set(mytext)
    choice=tk.StringVar()  # plot or abort
    
    tk.Label(root, text='Enter list of X, Y pixels for plot').grid(row=0, column=0)
    tk.Entry(root, textvariable=pixstr).grid(row=0, column=1)

    d=tk.Checkbutton(root, variable=backfitbool, text='Plot fitted backgrounds?')
    d.grid(row=1, column=0)
    d=tk.Checkbutton(root, variable=labelbool, text='Label peak centers?')
    d.grid(row=1, column=1)
    d=tk.Checkbutton(root, variable=shiftbool, text='Label shift from shift array?')
    d.grid(row=2, column=0)
    # option to reselect labeled elemental peaks 
        
    def abort(event):
        choice.set('abort')        
        root.destroy()  
    def plot(event):
        choice.set('plot')        
        root.destroy()  
    
    d=tk.Button(root, text='Abort')
    d.bind('<Button-1>', abort)
    d.grid(row=3, column=0)

    d=tk.Button(root, text='Plot')
    d.bind('<Button-1>', plot)
    d.grid(row=3, column=1)

    root.mainloop()
        
    mychoice=choice.get()
    
    if mychoice=='plot':
        # Set up kwargs for plot 
        kwargs={}
        if pixstr.get()!='': # extract list of pixels
            pixels=pixstr.get()
            pixels=pixels.split(';')
            # Convert pixels string to X, Y tuple list
            try:
                Xvals=[int(s.split("(")[1].split(',')[0]) for s in pixels]
                Yvals=[int(s.split(",")[1].split(')')[0]) for s in pixels]
                pixels=[list(i) for i in zip(Xvals,Yvals)] # zip together into single list of x,y pixels
                kwargs.update({'pixstr' : pixstr.get()})
            except:
                print('Failed conversion of pixel values from', pixstr.get())
            #return kwargs
        if backfitbool.get(): # include linear background fit
            kwargs.update({'backfit':True})
        if shiftbook.get(): # include linear background fit
            kwargs.update({'showshift':True})
        if labelbool.get(): # label peak centers 
            kwargs.update({'plotelems':Elements})
        plotpixels(specimage, energy, pixels, Elements, backarray, AESquantparams, **kwargs)
        return kwargs
    return
