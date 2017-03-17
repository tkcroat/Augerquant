# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 17:10:19 2016

@author: tkc
"""
import pandas as pd
import numpy as np
import shutil, sys, fileinput, os, math

#%%

# Can make basis and compositions for each area (normal quant process) 
# Superimpose on grayscale image?

# Sub-region selector for 10x10 arrays
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib import path
from PIL import Image, ImageDraw, ImageFont 
from collections import defaultdict
import numpy as np
# Plotting all numpy histograms from elementmaps
def plothistograms(elementmaps, Elements)

Mgmap=elementmaps[0][1]
plt.imshow(image)
plt.imshow(Mgmap, alpha=0.4) # transparent overlay

filenum=148167

x1=100
x2=150
y1=0
y2=50
val=20000
testmap=assignmaprange(Mgmap, x1,x2,y1,y2, val)
plt.figimage(Simap)

plt.hist(Mgmap.ravel(), bins=20)


ROIarray=selectROI()

myarray=makeROI('tr184.171.jpg')
myarray=makeROI('F01_3Mar17_ibeam_075.jpg')

def makeROI(imagename):
    '''Returns ROI from image either for QM data averaging or for setting active areas for mapping 
    has several inner functions '''
    fig = plt.figure()
    ax1 = fig.add_subplot(121) # 1 x 2 grid -first subplot
    ax1.set_title('Lasso selection:')
    ax1.plot()
    ax1.set_xlim([0, 512])
    ax1.set_ylim([0, 512])    
    ax1.set_aspect('equal')
    thisim=Image.open(imagename)
    ax1.imshow(thisim) # added image
    
    # Empty array to be filled with lasso selector
    array = np.zeros((512,512))
    ax2 = fig.add_subplot(122) # 1 x 2 grid -second subplot
    ax2.set_title('numpy array:')
    msk = ax2.imshow(array, origin='lower',vmax=1, interpolation='nearest')
    ax2.set_xlim([-1, 512])
    ax2.set_ylim([-1, 512])
    
    # Pixel coordinates
    pix = np.arange(512)
    xv, yv = np.meshgrid(pix,pix) # just makes rectangular array
    pix = np.vstack( (xv.flatten(), yv.flatten()) ).T # .T returns transpose of array or matrix
    return ax1, array
    
    def onselect(verts): # verts are vertices of selected path upon lasso release
        '''Lasso pixel selector for generating ROIs''' 
        global array, pix
        p = path.Path(verts)
        ind = p.contains_points(pix, radius=1)
        array = updateArray(array, ind)
        msk.set_data(array)
        fig.canvas.draw_idle() # similar to draw but renders figure only when idle

    def updateArray(array, indices):
        '''ROI display tool for quick region of interest selections (used by select ROI)'''
        lin = np.arange(array.size)
        newArray = array.flatten()
        newArray[lin[indices]] = 1
        return newArray.reshape(array.shape) 
        
    lasso = LassoSelector(ax1, onselect) # associated with subplot 1
    
    plt.show()
    
    return array


def createRGB(elementmaps, elemlist, savename=''):
    '''Pass 3 elements  '''    
    
def assignmaprange(nparr, x1,x2,y1,y2, val):
    ''' Assign single value to range in 512 x512 array '''
    for i in range(x1,x2):
        for j in range(y1,y2):
            nparr[i][j]=val
    return nparr
    
# some problem with wrapping this lasso in a function -- event handling problem and selection is thrown away 
# see http://stackoverflow.com/questions/17111859/matplotlib-spanselector-widget-how-to-use-inside-a-function

# can you return axis from first function then run 

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
    

def getmappedelements(filenum):
    ''' Returns list of element's mapped in given ... retrieved from multiplex  '''
    # direct read from multiplex spectral conditions
    if filenum>1000: # indicates combined file w/ firstfile-lastfile
        mystr=str(filenum)
        thislen=int(len(mystr)/2)
        filenum=int(mystr[0:thislen]) # convert name to first file number
    mult=pd.read_csv('multiplexspectralregionslog.csv', encoding='utf-8')
    match=mult[mult['Filenumber']==filenum]
    elements=np.ndarray.tolist(match.Element.unique())
    Multielem = defaultdict(list)
    for i, item in enumerate(elements):
        Multielem[item].append(i)
    Multielem = {k:v for k,v in Multielem.items() if len(v)>1} # dictionary with duplicated item and list with indices
   
    newelemlist, Multielem= parseelemlist(elements)
    # remove number suffix from most elements 
    newlist=[re.match('\D+',i).group(0) for i in elements] # D is non-digit
    # if duplicate peaks, return that element with attached numbers, else return as string only


    
plt.imshow(array)

plt.imshow(thisim)

# Testing 
ax1.autoscale(True)
ax1.imshow(thisim)

    return array
    



       
