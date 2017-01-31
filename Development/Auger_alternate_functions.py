# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 12:02:10 2016

@author: tkc
"""

#%% Load modules
import matplotlib.pyplot as plt
import re, os, glob, sys, csv # already run with functions 
import pandas as pd
import numpy as np
from pylab import *
import scipy
import scipy.stats
from matplotlib.backends.backend_pdf import PdfPages # Only needed for plotting
from collections import defaultdict
from math import factorial # used by Savgol matrix

#%%
# alt version using np.polyfit instead of scipy stats
def fitlinearregionalt(df, areanum, AugerFileName):
    '''Pass appropriate chunk from Auger spectral dataframe, perform linear fit
    return chunk with backfit column added '''
    colname='Counts'+str(areanum)
    backfitname='Backfit'+str(areanum)
    xcol=df['Energy']
    ycol=df[colname] # Counts1, Counts2 or whatever
    try:
        slope,intercept=np.polyfit(xcol, ycol, 1, cov=True) # numpy fitting version   
        # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xcol, ycol)
    except: # deal with common problems with linregress
        print('Fitting error for file ', AugerFileName)
        fitparams=('n/a','n/a','n/a','n/a','n/a') # return all n/a
        return df, fitparams
    
    fitparams=(slope, intercept, r_value, p_value, std_err) # tuple to return fitting results
    for index,row in df.iterrows():
        xval=df.loc[index]['Energy']
        yval=slope * xval + intercept
        df=df.set_value(index, backfitname, yval)
    return df, fitparams