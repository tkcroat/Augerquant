# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 07:45:29 2017
Auger plot/ shift determination 
@author: tkc
"""
import os, glob
from PyQt5 import QtGui, QtCore, QtWidgets
import numpy as np 
import pandas as pd

from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib as mpl
mpl.use('Qt5Agg') # choose option for interactive mpl 
from matplotlib.widgets import RectangleSelector


class AESdataset():
    ''' loads all dataframes with Auger parameters from current project folder '''
    def __init__(self, *args, **kwargs):
        # open files 
        self.Augerparamlog, self.spelist, self.Smdifpeakslog, self.Backfitlog, self.Integquantlog, self.AESquantparams=self.open_main_files()
        self.filelist=np.ndarray.tolist(self.Augerparamlog.Filenumber.unique())
        self.numfiles=len(self.Augerparamlog)
        print(str(self.numfiles),' loaded from EDXdataset.')
    
    def open_main_files(self):
        ''' Auto loads Auger param files from working directory including 
        EDXparalog- assorted params associated w/ each SEM-EDX or TEM-EDX emsa file 
        Backfitparamslog - ranges and parameters for EDX background fits
        Integquantlog - subtracted and corrected counts for chosen elements
        Peakfitlog - params of gaussian fits to each element (xc, width, peakarea, Y0, rsquared)'''
        if os.path.exists('Augerparamlog.csv'):
            Augerparamlog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
            spelist=Augerparamlog[pd.notnull(Augerparamlog['Areas'])]
        else:
            Augerparamlog=pd.DataFrame()
            spelist=pd.DataFrame()
        if os.path.exists('Smdifpeakslog.csv'):
            Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437')
        else:
            Smdifpeakslog=pd.DataFrame()
        if os.path.exists('Backfitlog.csv'):
            Backfitlog=pd.read_csv('Backfitlog.csv', encoding='cp437')
        else:
            Backfitlog=pd.DataFrame()
        if os.path.exists('Integquantlog.csv'):
            Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
        else:
            Integquantlog=pd.DataFrame()
        # Print TEM or SEM to console based on beam kV
        try:
            AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\AESquantparams.csv', encoding='utf-8')
        except:
            AESquantparams=pd.DataFrame()
        return Augerparamlog, spelist, Smdifpeakslog, Backfitlog, Integquantlog, AESquantparams

