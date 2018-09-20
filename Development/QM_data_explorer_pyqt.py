# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 07:45:29 2017
pixel by pixel spectral plotter for quant mapped arrays 
choose x and y pixels and underlying spectrum appears
optional choice of elemental subregion

Far from working 
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


class QMdataset():
    ''' loads all dataframes with Auger parameters from current project folder '''
    def __init__(self, specimage, energy, Elements, Elemdata, backarray, **kwargs):
        # open files
        self.energy=energy
        self.specimage=specimage
        self.Elements=Elements
        self.Elemdata=Elemdata
        self.backarray=backarray