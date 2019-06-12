# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:10:58 2017

@author: tkc
"""
import os
import pandas as pd
import numpy as np
from tkinter import filedialog
AESQUANTPARAMFILE='C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\AESquantparams.csv'

class AESspectrum():
    ''' Single instance of AES spectra file created from row of spelist (child of AESdataset)
    load file from AESdataset (pd dataframe row) 
    #TODO add direct file load? '''
    def __init__(self, AESdataset, rowindex):
        # can be opened with AESdataset parent and associated row from 
        # open files from directory arg
        self.AESdataset=AESdataset
        self.path=self.parent.path # same path as AESdataset parent
        # load params from batch processing of AESdataset
        row=AESdataset.Augerparamlog.iloc[rowindex]
        self.filename=row.Filename
        self.numareas=row.Numareas
        self.evbreaks=row.Evbreaks # TODO data type?
        
        self.spectype = row.Type.lower() # multiplex or survey 
        
        self.AESdf = None # entire AES dataframe (all areas)
        self.energy = None # same for all cols
        
        self.aesquantparams = None
        self.loadAESquantparams()
        
        self.elems_smdiff = None 
        self.get_elems_smdiff() # get quant from existing smdifpeakslog
        
        self.elems_integ = None # 
        
        self.elemdata = None
        self.getelemdata()
        print('Auger QM file', self.filename, 'loaded.')
    
    def open_csvfile(self):
        ''' Read Auger spectral file '''
        self.AESdf=pd.read_csv(self.filename.replace('.spe','.csv'))
        self.colset=self.AESdf.columns # Counts1, Counts2, S7D71, S7D72, etc.
        self.energy=self.AESdf['Energy']
        self.backfit=self.EDXdf['Backfit']
        self.subdata=self.EDXdf['Subdata']
        print('EDXfile ', self.filename,' loaded.')
        
    def loadAESquantparams(self):
        ''' Loads standard values of Auger quant parameters 
        TODO what about dealing with local shifts   '''
        # Checkbutton option for local (or standard) AESquantparams in file loader?
        print('AESquantparams loaded')
        self.aesquantparams=pd.read_csv(AESQUANTPARAMFILE, encoding='utf-8')
    
    def get_elems_smdiff(self):
        ''' Finds element quant already performed from smdifflog (within AESdataset) '''
        match=self.AESdataset.Smdifpeakslog.loc[ (self.AESdataset.Smdifpeakslog['Filename']==self.filename)]
        # should contain row for each element included in quant
        self.elems_smdiff=[] # elem/ peak name
        self.smdiff_shifts=[] # shifts from ideal position
        self.smdiff_ampl=[] # negintensity - posintensity
        self.smdiff_widths=[] # ev diff between negpeak and pospeak
        for index, row in match.iterrows():
            self.elems_smdiff.append(row.PeakID)
            self.smdiff_shifts.append(row.Shift)
            self.smdiff_ampl.append(row.Amplitude)
            self.smdiff_widths.append(row.Peakwidth)

    def get_elems_integ(self):
        ''' Pull existing quant results from integ log file (if present)  '''
        pass

            
    def savecsv():
        ''' Save any changes to underlying csv file '''

class AESdataset():
    ''' loads all dataframes with Auger parameters from current project folder '''
    def __init__(self, *args, **kwargs):
        self.path = filedialog.askdirectory()
        # open files 
        self.open_main_files()
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
            self.Augerparamlog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
            self.spelist=self.Augerparamlog[pd.notnull(self.Augerparamlog['Areas'])]
        else:
            self.Augerparamlog=pd.DataFrame()
            self.spelist=pd.DataFrame()
        if os.path.exists('Smdifpeakslog.csv'):
            self.Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437')
        else:
            self.Smdifpeakslog=pd.DataFrame()
        if os.path.exists('Backfitlog.csv'):
            self.Backfitlog=pd.read_csv('Backfitlog.csv', encoding='cp437')
        else:
            self.Backfitlog=pd.DataFrame()
        if os.path.exists('Integquantlog.csv'):
            self.Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
        else:
            self.Integquantlog=pd.DataFrame()
        # Print TEM or SEM to console based on beam kV
        try:
            self.AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\AESquantparams.csv', encoding='utf-8')
        except:
            self.AESquantparams=pd.DataFrame()
 
