# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:10:58 2017

@author: tkc
"""
import os
import pandas as pd
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
        self.path=self.AESdataset.path # same path as AESdataset parent
        # load params from spelist only (not AESlog which has images)
        row=AESdataset.spelist.iloc[rowindex]
        self.filename=row.Filename
        self.sample=str(row.Sample)
        self.numareas=int(row.Areas)
        self.evbreaks=row.Evbreaks # TODO data type?
        
        self.spectype = row.Type.lower() # multiplex or survey 
        
        self.AESdf = None # entire AES dataframe (all areas)
        self.energy = None # same for all cols
        self.open_csvfile()
        
        self.aesquantparams = None
        self.loadAESquantparams()
        
        # load peaks, shifts, ampls, widths
        self.smdifpeakinfo=None # dataframe w/ smdiff peak info
        self.get_peaks() # get quant info from smdifpeakslog
        self.integpeakinfo=None # dataframe w/ smdiff peak info
        self.get_integ_peaks() # get quant info from smdifpeakslog
        
        self.elems_integ = None # 
        
        print('Auger file', self.filename, 'loaded.')
    
    def open_csvfile(self):
        ''' Read Auger spectral file '''
        self.AESdf=pd.read_csv(self.filename.replace('.spe','.csv'))
        self.colset=self.AESdf.columns # Counts1, Counts2, S7D71, S7D72, etc.
        self.energy=self.AESdf['Energy']
        print('AESfile ', self.filename,' loaded.')
        
    def loadAESquantparams(self):
        ''' Loads standard values of Auger quant parameters 
        TODO what about dealing with local shifts   '''
        # Checkbutton option for local (or standard) AESquantparams in file loader?
        print('AESquantparams loaded')
        self.aesquantparams=pd.read_csv(AESQUANTPARAMFILE, encoding='utf-8')

    def get_peaks(self):
        ''' Finds element quant already performed from smdifflog (within AESdataset)
        needed for plots:  negpeak, pospeak (both indirect calc from shift, peakwidth)
        negint, posint (but usually indirectly

        ideal positions for peaks already loaded in AESdataset
        smdifpeakinfo contains: 0 peak, 1 negpeak energy, 2 pospeak energy, 3) negint 
        4) posint 5) ampl, 6) adjampl -- important stuff for graphical display
        of quant results
        returns dataframe for this filename
        '''
        mycols=['Areanumber', 'Peakenergy', 'Peakindex', 'PeakID', 'Shift',
       'Negintensity', 'Posintensity', 'Pospeak', 'Amplitude', 'Peakwidth','Adjamp']
        self.smdifpeakinfo=self.AESdataset.Smdifpeakslog[ (self.AESdataset.Smdifpeakslog['Filename']==self.filename)]
        self.smdifpeakinfo=self.smdifpeakinfo[mycols]
        
    def get_integ_peaks(self):
        ''' Pull existing quant results from integ log file (if present)  '''
        mycols=['Areanumber', 'Element', 'Integcounts', 'Backcounts', 
                'Significance', 'Adjcnts','Erradjcnts']
        self.integpeakinfo=self.AESdataset.Integquantlog[ (self.AESdataset.Integquantlog['Filename']==self.filename)]
        self.integpeakinfo=self.integpeakinfo[mycols]


    def savecsv():
        ''' Save any changes to underlying csv file '''

class AESdataset():
    ''' loads all dataframes with Auger parameters from current project folder '''
    def __init__(self, *args, **kwargs):
        self.path = filedialog.askdirectory()
        # open files 
        
        self.AESlog=None
        self.spelist=None
        self.Smdifpeakslog=None
        self.Integquantlog=None
        self.Backfitlog=None
        self.open_main_files() # loads above
        # self.filelist=np.ndarray.tolist(self.AESlog.Filenumber.unique())
        self.numfiles=len(self.AESlog)
        print(str(self.numfiles),' loaded from AESdataset.')

        self.peaks=None
        self.peakdata=None
        self.get_peakinfo() # load needed Auger peak params (Peaks and Peakdata)
    
    def get_peakinfo(self):
        ''' takes element strings and energies of background regs and returns tuple for each elem symbol containing all params necessary to find each Auger peak from given spe file 
        also returns 2-tuple with energy val and index of chosen background regions
        ''' 
        # elemental lines (incl Fe2, Fe1, etc.)
        self.peaks=self.Smdifpeakslog.PeakID.unique()
        self.peakdata=[]

        for peak in self.peaks:
            try:
                # find row in AESquantparams for this element
                thispeakdata=self.AESquantparams[(self.AESquantparams['element']==peak)]
                thispeakdata=thispeakdata.squeeze() # series with this elements params
                # return list of length numelements with 5-tuple for each containing 
                # 1) peak symbol, 2) ideal negpeak (eV)  3) ideal pospeak (in eV) 
                # 4)sensitivity kfactor.. and 5) error in kfactor
        
                peaktuple=(peak, thispeakdata.negpeak, thispeakdata.pospeak, 
                           thispeakdata.kfactor, thispeakdata.errkf1) # add tuple with info for this element     
                self.peakdata.append(peaktuple)
            except:
                print('AESquantparams not found for ', peak)
        print('Found', len(self.peakdata), 'quant peaks in smdifpeakslog' )

    def open_main_files(self):
        ''' Auto loads Auger param files from working directory including 
        AESparalog- assorted params associated w/ each SEM-AES or TEM-AES emsa file 
        Backfitparamslog - ranges and parameters for AES background fits
        Integquantlog - subtracted and corrected counts for chosen elements
        Peakfitlog - params of gaussian fits to each element (xc, width, peakarea, Y0, rsquared)'''
        if os.path.exists('Augerparamlog.csv'):
            self.AESlog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
            self.spelist=self.AESlog[pd.notnull(self.AESlog['Areas'])]
        else:
            self.AESlog=pd.DataFrame()
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
 
