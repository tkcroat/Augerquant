import numpy as np
import os, glob, re
import pandas as pd
import struct
from scipy.ndimage import filters  # uniform filter and others
from tkinter import filedialog
from operator import itemgetter
from itertools import groupby
from scipy.signal import argrelmin
from PIL import Image

AESQUANTPARAMFILE='C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\AESquantparams.csv'

class AESquantmap():
    ''' Data structure for auger quant maps (using multiplex mode) '''
    def __init__(self, directory):
        # os.chdir(directory) # use loader instead 
        self.directory=directory # single quantmap per directory
        self.pixarray = pd.DataFrame() # linkage between pixels in quant map and spe data files
        self.specregions = None # elem, evrange, etc. info from multiplex 
        self.uniquename = None # unique prefix string for all files
        # pixarray csv contains instructions to assemble spectral image from assorted files
        self.loadQMpixarray() 
        
        self.aesquantparams = pd.DataFrame()
        # all info on Auger peaks
        self.loadAESquantparams()
        
        self.QMfilename = None # name of existing npy file
        self.specimage = np.empty([0,0,0]) # 3d numpy array w/ entire multiplex
        
        self.dim = None #  (X, Y, Z) of spectral image
        self.shiftmaps = [] # peakshift 2D array of dimension X * Y (one per element) 
        self.amplmaps = [] # s7d7 amplitudes 2D array of dimension X * Y (one per element) 
        self.integmaps = []
        self.peakintmaps = []  # max intensity values at direct peak
        self.elemmaps = [] # for storage of normalized element quant maps
        self.specfilter = None # short term storage for filtered spectrum
        
        self.energy = None # multiplex energy values (sorted)
        self.evbreaks = None # multiplex energy breaks
        self.spectralregs=pd.DataFrame()        
        self.elements = [] # returned by get_elemdata
        self.elemdata = []
        self.loadspecimage() # attempt load/create of spectral image
        # calls makespecimage if it doesn't exist
        # also calls get_spectral_regs and get_elemdata        
        
        self.extracted = None # temp spectrum extracted from 
        self.extracts7d7 = None
        self.exttype = None # pixel or lasso
        self.quant =[] # list for quant results from extracted spectrum
        self.quant2=[] # quant results from direct integ (not deriv)
        self.derivparams=[] # pospeak and negpeak vals for deriv spectral plot
        self.integparams=[] # peak intensity and energy for direct spectral plot
        # Go ahead and attempt load of existing amplmaps, integmaps
        self.load_maps()
        
    def loadAESquantparams(self):
        ''' Loads standard values of Auger quant parameters 
        TODO what about dealing with local shifts   '''
        # Checkbutton option for local (or standard) AESquantparams in file loader?
        self.aesquantparams=pd.read_csv(AESQUANTPARAMFILE, encoding='utf-8')
        print('AESquantparams loaded')
        
    def loadQMpixarray(self):
        ''' Load of standard files from main Auger data directory
        button linked or auto-loaded?  '''
        os.chdir(self.directory)
        pixfiles=glob.glob('*QMpixarray.csv')
        if len(pixfiles)==1:
            self.pixarray=pd.read_csv(pixfiles[0])
            if 'Filename' not in self.pixarray:
                print('spe data file names not yet linked with quantmap pixel array.')
            # file naming conventions is 'uniquestring'+_QMpixarray.csv
            self.uniquename = pixfiles[0].split('_QMpixarray.csv')[0]
            print('QMpix array loaded')
        else:
            print("Couldn't locate single pixarray definition file")
        if 'Filename' not in self.pixarray.columns:
            print('QM pixels not yet linked to multiplex files')
            self.link_files()
    
    def loadspecimage(self):
        ''' Load existing numpy specimage array (if already created) '''
        os.chdir(self.directory)
        npyfiles=glob.glob('*specimage.npy')
        if len(npyfiles)==1:
            self.specimage=np.load(npyfiles[0])
            self.dim=self.specimage.shape
            self.QMfilename=npyfiles[0]
            self.get_spectral_regs() # get energy structure & eVbreaks
            self.get_elemdata() # extracts elements and elemdata
            print('Spectral image ', self.QMfilename,' loaded')
        else:
            print("Couldn't locate single spectral image numpy 3D array in current directory")
            # Attempt to create this from pix array and all spe files
            self.make_specimage()

    def save_specimage(self):
        ''' Called by menu to GUIroi  '''
        if self.specimage.shape[0]==0:
            return
        fname=self.directory+'/'+self.uniquename+'_specimage.npy'
        np.save(fname,self.specimage)
        print('Spectral image array saved.')

    def save_pixarray(self):
        ''' Called after link_files finds correct associated spe files '''
        if self.pixarray.empty:
            return
        fname=self.uniquename + '_QMpixarray.csv'
        self.pixarray.to_csv(fname, index=False)
        print('Altered pixarray csv file saved.')
        
    def make_specimage(self):
        ''' Create 3D numpy array from underlying multiplex files '''
        print('Making spectral image')
        # Check to ensure file linkage made in QMpixarray
        if self.pixarray.empty:
            print('QM pix array file is not loaded')
            return
        if 'Filename' not in self.pixarray.columns:
            print('QM pixels not yet linked to multiplex files.. link_files!!')
            self.link_files()
            # Get energy range for multiplex from first data file (same for all pixels)
        # get energy values, evbreaks, and multiplex details from first spe file
        self.get_spectral_regs()
        xmax=self.pixarray.Xindex.max()
        ymax=self.pixarray.Yindex.max()
        # make 3D blank numpy array of dimension xmax+1 * ymax+1 * len(energy)-duplicates
        # using set since inadvertent duplicate ev vals are removed (using keep_best_vals)
        self.specimage=np.empty([xmax+1,ymax+1,len(set(self.energy))])
    
        # Now extract counts from all areas in each multiplex (20 max per spe)
        spefilelist=np.ndarray.tolist(self.pixarray.Filename.unique())
        for i, file in enumerate(spefilelist):
            thisfile=self.pixarray[self.pixarray['Filename']==file]
            numareas=len(thisfile)
            # dataframe containing all spatial areas from this spe
            thismult=self.get_multiplex_data(file, numareas)
            thismult=self.keep_best_vals(thismult)
            if thismult.empty: 
                print('File missing ... no data for ', file)
                continue
            for index, row in thisfile.iterrows():
                xind=thisfile.loc[index]['Xindex']
                yind=thisfile.loc[index]['Yindex']
                thiscolnum=thisfile.loc[index]['Subnumber']+1
                if len(thismult)!=self.specimage.shape[2]:
                    print('spe spectrum has different length than specimage array!')
                    continue
                else:
                    self.specimage[xind, yind,:]=thismult['Counts'+str(thiscolnum)]
        # now remove duplicate eV vals
        self.energy=list(set(self.energy))
        self.energy.sort() # data vals are auto-sorted in get_multiplex_data
        print('Spectral image created.')
        self.dim=self.specimage.shape # Get spectral image dimensions
        self.get_elemdata() # now grab elemdata 
        
    def keep_best_vals(self, thismult):
        '''For duplicate energy values in multiplex scans, keep the one with 
        largest # of sweeps '''
        # make temp index ranges for spectralregs (can't use energy values)
        start=0
        self.spectralregs['Start']=0 # for index #s
        self.spectralregs['End']=0
        for index,row in self.spectralregs.iterrows():
            lower=self.spectralregs.loc[index]['Lower']
            upper=self.spectralregs.loc[index]['Upper']
            thisrange=int(upper-lower)
            self.spectralregs=self.spectralregs.set_value(index,'Start',start)
            self.spectralregs=self.spectralregs.set_value(index,'End',start+thisrange)
            start=start+thisrange+1 # adjust for next loop
        dupl=thismult.duplicated(['Energy'], keep=False) # only duplicate vals and keep both
        dupl=thismult.loc[dupl]
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
                specmatch1=self.spectralregs[(self.spectralregs['Start']<=index1)&(self.spectralregs['End']>=index1)]
                specmatch2=self.spectralregs[(self.spectralregs['Start']<=index2)&(self.spectralregs['End']>=index2)]
                
                try:
                    if specmatch1.iloc[0]['Sweeps']>=specmatch2.iloc[0]['Sweeps']:
                        # first is best value... remove it from dupl (which will be used as knockout df)
                        removelist.append(index1)
                    else:
                        removelist.append(index2)
                except:
                    print('Problem with duplicate removal of ', index1, index2)
        thismult=thismult[-thismult.index.isin(removelist)]
        print (len(removelist), ' duplicated energy values removed from multiplex')
        return thismult
    
    def get_spectral_regs(self):
        ''' Gets energy range and info about spectral regions that comprise the underlying multiplex spectra
        run once on first spectra?  '''
        print('Loading spectral regions.')
        if self.pixarray.empty:
            print('QM pixarray file must be loaded')
            return
        AugerFileName=self.pixarray.iloc[0]['Filename']
        filenumber=AugerFileName.split('.')[1].split('.')[0]
        # check for file in cwd or sub
        if not os.path.exists(self.directory+'/'+AugerFileName):
            if not os.path.exists(self.directory+'/sub/'+AugerFileName): # check sub directory
                print('Datafile ', AugerFileName,' not found in cwd or sub directories.')
            else:
                AugerFileName='/sub/'+AugerFileName
                if not os.path.exists(self.directory+'/'+AugerFileName):
                    print('Failed load of spectral regs... first data file not found')
                    return
        with open(self.directory+'/'+AugerFileName, 'rb') as file:
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
        # Spectral image data cube signal always sorted by energy (in get_multiplex_data)
        # need to auto-sort here as well
        tempstring=header.split('NoSpectralReg: ')[1] # unlike NoSpectralRegFull inactives are already removed
        match=re.search(r'\d+',tempstring)
        numdefregions=int(match.group(0)) # number of defined regions (can be larger than # active regions)
        self.spectralregs=pd.DataFrame(columns=['Filenumber','Filename','Numcycles','Timestep','Element','Sweeps','Evstep','Lower','Upper','Time'])
        numregions=0 # active regions (as opposed to nominal defined regions)
        timeperarea=0 #initialize as zero
        self.energy=[] # list for energy x values
        self.evbreaks=[0] # region boundaries needed for smoothdiff include first
        for i in range(numdefregions):
            tempstr=tempstring.split('SpectralRegDef: ')[i+1] # should work to split 
            element=(tempstr.split(' ')[2]) # name of elemental line
            numpts=int(tempstr.split(' ')[4])        
            evstep=float(tempstr.split(' ')[5]) #eV/step        
            startev=float(tempstr.split(' ')[6]) # starting eV
            endev=float(tempstr.split(' ')[7]) # ending eV
            for j in range(0,numpts): # generate energy values for survey
                self.energy.append(startev+evstep*j) # cannot remove energy duplicates yet
            self.evbreaks.append(len(self.energy)-1) # gives indices of discontinuities in energy scan (needed for smoothdiffS7D7)   
            tempstr=tempstr.split('SpectralRegDef2: ')[1]
            sweeps=int(tempstr.split(' ')[2]) # number of sweeps through element        
            time=numcycles*timestep*sweeps*(endev-startev)/1000 # acquisition time in seconds for this multiplex region        
            timeperarea+=time # add on acquisition time for this element eV range only
            # add ith new data row to spatialareas dataframe
            # make a list and then append as row? 
            self.spectralregs.loc[i]=[filenumber, AugerFileName, numcycles, timestep, element, sweeps, evstep, startev, endev, time]
            numregions+=1
            
        self.spectralregs=self.spectralregs.sort_values(['Lower'])
        self.spectralregs=self.spectralregs.reset_index(drop=True)
        # make indexrange showing Ev range of elements as index #s with specimage
        self.spectralregs['Start']=np.nan
        self.spectralregs['End']=np.nan
        # Need to remove overlapping vals from spectralregs (keep one w/ large # of sweeps)
        # These duplicates are removed from multiplex spectra and energy elsewhere
        for i in range(0, len(self.spectralregs)-1):
            if self.spectralregs.iloc[i+1]['Lower']<=self.spectralregs.iloc[i]['Upper']:
                redval=self.spectralregs.iloc[i]['Upper']-self.spectralregs.iloc[i+1]['Lower']+1
                # remove overlapping vals from inferior one
                if self.spectralregs.iloc[i+1]['Sweeps']>self.spectralregs.iloc[i]['Sweeps']:
                    self.spectralregs=self.spectralregs.set_value(i+1,'Lower',self.spectralregs.iloc[i+1]['Lower']+redval)
                else:
                    self.spectralregs=self.spectralregs.set_value(i,'Upper',self.spectralregs.iloc[i]['Upper']-redval)
        count=0
        for index, row in self.spectralregs.iterrows():
            thisrange=row.Upper-row.Lower+1
            # zero based indexing
            self.spectralregs=self.spectralregs.set_value(index, 'Start', count) 
            self.spectralregs=self.spectralregs.set_value(index, 'End', count+thisrange-1)
            count+=thisrange
        # Rename certain peaks/elem to align peak naming conventions (generally main peak w/o appended number)
        peakdict={'Mg2':'Mg','Si2':'Si','S1':'S','Fe3':'Fe'}
        for key, val in peakdict.items():
            self.spectralregs['Element']=self.spectralregs['Element'].str.replace(key, val)
        self.energy.sort() # extracted multiplex spectra always sorted by increasing eV
        ''' Get rid of accidental duplicated energy vals .. keep_best_values removes 
        these from multiplex spectra themselves
        '''
        print('Spectral regs loaded.',str(len(self.energy)), ' energy values.')

    def get_elemdata(self, **kwargs):
        ''' Return index ranges for each element's peak, low background, and high background
        these are sometimes needed separate from spectralregs
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
        10) ideal negpeak ev (position of negative peak in smooth-diff deriv s7d7).. 
           slightly higher than idealev (9 in list)
        11) typical integration peak width
        12) chargeshift -- # of eV of charge-compensation shift applied to scan  
		     figure out by comparison of AESquantparams w/ scan center
        13) peakwidth - element specific negpeak-pospeak
        14) searchwidth - region 
        kwarg: Elements - explicit element list (normally autogenerated)...
        only needed if multiplex setup is unusual
        '''
        print('Loading element data.')
        if self.spectralregs.empty or self.aesquantparams.empty:
            print('Cannot load element data... extract spectral regs from multiplex and load AESquantparams!')
            return

        self.elemdata=[] # list of lists with important element data
        self.elements=[]
        # determine if spectral regs are split into OL, O, OH or not
        Elements=np.ndarray.tolist(self.spectralregs.Element.unique())
        if self.spectralregs.Element[0].endswith('L'):
            # works only with OL, O, OH type setup (not continuous single O region)
            Elements=Elements[1::3]
            # if not this structure then above Elements list should be correct
        validelem=np.ndarray.tolist(self.aesquantparams['element'].unique())
        missing=[elem for elem in Elements if elem not in validelem]
        if len(missing)!=0:
            print(','.join(missing),' are not valid elements/peaks in AESquantparams')
        self.elements=Elements
        print('Elements are',','.join(Elements))
        # Ensure all are standard peak names
        for i, elem in enumerate(Elements):
            thispeak=self.spectralregs[self.spectralregs['Element']==elem]
            lowback=self.spectralregs[self.spectralregs['Element']==elem+'L']
            hiback=self.spectralregs[self.spectralregs['Element']==elem+'H']
            match=self.aesquantparams[self.aesquantparams['element']==elem]
            if len(match)!=1 | len(thispeak)!=1 :
                print('Element', elem,' not found in AESquantparams or associated peak missing')
                continue
            idealev=int(match.iloc[0]['negpeak']+match.iloc[0]['integpeak'])
            idealnegpeak=int(match.iloc[0]['negpeak'])
            width=int(match.iloc[0]['integwidth'])
            # typical distance between neg and pospeaks in smooth-diff data
            peakwidth=int(match.iloc[0]['peakwidth'])
            searchwidth=int(match.iloc[0]['searchwidth'])
            kfact=int(match.iloc[0]['kfactor']) # kfactor associated with smooth-diff method
            kfact2=int(match.iloc[0]['kfactor2'])
            mass=int(match.iloc[0]['mass']) 
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
                # Artificially split single peak region into low, peak, high 
                # Set index numbers
                fullrange=maxrange-minrange
                minlow=minrange
                maxlow=int(minrange+fullrange*0.25)
                minhigh=int(maxrange-fullrange*0.25)
                maxhigh=maxrange
                minrange=maxlow+1
                maxrange=minhigh-1
                # Set these values in eV
                fullrange=maxevrange-minevrange
                minevlow=minevrange
                maxevlow=int(minevrange+fullrange*0.25)
                minevhigh=int(maxevrange-fullrange*0.25)
                maxevhigh=maxevrange
                minevrange=maxevlow+1
                maxevrange=minevhigh-1                
            self.elemdata.append([elem, i, [minrange, maxrange], [minlow, maxlow], [minhigh, maxhigh], 
                [minevrange, maxevrange], [minevlow, maxevlow], [minevhigh, maxevhigh], idealev, idealind, 
                idealnegpeak, width, chargeshift, peakwidth, searchwidth, kfact, kfact2, mass])
        print('Elem data loaded with ',str(len(self.elemdata)), ' regions.')

    def get_multiplex_data(self, AugerFileName, numareas):
        ''' Extracts multiplex spectra from all spatial areas within single spe multiplex file
        also uses energy and evbreaks attribs '''
        if not os.path.exists(self.directory+'/'+AugerFileName):
            if os.path.exists(self.directory+'/sub/'+AugerFileName):
                AugerFileName='sub/'+AugerFileName # open from sub directory
        try:
            with open(self.directory+'/'+AugerFileName, 'rb') as file:
                filedata = file.read()
        except:
            print(AugerFileName," missing from data directory")
            return pd.DataFrame() # return empty frame?
        end=filedata.find(b'EOFH')
        bindata=filedata[end+6:] # binary data section of file (header removed)
        mystruct=struct.Struct('f')
        # binary header of variable length.. just count backward using known data length
        startbyte=len(bindata)-4*numareas*len(self.energy)
        
        # TEST for correct byte location... data regions begins right after last occurence of triple zero 4 byte float values.. for survey usually found in bytes 100:111 
        for i in range(startbyte-12,startbyte,4):
            byteval=bindata[i:i+4]
            unpackbyte=mystruct.unpack(byteval) # little endian encoding
            if unpackbyte[0]!=0:
                print('Possible binary read error... leading zeros are missing for ',AugerFileName)
    
        # create data frame for multiplex and assign energy values   
        multiplex=pd.DataFrame() # data frame for energy, and counts col for each area   
        multiplex['Energy']=self.energy # energy values found in SpectralRegions and passed as argument
        
        # Read in and convert all numeric values
        alldata=[] # single list for all Auger counts values in file
        for i in range(startbyte,len(bindata),4):
            byteval=bindata[i:i+4]
            unpackbyte=mystruct.unpack(byteval)
            alldata.append(unpackbyte[0])
        
        if len(alldata)!=len(self.energy)*numareas: # number of data values expected for each area
            print('Possible binary data reading error: Data string length does not match expected value for ',AugerFileName)
        ''' Multiplex file structure has same energy region of all areas bundled together 
        (split counts into spectral regions based on evbreaks)so organization is 
        multiplex spectral region 1 (all areas), spectral region 2 (all areas), etc.  
        '''
        datachunks=[]
        for i in range(0,len(self.evbreaks)-1):
            datachunks.append(self.evbreaks[i+1]-self.evbreaks[i])
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
        
    def link_files(self):
        ''' Create links between pix array and data structure (before'''
        # Prompt to get first spe file of this quantmap
        fullpath= filedialog.askopenfilename(title='Open first spe file of quantmap', 
                                             filetypes=[("spectrum","*.spe")])
        (directory, filename)=os.path.split(fullpath)
        basename=filename.split('.')[0]
        startnum=int(filename.split('.')[1])
        print('First file of quantmap is', basename, ' ', str(startnum))
        self.pixarray['Filename']='' # initialize new string column
        for index,row in self.pixarray.iterrows():
            self.pixarray=self.pixarray.set_value(index, 'Filename', 
                        basename+'.'+str(startnum+index//20)+'.spe')
        # Ensure that expected spe data files are present
        self.check_missing_spe(directory)
        
    def check_missing_spe(self, directory):
        ''' If link files is invoked, check to ensure all spe files exist in 
        proper directory '''
        spefiles=np.ndarray.tolist(self.pixarray.Filename.unique())
        datafiles=glob.glob(directory+'\\*.spe')+glob.glob(directory+'\\sub\\*.spe')
        datafiles=[s.replace('sub\\','') for s in datafiles]
        print(len(datafiles), 'files found')
        missing=[f for f in spefiles if f not in datafiles]  
        if len(missing)!=0:
            try:
                fnums=[int(i.split('.')[1].split('.')[0]) for i in missing]
                franges=[] # ranges of files for missing file output
                # TODO FIX this groupby can have int, str datatype problems
                for key, group in groupby(enumerate(fnums), lambda x: x[0]-x[1]):
                    thisgroup=list(map(itemgetter(1), group))
                    if len(thisgroup)>1:
                        # more than one consecutive so group as min-max in frange
                        franges.append(str(min(thisgroup))+'-'+ str(max(thisgroup)))
                    else:
                        franges.append(str(thisgroup[0])) # single non-consecutive filenumber
                print('Filenumbers ',', '.join(franges),' are missing from data directory.')
            except:
                print('Filenumbers', ', '.join(missing),' are missing from data directory.')
                
    def find_all_peaks(self):
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
        self.shiftmaps=[]
        ''' Parameters for s7d7 derivative peaks:  [0] ampl [1] negval [2] negpeakind
        [3] width (can get pospeak index and val indirectly for use in spectral plots)
        '''
        self.amplmaps=[]
        ''' Parameters for direct integrations: [0] integcnts [1] peak index (eV value
        available indirectly through shiftmaps) [2] countsmax (unsubtracted int at peak)
        [3] slope [4] intercept
        '''
        self.integmaps=[]
    
        for i, [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, 
                idealev, idealind, idealnegpeak, integwidth, chargeshift, peakwidth, 
                searchwidth, kfact, kfact2, mass] in enumerate(self.elemdata):
            print('Extracting maps for element', elem)
            if str(lowind[0])=='nan' and str(highind[0])=='nan':
                lowind, lowrange, highind, hirange, peakind, peakrange= self.fixoddsetup(peakind, peakrange)
            shiftmap=np.empty([self.specimage.shape[0],self.specimage.shape[1], 2]) 
            amplmap=np.empty([self.specimage.shape[0],self.specimage.shape[1], 4])
            integmap=np.empty([self.specimage.shape[0],self.specimage.shape[1], 5])
            for X in range(0,self.specimage.shape[0]):
                for Y in range(0,self.specimage.shape[1]):
                    # Get raw data associated w/ bigpeak (whole region not just peak subset)
                    rawdata=self.specimage[X,Y,lowind[0]:highind[1]+1]
                    if rawdata.max()<0.0001: # data missing.. only zeros for counts
                        shiftmap[X,Y]=np.nan
                        amplmap[X,Y]=np.nan
                        continue
                    s7d7=self.calcderiv(rawdata, peakind, lowind, highind, peakrange, lowrange, hirange)
                    # Find all relative minima
                    foundind=argrelmin(s7d7) # scipy.stat returned as tuple
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
                        try:
                            pospeakval=s7d7[lowlim:uplim].max()
                            # index in terms of original s7d7 dataset
                            pospeakind=np.unravel_index(s7d7[lowlim:uplim].argmax(), s7d7.shape)[0]+lowlim
                            # now calculate amplitude and s7d7 peak width
                            mypeaks=mypeaks.set_value(index,'Ampl',pospeakval-row.negpeakval)
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
                        slope, intercept = self.calcbackground(peakdata, lowind, highind)
                    except:
                        print('Slope fit problem for pixel', str(X), str(Y))
                        # print('Lowind/ highind are', print(lowind), print(highind))
                    # Add subtracted data
                    peakdata['Subdata']=np.nan
                    for index, row in peakdata.iterrows():
                        yval=row.Counts-(slope*row.Energy+intercept)
                        peakdata=peakdata.set_value(index,'Subdata', yval)
                    # Do direct integration over all suspected peaks
                    mypeaks=self.integpeaks(peakdata, mypeaks, slope, intercept, 
                                       idealnegpeak, idealev, integwidth)
                    # Get best peak data from available selections 
                    try:
                        mypeak=self.pickbestderivpeak(mypeaks, idealev, idealnegpeak)
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
            self.amplmaps.append(amplmap)
            self.shiftmaps.append(shiftmap)
            self.integmaps.append(integmap)
        print('Finished with map data extraction.')
                
    def calcderiv(self, rawdata, peakind, lowind, highind, peakrange, lowrange, hirange):
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
            thissect=self.smoothdiffS7D7(rawdata[start:stop+1])
            thiss7d7=np.concatenate((thiss7d7, thissect), axis=0)
        if len(rawdata)!=len(thiss7d7):
            print('Problem... data and deriv of data are different lengths', len(rawdata), len(thiss7d7))
        return thiss7d7

    def smoothdiffS7D7(self, cnts):
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
    
    def pickbestderivpeak(self, mypeaks, idealev, idealnegpeak):
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
    
    def calcbackground(self, peakdata, lowind, highind):
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
    
    def fixoddsetup(self, lowind, lowrange, highind, hirange, peakind, peakrange):
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
    
    def integpeaks(self, peakdata, mypeaks, slope, intercept, idealnegpeak, idealev, integwidth):
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
    
    def find_negpeaks(self, thispos):
        ''' Run negpeaks on active elements (called from GUIrois) 
        Find charging shift associated w/ some large peak (e.g. O) for all pixels 
        return charge shift and peak amplitude (measure of significance) 
        works if scan regions are sufficiently large to allow smooth-diff peak 
        calcshifts performs similar function on direct peaks (but needs mod for big shifts)
        thipos selec
        peakind- index # of this peak in np Z direction; peakrange - associated energy range
        lowind, highind, lowrange, hirange - above and below background region
        REPLACED BY find_all_peaks
        '''
        print('Finding negpeaks using sm-diff spectrum.')
        # TODO combine find_negpeaks (on deriv data) and calcshifts (on direct )
        # into more accurate combined method             
        [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, 
         idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth, kfact, kfact2, mass]=self.elemdata[thispos]
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
        chargeshiftmap=np.empty([self.dim[0],self.dim[1]])
        peakamplmap=np.empty([self.dim[0],self.dim[1]])
        peakintmap=np.empty([self.dim[0],self.dim[1]]) # counts value at maximum
        badcount=0
        for X in range(0,self.dim[0]):
            for Y in range(0,self.dim[1]):
                # Get raw data associated w/ bigpeak (whole region not just peak subset)
                rawdata=self.specimage[X,Y,lowind[0]:highind[1]+1]
                if rawdata.max()<0.0001: # data missing.. only zeros for counts
                    chargeshiftmap[X,Y]=np.nan
                    peakamplmap[X,Y]=np.nan
                    continue
                s7d7=self.calcderiv(rawdata, peakind, lowind, highind, peakrange, lowrange, hirange)
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
                    print('Negpeak location not found. Ampl:', str(negpeakamp),' neg index:', str(negpeak),'X,Y',str(X), str(Y))
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
                    print('Problem with ampl for pixel',str(X), str(Y))
                    print('Size of pospeak search region is',len(posreg))
                    print('Limits are ', str(lowlim), str(uplim))
                    badcount+=1
                # TODO fix problem with peakintmap acquisition
                try:   # also get direct peak intensity near peak
                    posreg=rawdata[lowlim+idealev-idealnegpeak:uplim+idealev-idealnegpeak]
                    peakintmap[X,Y]=posreg.max()
                except:
                    print('Problem with peakint for pixel',str(X), str(Y))
                    print('Size of pospeak search region is',len(posreg))
                    print('Limits are ', str(lowlim), str(uplim))
                    badcount+=1
        # save to list of shift and ampl maps corresponding to this element peak
        # thispos is this element peaks position in existing nan list
        self.shiftmaps[thispos]=chargeshiftmap
        self.amplmaps[thispos]=peakamplmap
        self.peakintmaps[thispos]=peakintmap


    def uniform_filter(self, size=1):
        ''' Uniform filter does adjacent averaging in spatial domain; returns entire
        spectral image but with adjacent pixel spatial averaging 
        other scipy options include: gaussian_filter; convolve2d '''
        self.specsmooth=np.empty([self.dim[0],self.dim[1],self.dim[2]])
        for z in range(0, self.dim[2]):
            slice2d=self.specimage[:,:, z]
            self.specsmooth[:,:,z]=filters.uniform_filter(slice2d, size)
    
    def save_maps(self):
        ''' Save shiftmaps (detected peak shift from ideal position for each element)
        , amplmaps (smooth-diff amplitude  and elemmaps if they exist
        np array is X by Y (spatial dimensions of scan) by number of elements
        (stored in single npy file) '''
        # Save as stack (one layer for each element)
        if len(self.shiftmaps)>0:
            temp=np.zeros((self.dim[0],self.dim[1], 2*len(self.shiftmaps)))
            for i, thismap in enumerate(self.shiftmaps):
                temp[:,:,2*i:2*i+2]=thismap
            np.save(self.directory+'/'+self.uniquename+'_shiftmaps.npy', temp)
            print('Shiftmaps saved')
        if len(self.amplmaps)>0:
            # each element's amplmap has 4 layers (ampl, negpeakval, negpeakind, peakwidth)
            temp=np.zeros((self.dim[0],self.dim[1], 4*len(self.amplmaps)))
            for i, thismap in enumerate(self.amplmaps):
                temp[:,:,4*i:4*i+4]=thismap
            np.save(self.directory+'/'+self.uniquename+'_amplmaps.npy', temp)
            print('Amplmaps saved')
        if len(self.integmaps)>0:
            # Each element's amplmap has 5 layers (integcounts, peakindex, countsmax, slope, intercept)
            temp=np.zeros((self.dim[0],self.dim[1], 5*len(self.integmaps)))
            for i, thismap in enumerate(self.integmaps):
                temp[:,:,5*i:5*i+5]=thismap
            np.save(self.directory+'/'+self.uniquename+'_integmaps.npy', temp)
            print('Integmaps saved')
        if len(self.elemmaps)>0:
            temp=np.zeros((self.dim[0],self.dim[1], len(self.elemmaps)))
            for i, thismap in enumerate(self.elemmaps):
                temp[:,:,i]=thismap
            np.save(self.directory+'/'+self.uniquename+'_elemmaps.npy', temp)
            print('Elemmaps saved')
    
    def save_ampl_images(self):
        ''' Save each element's extracted amplitude as image
        '''
        # Save as stack (one layer for each element)
        for i, amplmap in enumerate(self.amplmaps):
            elem=self.elements[i]
            print('amplmap data of type', type(amplmap))
            im=Image.fromarray(amplmap).convert('RGB') # convert np array to image
            fname=self.uniquename+'_'+elem+'map.jpg'
            im.save(fname)
            # TODO alpha channel problem with jpg... need to flatten or remove alpha
            # scipy.misc.imread  then flatten
            # scipy.misc.imsave(self.uniquename+'_'+elem+'map.png', amplmap)

    def load_maps(self):
        ''' Load saved shiftmaps, amplmaps and integmaps (and elemmaps) from file 
        , amplmaps (smooth-diff amplitude  and elemmaps if they exist
        np array is X by Y (spatial dimensions of scan) by number of elements
        (stored in single npy file) '''
        if len(self.elemdata)==0:
            return
        fname=self.directory+'/'+self.uniquename+'_shiftmaps.npy'
        # each elements shiftmap has two layers ([0]from deriv, [1] from direct)
        if os.path.exists(fname):
            self.shiftmaps=[]
            print('Loading peak shift maps from file')
            temp=np.load(fname) # contains one z layer per element
            for i, elem in enumerate(self.elements):
                self.shiftmaps.append(temp[:,:,2*i:2*i+2])
        else:
            print('Existing peak shift map file',fname,' not found')
        # each element contains [0]amplitude [1]negpeakval [2] negpeakindex
        # [3] peakwidth
        fname=self.directory+'/'+self.uniquename+'_amplmaps.npy'
        if os.path.exists(fname):
            self.amplmaps=[]
            print('Loading sm-diff amplitude maps from file')
            temp=np.load(fname) # contains one z layer per element
            for i, elem in enumerate(self.elements):
                self.amplmaps.append(temp[:,:,4*i:4*i+4])
        else:
            print('existing sm-diff amplitude map file', fname, ' not found')

        fname=self.directory+'/'+self.uniquename+'_integmaps.npy'
        if os.path.exists(fname):
            self.integmaps=[]
            print('Loading integcounts maps from file')
            temp=np.load(fname) # contains one z layer per element
            for i, elem in enumerate(self.elements):
                self.integmaps.append(temp[:,:,5*i:5*i+5])
        else:
            print('existing integcount maps file ', fname, ' not found')
            
        fname=self.directory+'/'+self.uniquename+'_elemmaps.npy'
        if os.path.exists(fname):
            print('Loading element quant maps from file')
            temp=np.load(fname) # contains one z layer per element
            for i, elem in enumerate(self.elements):
                self.elemmaps.append(temp[:,:,i])
        else:
            print('existing element quant maps file not found')

    def extract_spectrum(self, indices):
        ''' create extracted/averaged spectrum from passed selection
        X,Y '''
        # lasso passes back index #s from below flattened list
        print('Extracting spectrum from spectral image')
        # dim is # of rows then # of columns
        # X coord is column number, Y is row number (so reversed)
        xys=[[i,j] for i in range(0,self.dim[0]) for j in range(0,self.dim[1])]
        # Get subset of x,y vals using index passed by lasso 
        selectxys=[xys[i] for i in indices]
        if len(selectxys)==1:
            self.exttype='pixel' # single pixel extraction
        else:
            self.exttype='lasso'
        print('Selection is ', selectxys)
        self.extracted=np.zeros((self.dim[2]))
        speccount=0
        # numpy arrays are indexed as rownum, colnum which is reverse of typical
        # X, Y coords (basically column, row) 
        for col in range(0,self.dim[1]): # X coord is columns
            for row in range(0,self.dim[0]): # Y coord is rows
                if [row,col] in selectxys:
                    speccount+=1
                    self.extracted=np.add(self.extracted,self.specimage[row,col,:])
        if speccount>0:
            self.extracted=np.divide(self.extracted, speccount)
            print('Spectrum extracted from', str(speccount),' pixels.')
        # Make associated smooth-deriv spectrum, do deriv and integ-based quant
        self.find_peaks_extracted()
    
    def save_extracted(self):
        ''' Save average-extracted spectrum to file '''
        ftypes=[('CSV files', '.csv'), ('All files', '*')]
        filename=filedialog.asksaveasfilename(filetypes=ftypes, title='Save extracted spectrum', 
                        defaultextension='csv')
        if not filename:
            return # if cancelled
        extract=pd.DataFrame()
        extract['Energy']=self.energy
        extract['Counts']=self.extracted
        extract.to_csv(filename, index=False)
        
    def find_peaks_extracted(self):
        ''' Deriv and integ peak analysis on extracted spectrum (runs automatically)
        for single pixel these vals could instead be grabbed from existing maps
        but go ahead and redo for single extracted spectrum
        
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
        '''
        numpeaks=3 # max number of negpeaks in deriv spectrum to evaluate

        ''' Parameters for s7d7 derivative peaks:  [0] ampl [1] negval [2] negpeakind
        [3] width (can get pospeak index and val indirectly for use in spectral plots)
        '''
        ''' Parameters for direct integrations: [0] integcnts [1] peak index (eV value
        available indirectly through shiftmaps) [2] countsmax (unsubtracted int at peak)
        [3] slope [4] intercept
        '''
        self.quant = []
        self.quant2 = [] # quant2 (direct integ)
        self.derivparams= []
        self.integparams= []

        fulls7d7=np.zeros(len(self.extracted))
        for i, [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, 
                idealev, idealind, idealnegpeak, integwidth, chargeshift, peakwidth, 
                searchwidth, kfact, kfact2, mass] in enumerate(self.elemdata):
            if str(lowind[0])=='nan' and str(highind[0])=='nan':
                lowind, lowrange, highind, hirange, peakind, peakrange= self.fixoddsetup(peakind, peakrange)
            rawdata=self.extracted[lowind[0]:highind[1]+1]
            s7d7=self.calcderiv(rawdata, peakind, lowind, highind, peakrange, lowrange, hirange)
            fulls7d7[lowind[0]:highind[1]+1]=s7d7 # Write this chunk to full s7d7
            # Find all relative minima
            foundind=argrelmin(s7d7) # scipy.stat returned as tuple
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
                pospeakval=s7d7[lowlim:uplim].max()
                # index in terms of original s7d7 dataset
                pospeakind=np.unravel_index(s7d7[lowlim:uplim].argmax(), s7d7.shape)[0]+lowlim
                # now calculate amplitude and s7d7 peak width
                mypeaks=mypeaks.set_value(index,'Ampl',pospeakval-row.negpeakval)
                mypeaks=mypeaks.set_value(index,'Ampl',pospeakval-row.negpeakval)
                # peakwidth = negpeak- pospeak
                mypeaks=mypeaks.set_value(index,'Peakwidth',row.negpeakind-pospeakind)
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
                slope, intercept = self.calcbackground(peakdata, lowind, highind)
            except:
                print('Slope fit problem for extracted spectrum')
                # print('Lowind/ highind are', print(lowind), print(highind))
            # Add subtracted data
            peakdata['Subdata']=np.nan
            for index, row in peakdata.iterrows():
                yval=row.Counts-(slope*row.Energy+intercept)
                peakdata=peakdata.set_value(index,'Subdata', yval)
            # Do direct integration over all suspected peaks
            mypeaks=self.integpeaks(peakdata, mypeaks, slope, intercept, 
                               idealnegpeak, idealev, integwidth)
            # Get best peak data from available selections 
            try:
                mypeak=self.pickbestderivpeak(mypeaks, idealev, idealnegpeak)
            except:
                print('Problem with choosing best peak for extracted spectrum')
            # skip shift maps... shift will be apparent in spectra plot
            ''' Parameters for s7d7 derivative peaks:  [0] ampl [1] negval [2] negpeakind
            [3] width (can get pospeak index and val indirectly for use in spectral plots)
            '''
            # Save deriv-based quant (elem, ampl, k-factor corrected ampl )
            self.quant.append([elem, mypeak.Ampl, mypeak.Ampl*kfact/mass])
            # pospeak/negpeak energies
            print('Negpeak/pospeak at', int(mypeak.Energy)-int(mypeak.Peakwidth), int(mypeak.Energy))
            xvals=np.array([int(mypeak.Energy)-int(mypeak.Peakwidth), int(mypeak.Energy)])
            # pospeak and negpeak values
            yvals=np.array([ mypeak.Ampl+mypeak.negpeakval, mypeak.negpeakval])
            # Save element, xvals/yvals, and ampl value for deriv plot
            self.derivparams.append([elem, xvals, yvals, mypeak.Ampl])
            ''' Parameters saved in quant2: [0] elem 1 integcnts [2] energy [3] corrcnts
            [4] slope [5] intercept
            '''
            self.quant2.append([elem, int(mypeak.Integcounts), 
                                mypeak.Integcounts*kfact2/mass])
            # elem, peak energy (integration center), integcnts, slope/ intercept of backfit
            print('Direct peak at ', int(mypeak.Energy2),' for ', elem, 'extracted spectrum.')
            self.integparams.append([elem, int(mypeak.Energy2), int(mypeak.Integcounts), 
                                     slope, intercept])
        self.extracts7d7=fulls7d7
        print('Finished with deriv/integ quant on extracted spectrum')

    def smoothdiff_extract(self):
        ''' Piece-wise smooth-diff of extracted spectrum (S7D7)
        for quant calculate s7d7 amplitude and shift 
        REPLACED by find_peaks_extracted
        
        '''
        self.quant=[]  # reset quant list (deriv)
        self.quant2=[] # quant2 (direct integ)
        print('Calculating derivative for extracted spectrum')
        if self.extracted == None:
            print('No spectrum yet extracted from map')
            return
        fulls7d7=np.zeros(len(self.extracted))
        for i, [elem, order, peakind, lowind, highind, peakrange, lowrange, hirange, idealev, idealind, 
         idealnegpeak, integwidth, chargeshift, peakwidth, searchwidth, kfact, kfact2, mass] in enumerate(self.elemdata):
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
            # Get raw data associated w/ bigpeak (whole region not just peak subset)
            rawdata=self.extracted[lowind[0]:highind[1]+1]
            
            s7d7=self.calcderiv(rawdata, peakind, lowind, highind, peakrange, lowrange, hirange)
            fulls7d7[lowind[0]:highind[1]+1]=s7d7 # Write this chunk to full s7d7
            
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
                print('Negpeak location not found. Ampl:', str(negpeakamp),' neg index:', str(negpeak))
                thisshift=np.nan
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
                corrampl=ampl*kfact/mass # adjusted by known k-factor for quant
            except:
                ampl=np.nan
                corrampl=np.nan
                print('Size of pospeak search region is',len(posreg))
                print('Limits are ', str(lowlim), str(uplim))
            # now add amplitude, adjamp, and shift to active quant
            self.quant.append([elem, ampl, thisshift, corrampl])
        self.extracts7d7=fulls7d7
        
