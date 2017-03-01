# -*- coding: utf-8 -*-
"""
Spyder Editor
Auger_batch_header
Designed to read out pertinent header information from all Auger files within a folder.
Output into single log file for import into Excel or elsewhere
"""

#%% Load modules
import os, glob,sys # already run with functions 
import pandas as pd
#import re, struct (only used in sub-functions)
# import csv, fileinput
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Modules')
import Auger_batch_import_functions as Auger
import Auger_utility_functions as AESutils
import Auger_quantmap_functions as QM

#%% Set data directory with Auger data files (via cd in Ipython console, using tinker or using os.chdir)
os.chdir('C:\Temp\Auger')
# from tkinter import filedialog
# datapath = filedialog.askdirectory(initialdir="H:\\Research_data", title = "choose data directory")

#%% Load list of all Auger data files (images, spectral maps and spectra)

filelist=glob.glob('*.sem')+glob.glob('*.map')+glob.glob('*.spe')
filelist=glob.glob('*.spe')

# Find and read Auger logbook (xls file that contains phrase "Auger_logbook")
# this file is to associate sample/project names with Auger file numbers and can be used to combine/average multiple spe files
Augerlogbook=Auger.openorcreatelogbook(filelist)

# Check log file for consistency between xls log and data files in directory
# any errors output to iconsole... i.e. combinable files must all be spe, no missing entries in logbook or missing data files
Auger.checklogfile(filelist, Augerlogbook)

#%% Main file processing loop for spe, sem and map files (header param extraction/ create csv from binaries)
# If file in Augerparamlog.csv it won't be reprocessed.. for reprocessing delete filenumber from that log or delete entire log
AugerParamLog = Auger.Augerbatchimport(filelist, Augerlogbook) # Extract params and process SEM, SPE and MAP files

spelist=AugerParamLog[(AugerParamLog['Areas']>=1)]  # remove image & map files 

# Create jpg images annotated with spatial areas for spe files (assumes existence of .sem image taken just before .spe file)
# easiest to do this before moving combined files (into separate log)
SpatialAreas=pd.read_csv('spatialareaslog.csv') # open automatically created spatial areas log for spe files
Auger.makeannotatedjpg(AugerParamLog, SpatialAreas) # makes annotated jpgs for all spe files with prior sem image

# Manually create single annotated image with spatial areas for selected jpg file and spe file
Auger.annotateone('tr184.209.jpg','tr184.210.csv', SpatialAreas)

# log files are overwritten if they exist by default (but data binaries are not)
AugerParamLog.to_csv('Augerparamlog.csv',index=False) # save all params to new more complete log file (not autosaved)

# Determine stage drift/pixel shift from pre and post images by filenumber
shift, error = Auger.findshift(127, 129, AugerParamLog)

# Combine quantmap spe files (renumbers areas and moves all areas to single combined file with filenumber firstlast.csv)
AugerParamLog=QM.combineQMdata(AugerParamLog,'221-222',QMname='')
# Renumber/rename spatial areas for quant map files
SpatialAreas=pd.read_csv('spatialareaslog.csv') 
SpatialAreas=QM.renumberQMareas(SpatialAreas, '221-222',QMname='') # copy spatial areas for combined QM file, autosaved 

# log files are overwritten if they exist by default (but data binaries are not)
AugerParamLog.to_csv('Augerparamlog.csv',index=False) # save all params to new more complete log file (not autosaved)

#%% Section to average and combine multiple data passes
'''  If starting from here (combine multiple spe), just cd to data directory, reload "Auger_logbook" above,  and then load AugerParamLog
AugerParamLog2=pd.read_csv('Augerparamlog.csv', encoding='cp437')
AugerParamLogsubs=pd.read_csv('sub\\Augerparamlog_subs.csv', encoding='cp437')
AugerParamLog=pd.concat([AugerParamLog,AugerParamLogsubs])
AugerParamLog=AugerParamLog.drop_duplicates(['Filenumber'])
'''
#AVERAGE-COMBINE METHOD 1 (with unique filenumbers for entire project and combinations defined in xls logbook
combinelist=Augerlogbook[(Augerlogbook['Lastnumber']>0)] # gets file ranges to combine via averaging

# Combine consecutive files via averaging, move underlying files to /sub
AugerParamLog=Auger.combinespeloop(combinelist, AugerParamLog, movefiles=True)
AugerParamLog=combinespeloop(combinelist, AugerParamLog, movefiles=True)

AugerParamLog.to_csv('Augerparamlog.csv',index=False) # save all params to new more complete log file

# AVERAGE-COMBINE METHOD 2: Search for files with identical basename, X and Y, combine via averaging, save, create log entry
AugerParamLog=Auger.autocombinespe(AugerParamLog)
AugerParamLog=autocombinespe(AugerParamLog)

# Combine via averaging directly from list of files (different entire spectra with same defined areas)
fullfilelist=glob.glob('*.csv') # Select based on any naming rule (or use any other means to get list of csv filenames to combine)
filelist=fullfilelist[37:40]
# TODO select list of files with tkinter
AugerParamLog=Auger.combinespelist(filelist, AugerParamLog, csvname='', movefiles=False)  # if csvname blank, autonamed based on file#s

# TODO average multiple areas (select by number) within single file

# Check csv files against AugerParamLog
AugerParamLog=AESutils.checkparamlog(AugerParamLog, makeentry=False)

# Assembling larger parameter log files from subdirectories
paramloglist=glob.glob('**/Augerparamlog.csv', recursive=True)
integloglist=glob.glob('**/Integquantlog.csv', recursive=True)
smdifloglist=glob.glob('**/Smdifpeakslog.csv', recursive=True)

Masterparamlog, Masterinteglog, Mastersmdiflog=AESutils.assembledataset(paramloglist, integloglist, smdifloglist)
Masterparamlog, Masterinteglog, Mastersmdiflog=assembledataset(paramloglist, integloglist, smdifloglist)
Masterparamlog.to_csv('Augerparamlog.csv', index=False)
Masterinteglog.to_csv('Integquantlog.csv', index=False)
Mastersmdiflog.to_csv('Smdifpeakslog.csv', index=False)