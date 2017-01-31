# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 08:45:42 2016

@author: tkc
"""

from scipy.optimize import curve_fit
from scipy import optimize
from scipy import asarray as ar,exp
import numpy as np
from numpy import sqrt, log
import matplotlib.pyplot as plt
#%%

''' TESTING OF BELOW FITS 
plt.plot(xcol,ycol,'b-') # actual data in blue
plt.plot(xcol,gaussian(fitparams, xcol),'r-') # Gaussian fit in red
'''
 
def refitCamanual(df):
    
def fitgauss(df, areanum, width, elem, AugerFileName):
    ''' Gaussian fit of direct peaks (pass Augerfile just around peaks region
    no need to save Gaussian fit, just return width and other params 
    integwidth pass from AESquantparams value'''
    peakname='Peaks'+str(areanum)
    # Remove nan values from peak region
    df=df.dropna(subset=[peakname]) # remove nan entries from peak
    # estimate initial Gaussian parameters from data
    xc=df[peakname].idxmax() # estimate center based on peak max index
    xc=df.loc[xc]['Energy'] # associated energy value near center
    area=df[peakname].sum()  # decent area estimate
    y0=0 #
    params0=[xc,width,area,y0] # initial params list 
    
    xcol=df['Energy']
    ycol=df[peakname] # Counts1, Counts2 or whatever
    xcol=xcol.as_matrix() # convert both to numpy matrices
    ycol=ycol.as_matrix()
    
    # define standard gaussian funct (xc, width, area and yoffset are init params)
    gaussian=lambda params, x: params[3]+params[2]/(params[1]*np.sqrt(2*np.pi))*np.exp(-((x-params[0])**2/(2*params[1]**2)))
    
    # thisgauss= gaussian(params0,xcol) 
    errfunc=lambda p, xcol, ycol: ycol- gaussian(p,xcol) # lambda error funct definition
    # sigma2FWHM = lambda sigma: sigma * sqrt(2 * log(2)) * 2 / sqrt(2) # convert Gaussian widths to FWHM?
    
    try:
        fitparams,converged =optimize.leastsq(errfunc,params0,args=(xcol,ycol))
    except: # fitting problem 
        print('Gaussian fitting error for', elem, ' peak in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a') # return all n/a
        return fitparams
    return fitparams, converged

''' TESTING
For background fit testing
df=fitregion
Augerfile=pd.read_csv('C2010W_18Nov15_12231225.csv')
areanum=1
elem=Elemdata[0][0]
fittype=Elemdata[0][1]
integpeak=Elemdata[0][2]
lower1=Elemdata[0][3]
lower2=Elemdata[0][4]
upper1=Elemdata[0][5]
upper2=Elemdata[0][6]

'''
def integpeaks(Augerfile, areanum, Elemdata, Shifts, logmatch):
    ''' Background fit for each direct peak, shift is list of energy shifts of negpeak (same order as Eledata (opens source spectrum as Augerfile, 
    fits peak backgrounds above and below using Elemdata, saves background to source csv (overwrites existing fits), also saves linear fit params to logdataframe with position/amplitude/etc;
    desired elements out of data range are skipped (in prior findindices function)
    '''
    #create Smdifpeaks dataframe for temp storage of each peak's params
    dim=len(Elemdata)# can't write to non-existant df row so set # of rows as numareas*len(Elements)    
    dfrow=0 # keep track of row # for Integparams dataframe
    # Create temp df to hold and pass linear fit data
    AugerFileName=logmatch.Filename # 
    mycols=['Filenumber', 'Filename', 'Filepath', 'Sample', 'Comments', 'Area', 'Element', 'Integcounts', 
    'Xc', 'Width', 'Area', 'Y0','Numchannels']
    Integparams=pd.DataFrame(index=np.arange(0,dim),columns=mycols)
    peakname='Peaks'+str(areanum)
        # all fit regions modify fit region boundaries for this spectrum based on smooth-differentiated peak (2nd deriv, Savgol (poly=2, pts=11))
        # global shifts from smdifpeaks and local shift based on smoothed 2nd derivative 
        # already incorporated into Elemdata values (lower1,2 and upper1,2 fully adjusted)
                
        # loop through and fit all peaks for each element in this spatial area            
    for i, (elem, fittype, integpeak, lower1, lower2, upper1, upper2, kfact, integwidth) in enumerate(Elemdata):
        # linear fit below this elem's peak (shifts and adjustments already made)
        fitregion=Augerfile[lower1:upper2+1]
        fitparams, converged =fitgauss(fitregion, areanum, integwidth, elem, AugerFileName)
        # find center from Gaussian fit (xc)
        if converged!=1: # indication of failed Gaussian fit
            print('Failed gaussian fit for ', elem, ' in ', AugerFileName)
            continue
        xc=fitparams[0] # center of gaussian fit
        center=int(round(xc,0))
        tempdf=fitregion[fitregion['Energy']==center]ww
        centerindex=tempdf[peakname].idxmax() # corresponding index # of peak maximum
        # perform integration over peak center channel + integwidth on either side 
        Augerpeak=Augerfile[centerindex-integwidth:centerindex+integwidth+1]
        integcounts=Augerpeak[peakname].sum() # get counts sum 
        # integ width is half-width at 1.2*half max in eV (channels)

        # Write fit params from tuple over to Integparams df            
        Integparams.iloc[dfrow]['Integcounts']=integcounts
        Integparams.iloc[dfrow]['Element']=elem            
        Integparams.iloc[dfrow]['Xc']=fitparams[0]
        Integparams.iloc[dfrow]['Width']=fitparams[1]
        Integparams.iloc[dfrow]['Area']=fitparams[2]
        Integparams.iloc[dfrow]['Y0']=fitparams[3]
        Integparams.iloc[dfrow]['Numchannels']=integwidth
        dfrow=dfrow+1 # increment to next data row

    # assign params that are common to all areas/all peaks into rows of df (copied from original log)
    for i in range(0,dfrow):
        Integparams.iloc[i]['Filenumber']=logmatch.Filenumber   
        Integparams.iloc[i]['Filename']=logmatch.Filename
        Integparams.iloc[i]['Filepath']=logmatch.FilePath
        Integparams.iloc[i]['Sample']=logmatch.Sample
        Integparams.iloc[i]['Comments']=logmatch.Comments
    Integparams=Integparams[mycols] # put back in original order
    return Augerfile, Integparams # df with direct peak fitting info for all areas/ all elements

