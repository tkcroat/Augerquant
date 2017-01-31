# -*- coding: utf-8 -*-
"""
Created on Wed May 11 08:08:52 2016

@author: tkc
"""
import re
from collections import defaultdict
import pandas as pd
import numpy as np
import scipy
import scipy.stats
from scipy import optimize
from math import factorial # used by Savgol matrix
from scipy.optimize import curve_fit
#%% 

def organizecolumns(df1,mycols):
    ''' Pass df and template (list of desired columns in desired order) and return reorganized newdf
    '''
    cols1=df1.columns.tolist()
    newdf=df1 # avoids modification of passed df
    uniquelist=[i for i in cols1 if i not in mycols]
    for i,colname in enumerate(uniquelist): # remove cols from df1 that are absent from df2
        # newdf.drop(colname, axis=1, inplace=True) # this modifies both passed and returned dfs
        newdf=newdf.drop(colname, axis=1)
    newdf=newdf[mycols] # reorder columns based on template df
    return newdf
    
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
    
def parseelem2(elemlist, Multielem):
    ''' After multielement peaks removed, also move secondary peaks used as primary to dict (handle separately)
    e.g. splits "S Mg Fe2 Si" into "S Mg Si" and "{Fe,[Fe2]} dictionary; same structure and df output 
    for averaging of Fe, Fe2, or straight Fe2 or straight Fe'''
    # starting elemlist will only have single entries (i.e Ti2 but not Ti & Ti2)
    newelemlist=[] 
    for i, elem in enumerate(elemlist):
        if re.search(r'\d',elem): # has number
            match=re.search(r'\d',elem)
            newkey=elem[0:match.start()]
            # store alt quant (i.e. on Ti2) with same structure as multiple quant (Ti & Ti2)
            # Another entry in multielement list... makes things easier for later quant comparisons
            templist=[] # peakIDs added as list (of length 1)
            templist.append(elem) # list containing single string (keeps identical data structure)
            Multielem.update({newkey:templist}) # add to existing dictionary for separate handling
        else:
            newelemlist.append(elemlist[i]) # just copy over 
    return newelemlist, Multielem # return altered element list and multielem dictionary    
    
def getelemthresholds(elemlist, AESquantparams):
    '''get element-dependent significance thresholds for each peak from AESquantparams
    return dictionary with element and associated significance level''' 
    thresholds={} # returns list of element dependent thresholds for this element set
    for i, elem in enumerate(elemlist):
        # find row in AESquantparams for this element
        thiselemdata=AESquantparams[(AESquantparams['element']==elem)]
        thiselemdata=thiselemdata.squeeze() # series with this elements params
        thresholds.update({elem:thiselemdata.siglevel})        
    return thresholds
    
def cloneparamrows(df):
    ''' Make param log entry for for each areanum - used by calccomposition to correctly process spe files with multiple spatial areas
    passed df is usually list of spe files
    this solves problem that AugerParamLog only has one entry (despite possibly having multiple distinct areas with different spectra'''
    df['Areanumber']=1 # set existing entries as area 1
    mycols=df.dtypes.index
    newrows=pd.DataFrame(columns=mycols) # blank df for new entries
    for index, row in df.iterrows():
        numareas=int(df.loc[index]['Areas'])
        for i in range(2,numareas+1):
            newrow=df.loc[index] # clone this row as series
            newrow=newrow.set_value('Areanumber',i)
            newrows=newrows.append(newrow)
    df=pd.concat([df,newrows], ignore_index=True) # merge new rows with existing ones
    df=df.sort_values(['Filenumber','Areanumber'])
    return df
    
def calccomp(df, Integquantlog, elemlist, AESquantparams):
    '''Calculate elemental composition of given files based on input element list 
    threshold - ratio of element peak to noise peak (0 means no threshold applied
    load element-dependent significance level from AESquantparams'''
    thresholds=getelemthresholds(elemlist, AESquantparams) # Get list of sigma levels for significance/inclusion 
    # thresholds for both single and multipeak
    elemlist, multipeaklist = parseelemlist(elemlist) # list of single peak elements and dict with multipeaks
    # check if any of the single peaks are secondary (i.e. quant on Fe2 not main Fe)
    elemlist, multipeaklist= parseelem2(elemlist, multipeaklist)
    # two element lists needed (elements with one peak and elements with compositions averaged from two peaks i.e. Fe2, Fe3)
    # to process compositions from multiple areas, clone rows from spe log (one for each areanum)
    df=cloneparamrows(df) # splits single entry for 5 spatial area spe into 5 rows with Areanumber 1-5
    df=df.reset_index(drop=True)
    df['AESbasis']=0.0 # resets to zero if already present from calcamplitude
    mycols=['Filenumber', 'Project', 'Filename', 'FilePath', 'Sample', 'Comments','AESbasis','Areanumber']
    for i, elem in enumerate(elemlist):  # add columns for basis
        df[elem]=0.0 # add col for each element to spelist
        df['sig'+elem]=0.0 # copy peak significance (ratio of integrated counts over 1 sigma of background)
        df['err'+elem]=0.0  # another for total error in adjusted counts basis
        mycols.append(elem)
        mycols.append('sig'+elem)
        mycols.append('err'+elem)
    for i,elem in enumerate(list(multipeaklist.keys())): # get elements (keys) from dict
        df[elem]=0.0
        df['sig'+elem]=0.0 
        df['err'+elem]=0.0  
        mycols.append(elem)
        mycols.append('sig'+elem)
        mycols.append('err'+elem)
    for i, elem in enumerate(elemlist):  # now add at.% columns (e.g. %S, %Mg)
        colname='%'+elem # at % columns named %S, %Mg, etc.
        errname='err%'+elem
        mycols.append(colname)  # add to column list template
        mycols.append(errname)
        df[colname]=0.0
        df[errname]=0.0
    for i,elem in enumerate(list(multipeaklist.keys())): # add multipeak elements
        colname='%'+elem # at % columns named %S, %Mg, etc.
        errname='err%'+elem 
        mycols.append(colname)  # add to column list template
        mycols.append(errname)        
        df[colname]=0.0
        df[errname]=0.0
    for i in range(0,len(df)): # loop through all desired spectrum (multiarea ones already have duplicated rows)
        filenum=df.iloc[i]['Filenumber']
        areanum=df.iloc[i]['Areanumber']
        match=Integquantlog[Integquantlog['Filenumber']==filenum] # find integ data for this filenumber        
        match=match[match['Areanumber']==areanum]
        basis=0.0 #
        for j, elem in enumerate(elemlist): # handle the single peak elements            
            temp=match[match['Element']==elem] # finds entry for this element 
            if len(temp)==1:
                # thresholds is dict with required significance level for each element 
                thisthresh=thresholds.get(elem) # sig level for this element
                df=df.set_value(i, 'sig'+elem, temp.iloc[0]['Significance']) # always copy peak significance level
                if temp.iloc[0]['Significance']>thisthresh: # if above set threshold then calculate elem's value and add to basis
                    df=df.set_value(i, elem, temp.iloc[0]['Adjcnts']) # copy adjusted counts of this element
                    df=df.set_value(i, 'err'+elem, temp.iloc[0]['Erradjcnts'])
                    basis+=temp.iloc[0]['Adjcnts'] # add this element's value to AES basis
        # now handle the multipeak elements (get average value from both peaks)
        for key, value in multipeaklist.items(): # key is element (aka colname in df), value is list of peaks in Smdifpeakslog
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # this is number of lines that are average (i.e. 2 for Fe&Fe2)            
            avgval=0.0 # working value for averaged adjamplitude
            erravgval=0.0 # combined error from erradjcnts of each line
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                temp=match[match['Element']==peak] # finds integquantlog entry for this peak (match already trimmed to filenum and area)
                if len(temp)==1:
                    thisthresh=thresholds.get(peak) # sig level for this element/peak
                    df=df.set_value(i, 'sig'+elem, temp.iloc[0]['Significance']) # copy peak significance level
                    if temp.iloc[0]['Significance']>thisthresh:
                        avgval+=temp.iloc[0]['Adjcnts']
                        thiserrperc=temp.iloc[0]['Erradjcnts']/temp.iloc[0]['Adjcnts']**2
                        erravgval+=thiserrperc # sum of square of relative error
                    else:
                        numlines=numlines-1 # if peak is zeroed out and not added, this reduces # peaks in average
            if numlines>0: # avoid divbyzero if peak is too small
                avgval=avgval/numlines # this is now average basis for given element
                erravgval=np.sqrt(erravgval) # sqrt of sum of squares is relative error 
            df=df.set_value(i, key, avgval) # copy adjusted amplitude of this element
            df=df.set_value(i, key+'err', avgval*erravgval) # combined actual error of this elem (as detemined from mulitple lines)
            # add value from this element to AESbasis
            basis+=avgval
        # end of multipeak elements loop
        df=df.set_value(i, 'AESbasis', basis) # write total basis value to df
        # Now compute at.% for each listed element (incl errors)
        for j, elem in enumerate(elemlist):
            colname='%'+elem
            ratio=df.iloc[i][elem]/df.iloc[i]['AESbasis'] # initialized to zero in cases where peak is below significance threshold
            df.set_value(i, colname, ratio)
            temp=match[match['Element']==elem] # again find peak entry and get finds entry for this peak
            # TODO maybe check threshold again (although element's value will be zero)
            if len(temp)==1: 
                thiserr=temp.iloc[0]['Erradjcnts']
                atpercerr=thiserr/df.iloc[i]['AESbasis']
                errname='err%'+elem # error column
                df.set_value(i, errname, atpercerr) # Writes absolute error in at% 
        # Also calculate for elements w/ multiple peaks (if present)
        for key, value in multipeaklist.items(): 
            templist=value # dictionary value is list of elem peak index positions
            numlines=len(templist) # this is number of lines that are average (i.e. 2 for Fe&Fe2)
            colname='%'+key
            ratio=df.iloc[i][key]/df.iloc[i]['AESbasis']
            df.set_value(i, colname, ratio)
            # TODO need to propagate errors through Fe & Fe2
            errlist=[] # list of errors in % (usually max of two)
            for k, peak in enumerate(templist): # create new list with original elem peak from index positions
                temp=match[match['Element']==peak] # finds entry for this peak
                if len(temp)==1:
                    if temp.iloc[0]['Adjcnts']>0: # skip negative values
                        err=temp.iloc[0]['Erradjcnts']/temp.iloc[0]['Adjcnts']
                        errlist.append(err) # add this to list 
            # combine errors in quadrature
            totalerr=0.0
            for j, err in enumerate(errlist):
                totalerr+=err**2
            totalerr=np.sqrt(totalerr) # percent error in at % 
            # now get  actual error
            thisval=df.iloc[i][key] # this is averaged value computed above (possibly zero if below thresholds )
            thiserr=thisval*totalerr # error (in Fe) as actual value based on average of multiple peaks
            atpercerr=thiserr/df.iloc[i]['AESbasis']
            errname='err%'+ key  # error column
            df.set_value(i, errname, atpercerr) # Writes absolute error in at% 
        # end of loop calculation for each spectrum 
                
    # organize data based on mycols template
    df=organizecolumns(df,mycols)
    return df

def calcadjcounts(df, AESquantparams, sig=2, kerrors=True):
    '''For each elemental peak in interquantlog, calculate or recalcuated adjusted counts using k-factor2 and mass
    result stored in adjcnts column and used for subsequent compositional determinations
    can change AESquantresults and recalc at any time; sig (aka 2 sigma errors) is default setting
    kerrors -- include error associated with kfactor (along with Poisson errors)''' 
    if 'Adjcnts' not in df:
        df['Adjcnts']=0.0 # new column for adjusted amplitude (if not already present)
    if 'Erradjcnts' not in df:
        df['Erradjcnts']=0.0 # new column for associated error
    if 'err%cnts' not in df:
        df['err%cnts']=0.0 # percentage error only from counting statistics (not including kfactor err)
    if 'err%total' not in df:
        df['err%total']=0.0 # percentage error only from counting statistics (not including kfactor err)
    # loop for each element, mask df, get appropriate k-factor & mass
    df=df.reset_index(drop=True) # go ahead and reset index
    elemlist=np.ndarray.tolist(df.Element.unique()) # list of unique elements from df
    for i,elem in enumerate(elemlist):
        match=AESquantparams[(AESquantparams['element']==elem)]
        match=match.reset_index(drop=True)
        kfactor2=match.iloc[0]['kfactor2'] # kfactor and mass for this element/peak
        errkf2=match.iloc[0]['errkf2'] # percent error in above for integ method
        mass=match.iloc[0]['mass']      
        elemmask=(df['Element']==elem) # mask for this element in loop 
        
        for j in range(0,len(df)): # loop and set adjamplitude to amp*kfact/mass
            if elemmask[j]==True: # row has this element
                newval=df.iloc[j]['Integcounts']*kfactor2/mass
                percerr=sig/np.sqrt(df.iloc[j]['Integcounts']) # 2/sqrt(N) is percent error
                totalerr=np.sqrt(errkf2**2+percerr**2) # combine in quadrature
                err=newval*totalerr # error value is adjusted counts * 2 sig error percentage
                df=df.set_value(j,'Adjcnts',newval)
                df=df.set_value(j,'err%cnts',percerr)
                df=df.set_value(j,'err%total',totalerr)
                df=df.set_value(j,'Erradjcnts',err)
    return df

''' TESTING
df=lowerfitpeak
'''
    
def makelinebackground(df, areanum, fitparams):
    '''Create linear background under peak region
    passed small slice of Augerfile df just peak region and small adjacent background ''' 
    if fitparams[0]=='n/a': # prior linregresss problem
        return df # return unmodified file
    slope=fitparams[0]
    intercept=fitparams[1]
    backfitname='Backfit'+str(areanum) 
    for index,row in df.iterrows(): # blend between lower line and upper line
        xval=df.loc[index]['Energy']
        yval=slope*xval+intercept
        df=df.set_value(index,backfitname,yval)     
    return df # return same df with interpolated background region added 
    
def makeinterplinebackground(df, areanum, lowerfitparams, upperfitparams):
    '''Create interpolated background from lower and upper peak fits 
    passed small slice of Augerfile df just peak region and small adjacent background ''' 
    # check for n/a values 
    if lowerfitparams[0]=='n/a' or upperfitparams[0]=='n/a': # prior linregresss problem
        return df # return unmodified file
    lowslope=lowerfitparams[0]
    lowintercept=lowerfitparams[1]
    upslope=upperfitparams[0]
    upintercept=upperfitparams[1]
    backfitname='Backfit'+str(areanum) 
    if len(df)>0: # entire region passed should have no vals in backfit (only interpolated region)
        evstep=1/(len(df)+1)   
    else:
        print('Unspecified error in creating background') 
        return
    startrow=df.iloc[0].name # index of first value
    for index,row in df.iterrows(): # blend between lower line and upper line
        xval=df.loc[index]['Energy']
        yval=(1-evstep*(index-startrow))*(lowslope*xval+lowintercept)+evstep*(index-startrow)*(upslope*xval+upintercept)
        df=df.set_value(index,backfitname,yval)     
    return df # return same df with interpolated background region added 

def fitCapeak(df, areanum, elem, AugerFileName):
    '''Pass appropriate chunk from Auger spectral dataframe, perform linear fit
    return chunk with backfit column added '''
    colname='Smcounts'+str(areanum)
    backfitname='Backfit'+str(areanum)
    xcol=df['Energy']
    ycol=df[colname] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
        parabfunc=lambda x, a, b, c: a*x**2 + b*x + c # lambda definition of cubic poly
        fitparams, cov =curve_fit(parabfunc, xcol, ycol) # scipy optimize        
        ss_res=np.dot((ycol-parabfunc(xcol,*fitparams)), (ycol-parabfunc(xcol,*fitparams))) # dot product of diff between data and function
        ymean=np.mean(ycol) # mean of dataset
        ss_tot=np.dot((ycol-ymean),(ycol-ymean)) 
        R2=1-(ss_res/ss_tot) # coeff of determination
        # diagonal of covariance matrix contains variances for fit params
    except: # deal with common problems with linregress
        print('Fitting error for', elem, ' in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a') # return all n/a
        R2='n/a'
        return df, fitparams, R2
    for index,row in df.iterrows():
        xval=df.loc[index]['Energy']
        yval= fitparams[0] * xval**2+ fitparams[1] * xval + fitparams[2]
        df=df.set_value(index, backfitname, yval)
    return df, fitparams, R2

def makeCabackground(df, areanum, fitparams):
    ''' Fill background col of auger spe file with values derived from 2nd order poly fit (pass region under peak
    not fitted by fit Ca peak (which only grabs adjacent background)''' 
    backfitname='Backfit'+str(areanum)
    if len(fitparams)!=3: # prior fitting error already reported via print
        return df
    A=fitparams[0]
    B=fitparams[1]
    C=fitparams[2]
    for index,row in df.iterrows(): # blend between lower line and upper line
        xval=df.loc[index]['Energy']
        yval=A*xval**2+ B* xval +C
        df=df.set_value(index,backfitname,yval)     
    return df
    
def fitcubic(df, areanum, elem, AugerFileName):
    '''Pass appropriate chunk from Auger spectral dataframe, perform cubic fit
    return chunk with backfit column added '''
    colname='Smcounts'+str(areanum) # use smoothed data for background fits
    backfitname='Backfit'+str(areanum)
    xcol=df['Energy']
    ycol=df[colname] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
        cubicfunc=lambda x, a, b, c, d: a*x**3 + b*x**2 + c*x + d # lambda definition of cubic poly
        fitparams, cov =curve_fit(cubicfunc, xcol, ycol) # scipy optimize        
        ss_res=np.dot((ycol-cubicfunc(xcol,*fitparams)), (ycol-cubicfunc(xcol,*fitparams))) # dot product of diff between data and function
        ymean=np.mean(ycol) # mean of dataset
        ss_tot=np.dot((ycol-ymean),(ycol-ymean)) 
        R2=1-(ss_res/ss_tot) # coeff of determination
    # TODO insert special handling for failed fits (some R2 threshold)
    #  Maybe restrictions on curvature 
    except: # deal with failed fit
        print('Fitting error for', elem, ' in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a','n/a') # return all n/a
        return df, fitparams
    for index,row in df.iterrows():
        xval=df.loc[index]['Energy']
        yval= fitparams[0] * xval**3+ fitparams[1] * xval**2 + fitparams[2] * xval + fitparams[3]
        df=df.set_value(index, backfitname, yval)
    return df, fitparams, R2
    
def makecubicbackground(df, areanum, fitparams):
    ''' Fill background col of auger spe file with values derived from 2nd order poly fit (pass region under peak
    not fitted by fit Ca peak (which only grabs adjacent background)''' 
    backfitname='Backfit'+str(areanum)
    if len(fitparams)!=4: # prior fitting error already reported via print
        return df
    A=fitparams[0]
    B=fitparams[1]
    C=fitparams[2]
    D=fitparams[3]
    for index,row in df.iterrows(): # blend between lower line and upper line
        xval=df.loc[index]['Energy']
        yval= A * xval**3+ B * xval**2 + C * xval + D
        df=df.set_value(index,backfitname,yval)     
    return df  


'''
For background fit testing
Augerfile=pd.read_csv('C2010W_18Nov15_12231225.csv')
areanum=1
elem=Elemdata[0][0]
fittype=Elemdata[0][1]
integpeak=Elemdata[0][2]
lower1=Elemdata[0][3]
lower2=Elemdata[0][4]
upper1=Elemdata[0][5]
upper2=Elemdata[0][6]
df=fitregion
Augerfile.to_csv('C2010W_18Nov15_12231225.csv', index=False)
'''
''' TESTING OF BELOW FITS 
plt.plot(xcol,ycol,'b-') # actual data in blue
plt.plot(xcol,gaussian(fitparams, xcol),'r-') # Gaussian fit in red
'''
 
def fitgauss(df, areanum, width, elem, AugerFileName, addgauss=True):
    ''' Gaussian fit of direct peaks (pass Augerfile just around peaks region
    no need to save Gaussian fit, just return width and other params 
    integwidth pass from AESquantparams value'''
    peakname='Peaks'+str(areanum)
    # Remove nan values from peak region
    df=df.dropna(subset=[peakname]) # remove nan entries from peak
    # estimate initial Gaussian parameters from data
    if df.empty: # deal with prior failed background fits (no data in this region after dropna
        print('Gaussian fitting error for', elem, ' peak in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a','n/a') # return all n/a
        rsquared='n/a'
        ier='n/a'
        return df, fitparams, rsquared, ier
    xc=df[peakname].idxmax() # estimate center based on peak max index
    xc=df.loc[xc]['Energy'] # associated energy value near center
    peakarea=df[peakname].sum()  # decent area estimate
    y0=0 #
    params0=[xc,width,peakarea,y0] # initial params list (first guess at gaussian params)
    
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
        fitparams, cov, infodict, mesg, ier =optimize.leastsq(errfunc,params0,args=(xcol,ycol),full_output=True)    
        ss_err=(infodict['fvec']**2).sum()
        ss_tot=((ycol-ycol.mean())**2).sum()
        rsquared=1-(ss_err/ss_tot)
        
    except: # fitting problem 
        print('Gaussian fitting error for', elem, ' peak in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a','n/a') # return all n/a
        rsquared='n/a'
        ier='n/a'
        return df, fitparams, rsquared, ier
    if addgauss==True:
        gaussname="Gauss"+str(areanum)
        df[gaussname]='' # add col for gaussian fit
        for index,row in df.iterrows():
            xval=df.loc[index]['Energy']
            yval=fitparams[3]+fitparams[2]/(fitparams[1]*np.sqrt(2*np.pi))*np.exp(-((xval-fitparams[0])**2/(2*fitparams[1]**2)))
            df.set_value(index,gaussname,yval)
    return df, fitparams, rsquared, ier

''' TESTING
For background fit testing
df=fitregion
Augerfile=pd.read_csv('C2010W_18Nov15_12231225.csv')
areanum=1
elem=Elemdata[1][0]
fittype=Elemdata[1][1]
integpeak=Elemdata[1][2]
lower1=Elemdata[1][3]
lower2=Elemdata[1][4]
upper1=Elemdata[1][5]
upper2=Elemdata[1][6]
integwidth=Elemdata[0][8]
if ier in [1,2,3,4]: print ('true')
'''
def findintegparams(Augerfile, Elements, AESquantparams, Shifts):
    '''Grab integration width and expected counts peak position (also incorporates shift from deriv method)''' 
    halfwidths=[]
    peakcenters=[]
    Energyvals = Augerfile.Energy # for finding index #s corresponding to energy vals for this spectrum
    for i, elem in enumerate(Elements):    
        thiselem=AESquantparams[AESquantparams['Element']==elem]
        if len(thiselem)!=1:
            print('WARNING ... AES quant parameters not found for ', elem)
            halfwidths.append(4) # default integration width
            peakcenters.append('n/a') # 
            return halfwidths, peakcenters
        halfwidths.append(int((thiselem.iloc[0]['integwidth']-1)/2)) # integration uses half-width on either side of center 
        integpeakeV=thiselem.iloc[0]['negpeak']-thiselem.iloc[0]['integpeak']+Shifts[i] # shift of direct peak (defined relative to deriv peak)
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-integpeakeV)) # tuple with index of closest and closest value
        peakcenters.append(temptuple[0]) # first of tuple is closest index # 
    return halfwidths, peakcenters

def integpeaks(Augerfile, Backfitparams, areanum, Elements, Shifts, logmatch, AESquantparams):
    ''' Background fit for each direct peak, shift is list of energy shifts of negpeak (same order as Eledata (opens source spectrum as Augerfile, 
    fits peak backgrounds above and below using Elemdata, saves background to source csv (overwrites existing fits), also saves linear fit params to logdataframe with position/amplitude/etc;
    desired elements out of data range are skipped (in prior findindices function)
    backfitparams is all elements but only this Augerfile
    '''
    #create Smdifpeaks dataframe for temp storage of each peak's params
    Backfitparams=Backfitparams.dropna(subset=['Rval1']) # skip integration/Gaussian fit if background fit failed
        
    AugerFileName=logmatch.Filename # 
    # Create temp df to hold and pass linear fit data
    mycols=['Filenumber', 'Filename', 'Filepath', 'Sample', 'Comments', 'Areanumber', 'Element', 'Integcounts', 
    'Backcounts', 'Significance', 'Xc', 'Width', 'Peakarea', 'Y0','Rsquared','Numchannels']
    Integresults=pd.DataFrame(columns=mycols) # empty df for all integ results for elems in this spe file
    peakname='Peaks'+str(areanum) # this is counts - background (but only calculated in vicinity of known elemental peaks)
    backfitname='Backfit'+str(areanum)
    # global shifts from smdifpeaks and local shift based on smoothed 2nd derivative 
    halfwidths, peakcenters=findintegparams(Augerfile, Elements, AESquantparams, Shifts)
    # loop through and fit all peaks for each element in this spatial area            
    for i, elem in enumerate(Elements):
        if i not in Backfitparams.index: # skips integ calc if backfit is n/a
            continue
        thisbackfit=Backfitparams[Backfitparams['Element']==elem]
        if len(thisbackfit)!=1:
            print('Problem retrieving fit boundaries for ',elem, ' in ', AugerFileName)
            continue
        lower1=thisbackfit.iloc[0]['Lower1']
        upper2=thisbackfit.iloc[0]['Upper2']                
        fitregion=Augerfile[lower1:upper2+1]
        if fitregion.empty==True: # skip if no data present (already should be skipped in Elemdata)
            print('No data present for ', elem, ' in ', AugerFileName)
            continue
        # also need accurate lower/upper bounds ... available from backfitparams
        Integresult=pd.DataFrame(index=np.arange(0,1),columns=mycols) # blank df row for this element
        # get integpeak, kfact, integwidth, siglevel
    
        # addgauss if save of gaussian peak fit in Augerfile is desired
        # Probably could skip Gaussian fitting entirely if peak is weak (check smdiff)
        fitregion, fitparams, rsquared, ier = fitgauss(fitregion, areanum, halfwidths[i], elem, AugerFileName, addgauss=True)
        addgauss=True # maybe pass this arg from elsewhere
        if addgauss==True and ier in [1,2,3,4]: # copy gaussian fit over to csv file if successful
            gaussname="Gauss"+str(areanum)
            if gaussname not in Augerfile.dtypes.index: # add col if not already present                
                Augerfile[gaussname]='' # add col for gaussian fit
            # Copy gaussian fit to Augerfile... fitregion only modified in new Gauss peak fit column
            Augerfile.loc[fitregion.index,fitregion.columns]=fitregion 
        # if gaussian fit is successful set center integration channel to index nearest xc
        # ier flag of 1,2,3,4 if fit succeeds but rsquared threshold is better
        if rsquared!='n/a': # skip integcounts calc but do put 'n/a' entries in df
            if rsquared>0.4:
                xc=fitparams[0] # center of gaussian fit
                center=int(round(xc,0))
                tempdf=fitregion[fitregion['Energy']==center]
                try:
                    centerindex=tempdf[peakname].idxmax() # corresponding index # of peak maximum
                except:
                    print('Gaussian fit center out of data range for ', elem, ' in ', AugerFileName)
                    # use center based on deriv shift and relative offset (index corresponding to integpeakeV)
                    centerindex=peakcenters[i] # backup method of finding center of integration region
            else: # indication of poor Gaussian fit R2<0.4 (use prior knowledge of peak position)
                print('Failed gaussian fit for ', elem, ' in ', AugerFileName)
                # set center integration channel to value passed by integpeak 
                # this is ideal energy value but adjusted by shift found using smooth-diff quant method
                centerindex=peakcenters[i] # already stores index number of central peak (ideal - sm-diff shift value)
                # Still do the counts integration for poor gaussian fits 
            # perform integration over peak center channel + integwidth on either side 
            Augerpeak=Augerfile[centerindex-halfwidths[i]:centerindex+halfwidths[i]+1]
            integcounts=Augerpeak[peakname].sum() # get counts sum 
            backgroundcnts=Augerpeak[backfitname].sum() # sum counts over identical width in background fit
            # Used for peak significance i.e. typically 2 sigma of background integration over identical width
            # full integ width is 1.2*FWHM but integwidth here is closest integer half-width
            # Write fit params from tuple over to Integresult df            
            Integresult.iloc[0]['Integcounts']=integcounts
            Integresult.iloc[0]['Backcounts']=backgroundcnts
            Integresult.iloc[0]['Significance']=round(integcounts/(np.sqrt(backgroundcnts)),3)
        # TODO add 2/sqrt(n) calc of associated percent error (also can calculate later)
        Integresult.iloc[0]['Numchannels']=halfwidths[i]*2+1
        Integresult.iloc[0]['Rsquared']=rsquared
        Integresult.iloc[0]['Element']=elem         
        # These will be n/a if fit fails
        Integresult.iloc[0]['Xc']=fitparams[0]
        Integresult.iloc[0]['Width']=fitparams[1]
        Integresult.iloc[0]['Peakarea']=fitparams[2]
        Integresult.iloc[0]['Y0']=fitparams[3]
        Integresults=Integresults.append(Integresult, ignore_index=True) # add row to list with valid 
        # end of loop through each element 

    # assign params that are common to all areas/all peaks into rows of df (copied from original log)
    for index,row in Integresults.iterrows():
        Integresults.loc[index]['Filenumber']=logmatch.Filenumber   
        Integresults.iloc[index]['Filename']=logmatch.Filename
        Integresults.iloc[index]['Filepath']=logmatch.FilePath
        Integresults.iloc[index]['Sample']=logmatch.Sample
        Integresults.iloc[index]['Comments']=logmatch.Comments
        Integresults.loc[index]['Areanumber']=areanum 
    Integresults=Integresults[mycols] # put back in original order
    return Augerfile, Integresults # df with direct peak fitting info for all areas/ all elements
    
''' TESTING BACKGROUNDS
elem, fittype, integpeak, lower1, lower2, upper1, upper2, kfact, integwidth, siglevel=Elemdata[5]
'''

def fitlinregions(Augerfile, fitbounds, areanum, maxshift, elem, AugerFileName):
    '''Simultaneously deal with linear fits below and above peak of interest
    use residual/outlier detection to trim boundaries of linear fit regions (more robust than deriv threshold style
    return lower and upper slopes/intercepts
    background is either single linear fit or interpolated between the two
    '''
    cntname='Smcounts'+str(areanum) # can use Counts or Smcounts ... maybe smcounts is better for background fit
    lowfitreg=Augerfile[fitbounds[0]:fitbounds[1]+1] # already tested for out-of-range in definefitreg
    upfitreg=Augerfile[fitbounds[2]:fitbounds[3]+1]
    
    # combine above and below regions (without possible added adjacent points)     
    xdata=np.concatenate((lowfitreg['Energy'][maxshift:-maxshift].as_matrix(),upfitreg['Energy'][maxshift:-maxshift].as_matrix()),axis=0)
    ydata=np.concatenate((lowfitreg[cntname][maxshift:-maxshift].as_matrix(),upfitreg[cntname][maxshift:-maxshift].as_matrix()),axis=0)
    
    # linear fit over both background regions
    slope,intercept=np.polyfit(xdata, ydata, 1)
    # manually compute residuals over entire range
    
    xdata=np.concatenate((lowfitreg['Energy'].as_matrix(),upfitreg['Energy'].as_matrix()),axis=0)
    ydata=np.concatenate((lowfitreg[cntname].as_matrix(),upfitreg[cntname].as_matrix()),axis=0)
    
    yfit=slope*xdata+intercept
    resid=np.subtract(ydata,yfit)
    thresh=2*resid.std() # set threshold above which point is removed for having high residual
    backregs=pd.concat([lowfitreg,upfitreg]) # combine both regions 
    backregs['Resid']=resid # paste in calculated residuals
    
    excludelist=[] # clunky way to remove extra backfit points from fit if they exceed threshold (indicating beginning of a peak)
    for index in range(int(lowfitreg.index.min()),int(lowfitreg.index.min())+maxshift): 
        if backregs.loc[index]['Resid']>thresh:
            excludelist.append(index)
        else:
            break
    for index in range(int(lowfitreg.index.max()),int(lowfitreg.index.max())-maxshift,-1):
        if backregs.loc[index]['Resid']>thresh:
            excludelist.append(index)
        else:
            break
    for index in range(int(upfitreg.index.min()),int(upfitreg.index.min())+maxshift): 
        if backregs.loc[index]['Resid']>thresh:
            excludelist.append(index)
        else:
            break
    for index in range(int(lowfitreg.index.max()),int(lowfitreg.index.max())-maxshift,-1):
        if backregs.loc[index]['Resid']>thresh:
            excludelist.append(index)
        else:
            break
    # Need to return modified lower1, lower2 ,upper1, upper2, lowrange, highrange make text string showing extent 
    # now filter dataframe to remove any points indicating presence of peak   
    thismask=lowfitreg.index.isin(excludelist)
    lowfitreg=lowfitreg.loc[~thismask]
    thismask=upfitreg.index.isin(excludelist)
    upfitreg=upfitreg.loc[~thismask]
    fitbounds=[lowfitreg.index.min(), lowfitreg.index.max(),upfitreg.index.min(), upfitreg.index.max()] # return list of refined boundaries
    energybounds=[lowfitreg.Energy.min(), lowfitreg.Energy.max(),upfitreg.Energy.min(), upfitreg.Energy.max()] 
    thismask=backregs.index.isin(excludelist) # 
    backregs=backregs.loc[~thismask]
    # now remove temporary residuals column
    backregs=backregs.drop('Resid', axis=1, inplace=False)

    # Now can refit using expanded data range
    backfitname='Backfit'+str(areanum)
    xcol=backregs['Energy']
    ycol=backregs[cntname] # Counts1, Counts2 or Smcounts1, 2 whatever
    try:
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xcol, ycol)
    except: # deal with common problems with linregress
        print('Fitting error for', elem, ' in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a','n/a','n/a') # return all n/a
        return backregs, fitparams, fitbounds, energybounds
        # return the dataframe unmodified
    fitparams=(slope, intercept, r_value, p_value, std_err) # tuple to return fitting results
    for index,row in backregs.iterrows():
        xval=backregs.loc[index]['Energy']
        yval=slope * xval + intercept
        backregs=backregs.set_value(index, backfitname, yval)
    # fitbounds is index # boundaries of lower and upper fits (list of 4) 
    # energy bounds is energy range over which background was fitted
    return backregs, fitparams, fitbounds, energybounds

def definefitreg(bound1, bound2, maxshift, Augerfile, evbreaks):
    ''' Widen fit region from standard size (based on allowed maxshift) and ensure that multiplex evbreaks are not included in the region
    also make sure one doesn't go into region with no data'''
    lowbound=bound1
    for i in range(bound1, bound1-maxshift-1, -1): # lower by allowed shift but ensure not at data boundary
        if i not in evbreaks and i in Augerfile.index: # also ensure we don't exit boundaries of dataset here 
            lowbound=i
        else:
            break
    upbound=bound2
    for i in range(bound2, bound2+maxshift+1): # lower by allowed shift but ensure not at data boundary
        if i not in evbreaks and i in Augerfile.index:
            upbound=i
        else:
            break
    return lowbound, upbound # this is range of Auger slice that'll be used
    

def fitbackgrounds(Augerfile, areanum, Elements, Shifts, AESquantparams, logmatch):
    ''' takes element strings and element list and returns tuple for each elem symbol containing all params necessary to find each Auger peak from given spe file 
    tuple for integ peak is symbol, ideal peak index #, and integ kfactor''' 
    # integpeak is position of direct peak relative to ideal negative peak in smooth-diff S7D7
    # lower1,lower2 and upper1,upper2 are boundaries of lower and higher energy linear backgroundfit (again energies relative to negpeak)
    # Incorporate shifting of background fit regions into this section from ideal position based on savgol deriv
    
    Energyvals = Augerfile.Energy # 
    evbreaks=logmatch.Evbreaks # needed to ensure fit boundaries don't cross into adjacent element
    tempstring=evbreaks.split('[')[1] # remove brackets from list
    tempstring=tempstring.split(']')[0]
    evbreaks=[int(s) for s in tempstring.split(',')] # convert string to list of break index values
    AugerFileName=logmatch.Filename # 
    mycols=['Filenumber', 'Filename', 'Filepath', 'Sample', 'Comments', 'Date', 'Areanumber', 'Element', 'Lower1', 'Lower2', 'Upper1', 
    'Upper2', 'Lowrange','Highrange','Peakshift', 'Fittype', 'P1','P2','P3','P4','Rval1', 'Pval1', 'Stderr1','Rval2', 'Pval2', 'Stderr2']
    Backfitparams=pd.DataFrame(columns=mycols) # empty df 
    
    for i, elem in enumerate(Elements):
        # find row in AESquantparams for this element
        thiselemdata=AESquantparams[(AESquantparams['element']==elem)]
        thiselemdata=thiselemdata.squeeze() # series with this elements params
        thisshift=Shifts[i] # shift in eV/index # corresponding to this peak from prior smdif quant
        if thisshift=='n/a': # peak not in smdifpeakslog ... usually data out of range
            thisshift=0 # just set shift to zero to avoid problems
        # integ peak position value is relative to negpeak in smooth-diff (i.e. -5 is 5 eV below ideal negpeak)
        integpeakev=thiselemdata.negpeak + thiselemdata.integpeak # ideal energy value of negative Auger peak in smooth-diff spectrum
        lower1ev=thiselemdata.negpeak + thiselemdata.lower1 + thisshift # lower bound of lower energy fit region
        lower2ev=thiselemdata.negpeak + thiselemdata.lower2 + thisshift # upper  bound of lower energy fit region
        upper1ev=thiselemdata.negpeak + thiselemdata.upper1 + thisshift # lower bound of higher energy fit region
        upper2ev=thiselemdata.negpeak + thiselemdata.upper2 + thisshift # upper bound of higher energy fit region
        # width=int(thiselemdata.searchwidth) # search width used to find actual peak in real data
        
        # find index # for ideal neg and pos peaks... use lambda funct.
        #   min(Energyvals, key=lambda x:abs(x-negpeakev)) gives value but not index #
        
        # convert each energy value into index # (global shift already applied)
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-integpeakev)) # tuple with index of closest and closest value
        integpeak=temptuple[0] # first of tuple is index #
        peakinrange=temptuple[1]-integpeakev # should be ~0 if desired peak is in data range
        if abs(peakinrange)>0.5: # Must skip entire desired element here if it's out of range of the data in this particular spe        
            print(elem,' is out of data range for ', AugerFileName)
            continue
        fitbounds=[]
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-lower1ev)) # tuple with index of closest and closest value
        fitbounds.append(temptuple[0]) # first of tuple is index #
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-lower2ev)) # tuple with index of closest and closest value
        fitbounds.append(temptuple[0]) # first of tuple is index #
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-upper1ev)) # tuple with index of closest and closest value
        fitbounds.append(temptuple[0]) # first of tuple is index #
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-upper2ev)) # tuple with index of closest and closest value
        fitbounds.append(temptuple[0]) # first of tuple is index #
    
        maxshift=int(thiselemdata.peakshift) # get allowed max energy shift in channels (normally 1eV/chan)
        
        fittype=thiselemdata.fittype # default type of peak fit for given element
        if fittype=='line': 
            fitbounds[0], fitbounds[1]= definefitreg(fitbounds[0], fitbounds[1], maxshift, Augerfile, evbreaks) # bounds for lower fit region
            fitbounds[2], fitbounds[3]= definefitreg(fitbounds[2], fitbounds[3], maxshift, Augerfile, evbreaks) # bounds for upper fit region
        
            # return fitpeakdf (new background fits), fitparams (slope,intercept, point fit range), R2 val (for tossing vals)
            # Since linear fit may span both, pass both regions and deal with  them simultaneously
            fitpeak, fitparams, fitbounds, energybounds=fitlinregions(Augerfile, fitbounds, areanum, maxshift, elem, AugerFileName)
            if fitparams[2]!='n/a': # holds R2 value and skip for failed fits 
                Augerfile.loc[fitpeak.index,fitpeak.columns]=fitpeak # Copy/save to original file
                # Need to generate values for actual peak region from single linear fit
                thispeak=Augerfile[fitbounds[1]:fitbounds[2]] # parse to get actual determined peak region 
                thispeak=makelinebackground(thispeak, areanum, fitparams)
                Augerfile.loc[thispeak.index,thispeak.columns]=thispeak # copy peak region to source data file

        elif fittype=='Ca': # special treatment
            # find relative minimum if present between C falling edge and Ca peak
            smcountname='Smcounts'+str(areanum)
            minindex=Augerfile[fitbounds[0]:fitbounds[0]+10][smcountname].idxmin() # index value of min left of Ca peak (counts or smoothed counts)
            # minval=Augerfile[lower1:lower1+10][countname].min()
            # maxindex=Augerfile[integpeak-5:integpeak+5][countname].idxmax() # Ca peak index if present
            # maxval=Augerfile[integpeak-5:integpeak+5][countname].max()
            # polynomial fit over two pts at relative min left of peak and small region right of peak
            thispeak=pd.concat([Augerfile[minindex-1:minindex+1],Augerfile[integpeak+10:integpeak+15]])
            # Grab several points on low energy side, 2-3 pts  
            lowevrange=str(round(Augerfile[minindex-1:minindex+1]['Energy'].min(),0))+'-'+ str(round(Augerfile[minindex-1:minindex+1]['Energy'].min(),0)) 
            # Get a few more at upper energy end
            upperevrange=str(round(Augerfile[integpeak+10:integpeak+15]['Energy'].min(),0))+'-'+ str(round(Augerfile[integpeak+10:integpeak+15]['Energy'].max(),0)) 
            thispeak, fitparams, R2 =fitCapeak(thispeak, areanum, elem, AugerFileName) # polynomial fit 
            if R2!='n/a': # only copy successful fits (skip n/a)
                Augerfile.loc[thispeak.index,thispeak.columns]=thispeak # copy over to full spe file
                thispeak=Augerfile[minindex+1:integpeak+11] # actual peak region
                thispeak = makeCabackground(thispeak, areanum, fitparams) # now fill peak region with 2nd order poly background
                Augerfile.loc[thispeak.index,thispeak.columns]=thispeak # copy peak region to source data file
            # Make subtracted peak 
            countname='Counts'+str(areanum)
            peakname='Peaks'+str(areanum)
            backfitname='Backfit'+str(areanum)
            for index in range(fitbounds[1],fitbounds[2]):
                Augerfile.set_value(index, peakname, Augerfile.loc[index][countname]-Augerfile.loc[index][backfitname])   
        else:
            print('Need to write fitting functions for fittype', fittype)
            continue # next in loop to avoid errors below
            
        # Make subtracted peak 
        countname='Counts'+str(areanum)
        peakname='Peaks'+str(areanum)
        backfitname='Backfit'+str(areanum)
        for index in range(fitbounds[1],fitbounds[2]):
            Augerfile.set_value(index, peakname, Augerfile.loc[index][countname]-Augerfile.loc[index][backfitname])
        
        # Integration 
        # create single-rowed dataframe for backfitparams of this element (out-of-range data already skipped)
        Backfitparamrow=pd.DataFrame(index=np.arange(0,1),columns=mycols) 
        # transfer common parameters
        Backfitparamrow.iloc[0]['Areanumber']=areanum
        Backfitparamrow.iloc[0]['Element']=elem            
        Backfitparamrow.iloc[0]['Peakshift']=Shifts[i] # shift of this elem's peak based on derivative method        
        Backfitparamrow.iloc[0]['Filenumber']=logmatch.Filenumber   
        Backfitparamrow.iloc[0]['Filename']=logmatch.Filename
        Backfitparamrow.iloc[0]['Filepath']=logmatch.FilePath
        Backfitparamrow.iloc[0]['Sample']=logmatch.Sample
        Backfitparamrow.iloc[0]['Comments']=logmatch.Comments
        Backfitparamrow.iloc[0]['Date']=logmatch.Date
        Backfitparamrow.iloc[0]['Fittype']=fittype # string with type of background fit to attempt
        
        if fittype=='line':        
            Backfitparamrow.iloc[0]['Lower1']=fitbounds[0] # save boundaries of fit regions
            Backfitparamrow.iloc[0]['Lower2']=fitbounds[1]
            Backfitparamrow.iloc[0]['Upper1']=fitbounds[2]
            Backfitparamrow.iloc[0]['Upper2']=fitbounds[3]
            Backfitparamrow.iloc[0]['Lowrange']=str(energybounds[0])+'-'+str(energybounds[1]) # string with lower fitted eV range
            Backfitparamrow.iloc[0]['Highrange']=str(energybounds[2])+'-'+str(energybounds[3])# string with upper fitted eV range
            Backfitparamrow.iloc[0]['P1']=fitparams[0] # slope for single fit
            Backfitparamrow.iloc[0]['P2']=fitparams[1] # intercept for single fit
            Backfitparamrow.iloc[0]['Rval1']=fitparams[2]
            Backfitparamrow.iloc[0]['Pval1']=fitparams[3]
            Backfitparamrow.iloc[0]['Stderr1']=fitparams[4]
        
        if fittype=='Ca':        
            # copy from lowerfitparams to log df
            Backfitparamrow.iloc[0]['Lowrange']=lowevrange
            Backfitparamrow.iloc[0]['Highrange']=upperevrange
            Backfitparamrow.iloc[0]['P1']=fitparams[0] # A*x2 coeff
            Backfitparamrow.iloc[0]['P2']=fitparams[1] # B*x coeff
            Backfitparamrow.iloc[0]['P3']=fitparams[2] # C coeff
            Backfitparamrow.iloc[0]['Rval1']=R2
        Backfitparams=Backfitparams.append(Backfitparamrow)
    Backfitparams=Backfitparams[mycols]
    return Augerfile, Backfitparams 

def findfitregions(Augerfile, areanum, Elements, Shifts, AESquantparams, logmatch):
    ''' takes element strings and element list and returns tuple for each elem symbol containing all params necessary to find each Auger peak from given spe file 
    tuple for integ peak is symbol, ideal peak index #, and integ kfactor''' 
    # integpeak is position of direct peak relative to ideal negative peak in smooth-diff S7D7
    # lower1,lower2 and upper1,upper2 are boundaries of lower and higher energy linear backgroundfit (again energies relative to negpeak)
    # Incorporate shifting of background fit regions into this section from ideal position based on savgol deriv
    
    Elemdata=[] # returns list of length5 tuples for all elements
    Energyvals = Augerfile.Energy # 
    evbreaks=logmatch.Evbreaks # needed to ensure fit boundaries don't cross into adjacent element
    tempstring=evbreaks.split('[')[1] # remove brackets from list
    tempstring=tempstring.split(']')[0]
    evbreaks=[int(s) for s in tempstring.split(',')] # convert string to list of break index values
    for i, elem in enumerate(Elements):
        # find row in AESquantparams for this element
        thiselemdata=AESquantparams[(AESquantparams['element']==elem)]
        thiselemdata=thiselemdata.squeeze() # series with this elements params
        thisshift=Shifts[i] # shift in eV/index # corresponding to this peak from prior smdif quant
        if thisshift=='n/a': # peak not in smdifpeakslog ... usually data out of range
            thisshift=0 # just set shift to zero to avoid problems
        # integ peak position value is relative to negpeak in smooth-diff (i.e. -5 is 5 eV below ideal negpeak)
        integpeakev=thiselemdata.negpeak + thiselemdata.integpeak # ideal energy value of negative Auger peak in smooth-diff spectrum
        lower1ev=thiselemdata.negpeak + thiselemdata.lower1 + thisshift # lower bound of lower energy fit region
        lower2ev=thiselemdata.negpeak + thiselemdata.lower2 + thisshift # upper  bound of lower energy fit region
        upper1ev=thiselemdata.negpeak + thiselemdata.upper1 + thisshift # lower bound of higher energy fit region
        upper2ev=thiselemdata.negpeak + thiselemdata.upper2 + thisshift # upper bound of higher energy fit region
        # width=int(thiselemdata.searchwidth) # search width used to find actual peak in real data
        
        # find index # for ideal neg and pos peaks... use lambda funct.
        #   min(Energyvals, key=lambda x:abs(x-negpeakev)) gives value but not index #
        
        # convert each energy value into index # (global shift already applied)
        temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-integpeakev)) # tuple with index of closest and closest value
        integpeak=temptuple[0] # first of tuple is index #
        peakinrange=temptuple[1]-integpeakev # should be ~0 if desired peak is in data range
        if abs(peakinrange)<1: # Must skip entire desired element here if it's out of range of the data in this particular spe        
            temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-lower1ev)) # tuple with index of closest and closest value
            lower1=temptuple[0] # first of tuple is index #
            temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-lower2ev)) # tuple with index of closest and closest value
            lower2=temptuple[0] # first of tuple is index #
            temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-upper1ev)) # tuple with index of closest and closest value
            upper1=temptuple[0] # first of tuple is index #
            temptuple=min(enumerate(Energyvals), key=lambda x: abs(x[1]-upper2ev)) # tuple with index of closest and closest value
            upper2=temptuple[0] # first of tuple is index #
        
            shift=int(thiselemdata.peakshift) # get allowed max energy shift in channels (normally 1eV/chan)
            
            tempparams =(lower1, lower2, upper1, upper2, shift) #             
            # Now call for adjustment of lower1,2 and upper 1,2 based on savgol column (single spe, single area, single elem region)
            fitlimits = modfitreg(Augerfile, areanum, tempparams, evbreaks)
            kfact=thiselemdata.kfactor2 # typical sensitivity k-factor associated with element for integration
            siglevel=thiselemdata.siglevel # element dependent threshold for significance (# sigmas above background)
            integwidth=int((thiselemdata.integwidth-1)/2) # integration width in channels for direct integration for this element
            # total # of channels in AESquantparams but include n-1/2 channels on either side of peak center
            fittype=thiselemdata.fittype # default type of peak fit for given element
            #Elemdata is a list (of length number of elements) containing length5 tuples
            elemtuple=(elem, fittype, integpeak, fitlimits[0], fitlimits[1], fitlimits[2], fitlimits[3], kfact, integwidth, siglevel) # add tuple with info for this element     
            Elemdata.append(elemtuple) # now contains proper limits on fitting regions 
        else:
            AugerFileName=logmatch.Filename # logmatch is series
            print('Warning: No quant for ',elem,' for ',AugerFileName, 'data not collected in this energy range.')
    return Elemdata

def findpeakshifts(logmatch, areanum, Smdifpeakslog, Elements):
    ''' Find shifts of negpeak positions for each element in list for single spe file, return as list of floats
    pass series with filename and given area
    ''' 
    # TODO problem if len(Elements)!=len(Shifts) due to couldn't find peak error
    filename=logmatch.Filename # get number from Series
    thispeakslog= Smdifpeakslog[(Smdifpeakslog['Filename']==filename)&(Smdifpeakslog['Areanumber']==areanum)]
    # need to match area number and file number for finding unique shift for this elem
    Shifts=[] # shift in peak position suggested by smdiff quant method
    for i, elem in enumerate(Elements):
        thiselem= thispeakslog[(thispeakslog['PeakID']==elem)]
        if len(thiselem)!=1: # peaks not present should have already been removed 
            print ("Couldn't find ", elem, " peak for area", str(areanum),"of spectrum ", filename)
            Shifts.append('n/a') # keeps len(Elements)== len(Shifts)
        if len(thiselem)==1: # should be match for all peaks that are present
            shift=thiselem.iloc[0]['Shift']
            Shifts.append(shift)
    return Shifts # list of energy shift relative to ideal negpeak for each elemental peak

def makesavgol(df, areanum, evbreaks):
    '''Perform python smooth-diff used to guide selection of background regions
    perform this in chunks between evbreaks (list), works for survey or multiplex, adds col to Augerfile and returns
    evbreaks is list of index #s
    ''' 
    countsname='Counts'+str(areanum)
    # add savgol column (only called if not present)
    savgolname='Savgol'+str(areanum)
    df[savgolname]=0.0  # add/initialize col for 2nd deriv Sav-gol
    # Add 1 to last region boundary to avoid data truncation problem 
    evbreaks[-1]=evbreaks[-1]+1
    for i in range(1,len(evbreaks)): # region 1 to nth region
        thisreg=df.loc[evbreaks[i-1]:evbreaks[i]-1] # slice into separate multiplex regions and process separately        
        thisreg=thisreg[countsname] # convert to Series (keep these index)
        myarr=np.asarray(thisreg) # convert to numpy array
        window_size=11
        deriv=2 
        order=2 # order of savgol fit 
        rate=1
        order_range = range(order+1) # range object
        half_window = (window_size -1) // 2 # type int
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        # b is matrix 3 by window size
        m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv) # series as long as array
        # linalg.pinv gets pseudo-inverse of a matrix (window-sized series)
        # .A of any matrix returns it as ndarray object 
        
        # Pad the signal at the extremes with values taken from the signal itself
        firstvals = myarr[0] - np.abs(myarr[1:half_window+1][::-1] - myarr[0] )
        lastvals = myarr[-1] + np.abs(myarr[-half_window-1:-1][::-1] - myarr[-1])
        myarr= np.concatenate((firstvals, myarr, lastvals))
        # Now convolve input signal and sav-gol processing 1D array .. thisreg is numpy array w/ savgol results
        myarr=np.convolve( myarr, m[::-1], mode='valid')
        
        thisreg.loc[evbreaks[i-1]:evbreaks[i]-1]=myarr # copies numpy array but keeps same indices
        # for loop endpoint is 1 off from df indexing (due to different inclusion rules for last point of range)
        for index in range(evbreaks[i-1],evbreaks[i]): 
            df.set_value(index,savgolname,thisreg.loc[index]) # copy vals from series into entire spe df        
    return df #  returns savitsky-golay smooth diff over same full region 

def integbatchquant(spelist, Smdifpeakslog, AESquantparams, Elements, overwrite=True):
    ''' Batch quantification of all peaks in Elements list and noise amplitude at all chosen background regions (Backregs) 
    returns df with peak positions, amplitudes, width, energy shift, etc. '''   
    # create empty dataframe for storing/passing linear fit params (same structure as in fitbackgrounds)
    mycols=['Filenumber', 'Filename', 'Filepath', 'Sample', 'Comments', 'Date', 'Areanumber', 'Element', 'Lower1', 'Lower2', 'Upper1', 
    'Upper2', 'Lowrange','Highrange','Peakshift', 'Fittype', 'P1','P2','P3','P4','Rval1', 'Pval1', 'Stderr1','Rval2', 'Pval2', 'Stderr2']
    Linearfitlog=pd.DataFrame(columns=mycols)
    # TODO set up log for integration results    
    Integquantlog=pd.DataFrame(columns=['Filenumber', 'Filename', 'Filepath', 'Sample', 'Comments', 'Areanumber', 'Element', 'Integcounts', 
    'Backcounts', 'Significance', 'Xc', 'Width', 'Peakarea', 'Y0','Rsquared','Numchannels'])
    for i in range(0,len(spelist)):
        # get ith row from parameters log for subset of selected spe files (i.e. from spelist)
        logmatch=spelist.iloc[i] #contains row with filename and all other parameters from a given spectra 
        logmatch=logmatch.squeeze() # convert/flatten to Series
        numareas=int(logmatch.Areas) # get # of spatial areas for this spe
        # load Auger spe file of interest here
        AugerFileName=logmatch.Filename  # get Auger filename from Series
        Augerfile=pd.read_csv(AugerFileName) # read entire spectra into df
            
        # now loop through any areas within this spectrum (typically only 1 area)
        for areanum in range(1,numareas+1): # loop over each separate area in spe
            # Now check to ensure this Augerfile has all necessary columns for this area
            # print('Processing area ', areanum)  TESTING
            colname='Counts'+str(areanum)
            if colname not in Augerfile:
                print(colname, ' not present in file ', AugerFileName)
                continue # skip to next area
            backfitname='Backfit'+str(areanum)
            if backfitname not in Augerfile: # add this background fit column if not present
                Augerfile[backfitname]=np.nan
            if overwrite==True: # clear all prior background regions
                Augerfile[backfitname]=np.nan
            savgolname='Savgol'+str(areanum) # Sav-gol 2nd deriv column used to guide selection of fitting regions
            if savgolname not in Augerfile: # returns df with this Savgol column added
                evbreaks=logmatch.Evbreaks # needed for possible savgol smooth-diff
                tempstring=evbreaks.split('[')[1] # remove brackets from list
                tempstring=tempstring.split(']')[0]
                evbreaks=[int(s) for s in tempstring.split(',')] # convert string to list of break index values
                Augerfile=makesavgol(Augerfile, areanum, evbreaks)  # FUNCT pass full spectrum for given area (saved below)   
                
            peakname='Peaks'+str(areanum)
            if peakname not in Augerfile: # add col for subtracted peak data
                Augerfile[peakname]=np.nan
            
            # Get list of negpeak shift for these elements (from Shift column of Smdifpeakslog)
            Shifts=findpeakshifts(logmatch, areanum, Smdifpeakslog, Elements) # single shift val in eV for each elem 
            # Each area has its own Elemdata (selected background fit regions)
            # Elemdata=findfitregions(Augerfile, areanum, Elements, Shifts, AESquantparams, logmatch)
            
            Augerfile, Backfitparams=fitbackgrounds(Augerfile, areanum, Elements, Shifts, AESquantparams, logmatch)
            # Linear background fits above and below plus interpolation between
            # All params from linear fits of pre-peak and post-peak background stored in Backfitparams 
            
            # Peak gaussian fitting and integration subroutine
            Augerfile, Integresults=integpeaks(Augerfile, Backfitparams, areanum, Elements, Shifts, logmatch, AESquantparams) 
            # append linear fit result from this spe/this area to longer master list
            Linearfitlog=Linearfitlog.append(Backfitparams, ignore_index=True)
            Integquantlog=Integquantlog.append(Integresults, ignore_index=True)
        # direct save of modified auger csv with new linear background fits (after all areas processed)
        Augerfile.to_csv(AugerFileName, index=False)
        Linearfitlog=Linearfitlog[mycols] # put back in original order
    return Linearfitlog, Integquantlog

