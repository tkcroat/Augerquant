# Auger legacy functions

# old version of scattercompplot before outlier determination/plotting

def plotcntsmajor(Params, areanum):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    calls findevbreaks
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    #numareas=int(Params.Areas)
    #myareas=parseareas(areas, numareas) # set of areas for plotting
    
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max())
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
    colname='Counts'+str(areanum)
    backname='Backfit'+str(areanum)
    
    # find multiplex evbreaks 
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float

    # S region
    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
    if Augerslice.empty==False: # only plottable if data exists in this range          
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region        
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,0])     
        axes[0,0].axvline(x=150.3, color='b') # at ideal S position
               
        for j, val in enumerate(energyvals): # plot multiplex ev breaks in red
            if val > 115 and val < 200: 
                axes[0,0].axvline(x=val, color='r') # on counts plot
                
    # C/Ca region
    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
    if Augerslice.empty==False: # only plottable if data exists in this range       
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,0])         
        axes[1,0].axvline(x=271, color='b') # at ideal C position
        axes[1,0].axvline(x=291, color='b') # at ideal Ca position
        
        # add red vert line at multiplex energy break if present
        for j, val in enumerate(energyvals):
            if val > 225 and val < 320: 
                axes[1,0].axvline(x=val, color='r') # on counts plot
            
    # O regions
    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]     
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region        
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,1])
        axes[0,1].axvline(x=508, color='b') # at ideal O position
        
        for j, val in enumerate(energyvals):
            if val > 475 and val < 535: 
                axes[0,1].axvline(x=val, color='r') # on counts plot
    # Fe region
    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region         
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,1])       
        axes[1,1].axvline(x=595, color='b') # at ideal Fe1 position
        axes[1,1].axvline(x=647.2, color='b') # at ideal Fe2 position
        axes[1,1].axvline(x=702.2, color='b') # at ideal Fe3 position
       
        for j, val in enumerate(energyvals):
            if val > 560 and val < 750: 
                axes[1,1].axvline(x=val, color='r') # on counts plot
    # Mg region
    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region 
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,2])
        axes[0,2].axvline(x=1181.7, color='b') # at ideal Mg position
        
        for j, val in enumerate(energyvals):
            if val > 1140 and val < 1220: 
                axes[0,2].axvline(x=val, color='r') # on counts plot
                 
    # Si region
    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
    if Augerslice.empty==False: # only plottable if data exists in this range 
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,2]) # Si2 region                          
        axes[1,2].axvline(x=1609, color='b') # at ideal Si position  
        
        for j, val in enumerate(energyvals):
            if val > 1550 and val < 1650: 
                axes[1,2].axvline(x=val, color='r') # on counts plot
    return

	# Legacy versions replaced by above generalized functions
def reportpeaksmajor(paramlog, addgauss=True):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass list of files and selected background regions from automated fitting
    '''
    with PdfPages('Peak_report.pdf') as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.iloc[index]['Filename']
            numareas=int(paramlog.iloc[index]['Areas'])            
            Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.iloc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 # 
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
                peakname='Peaks'+str(areanum)
                gaussname='Gauss'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                # S region
                if myplotrange[0] < 135 and  myplotrange[1] > 165:        
                    Augerslice=Augerfile[(Augerfile['Energy']>135) & (Augerfile['Energy']<165)]
                    axes[0,0].axvline(x=150.3, color='b') # at ideal direct peak S position 
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[0,0]) # S region
                        if addgauss==True and gaussname in Augerslice.dtypes.index: # ensure Gaussian fit col exists
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[0,0]) # S region
                # C/Ca region
                if myplotrange[0] < 275 and  myplotrange[1] > 305:         
                    Augerslice=Augerfile[(Augerfile['Energy']>275) & (Augerfile['Energy']<305)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[1,0]) # C/Ca region               
                        if addgauss==True and gaussname in Augerslice.dtypes.index: 
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[1,0]) # C/Ca region                    
                    axes[1,0].axvline(x=291, color='b') # at ideal direct peak Ca position
                # Fe region
                if myplotrange[0] < 625 and  myplotrange[1] > 725:
                    Augerslice=Augerfile[(Augerfile['Energy']>625) & (Augerfile['Energy']<725)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[1,1]) # Fe region
                        if addgauss==True and gaussname in Augerslice.dtypes.index: 
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[1,1])
                    axes[1,1].axvline(x=648, color='b') # at ideal direct peak Fe2 position
                    axes[1,1].axvline(x=702, color='b') # at ideal direct peak Fe3 position
                # Mg region
                if myplotrange[0] < 1165 and  myplotrange[1] > 1200:
                    Augerslice=Augerfile[(Augerfile['Energy']>1165) & (Augerfile['Energy']<1200)]    
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[0,2]) # Mg region
                        if addgauss==True and gaussname in Augerslice.dtypes.index: 
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[0,2])
                    axes[0,2].axvline(x=1185, color='b') # at ideal Mg position
                # Si region
                if myplotrange[0] < 1590 and  myplotrange[1] > 1625:   
                    Augerslice=Augerfile[(Augerfile['Energy']>1590) & (Augerfile['Energy']<1625)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=peakname, ax=axes[1,2]) # Si2 region  
                        if addgauss==True and gaussname in Augerslice.dtypes.index: 
                            Augerslice.plot(x='Energy', y=gaussname, ax=axes[1,2])
                    axes[1,2].axvline(x=1614, color='b') # at ideal Si position  
                pdf.savefig(fig)
                plt.close('all') # close all open figures
    return
    
def reportSDmajor(paramlog, Smdifpeakslog, PDFname='SDplots_report.pdf'):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    '''
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])            
            Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            mypeaks=Smdifpeakslog[(Smdifpeakslog['Filename']==AugerFileName)] # retrieve assoc. subset of peaks data
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 # 
                Peaks=mypeaks[(mypeaks['Areanumber']==areanum)] # then only this areanum
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
                colname='S7D7'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                # S region
                Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]   
                if not Augerslice.empty: # skip entire plot if no data
                    Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>115) & (Peaks['Peakenergy']<200)]
                    if not plotpts.empty: 
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,0], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,0], color='r')                                        
                        titlestring=maketitlestring(plotpts)
                        axes[0,0].set_title(titlestring, fontsize=10)
                        axes[0,0].axvline(x=154, color='b') # at ideal S position
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[0,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[0,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
                    
                    # add red vert line at multiplex energy break if present                    
                    for j, val in enumerate(energyvals):
                        if val > 115 and val < 200: 
                            axes[0,0].axvline(x=val, color='r') # on counts plot
                # C/Ca region
                Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]   
                if not Augerslice.empty: 
                    Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
                    axes[1,0].axvline(x=276, color='b') # at ideal C position
                    axes[1,0].axvline(x=296, color='b') # at ideal Ca position
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>225) & (Peaks['Peakenergy']<320)]
                    if not plotpts.empty:
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,0], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,0], color='r')                                                                                       
                        titlestring=maketitlestring(plotpts)
                        axes[1,0].set_title(titlestring, fontsize=10)        
                       
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[1,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[1,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

                    # add red vert line at multiplex energy break if present
                    for j, val in enumerate(energyvals): 
                        if val > 225 and val < 320: 
                            axes[1,0].axvline(x=val, color='r') # on counts plot
                        
                # O regions
                Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]
                if not Augerslice.empty:
                    Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region
                    axes[0,1].axvline(x=513, color='b') # at ideal O position
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>475) & (Peaks['Peakenergy']<535)]
                    if not plotpts.empty:
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,1], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,1], color='r')                                                              
                        titlestring=maketitlestring(plotpts)
                        axes[0,1].set_title(titlestring, fontsize=10)                
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[0,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[0,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
                    
                    for j, val in enumerate(energyvals):
                        if val > 475 and val < 535: 
                            axes[0,1].axvline(x=val, color='r') # on counts plot
                # Fe region
                Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
                if not Augerslice.empty:
                    Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region   
                    axes[1,1].axvline(x=600, color='b') # at ideal Fe1 position
                    axes[1,1].axvline(x=654, color='b') # at ideal Fe2 position
                    axes[1,1].axvline(x=707, color='b') # at ideal Fe3 position
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>560) & (Peaks['Peakenergy']<750)]
                    if not plotpts.empty:
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,1], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,1], color='r') 
                        titlestring=maketitlestring(plotpts)
                        axes[1,1].set_title(titlestring, fontsize=10)                
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[1,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[1,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

                    for j, val in enumerate(energyvals):
                        if val > 560 and val < 750: 
                            axes[1,1].axvline(x=val, color='r') # on counts plot
                # Mg region
                Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]    
                if not Augerslice.empty:                    
                    Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region
                    axes[0,2].axvline(x=1185, color='b') # at ideal Mg position
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>1140) & (Peaks['Peakenergy']<1220)]
                    if not plotpts.empty: 
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,2], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,2], color='r')                     
                        titlestring=maketitlestring(plotpts)
                        axes[0,2].set_title(titlestring, fontsize=10)        
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[0,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[0,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

                    for j, val in enumerate(energyvals):
                        if val > 1140 and val < 1220: 
                            axes[0,2].axvline(x=val, color='r') # on counts plot
                # Si region
                Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
                if not Augerslice.empty:                    
                    Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region         
                    axes[1,2].axvline(x=1614, color='b') # at ideal Si position 
                    
                    plotpts=Peaks[(Peaks['Peakenergy']>1550) & (Peaks['Peakenergy']<1650)]
                    if not plotpts.empty: 
                        plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,2], color='r')
                        plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,2], color='r') 
                        titlestring=maketitlestring(plotpts)
                        axes[1,2].set_title(titlestring, fontsize=10)
                        noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
                        axes[1,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
                        axes[1,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
                    for j, val in enumerate(energyvals):
                        if val > 1550 and val < 1650: 
                            axes[1,2].axvline(x=val, color='r') # on counts plot
                    pdf.savefig(fig)
                    plt.close(fig) # closes recently displayed figure (since it's saved in PDF report)
    return

def plotcntsmajor(Params, areanum):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    calls findevbreaks
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    #numareas=int(Params.Areas)
    #myareas=parseareas(areas, numareas) # set of areas for plotting
    
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max())
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
    colname='Counts'+str(areanum)
    backname='Backfit'+str(areanum)
    
    # find multiplex evbreaks 
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float

    # S region
    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
    if Augerslice.empty==False: # only plottable if data exists in this range          
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region        
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,0])     
        axes[0,0].axvline(x=150.3, color='b') # at ideal S position
               
        for j, val in enumerate(energyvals): # plot multiplex ev breaks in red
            if val > 115 and val < 200: 
                axes[0,0].axvline(x=val, color='r') # on counts plot
                
    # C/Ca region
    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
    if Augerslice.empty==False: # only plottable if data exists in this range       
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,0])         
        axes[1,0].axvline(x=271, color='b') # at ideal C position
        axes[1,0].axvline(x=291, color='b') # at ideal Ca position
        
        # add red vert line at multiplex energy break if present
        for j, val in enumerate(energyvals):
            if val > 225 and val < 320: 
                axes[1,0].axvline(x=val, color='r') # on counts plot
            
    # O regions
    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]     
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region        
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,1])
        axes[0,1].axvline(x=508, color='b') # at ideal O position
        
        for j, val in enumerate(energyvals):
            if val > 475 and val < 535: 
                axes[0,1].axvline(x=val, color='r') # on counts plot
    # Fe region
    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region         
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,1])       
        axes[1,1].axvline(x=595, color='b') # at ideal Fe1 position
        axes[1,1].axvline(x=647.2, color='b') # at ideal Fe2 position
        axes[1,1].axvline(x=702.2, color='b') # at ideal Fe3 position
       
        for j, val in enumerate(energyvals):
            if val > 560 and val < 750: 
                axes[1,1].axvline(x=val, color='r') # on counts plot
    # Mg region
    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region 
        Augerslice.plot(x='Energy', y=backname, ax=axes[0,2])
        axes[0,2].axvline(x=1181.7, color='b') # at ideal Mg position
        
        for j, val in enumerate(energyvals):
            if val > 1140 and val < 1220: 
                axes[0,2].axvline(x=val, color='r') # on counts plot
                 
    # Si region
    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
    if Augerslice.empty==False: # only plottable if data exists in this range 
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region
        Augerslice.plot(x='Energy', y=backname, ax=axes[1,2]) # Si2 region                          
        axes[1,2].axvline(x=1609, color='b') # at ideal Si position  
        
        for j, val in enumerate(energyvals):
            if val > 1550 and val < 1650: 
                axes[1,2].axvline(x=val, color='r') # on counts plot
    return
  
def plotSDmajor(Params, Peaks, areanum):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    #numareas=int(Params.Areas)
    #myareas=parseareas(areas, numareas) # set of areas for plotting
    
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max())
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
    colname='S7D7'+str(areanum)
    
    # find multiplex evbreaks 
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float

    # S region
    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
    if Augerslice.empty==False: # only plottable if data exists in this range          
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region        
        axes[0,0].axvline(x=154, color='b') # at ideal S position
        
        plotpts=Peaks[(Peaks['Peakenergy']>115) & (Peaks['Peakenergy']<200)] 
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,0], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,0], color='r')       
            titlestring=maketitlestring(plotpts)
            axes[0,0].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals): # plot multiplex ev breaks in red
            if val > 115 and val < 200: 
                axes[0,0].axvline(x=val, color='r') # on counts plot
                
    # C/Ca region
    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
    if Augerslice.empty==False: # only plottable if data exists in this range       
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
        axes[1,0].axvline(x=276, color='b') # at ideal C position
        axes[1,0].axvline(x=296, color='b') # at ideal Ca position
        
        plotpts=Peaks[(Peaks['Peakenergy']>225) & (Peaks['Peakenergy']<320)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,0], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,0], color='r')
            titlestring=maketitlestring(plotpts)
            axes[1,0].set_title(titlestring, fontsize=10)        
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
        # add red vert line at multiplex energy break if present
        for j, val in enumerate(energyvals):
            if val > 225 and val < 320: 
                axes[1,0].axvline(x=val, color='r') # on counts plot
            
    # O regions
    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]     
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region        
        axes[0,1].axvline(x=513, color='b') # at ideal O position
        
        plotpts=Peaks[(Peaks['Peakenergy']>475) & (Peaks['Peakenergy']<535)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,1], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,1], color='r')                    
            titlestring=maketitlestring(plotpts)
            axes[0,1].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals):
            if val > 475 and val < 535: 
                axes[0,1].axvline(x=val, color='r') # on counts plot
    # Fe region
    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region         
        axes[1,1].axvline(x=600, color='b') # at ideal Fe1 position
        axes[1,1].axvline(x=654, color='b') # at ideal Fe2 position
        axes[1,1].axvline(x=707, color='b') # at ideal Fe3 position

        plotpts=Peaks[(Peaks['Peakenergy']>560) & (Peaks['Peakenergy']<750)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,1], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,1], color='r')                   
            titlestring=maketitlestring(plotpts)
            axes[1,1].set_title(titlestring, fontsize=10)                
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
        
        for j, val in enumerate(energyvals):
            if val > 560 and val < 750: 
                axes[1,1].axvline(x=val, color='r') # on counts plot
    # Mg region
    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region        
        axes[0,2].axvline(x=1185, color='b') # at ideal Mg position
        
        plotpts=Peaks[(Peaks['Peakenergy']>1140) & (Peaks['Peakenergy']<1220)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,2], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,2], color='r') 
            titlestring=maketitlestring(plotpts)
            axes[0,2].set_title(titlestring, fontsize=10)               
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
 
        for j, val in enumerate(energyvals):
            if val > 1140 and val < 1220: 
                axes[0,2].axvline(x=val, color='r') # on counts plot
                 
    # Si region
    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
    if Augerslice.empty==False: # only plottable if data exists in this range 
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region                          
        axes[1,2].axvline(x=1614, color='b') # at ideal Si position  
        
        plotpts=Peaks[(Peaks['Peakenergy']>1550) & (Peaks['Peakenergy']<1650)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,2], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,2], color='r')         
            titlestring=maketitlestring(plotpts)
            axes[1,2].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals):
            if val > 1550 and val < 1650: 
                axes[1,2].axvline(x=val, color='r') # on counts plot
    return

def reportcountsmajor(paramlog, Smdifpeakslog, PDFname='countsplot_report.pdf'):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    REDUNDANT... just use reportcountsback with background switch
    '''
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])            
            Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 # 
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
                colname='Counts'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                # S region
                if myplotrange[0] < 115 and  myplotrange[1] > 200:        
                    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region
                    # TODO look up these constants from AESquantparams
                    axes[0,0].axvline(x=154, color='r') # at ideal S position
                    for j, val in enumerate(energyvals):
                        if val > 115 and val < 200: 
                            axes[0,0].axvline(x=val, color='r') # on counts plot
                # C/Ca region
                if myplotrange[0] < 225 and  myplotrange[1] > 320:         
                    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region               
                    axes[1,0].axvline(x=276, color='r') # at ideal C position
                    axes[1,0].axvline(x=296, color='r') # at ideal Ca position
                    # add red vert line at multiplex energy break if present
                    for j, val in enumerate(energyvals):
                        if val > 225 and val < 320: 
                            axes[1,0].axvline(x=val, color='r') # on counts plot
                        
                # O regions
                if myplotrange[0] < 475 and  myplotrange[1] > 535:
                    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]    
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region
                    axes[0,1].axvline(x=513, color='r') # at ideal O position
                    for j, val in enumerate(energyvals):
                        if val > 475 and val < 535: 
                            axes[0,1].axvline(x=val, color='r') # on counts plot
                # Fe region
                if myplotrange[0] < 560 and  myplotrange[1] > 750:
                    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region                              
                    axes[1,1].axvline(x=600, color='r') # at ideal Fe1 position
                    axes[1,1].axvline(x=654, color='r') # at ideal Fe2 position
                    axes[1,1].axvline(x=707, color='r') # at ideal Fe3 position
                    for j, val in enumerate(energyvals):
                        if val > 560 and val < 750: 
                            axes[1,1].axvline(x=val, color='r') # on counts plot
                # Mg region
                if myplotrange[0] < 1140 and  myplotrange[1] > 1220:
                    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]    
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region
                    axes[0,2].axvline(x=1185, color='r') # at ideal Mg position
                    for j, val in enumerate(energyvals):
                        if val > 1140 and val < 1220: 
                            axes[0,2].axvline(x=val, color='r') # on counts plot
                # Si region
                if myplotrange[0] < 1550 and  myplotrange[1] > 1650:   
                    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
                    if not Augerslice.empty:                    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region                  
                    axes[1,2].axvline(x=1614, color='b') # at ideal Si position  
                    for j, val in enumerate(energyvals):
                        if val > 1550 and val < 1650: 
                            axes[1,2].axvline(x=val, color='r') # on counts plot
                pdf.savefig(fig)
                plt.close('all') # close all open figures
       
def plotSDmajor(Params, Peaks, areanum):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass pre-sliced params and peaks
    '''
    AugerFileName=Params.Filename # filename if logmatch is series
    #numareas=int(Params.Areas)
    #myareas=parseareas(areas, numareas) # set of areas for plotting
    
    Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df
    # myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max())
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
    colname='S7D7'+str(areanum)
    
    # find multiplex evbreaks 
    energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float

    # S region
    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
    if Augerslice.empty==False: # only plottable if data exists in this range          
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region        
        axes[0,0].axvline(x=154, color='b') # at ideal S position
        
        plotpts=Peaks[(Peaks['Peakenergy']>115) & (Peaks['Peakenergy']<200)] 
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,0], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,0], color='r')       
            titlestring=maketitlestring(plotpts)
            axes[0,0].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals): # plot multiplex ev breaks in red
            if val > 115 and val < 200: 
                axes[0,0].axvline(x=val, color='r') # on counts plot
                
    # C/Ca region
    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<320)]
    if Augerslice.empty==False: # only plottable if data exists in this range       
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region           
        axes[1,0].axvline(x=276, color='b') # at ideal C position
        axes[1,0].axvline(x=296, color='b') # at ideal Ca position
        
        plotpts=Peaks[(Peaks['Peakenergy']>225) & (Peaks['Peakenergy']<320)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,0], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,0], color='r')
            titlestring=maketitlestring(plotpts)
            axes[1,0].set_title(titlestring, fontsize=10)        
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,0].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,0].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
        # add red vert line at multiplex energy break if present
        for j, val in enumerate(energyvals):
            if val > 225 and val < 320: 
                axes[1,0].axvline(x=val, color='r') # on counts plot
            
    # O regions
    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]     
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region        
        axes[0,1].axvline(x=513, color='b') # at ideal O position
        
        plotpts=Peaks[(Peaks['Peakenergy']>475) & (Peaks['Peakenergy']<535)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,1], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,1], color='r')                    
            titlestring=maketitlestring(plotpts)
            axes[0,1].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals):
            if val > 475 and val < 535: 
                axes[0,1].axvline(x=val, color='r') # on counts plot
    # Fe region
    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region         
        axes[1,1].axvline(x=600, color='b') # at ideal Fe1 position
        axes[1,1].axvline(x=654, color='b') # at ideal Fe2 position
        axes[1,1].axvline(x=707, color='b') # at ideal Fe3 position

        plotpts=Peaks[(Peaks['Peakenergy']>560) & (Peaks['Peakenergy']<750)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,1], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,1], color='r')                   
            titlestring=maketitlestring(plotpts)
            axes[1,1].set_title(titlestring, fontsize=10)                
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,1].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,1].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
        
        for j, val in enumerate(energyvals):
            if val > 560 and val < 750: 
                axes[1,1].axvline(x=val, color='r') # on counts plot
    # Mg region
    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]
    if Augerslice.empty==False: # only plottable if data exists in this range
        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region        
        axes[0,2].axvline(x=1185, color='b') # at ideal Mg position
        
        plotpts=Peaks[(Peaks['Peakenergy']>1140) & (Peaks['Peakenergy']<1220)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[0,2], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[0,2], color='r') 
            titlestring=maketitlestring(plotpts)
            axes[0,2].set_title(titlestring, fontsize=10)               
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[0,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[0,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')
 
        for j, val in enumerate(energyvals):
            if val > 1140 and val < 1220: 
                axes[0,2].axvline(x=val, color='r') # on counts plot
                 
    # Si region
    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
    if Augerslice.empty==False: # only plottable if data exists in this range 
        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region                          
        axes[1,2].axvline(x=1614, color='b') # at ideal Si position  
        
        plotpts=Peaks[(Peaks['Peakenergy']>1550) & (Peaks['Peakenergy']<1650)]
        if plotpts.empty==False: # plot smdif quant points if present
            plotpts.plot.scatter(x='Peakenergy', y='Negintensity', ax=axes[1,2], color='r')
            plotpts.plot.scatter(x='Pospeak', y='Posintensity', ax=axes[1,2], color='r')         
            titlestring=maketitlestring(plotpts)
            axes[1,2].set_title(titlestring, fontsize=10)
            noisebars=makenoisebars(plotpts) # vert lines showing noise ampl at low and high energies (list of two 3-tuples)         
            axes[1,2].vlines(x=noisebars[0][0], ymin=noisebars[0][1], ymax = noisebars[0][2], linewidth=2, color='r')
            axes[1,2].vlines(x=noisebars[1][0], ymin=noisebars[1][1], ymax = noisebars[1][2], linewidth=2, color='r')

        for j, val in enumerate(energyvals):
            if val > 1550 and val < 1650: 
                axes[1,2].axvline(x=val, color='r') # on counts plot
    return

def fitcubic(df, areanum, elem, AugerFileName):
    '''Pass appropriate chunk from Auger spectral dataframe, perform cubic fit
    return chunk with backfit column added '''
    colname='Counts'+str(areanum)
    backfitname='Backfit'+str(areanum)
    xcol=df['Energy']
    ycol=df[colname] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
		cubicfunc=lambda x, a, b, c, d: a*x**3 + b*x**2 + c*x + d # lambda definition of cubic poly
        A,B,C, D=np.polyfit(xcol, ycol, 3)
    except: # deal with common problems with linregress
        print('Fitting error for', elem, ' in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a','n/a') # return all n/a
        return df, fitparams
    fitparams=(A, B, C, D) # tuple to return coeffs of 2nd order poly fit
    for index,row in df.iterrows():
        xval=df.loc[index]['Energy']
        yval= A * xval**3+ B * xval**2 + C * xval + D
        df=df.set_value(index, backfitname, yval)
    return df, fitparams

def fitCapeak(df, areanum, elem, AugerFileName):
    '''Pass appropriate chunk from Auger spectral dataframe, perform linear fit
    return chunk with backfit column added 
	replaced with curve_fit version'''
    colname='Counts'+str(areanum)
    backfitname='Backfit'+str(areanum)
    xcol=df['Energy']
    ycol=df[colname] # Counts1, Counts2 or whatever
    # find relative minimum 
    try:
        A,B,C=np.polyfit(xcol, ycol, 2)
        # diagonal of covariance matrix contains variances for fit params
    except: # deal with common problems with linregress
        print('Fitting error for', elem, ' in file ', AugerFileName)
        fitparams=('n/a','n/a','n/a') # return all n/a
        return df, fitparams
    fitparams=(A, B, C) # tuple to return coeffs of 2nd order poly fit
    for index,row in df.iterrows():
        xval=df.loc[index]['Energy']
        yval= A * xval**2+ B * xval + C
        df=df.set_value(index, backfitname, yval)
    return df, fitparams
	
def scattercompplot(comp1, comp2, elemlist, joinlist=['Sample','Areanum'], basis=False):
    '''Pass two versions of composition calculation (using different lines or whatever) and compare 
    major elements using scatter graphs .. single point for each sample
    uses inner merge to select only subset with values from each df
	use either sample or filenumber'''
    
    elemlist=[re.match('\D+',i).group(0) for i in elemlist] 
    # strip number from peaks like Fe2 if present; columns will be element names (Fe) not peak names (Fe2)
    if basis==False: # use atomic % (which is the default), not basis for each element
        elemlist=['%'+s for s in elemlist]
    numareas=len(elemlist)
    
    # set nrows and ncols for figure of proper size
    cols=divmod(numareas,2)[0]+ divmod(numareas,2)[1]
    if numareas>1:
        rows=2
    else:
        rows=1        
    fig, axes = plt.subplots(nrows=rows, ncols=cols) # axes is array
    # merge dfs with comp1 and comp2 using inner join    
    df=pd.merge(comp1, comp2, how='inner', on=joinlist, suffixes=('','b'))
    mycols=df.dtypes.index # same columns
    outliers=pd.DataFrame(columns=mycols) # empty dataframe for outlying points
    outliers['Element'] # Find outliers from linear fit for each listed element
    for i,elem in enumerate(elemlist):
         # determine which subplot to use
        if (i+1)%2==1:
            rownum=0
        else:
            rownum=1
        colnum=int((i+1)/2.1)       
        xcol=elem
        ycol=elem+'b' # same element from second dataset
        if numareas==1: # deal with single subplot separately
            df.plot.scatter(x=xcol, y=ycol, ax=axes) # single plot axes has no [#,#]
        else:
            df.plot.scatter(x=xcol, y=ycol, ax=axes[rownum,colnum])
            # linear regression: fitting, plot and add labels
        data1=df[elem] # just this data column
        colname=elem+'b'
        data2=df[colname]
        # slope,intercept=np.polyfit(data1, data2, 1)  numpy version
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(data1, data2) # imported from scipy.stats
        # set x range for linear plot 
        text1=str(round(slope,2))+' *x +' + str(round(intercept,2))
        text2='R = ' + str(round(r_value,3)) + ' p = '+str(round(p_value,3))
        xmax=max(max(data1),max(data2))*1.1 # set to slightly larger than max of dataset
        x=np.linspace(0,xmax,100) # setting range for 
        if numareas==1: # deal with single subplot separately
            axes.text(0.025,0.9, text1, fontsize=12, transform=axes.transAxes)
            axes.text(0.025,0.8, text2, fontsize=12, transform=axes.transAxes)
            plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        else: # typical multiarea plot
            axes[rownum,colnum].text(0.025,0.9, text1, fontsize=12, transform=axes[rownum,colnum].transAxes)
            axes[rownum,colnum].text(0.025,0.8, text2, fontsize=12, transform=axes[rownum,colnum].transAxes)
            plt.axes(axes[rownum,colnum]) # set correct axes as active
            plt.plot(x, x*slope+intercept, color='r') # plot appropriate line
        outliers=returnoutliers(df,slope, intercept, elem, basis=False)
    return df


# old version when entire datachunk including evbreaks was passed
def savgol(counts, evbreaks):
    '''Perform python smooth-diff used to guide selection of background regions
    perform this in chunks between evbreaks, works for survey or multiplex
    returns list with smooth-diff columns
    '''
    savgollist=[] # Empty list to hold savgol data
    for i in range(0,len(evbreaks)-1):
        thisreg=counts[evbreaks[i]:evbreaks[i+1]] # slice into separate multiplex regions and process separately
        thisreg=np.asarray(thisreg) # convert list to array
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
        firstvals = thisreg[0] - np.abs( thisreg[1:half_window+1][::-1] - thisreg[0] )
        lastvals = thisreg[-1] + np.abs(thisreg[-half_window-1:-1][::-1] - thisreg[-1])
        thisreg = np.concatenate((firstvals, thisreg, lastvals))
        # Now convolve input signal and sav-gol processing 1D array)
        thisreg=np.convolve( thisreg, m[::-1], mode='valid')
        thisreg=thisreg.tolist() # convert to list 
        savgollist.extend(thisreg) # append to full spectrum list
    while len(savgollist)<len(counts): # can be one element too short 
        savgollist.append(0)
    return savgollist #  returns savitsky-golay smooth diff over same full region 

# old (not generalized) version of counts report iwth backgrounds 
def reportcountsback(paramlog, plotelems, AESquantparams, plotback=True, PDFname='countsback_report.pdf'):
    ''' 2x3 plot of Si, Mg, S, Fe, C/Ca, and O
    pass list of files and selected background regions from automated fitting
    background fits themselves stored with auger csv files
    optional pass of backfitlog (w/ points defining region boundary for background fitting useful for troubleshooting fits)
    evbreaks from multiplex plotted as red lines (often different dwell times for different elements)
    plotback switch -- whether or not to plot 
    '''
    with PdfPages(PDFname) as pdf:
        for index,row in paramlog.iterrows(): # iterrows avoids reindexing problems
            AugerFileName=paramlog.loc[index]['Filename']
            numareas=int(paramlog.loc[index]['Areas'])            
            Augerfile=pd.read_csv(AugerFileName) # reads entire spectra into df (all areas)
            myplotrange=(Augerfile['Energy'].min(),Augerfile['Energy'].max()) # same range for all areas in spe            
            Params=paramlog.loc[index] # grab row for this spe file as Series
            filenumber=Params.Filenumber # retrieve filenumber          
            energyvals=findevbreaks(Params, Augerfile) # get x energy vals for evbreaks in this multiplex (if they exist) as float
            for i in range(0,numareas): # create plot for each area 
                areanum=i+1 # 
                Peaks=mypeaks[(mypeaks['Area']==areanum)] # get subset of peaks for this area number 
                fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9)) # 2 by 3 axes array
                colname='Counts'+str(areanum)
                backname='Backfit'+str(areanum)
                mytitle=str(filenumber)+' area #'+str(areanum)
                plt.suptitle(mytitle)
                # S region
                if myplotrange[0] < 115 and  myplotrange[1] > 200:        
                    Augerslice=Augerfile[(Augerfile['Energy']>115) & (Augerfile['Energy']<200)]
                    if not Augerslice.empty:      
                        # find elemental lines in this keV range
                        thisrange='115-200'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[0,0].axvline(x=elemtuple[0], color='b') # O line
                            axes[0,0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,0]) # S region                        
                        if bflog!=False: # plotting of Auger background fits from integral quant method
                            Augerslice.plot(x='Energy', y=backname, ax=axes[0,0])   
                            plotpts=Peaks[(Peaks['Element']=='S')]
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # fit boundaries for S
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[0,0])
                        for j, val in enumerate(energyvals): # any possible evbreaks in multiplex
                            if val > 115 and val < 200: 
                                axes[0,0].axvline(x=val, color='r') # on counts plot
                # C/Ca region
                if myplotrange[0] < 225 and  myplotrange[1] > 330:         
                    Augerslice=Augerfile[(Augerfile['Energy']>225) & (Augerfile['Energy']<330)]
                    if not Augerslice.empty: 
                        # find elemental lines in this keV range
                        thisrange='225-330'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[1,0].axvline(x=elemtuple[0], color='b') # O line
                            axes[1,0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,0]) # C/Ca region
                        if bflog!=False:
                            Augerslice.plot(x='Energy', y=backname, ax=axes[1,0])
                            plotpts=Peaks[(Peaks['Element']=='Ca')]                        
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # fit boundaries for Ca
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[1,0])
                    # add red vert line at multiplex energy break if present
                    for j, val in enumerate(energyvals):
                        if val > 225 and val < 330: 
                            axes[1,0].axvline(x=val, color='r') # on counts plot
                        
                # O regions (Currently skip background fitting for O )
                if myplotrange[0] < 475 and  myplotrange[1] > 535:
                    Augerslice=Augerfile[(Augerfile['Energy']>475) & (Augerfile['Energy']<535)]    
                    if not Augerslice.empty: 
                        # find elemental lines in this keV range
                        thisrange='475-535'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[1,0].axvline(x=elemtuple[0], color='b') # O line
                            axes[1,0].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,1]) # O region
                        for j, val in enumerate(energyvals):
                            if val > 475 and val < 535: 
                                axes[0,1].axvline(x=val, color='r') # on counts plot
                # Fe region
                if myplotrange[0] < 560 and  myplotrange[1] > 750:
                    Augerslice=Augerfile[(Augerfile['Energy']>560) & (Augerfile['Energy']<750)]
                    if not Augerslice.empty:
                        # find elemental lines in this keV range
                        thisrange='560-750'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[1,1].axvline(x=elemtuple[0], color='b') # O line
                            axes[1,1].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,1]) # Fe region
                        if bflog!=False:
                            Augerslice.plot(x='Energy', y=backname, ax=axes[1,1])
                            plotpts=Peaks[Peaks['Element'].str.contains('Fe', na=False, case=False)] # match all Fe peaks
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # should get values from Fe2 and Fe3                
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[1,1])
                        for j, val in enumerate(energyvals):
                            if val > 560 and val < 750: 
                                axes[1,1].axvline(x=val, color='r') # on counts plot
                # Mg region
                if myplotrange[0] < 1140 and  myplotrange[1] > 1220:
                    Augerslice=Augerfile[(Augerfile['Energy']>1140) & (Augerfile['Energy']<1220)]    
                    if not Augerslice.empty:
                        # find elemental lines in this keV range
                        thisrange='1140-1220'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[0,2].axvline(x=elemtuple[0], color='b') # O line
                            axes[0,2].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[0,2]) # Mg region
                        if bflog!=False:
                            Augerslice.plot(x='Energy', y=backname, ax=axes[0,2])
                            plotpts=Peaks[(Peaks['Element']=='Mg')]
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # should get values from Fe2 and Fe3                
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[0,2])
                        for j, val in enumerate(energyvals):
                            if val > 1140 and val < 1220: 
                                axes[0,2].axvline(x=val, color='r') # on counts plot
                # Si region
                if myplotrange[0] < 1550 and  myplotrange[1] > 1650:   
                    Augerslice=Augerfile[(Augerfile['Energy']>1550) & (Augerfile['Energy']<1650)]
                    if not Augerslice.empty:                        # find elemental lines in this keV range
                        thisrange='1550-1650'
                        elemlines=getelemenergy(plotelems, thisrange, AESquantparams)
                        # list of tuples with energy,elemname
                        for i, elemtuple in enumerate(elemlines):
                            # elemtuple[0] is energy and [1] is element symbol
                            axes[1,2].axvline(x=elemtuple[0], color='b') # O line
                            axes[1,2].text(elemtuple[0],-250, elemtuple[1],rotation=90) # use standard -250 y val    
                        Augerslice.plot(x='Energy', y=colname, ax=axes[1,2]) # Si2 region
                        if bflog!=False:
                            Augerslice.plot(x='Energy', y=backname, ax=axes[1,2])
                            plotpts=Peaks[(Peaks['Element']=='Si')]
                            indexvals=[] # list of single points from Augerfile to plot as different colored scatter
                            for i in range(0,len(plotpts)): # should get values from Fe2 and Fe3                
                                indexvals.append(int(plotpts.iloc[i]['Lower1']))
                                indexvals.append(int(plotpts.iloc[i]['Lower2']))
                                indexvals.append(int(plotpts.iloc[i]['Upper1']))
                                indexvals.append(int(plotpts.iloc[i]['Upper2']))
                            Augerpts=Augerfile.iloc[indexvals] # boundaries of fitted regions 
                            Augerpts.plot.scatter(x='Energy', y=colname, ax=axes[1,2])
                        for j, val in enumerate(energyvals):
                            if val > 1550 and val < 1650: 
                                axes[1,2].axvline(x=val, color='r') # on counts plot
                pdf.savefig(fig)
                plt.close('all') # close all open figures
    return