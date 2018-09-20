# -*- coding: utf-8 -*-
"""
Interactive AES background refitter 
Created on Wed Oct 11 00:44:29 2017

@author: tkc
"""

import tkinter as tk
import tkinter.messagebox as tkmess
from tkinter import filedialog
import matplotlib as mpl # using path, figure, rcParams
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import sys
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant//Modules' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant//Modules')
from AES_data_classes import AESspectrum, AESdataset

PLOT_SIZE = (10,6) # 8, 5 or 
MPL_STYLE = {
    "text.color":"k",
    "axes.labelcolor":"black",
    "axes.edgecolor":"0.4",
    "axes.facecolor":"white",   
    "xtick.color": "lightblue",
    "ytick.color": "lightblue",
    "figure.facecolor":"white",
    "figure.edgecolor":"white",
    "text.usetex":False
}
mpl.rcParams.update(MPL_STYLE)

# TESTING launch_plotter(os.getcwd())
#-----------Misc.-------------#

def launch_plotter(currdir):
    ''' Launcher function for tk quantplotter AES GUI '''
    root = tk.Tk()
    root.wm_title("AES quantplotter")
    GUIMain(root, currdir)
    root.mainloop()
    return

class GUIMain():
    ''' Main container for plotter, options (at right), and fileloader (bottom) 
    pass current working directory as default directory'''
    def __init__(self,root, currdir):
        self.root = root
        self.root.wm_title("AES quant plotter")
        self.top_frame = tk.Frame(self.root)
        self.top_frame.pack(side=tk.TOP) 
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM)
        self.plot_frame = tk.Frame(self.top_frame)
        self.plot_frame.pack(side=tk.LEFT)
        self.refit_frame = tk.Frame(self.top_frame) 
        self.refit_frame .pack(side=tk.LEFT)
        self.loader_frame = tk.Frame(self.bottom_frame)
        self.loader_frame.pack(side=tk.LEFT,fill=tk.BOTH)
        self.plotter  = GUIPlotter(self.plot_frame,self)
        self.opts  = GUIOpts(self.refit_frame,self)
        # Direct call to project loader (asks for directory)
        self.loader  = GUIprojectloader(self.loader_frame,self, currdir)

class GUIPlotter():
    ''' Default plotter should be both direct and deriv (above and below) with 
    linked x scale
    
    '''
    def __init__(self,root, parent):
        self.root = root
        self.parent = parent
        self.figure = mpl.figure.Figure(figsize=PLOT_SIZE, dpi=100)
        self.ax_top = self.figure.add_subplot(211) # two rows, one column
        ''' TODO need to get sharex working
        pkwargs={'sharex':True}
        self.ax_bottom = self.figure.add_subplot(211, **pkwargs) # two rows, one column
        '''
        self.ax_bottom = self.figure.add_subplot(212) # two rows, one column
        self.figure.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
        self.canvas = FigureCanvasTkAgg(self.figure,self.root)
        # just use standard toolbar
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,self.root)
        self.toolbar.update()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.show()
        # instance of AESspectrum now held only in GUIopts

    def plot_both(self,**kwargs):
        ''' plot AESspectrum instance from GUIopts (only active areas)
        direct spectrum in top window and sm-diff (deriv) in bottom
        '''
        # ensure that AESspectrum exists 
        print('Running plot_both from GUIplotter')
        if self.parent.opts.AESspectrum is None:return
        print('Numareas is', self.parent.opts.AESspectrum.numareas)
        self.ax_top.cla()
        self.ax_bottom.cla()
        print('Number of tkareas bools is', len(self.parent.opts.tkareas))
        # Get active areas
        activeareas=[]
        for i in range(0, self.parent.opts.AESspectrum.numareas):
            if self.parent.opts.tkareas[i]:
                activeareas.append(i+1) # areanumber is not zero-based
        print('Active areas are', ','.join([str(i) for i in activeareas]))
        # Plot direct spectra in top
        for i, areanum in enumerate(activeareas):
            colname='Counts'+str(areanum)
            self.ax_top.plot(self.parent.opts.AESspectrum.energy, 
                self.parent.opts.AESspectrum.AESdf[colname], color='k', 
                    picker=True,**kwargs)
        # Plot deriv(s) in top
        for i, areanum in enumerate(activeareas):
            colname='S7D7'+str(areanum)
            self.ax_bottom.plot(self.parent.opts.AESspectrum.energy, 
                self.parent.opts.AESspectrum.AESdf[colname], 
                         color='k', picker=True,**kwargs)
        # Plot direct peak background (if selected in opts)
        if self.parent.opts.showback:
            for i, areanum in enumerate(activeareas):
                colname='Backfit'+str(areanum)
                self.ax_bottom.plot(self.parent.opts.AESspectrum.energy, 
                    self.parent.opts.AESspectrum.AESdf[colname], 
                    color='k', picker=True,**kwargs)
        # Plot smdiff negpeak, pospeak, etc.   
        if self.parent.opts.showsmdifpeaks:
            # filter 
            activepeaks=self.parent.opts.AESspectrum.smdifpeakinfo
            # Filter df for only active areas
            activepeaks=activepeaks['Areanumber'].isin(activeareas)
        self.canvas.show()

    def plot_elems(self, elemparams):
        ''' Add vertical lines showing element positions (toggle?)
        # TODO make this work via toggle
        '''
        if self.AESspectrum is None:return
        self.ax.cla() # clear
        self.plot_both() # Add direct and deriv plots again
        # Add vertical lines at known element energies 
        # (empty elemparams passed to toggle off)
        for i, [elem, energyval] in enumerate(elemparams):
            self.ax_bottom.axvline(x=energyval, color='b', linestyle='dashed', label=elem)
        if len(elemparams)!=0:
            xmax=int(max(map(lambda x:x[1],elemparams))*1.2) # max of 2nd element in list of lists 
            self.ax_bottom.set_xlim([0.1,xmax])
        # Add text labels at approx position
        self.canvas.show()

    def label_quant(self, elems, vals):
        ''' Add quant text label with active elems and at. % values '''
        if self.AESspectrum is None:return
        # Add vertical lines at known element energies
        fullstr=''
        for i, (elem,val) in enumerate(zip(elems, vals)):
            tempstr=r'$%s_{%.0f}$' %(elem, float(val))
            fullstr+=tempstr            
        # transform=ax.transAxes uses coords from 0 to 1 (instead of true x and y vals)
        self.ax_bottom.text(0.05,0.95, fullstr, fontsize=30, verticalalignment='top', transform=self.ax.transAxes)
        self.canvas.show()
        
    # maybe unnecessary below this point    
    def plot_elemranges(self):
        ''' TODO create multipanel element plot as done for AESreports?
        maybe unnecessary 
        '''
        pass
        
class GUIprojectloader():
    ''' Picks directory and loads main Auger param files
    needs current path (which should be set to working data directory) '''
    def __init__(self,root,parent, currdir):
        self.root = root
        self.parent = parent # GUImain is parent
        self.top_frame = tk.Frame(self.root)
        self.top_frame.pack(side=tk.TOP,anchor=tk.W)
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM,fill=tk.BOTH,expand=1)
        tk.Label(self.top_frame,text="Directory:",padx=8,pady=2,height=1).pack(side=tk.LEFT,anchor=tk.W)
        self.directory_entry = tk.Entry(self.top_frame,width=90,bg="lightblue",
                                    fg="black",highlightcolor="lightblue",insertbackground="black",
                                    highlightthickness=2)
        self.directory_entry.pack(side=tk.LEFT,fill=tk.BOTH,expand=1,anchor=tk.W)
        self.directory_entry.insert(0,currdir)
        tk.Button(self.top_frame,text="Browse",command=self.launch_dir_finder
                  ).pack(side=tk.LEFT,fill=tk.BOTH,expand=1,anchor=tk.W)
        self.load_button = tk.Button(self.bottom_frame,text="Load AES project folder",width=60,
                                       command=self.load_AESdataset)
        self.load_button.pack(side=tk.BOTTOM,expand=1,anchor=tk.CENTER)
        # auto call to directory selection and then autoload
        self.load_AESdataset()
        
    def load_AESdataset(self):
        ''' Load standard AES files (paramlogwith data returned to a DataManager '''
        directory = self.directory_entry.get()
        print('Chosen directory is', directory)
        AESdata = AESdataset(directory)
        # pass to GUIoptions and set spectral selector spinbox values
        self.parent.opts.associate_AESdataset(AESdata)        
        # TODO associate AESdataset with GUIplotter (or just selected AESspectrum)
    
    def launch_dir_finder(self):
        directory = filedialog.askdirectory()
        self.directory_entry.delete(0,tk.END)
        self.directory_entry.insert(0,directory)

class GUIOpts():
    ''' Parent is GUImain, manages AESspectrum displayed in GUIplotter
    handles addition/removal of points for background (re)fitting'''
    def __init__(self,root,parent):
        self.root = root
        self.parent = parent
        self.AESdataset = None  # created in GUIprojectloader but associated here
        # Instance of AESspectrum in opts 
        self.AESspectrum = None
        # spinbox 
        self.specselect_frame = tk.Frame(self.root,pady=10)
        self.specselect_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # for display of currently plotted file
        self.currfile_frame = tk.Frame(self.root,pady=10)
        self.currfile_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Spatial area selector checkboxes
        self.areas_frame = tk.Frame(self.root,pady=10)
        self.areas_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Element selector checkboxes
        self.elems_frame = tk.Frame(self.root,pady=10)
        self.elems_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Frame for display of counts/quant results
        self.quant_frame = tk.Frame(self.root,pady=10)
        self.quant_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        self.misc_frame = tk.Frame(self.root,pady=10)
        self.misc_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Simple spinbox for file selection in specselect frame
        self.specspin=tk.Spinbox(self.specselect_frame, command=self.on_specspin_change)
        # TODO does this need config before AESdataset is loaded??
        self.specspin.pack(side=tk.TOP) # throw into specselect sub-frame
        self.tkareas=[] # bools for active/inactive areas (checkboxes)
        self.tkelems=[] # bools list for elem display or quant
        self.activequant=[] # for at % results 
        self.showelems=False # toggle for showing elemental lines on plot
        self.showback=False # add background to direct plot
        self.showsmdifpeaks=False # add smdiff peak params to deriv plot
        
        # Element presets (top of misc frame)
        rowframe=tk.Frame(self.misc_frame)
        tk.Button(rowframe, text='Clear all', command=self.clearall).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Button(rowframe, text='Select all', command=self.selectall).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Button(rowframe, text='Si Fe Mg S Ca O', command=self.elempreset).pack(side=tk.LEFT,fill=tk.X,expand=1)
        rowframe.pack(fill=tk.X, expand=1)
        
        # Permanent buttons in misc_frame
        self.showback_check = tk.Checkbutton(
            self.misc_frame,variable=self.showback, text="Plot background fit?")
        self.showback_check.pack(side=tk.TOP,fill=tk.X,expand=1)

        self.showinfo_check = tk.Checkbutton(
            self.misc_frame,variable=self.showsmdifpeaks, text="Show smdiff peak info?")
        self.showinfo_check.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.update_plot_button = tk.Button(
            self.misc_frame,text="Update plot",command=self.update_plot)
        self.update_plot_button.pack(side=tk.TOP,fill=tk.X,expand=1)

        self.label_button = tk.Button(
            self.misc_frame,text="Label elements",command=self.label_peaks)
        self.label_button.pack(side=tk.TOP,fill=tk.X,expand=1)
        # button to calculate at.% (deriv and integ) and display results
        self.quant_button = tk.Button(
            self.misc_frame,text="Update quant", command=self.calc_atperc)
        self.quant_button.pack(side=tk.TOP,fill=tk.X,expand=1)

        self.quant_button2 = tk.Button(
            self.misc_frame,text="Add quant label", command=self.label_quant)
        self.quant_button2.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.quit_button = tk.Button(
            self.misc_frame, text="Quit", command=self.on_quitapp)
        self.quit_button.pack(side=tk.TOP,fill=tk.X,expand=1)
    
    def update_plot(self):
        ''' Change active areas or selected elements (readback)
        linked to above button 
        plotter uses GUIopts instance of AESspectrum
        includes activeareas lookup
        '''
        self.parent.plotter.plot_both()
        
    def elempreset(self):
        ''' Clear selected elements '''
        if not self.AESspectrum: return
        presets=['S','Fe','Mg','Si','Ca','O']
        for i, elem in enumerate(self.AESdataset.peaks):
            if elem in presets:
                self.tkelems[i].set(1)
            else:
                self.tkelems[i].set(0)

    def display_areas(self):
        ''' Checkbox selectors for spatial areas display in areas_frame (number 
        of areas in spe file from AESspectrum.numareas)
        '''
        # clear previous
        for child in self.areas_frame.winfo_children():
            child.destroy()
        # Write header row into backregs
        tk.Label(self.areas_frame, text='Areas').pack(side=tk.TOP,fill=tk.X,expand=1)

        # tkareas bool variables for show/hide of each area
        self.tkareas=[]
        for i in range(0, self.AESspectrum.numareas):
            self.tkareas.append(tk.BooleanVar())
            self.tkareas[i].set(False)
        self.tkareas[0].set(True) # default show of first only
        # Unfortunately tk/mpl combo requires use of pack (not grid) 
        # display areas in rows of 5
        print(self.AESspectrum.numareas,' areas in spectrum')
        for i in range(0, self.AESspectrum.numareas, 5):
            # associated checkbutton for each quantelem
            arealistframe=tk.Frame(self.areas_frame)
            tk.Checkbutton(arealistframe, variable=self.tkareas[i]).pack(side=tk.LEFT)
            tk.Label(arealistframe, text=str(i+1)).pack(side=tk.LEFT)
            try:
                tk.Checkbutton(arealistframe, variable=self.tkareas[i+1]).pack(side=tk.LEFT)
                tk.Label(arealistframe, text=str(i+2)).pack(side=tk.LEFT)
            except: # out of list range problem
                pass
            try:
                tk.Checkbutton(arealistframe, variable=self.tkareas[i+2]).pack(side=tk.LEFT)
                tk.Label(arealistframe, text=str(i+3)).pack(side=tk.LEFT) # actual area number (not zero based)
            except: # out of list range problem
                pass
            try:
                tk.Checkbutton(arealistframe, variable=self.tkareas[i+3]).pack(side=tk.LEFT)
                tk.Label(arealistframe, text=str(i+4)).pack(side=tk.LEFT)
            except: # out of list range problem
                pass
            try:
                tk.Checkbutton(arealistframe, variable=self.tkareas[i+4]).pack(side=tk.LEFT)
                tk.Label(arealistframe, text=str(i+5)).pack(side=tk.LEFT)
            except: # out of list range problem
                pass
            arealistframe.pack(fill=tk.X, expand=1) # pack after each row
        # area presets
        rowframe=tk.Frame(self.areas_frame)
        tk.Button(rowframe, text='Clear all', command=self.clearallareas).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Button(rowframe, text='Select all', command=self.selectallareas).pack(side=tk.LEFT,fill=tk.X,expand=1)
        rowframe.pack(fill=tk.X, expand=1)
    
    def selectallareas(self):
        ''' Select all areas '''
        for i, tkbool in enumerate(self.tkareas):
            self.tkareas[i].set(1)

    def clearallareas(self):
        ''' Clear all areas '''
        for i, tkbool in enumerate(self.tkareas):
            self.tkareas[i].set(0)
    
    def update_activeareas(self):
        ''' Readback currently selected areas (checkbox vals)
        makes list of active areas by Areanumber 
        # TODO is this necessary .. probably happens automatically
        '''
        self.activeareas=[]
        for i in range(0, self.AESspectrum.numareas):
            if self.tkareas[i]:
                self.activeareas.append(i+1)

    def selectall(self):
        ''' Select all elemental peaks'''
        for i, tkbool in enumerate(self.tkelems):
            self.tkelems[i].set(1)

    def clearall(self):
        ''' Clear all elemental peaks'''
        for i, tkbool in enumerate(self.tkelems):
            self.tkelems[i].set(0)
            
    def display_peaks(self):
        ''' Display available quant elements in elems_frame (self.AESdataset.peaks only 
        contains elements that already have quant results (not other random peaks)
        '''
        for child in self.elems_frame.winfo_children():
            child.destroy()
        # Need to call frame header Label and 

        # Write header row into backregs
        tk.Label(self.elems_frame, text='Available Elements').pack(side=tk.TOP,fill=tk.X,expand=1)
        
        print('Peaks attrib of type', type(self.AESdataset.peaks))
        print('Peaks value is:', self.AESdataset.peaks)

        # tkelems bool variables for active/inactive for each element
        self.tkelems=[]
        for i, elem in enumerate(self.AESdataset.peaks):
            self.tkelems.append(tk.BooleanVar())
            self.tkelems[i].set(True)
        # Unfortunately tk/mpl combo requires use of pack (not grid) 
        for i in range(0, len(self.AESdataset.peaks), 3):
            # associated checkbutton for each quantelem
            elemlistframe=tk.Frame(self.elems_frame)
            tk.Checkbutton(elemlistframe, variable=self.tkelems[i]).pack(side=tk.LEFT)
            tk.Label(elemlistframe, text=self.AESdataset.peaks[i]).pack(side=tk.LEFT)
            try:
                tk.Checkbutton(elemlistframe, variable=self.tkelems[i+1]).pack(side=tk.LEFT)
                tk.Label(elemlistframe, text=self.AESdataset.peaks[i+1]).pack(side=tk.LEFT)
            except: # out of list range problem
                pass
            try:
                tk.Checkbutton(elemlistframe, variable=self.tkelems[i+2]).pack(side=tk.LEFT)
                tk.Label(elemlistframe, text=self.AESdataset.peaks[i+2]).pack(side=tk.LEFT)
            except:
                pass
            elemlistframe.pack(fill=tk.X, expand=1)
            
    def associate_AESdataset(self, AESdataset):
        ''' associate loaded AESdataset with GUIoptions 
        (passed as arg from GUIprojectloader) 
        called by GUIprojectloader '''
        self.AESdataset= AESdataset
        print('AESdataset associated w/ GUIopts has ', len(AESdataset.AESlog),' files.')
        # Set specspin range (using zero-based indexing)
        self.specspin.config(from_=0, to=len(AESdataset.AESlog)-1)
        # clear any existing widgets in backreg frame
        for child in self.elems_frame.winfo_children():
            child.destroy()
        # load first AESspectrum into GUIplotter?
        self.load_AESspectrum(0) # zero-based indexing so row zero
        # auto-plot first spectrum in stack 
        self.parent.plotter.plot_both()
        # display spatial areas
        self.display_areas()
        # loads quant elements into elems frame
        self.display_peaks()
            
    def load_AESspectrum(self, rowindex):
        ''' Load an AESspectrum out of AESdataset using dataframe row (.iloc) 
        first spectrum auto-loads so need to display areas, peaks, etc (and any 
        other attributes looked up by plotter'''
        # Make new instance of AESspectrum class using parent (not AESdataset itself) 
        # and rowindex 
        self.AESspectrum=AESspectrum(self.AESdataset, rowindex)
        # Update displayed filename
        self.display_filename()
        # Need to create tkareas before auto-plot runs
        self.display_areas()
        self.display_peaks()

    def label_quant(self):
        ''' Add a text label with current quant results to plotter 
        launched via button  '''
        elems=[i[0] for i in self.activequant]
        vals=[f[1] for f in self.activequant]
        vals=[int(i) if i>1 else "%0.1f" % i for i in vals]
        self.parent.plotter.label_quant(elems, vals)
        
    def display_filename(self):
        ''' Displays csv name and sample name of currently-active emsa/csv file 
        called after every new load '''
        # clear filename display
        for child in self.currfile_frame.winfo_children():
            child.destroy()
        tk.Label(self.currfile_frame, text=self.AESspectrum.filename).pack()
        tempstr='Sample name: '+self.AESspectrum.sample
        tk.Label(self.currfile_frame, text=tempstr).pack()
        
    def on_specspin_change(self):
        ''' Load and plot chosen file, update smdifpeakinfo and integpeakinfo, etc. 
        '''
        # clear old entries from any prior file 
        for child in self.currfile_frame.winfo_children():
            child.destroy()
        for child in self.elems_frame.winfo_children():
            child.destroy()
        for child in self.quant_frame.winfo_children():
            child.destroy()
        # AESproject file must be loaded or no effect
        self.load_AESspectrum(int(self.specspin.get()))
        # Update displayed fitregions, backfitpts
        self.display_filename()
        persistlist=[] 
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                persistlist.append(True)
            else:
                persistlist.append(False)
        self.display_areas() # recreates area checkboxes
        self.display_peaks() # recreates tkelems bools list
        # reselect the same elements
        if len(persistlist)==len(self.tkelems):
            for i, thisbool in enumerate(persistlist):
                self.tkelems[i].set(thisbool)
        # update activeareas 

        # update plot in GUIplotter
        self.parent.plotter.plot_both()
        # Also rerun label_peaks (lines stay on if set to on)
        if self.showelems: # needs double toggle to stay the same
            self.showelems=False
        else:
            self.showelems=True
        self.label_peaks()
    
    def on_quitapp(self):
        msg = "Quitting:\nUnsaved progress will be lost.\nDo you wish to Continue?"
        if tkmess.askokcancel("AES refitter",msg):
            self.parent.root.destroy()
    
    def label_peaks(self):
        ''' Get currently selected elements from elems frame and associated energy
        pass to plotter.label_peaks
        toggle style button
        '''
        if self.showelems:
            self.showelems=False
        else:
            self.showelems=True
        peakparams=[] # list of [peak, energyval]
        # active (checked) tk elements will automatically be updated on check, right?
        if self.showelems:
            for i, tkbool in enumerate(self.tkelems):
                if tkbool.get():
                    match=self.AESdataset.AESquantparams[self.AESdataset.AESquantparams['element']==self.AESdataset.peaks[i]]
                    ev=match.iloc[0]['energy']
                    peakparams.append([self.AESdataset.peaks[i],ev])
        # now pass info to plotter (can be empty list)
        self.parent.plotter.plot_elems(peakparams)
    
    def calc_atperc(self):
        ''' Generate at % (without error) for subset of selected elements
        then display quant ... linked to quant button 
        display quant only for lowest areanumber (if many are selected)'''
        print('Running calc atperc')
        # Find lowest active area number
        for i, mybool in enumerate(self.tkareas):
            if mybool:
                lowareanum=i+1
                break
        # Data subset for lowest areanumber
        smpeaks=self.AESspectrum.smdifpeakinfo[self.AESspectrum.smdifpeakinfo['Areanumber']==lowareanum]
        # clear current values for smdiff deriv method
        self.activequant=[] # list of lists w/ [peak symbol, adjampl, atperc] (no decent error)
        AESbasis=0.0
        # Now filter by active/selected peaks (bools for )
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                thiselem=smpeaks[smpeaks['PeakID']==self.AESdataset.peaks[i]]
                # Check but this shouldn't happen
                if len(thiselem)!=1:
                    print('Problem finding unique smdiff entry for area', lowareanum,
                          ' and ', self.AESdataset.peaks[i])
                    return
                row=thiselem.iloc[0]
                self.activequant.append([self.AESdataset.peaks[i], row.Adjamp, 0.0])
                if row.Adjamp>0:
                   AESbasis+=row.Adjamp
                else: # Exclude negative values but list entry needed
                    pass
        # Now calculate at. % 
        for i, [symb, adjampl, atperc] in enumerate(self.activequant):
            self.activequant[i][2]=100*adjampl/AESbasis
        # Now update at. % calculation for integ method
        '''
        contains: adjcnts, erradjcnt, atperc, erratperc
        erradjcnts is actual (not fractional) error and includes k-factor error (which
        is almost always the dominant error source)
        '''
        # Data subset for lowest areanumber
        integpeaks=self.AESspectrum.integpeakinfo[self.AESspectrum.integpeakinfo['Areanumber']==lowareanum]
         
        # clear current values
        self.activequant2=[] # list w/ [elem symb, adjcnts, erradjcnt, atperc, erratperc]

        # Error % calc by Errcorrcnts/corrcnts (integlog has actual error in adjcnts)
        temperrperc=[] # temporary fractional error percent (for calc
        AESbasis=0.0
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                thiselem=integpeaks[integpeaks['Element']==self.AESdataset.peaks[i]]
                if len(thiselem)!=1:
                    print('Problem finding unique integquant entry for area', lowareanum,
                          ' and ', self.AESdataset.peaks[i])
                    return
                row=thiselem.iloc[0]
                self.activequant2.append([self.AESdataset.peaks[i], 
                    row.Adjcnts, row.Erradjcnts, 0.0, 0.0])
                if row.Adjcnts>0:
                   temperrperc.append(row.Erradjcnts/row.Adjcnts)
                   AESbasis+=row.Adjcnts
                else: # Exclude negative values but list entry needed
                    temperrperc.append(0.0) # percent error (intermediate calculation)
        for i, [elem, adjcnts, erradj, atper, erratper] in enumerate(self.activequant2):
            self.activequant2[i][3]=100*row.Adjcnts/AESbasis
            # Fractional error in at.% needed for display
            self.activequant2[i][4]=(100*row.Adjcnts/AESbasis)*temperrperc[i]
        self.display_quant() # run to display this
        
    def display_quant(self):
        ''' Need to display sm-diff and integ compositions for active elements
        entries in activequant (for smdiff) and activequant_integ should be 
        corresponding lists (recalc from above)
        '''
        # Clear any existing widgets in backreg frame
        for child in self.quant_frame.winfo_children():
            child.destroy()
        # Sort active quant elements by at percent (using deriv not direct)
        # Deriv and integ lists now in different order
        self.activequant.sort(key=lambda x: float(x[1]), reverse=True)
        # Headers for quant display (deriv then integ)
        rowframe=tk.Frame(self.quant_frame)
        tk.Label(rowframe, text='Element').pack(side=tk.LEFT)
        tk.Label(rowframe, text='At%_deriv').pack(side=tk.LEFT)
        # adjampl is adjusted amplitude of sm-diff (signal strength measure)
        tk.Label(rowframe, text='Adjampl').pack(side=tk.LEFT)
        tk.Label(rowframe, text='At%_integ').pack(side=tk.LEFT)
        tk.Label(rowframe, text='err at%_integ').pack(side=tk.LEFT)
        # adjcnts is adjusted counts of peak integration (signal strength measure)
        tk.Label(rowframe, text='Adjcnts').pack(side=tk.LEFT)
        rowframe.pack(fill=tk.X, expand=1)
        # For values display len(peaks)!=len(activeelems)
        for i, [elem, adjampl, atperc] in enumerate(self.activequant):
            rowframe=tk.Frame(self.quant_frame)
            tempstr=elem+'   '
            tk.Label(rowframe, text=tempstr).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.1f" % atperc).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.0f" % adjampl).pack(side=tk.LEFT)
            # Now find integvals for this peak/element
            # values are [elem symb, adjcnts, erradjcnt, atperc, erratperc]
            thisel=[e for e in self.activequant2 if e[0]==elem]
            print('Length of this element')
            if len(thisel)!=1:
                print('Error finding integ results for', elem)
            else:
                tk.Label(rowframe, text="%.1f" % thisel[3]).pack(side=tk.LEFT)
                tk.Label(rowframe, text="%.1f" % thisel[4]).pack(side=tk.LEFT)
                tk.Label(rowframe, text="%.1f" % thisel[1]).pack(side=tk.LEFT)
            rowframe.pack(fill=tk.X, expand=1)