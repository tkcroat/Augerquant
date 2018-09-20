'''
This is intended to be a plot interface w/ 
1)filenumber/filename text string filter (allowing update of current Augerfile)
2) adjustable charging energy value, 
3) comparison of data with plotted line energies (including shift), 

entry for element list? 


'''
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import random

# need to load AESdataset class from Auger_shift_pyqt
os.chdir('H:\Research_data\Miscellaneous\Auger\PH_NWA') 
        
def AESgui():
    ''' Launcher function for Auger quantmap GUI '''
    root = tk.Tk()
    root.wm_title("Auger spectra explorer")
    GUIMain(root)
    root.mainloop()
    return

class GUIMain():    
    ''' Main container for plotter, options (at right), and fileloader (bottom) 
    pass current working directory as default directory '''
    def __init__(self,root):
        self.root = root
        self.root.wm_title("Auger quant map")
        self.plots_frame = tk.Frame(self.root)
        self.plots_frame.pack(side=tk.TOP)
        self.spec_frame = tk.Frame(self.plots_frame)
        self.spec_frame.pack(side=tk.LEFT)
        self.roi_frame = tk.Frame(self.root) 
        self.roi_frame.pack(side=tk.TOP)
        self.specviewer= SpectraViewer(self.spec_frame,self)
        self.opts = GUIopts(self.roi_frame,self)
        
        # Menubars
        self.menubar=tk.Menu(self.root)
        filemenu = tk.Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="Load AES files", command=self.load_AESfiles)
        filemenu.add_command(label="Save", command=lambda: 
            self.args_popup_menu({'command':'save','radio':[]}))
        filemenu.add_command(label="Exit", command=self.on_quitapp)
        self.menubar.add_cascade(label="File", menu=filemenu)
        self.root.config(menu=self.menubar)

        specmenu = tk.Menu(self.menubar, tearoff=0)
        # finds negative peaks in sm-diff spectrum for active (checked) elements
        specmenu.add_command(label="Find all peaks (deriv + integ) ", command=self.rois.find_all_peaks)
        specmenu.add_command(label="Uniform filter", command=lambda: 
            self.args_popup_menu({'command':'uniform filter', 'entries':[['filter size',1]]}))
        specmenu.add_command(label="Mask ", command=lambda: 
            self.args_popup_menu({'command':'uniform filter', 'entries':[['filter size',1]]}))
        self.menubar.add_cascade(label="Spectral Commands", menu=specmenu)
        self.root.config(menu=self.menubar)
    
    def on_quitapp(self):
        msg = "Quitting:\nUnsaved progress will be lost.\nDo you wish to Continue?"
        if tkmess.askokcancel("Quantmap",msg):
            self.root.destroy()
       
    def load_AESfile(self):
        ''' Load standard Auger dataset file using directory 
        '''
        fullpath= filedialog.askopenfilename(title='Select pixarray file', 
                                             filetypes=[("pix array","*.csv")])
        (directory, filename)=os.path.split(fullpath)
        print('Directory is', directory)
        self.opts.open_AESfile(directory)

    def args_popup_menu(self, kwargs):
        ''' Menu launched top-level window to get entered args/kwargs entry and 
        then call GUIrois method (home of QMfile data object and other assorted methods)
        kwargs: command -  name of desired method
            param name & associated value (e.g.  kwargs={'filter size':1})
        implement radiobuttons using 'radio':['plottype',['scatter','line']]
        implement checkbuttons using 'check':['add backfit',True]
        implement entries using 'entry':['filter size',1]
        '''
        def abort():
            t.destroy()
        def runcmd():
            ''' run command w/ entered args/kwargs '''
            # Construct args, kwargs for method call
            myargs={}
            for i, (key, val) in enumerate(kwargs.items()):
                if key!='command':
                    myargs.update({val[0], tkvars[i].get()})
                else:
                    myargs.update({'command', kwargs.get('command')})
            self.rois.runcmd(myargs)
            t.destroy()
        t = tk.Toplevel(self.root) # open new toplevel window
        tkvars=[] # Display and alter params passed in kwargs
        
        bottomframe=tk.Frame(t)
        # Key gives type of tkinter object
        for i, (key, val) in enumerate(kwargs.items()):
            if 'rad' in key: # Make radiobutton w/ choices list 
                prow=tk.Frame(bottomframe)
                [param, choices]=kwargs.get(key,[])
                tk.Label(prow, text=param).pack(side=tk.LEFT)
                tkvars.append(tk.StringVar()) # single common variable for chosen radiobutton
                for j, val in enumerate(choices): # list of opts for radiobutton
                    tk.Radiobutton(prow, text=val, value=val, variable=tkvars[i]).pack(side=tk.LEFT)
                prow.pack(side="top", fill="both", expand=True, padx=100, pady=100)
            elif 'chk' in key: # each dict val has param name, default bool val as 2 item list
                prow=tk.Frame(bottomframe)
                [param, val]=kwargs.get(key,['',''])
                tkvars.append(tk.BooleanVar())
                tkvars[i].set(val)
                tk.Checkbutton(prow, text=param, variable=tkvars[i]).pack(side=tk.LEFT)
                prow.pack(side="top", fill="both", expand=True, padx=100, pady=100)
            elif 'ent' in key:
                prow=tk.Frame(bottomframe)
                [param, val]=kwargs.get(key,[])
                tk.Label(prow, text=param).pack(side=tk.LEFT)
                tkvars.append(tk.StringVar())
                tk.Entry(prow, textvariable=tkvars[i]).pack(side=tk.LEFT)
                prow.pack(side="top", fill="both", expand=True, padx=100, pady=100)
            elif key=='command': # put command name at top? 
                topframe=tk.Frame(t)
                tkvars.append(tk.StringVar()) # unused dummy
                tk.Label(topframe, text=key).pack(side=tk.LEFT)
                topframe.pack(side="top", fill="both", expand=True, padx=100, pady=100)
        bottomframe.pack(side="bottom", fill="both", expand=True, padx=100, pady=100)
        # Row for abort & run buttons
        prow=tk.Frame(t)
        tk.Button(prow, text='Abort', command=abort).pack(side=tk.LEFT)
        mystr='Run '+kwargs.get('command','')
        tk.Button(prow, text=mystr, command=runcmd).pack(side=tk.LEFT)
        prow.pack(side="top", fill="both", expand=True, padx=100, pady=100)
        
class App(tk.Tk):
    def __init__(self, parent=None):
        tk.Tk.__init__(self, parent)
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.title("Auger spectral plot")
        # load/instantiate Auger data set 
        AESdata=AESdataset() 
        self.Augerparamlog, self.spelist, self.Smdifpeakslog, self.Backfitlog, self.Integquantlog, self.AESquantparams=AESdataset() 
        self.filename=self.spelist.iloc[0]['Filename'] # for autoloading of first spe file
        self.Augerfile=pd.DataFrame()
        self.filterstr = tk.StringVar()# filenums, name filters
        
        tk.Label(self, text='Filename/number filter string').grid(row=0, column=0)
        filentry=tk.Entry(self, textvariable=self.filterstr)
        filentry.grid(row=0, column=1)
        
        button = tk.Button(self, text="Quit", command=self.on_click)
        button.grid(row=1, column=0)
        button = tk.Button(self, text="Update plot", command=self.on_click)
        button.grid(row=1, column=0)
        
        fig = Figure(figsize=(6, 4), dpi=96)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        x, y = self.data(self.n.get(), self.mu.get())
        self.line1, = ax.plot(x, y)

        self.graph = FigureCanvasTkAgg(fig, master=self)
        canvas = self.graph.get_tk_widget()
        canvas.grid(row=0, column=2)

    def on_click(self):
        self.quit()

    def filterfiles(self):
        ''' grab filter string and apply to  '''
        self.quit()

    def on_change(self, value):
        x, y = self.data(self.n.get(), self.mu.get())
        self.line1.set_data(x, y)  # update data

        # set plot limit
        # ax = self.graph.figure.axes[0]
        # ax.set_xlim(min(x), max(x))
        # ax.set_ylim(min(y), max(y))

        # update graph
        self.graph.draw()

    def data(self, n, mu):
        lst_y = []
        for i in range(n):
            lst_y.append(mu * random.random())
        return range(n), lst_y

class GUIrois():
    ''' Parent is GUImain, manages QMfile displayed in GUIplotter
    handles element and plot selections '''
    def __init__(self,root,parent):
        self.root = root
        self.parent = parent
        # instance of QMfile local to the roi/opts window
        self.QMfile = None  

        self.tkelems=[] # bools list for elem display or quant
        self.activequant=[] # for at % results (on extracted spectrum)
        self.showelems=False # toggle for showing elemental lines on plot
        self.currxy = None # x,y of extracted spectrum (or avg x,y if mult pixels)
        self.togglederiv =False # plot quant

        # Element selector checkboxes
        self.left_frame = tk.Frame(self.root)
        self.elems_frame = tk.Frame(self.left_frame, pady=10)
        self.elems_frame.pack(side=tk.TOP,fill=tk.X,expand=1)

        self.misc_frame = tk.Frame(self.left_frame, pady=10)
        self.misc_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        self.left_frame.pack(side=tk.LEFT)
        # Frame for display of counts/quant results (at right)
        self.quant_frame = tk.Frame(self.root, pady=10)
        self.quant_frame.pack(side=tk.LEFT,fill=tk.X,expand=1)        
        
        # Element presets (top of misc frame)
        rowframe=tk.Frame(self.misc_frame)
        tk.Button(rowframe, text='Clear all', command=self.clearall).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Button(rowframe, text='Select all', command=self.selectall).pack(side=tk.LEFT,fill=tk.X,expand=1)
        rowframe.pack(fill=tk.X, expand=1)
        
        # permanent buttons in misc_frame
        rowframe=tk.Frame(self.misc_frame)
        self.plottype=tk.StringVar()
        tk.Radiobutton(rowframe, text='Shiftmap',value='Shiftmap', 
                    variable=self.plottype).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Radiobutton(rowframe, text='Amplmap',value='Amplmap', 
                    variable=self.plottype).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Radiobutton(rowframe, text='Integmap',value='Integmap', 
                    variable=self.plottype).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Radiobutton(rowframe, text='Countsmax',value='Countsmax', 
                    variable=self.plottype).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Radiobutton(rowframe, text='Elemmap',value='Elemmap', 
                    variable=self.plottype).pack(side=tk.LEFT,fill=tk.X,expand=1)
        tk.Button(rowframe, text='Plot', command=self.plot_maps).pack(side=tk.LEFT,fill=tk.X,expand=1)
        rowframe.pack(fill=tk.X, expand=1)
        self.toggle_button = tk.Button(
            self.misc_frame,text="Toggle deriv/direct",command=self.toggle_deriv)
        self.toggle_button.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.label_button = tk.Button(
            self.misc_frame,text="Label elements",command=self.label_elems)
        self.label_button.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.quant_button = tk.Button(
            self.misc_frame,text="Update quant", command=self.do_quant)
        self.quant_button.pack(side=tk.TOP,fill=tk.X,expand=1)
    
    def create_QMfile(self, directory):
        ''' Creates QM file instance (called from menu)
        automatically finds pixarray and works from there '''
        print('Creating QM file.')
        self.QMfile = AESquantmap(directory)
        print("QMfile created with name ", self.QMfile.uniquename)
        #print("QMfile", QMfile.uniquename, "created.")
        for child in self.elems_frame.winfo_children():
            child.destroy()
        # no direct pass to GUIplotter (only 2D projections)
        # loads quant elements into elems frame
        self.display_elems()
        
    def save_specimage(self):
        ''' Call AESquantmap save_specimage  '''
        if not self.QMfile:
            return
        self.QMfile.save_specimage()
        
    def toggle_deriv(self):
        ''' Toggle plotting from direct counts plot to s7d7 smooth-deriv '''
        if self.togglederiv==False:
            self.togglederiv=True
        else:
            self.togglederiv=False
    def selectall(self):
        ''' Clear selected elements '''
        for i, tkbool in enumerate(self.tkelems):
            self.tkelems[i].set(1)

    def clearall(self):
        ''' Clear selected elements '''
        for i, tkbool in enumerate(self.tkelems):
            self.tkelems[i].set(0)

    def load_maps(self):
        '''Menu/main lauched '''
        if self.QMfile is not None:
            self.QMfile.load_maps()
    
    def save_maps(self):
        '''Menu/main lauched save of amplmaps, integmaps and shiftmaps '''
        if self.QMfile is not None:
            self.QMfile.save_maps()

    def save_pixarray(self):
        '''Menu/main lauched save of pixarray file (after linking with
        underlying data files '''
        if self.QMfile is not None:
            self.QMfile.save_pixarray()
        
    def save_ampl_images(self):
        ''' Save all extracted amplitude images as separate jpgs '''
        if self.QMfile is not None:
            self.QMfile.save_ampl_images()
    
    def find_all_peaks(self):
        ''' Normal and best method for data extraction from specimage '''
        if self.QMfile is not None:
            self.QMfile.find_all_peaks()
            
    def find_peaks(self):
        ''' For selected element(s), find peak center '''
        # check if shiftmaps and amplitudes maps have already been saved
        print('Running find_peaks in GUIroi')
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                self.QMfile.find_negpeaks(i)
                print('Negpeak positions found for element', str(i))
    
    def plot_multiplex(self):
        ''' Display current extracted spectrum in specviewer
        shows only active elems '''
        if self.QMfile.extracted is None: return
        actelemdata=[]
        vals=[] # for scatter points on spectral plots (active elems only)
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                actelemdata.append(self.QMfile.elemdata[i])
                if self.togglederiv==False:
                    vals.append(self.QMfile.integparams[i])
                else:
                    vals.append(self.QMfile.derivparams[i])
                
        print("Plotting current extracted spectrum")
        # pass current extracted spectrum or its deriv (1d np arr) and subset of active elem data info
        if self.togglederiv==False: # plot direct
            # pass integration center/ background fit for plotting 
            pkwargs={'type':'integ', 'vals':vals}
            self.parent.specviewer.plot_multiplex(self.QMfile.extracted, self.QMfile.energy, 
                                                  actelemdata, self.currxy, **pkwargs)
        else: # plot deriv
            # pass list of deriv params (xvals/yvals) for plot for each peak
            pkwargs={'type':'deriv', 'vals':vals}
            self.parent.specviewer.plot_multiplex(self.QMfile.extracts7d7, self.QMfile.energy, 
                                                  actelemdata, self.currxy, **pkwargs)
        
    def plot_maps(self):
        ''' Display 2D arrays of various types in mapviewer '''
        activeelems=[]
        plotmaps=[]
        title=''
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                if self.plottype.get()=='Shiftmap':
                    # Use togglederiv to decide between deriv shift and integ shift
                    if self.QMfile.shiftmaps[i] is not None:
                        activeelems.append(self.QMfile.elements[i])
                        if self.togglederiv: # use deriv based shift
                            plotmaps.append(self.QMfile.shiftmaps[i][:,:,0])
                            title='Peak shift deriv'
                        else: # use direct peak shift
                            plotmaps.append(self.QMfile.shiftmaps[i][:,:,1])
                            title='Peak shift direct'
                elif self.plottype.get()=='Amplmap':
                    if self.QMfile.amplmaps[i] is not None:
                        activeelems.append(self.QMfile.elements[i])
                        # 0th layer is amplitude
                        plotmaps.append(self.QMfile.amplmaps[i][:,:, 0]) 
                        title='Peak amplitude'
                elif self.plottype.get()=='Integmap':
                    if self.QMfile.integmaps[i] is not None:
                        activeelems.append(self.QMfile.elements[i])
                        plotmaps.append(self.QMfile.integmaps[i][:,:,0])                                
                        title='Integcounts map'
                elif self.plottype.get()=='Countsmax':
                    if self.QMfile.integmaps[i] is not None:
                        activeelems.append(self.QMfile.elements[i])
						# intensity before background subtraction at peak
                        plotmaps.append(self.QMfile.integmaps[i][:,:,2])                                
                        title='Integcounts map'
						# TODO finish element map creation
                elif self.plottype.get()=='Elemmap':
                    if self.QMfile.amplmaps[i] is not None:
                        activeelems.append(self.QMfile.elements[i])
                        plotmaps.append(self.QMfile.elemmaps[i][:,:,0]) 
                        title='Element map'
        print("Plotting ", self.plottype.get(), "for elements", ".".join(activeelems))
        self.parent.mapviewer.replot_maps(activeelems, plotmaps, title)
        
    def display_elems(self):
        ''' Display available quant elements in elems_frame (self.QMfile.elements only 
        contains elements that already have quant results (not other random peaks) '''
        for child in self.elems_frame.winfo_children():
            child.destroy()
        # Need to call frame header Label and 
        # Write header row into backregs
        tk.Label(self.elems_frame, text='Available Elements').pack(side=tk.TOP,fill=tk.X,expand=1)

        # tkelems bool variables for active/inactive for each element
        self.tkelems=[]
        for i, elem in enumerate(self.QMfile.elements):
            self.tkelems.append(tk.BooleanVar())
            self.tkelems[i].set(True)
        # Unfortunately tk/mpl combo requires use of pack (not grid) 
        for i in range(0, len(self.QMfile.elements), 4):
            # associated checkbutton for each quantelem
            elemlistframe=tk.Frame(self.elems_frame)
            tk.Checkbutton(elemlistframe, variable=self.tkelems[i]).pack(side=tk.LEFT)
            tk.Label(elemlistframe, text=self.QMfile.elements[i]).pack(side=tk.LEFT)
            try:
                tk.Checkbutton(elemlistframe, variable=self.tkelems[i+1]).pack(side=tk.LEFT)
                tk.Label(elemlistframe, text=self.QMfile.elements[i+1]).pack(side=tk.LEFT)
            except: # out of list range problem
                pass
            try:
                tk.Checkbutton(elemlistframe, variable=self.tkelems[i+2]).pack(side=tk.LEFT)
                tk.Label(elemlistframe, text=self.QMfile.elements[i+2]).pack(side=tk.LEFT)
            except:
                pass
            try:
                tk.Checkbutton(elemlistframe, variable=self.tkelems[i+3]).pack(side=tk.LEFT)
                tk.Label(elemlistframe, text=self.QMfile.elements[i+3]).pack(side=tk.LEFT)
            except:
                pass
            elemlistframe.pack(fill=tk.X, expand=1)
            
            
    def label_quant(self):
        ''' Add a text label with current quant results to map viewer?
        launched via button  '''
        #TODO adjust position of text quant label 
        elems=[i[0] for i in self.activequant]
        vals=[f[1] for f in self.activequant]
        vals=[int(i) if i>1 else "%0.1f" % i for i in vals]
        self.parent.mapviewer.label_quant(elems, vals)

    def label_elems(self):
        ''' Get currently selected elements from elems frame and associated energy
        pass to plotter.label_elems
        toggle style button
        TODO  not implemented ... is this even necessary for multiplex?
        '''
        if self.showelems:
            self.showelems=False
        else:
            self.showelems=True
        elemparams=[] # list of [elem, energyval]
        # active (checked) tk elements will automatically be updated on check, right?
        if self.showelems:
            for i, tkbool in enumerate(self.tkelems):
                if tkbool.get():
                    match=self.EDXdataset.EDXquantparams[self.EDXdataset.EDXquantparams['element']==self.QMfile.elements[i]]
                    ev=match.iloc[0]['energy']
                    elemparams.append([self.QMfile.elements[i],ev])
        # now pass info to plotter (can be empty list)
        self.parent.mapviewer.plot_elems(elemparams)
    
    def do_quant(self):
        ''' Generate at % and error for subset of selected elements
        then display quant ... linked to quant button 
        quant results for current extract spectrum in QMfile.quant (e.g.[elem, 
        ampl, shift, corrampl (kfactor adjusted)] 
         '''
        # Clear current values
        self.activequant=[] # list w/ [elem symb, atperc, ampl]
        self.activequant2=[]  # list w/ [elem, atperc, ]
        '''Parameters saved in quant2: [0] elem 1 integcnts [2] energy [3] corrcnts
            [4] slope [5] intercept
        '''
        # note erratperc (3rd in list) is actual error not % err
        corrsum=0.0
        corrsum2=0.0
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                # Add elem, val for at % and corrected amplitude (deriv-based)
                if self.QMfile.quant[i][2]>0:
                   corrsum+=self.QMfile.quant[i][2]
                   self.activequant.append([self.QMfile.quant[i][0], 0.0, 
                                            self.QMfile.quant[i][2]])
                else:
                    self.activequant.append([self.QMfile.quant[i][0], 0.0, 0.0])
                # Add element, at% value, corrcnts and integcounts (integ-based)
                # quant2 is elem, integcnts, corrcnts
                if self.QMfile.quant2[i][2]>0:
                   corrsum2+=self.QMfile.quant2[i][2]
                   # elem, at. %, corrcnts, integcnts
                   self.activequant2.append([self.QMfile.quant2[i][0], 0.0, 
                            self.QMfile.quant2[i][2], self.QMfile.quant2[i][1]])
                else:
                    self.activequant2.append([self.QMfile.quant[i][0], 0.0, 0.0, 
                                              self.QMfile.quant2[i][1]])
                
        for i, vals in enumerate(self.activequant): # calc at %
            self.activequant[i][1]=100*self.activequant[i][2]/corrsum
            self.activequant2[i][1]=100*self.activequant2[i][2]/corrsum2
        self.display_quant()
        
    def display_quant(self):
        ''' Display elemental quant (deriv and integ) for extracted spectrum
        '''
        # Clear any existing widgets in backreg frame
        for child in self.quant_frame.winfo_children():
            child.destroy()
        # sort active quant elements by at percent
        self.activequant.sort(key=lambda x: float(x[1]), reverse=True)
        self.activequant2.sort(key=lambda x: float(x[1]), reverse=True)
        # Write header row into backregs 
        rowframe=tk.Frame(self.quant_frame)
        tk.Label(rowframe, text='Element').pack(side=tk.LEFT)
        tk.Label(rowframe, text='At%').pack(side=tk.LEFT)
        tk.Label(rowframe, text='Corrampl').pack(side=tk.LEFT)
        #tk.Label(rowframe, text='Corrcnts').pack(side=tk.LEFT)
        #tk.Label(rowframe, text='Err corrcnts').pack(side=tk.LEFT)
        print('Element    At%    Corrampl')
        rowframe.pack(fill=tk.X, expand=1)
        # For values display len(quantelems)!=len(activeelems)
        for i, [elem, atper, corrampl] in enumerate(self.activequant):
            rowframe=tk.Frame(self.quant_frame)
            tempstr=elem+'   '
            tk.Label(rowframe, text=tempstr).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.1f" % atper).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.0f" % corrampl).pack(side=tk.LEFT)
            rowframe.pack(fill=tk.X, expand=1)
            print(elem+'    '+"%.1f" % atper+'    '"%.0f" % corrampl)
        # also print sorted quant to spyder console
        print('Element    At%    Corrcnts   Integcounts')
        rowframe.pack(fill=tk.X, expand=1)
        # For values display len(quantelems)!=len(activeelems)
        for i, [elem, atper, corrcnts, integcnts] in enumerate(self.activequant2):
            rowframe=tk.Frame(self.quant_frame)
            tempstr=elem+'   '
            tk.Label(rowframe, text=tempstr).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.1f" % atper).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.0f" % corrcnts).pack(side=tk.LEFT)
            tk.Label(rowframe, text="%.0f" % integcnts).pack(side=tk.LEFT)
            rowframe.pack(fill=tk.X, expand=1)
            print(elem+'    '+"%.1f" % atper+'    '+"%.0f" % corrcnts +'    '+"%.0f" % integcnts)
    
    def runcmd(self, ckwargs):
        ''' for menu launched commands with params '''
        cmd=ckwargs.get('command')
        if cmd=='uniform filter':
            self.QMfile.uniform_filter(self, size=ckwargs.get('filter size', 1))
    
    def getcurrxy(self, indices):
        ''' Find X, Y values from passed selection (list of indices) '''
        xys=[[i,j] for i in range(0,self.QMfile.dim[0]) for j in range(0,self.QMfile.dim[1])]
        # get subset of x,y vals using index passed by lasso 
        selectxys=[xys[i] for i in indices]
        if len(selectxys)==1:
            self.currxy=selectxys[0] # X, Y map location of current extracted spectrum
            print('Current xy is:', self.currxy)
        else: # Get average X,Y of multiple pixel selection 
            print(len(selectxys),' pixels selected')
            xvals=[i[0] for i in selectxys]
            yvals=[i[1] for i in selectxys]
            avgx=int(sum(xvals)/len(xvals))
            avgy=int(sum(yvals)/len(yvals))
            self.currxy=[avgx, avgy]
    
    def extract_spectrum(self, selected):
        ''' Pass back lasso-selected path or single point, create/extract average
        spectrum over this X-Y set '''
        if self.QMfile is not None:
            self.QMfile.extract_spectrum(selected)
            self.getcurrxy(selected)
        # Now plot multiplex (w/ active elements)
        self.plot_multiplex()
    
    def save_extracted(self):
        ''' Save currently extracted spectrum (calling method from open QMfile'''
        if self.QMfile is not None:
            self.QMfile.save_extracted()
    
    def print_spectralregs(self):
        ''' Pop up to show spectral regions, evbreaks, etc.   '''
        if self.QMfile is not None:
            print(self.QMfile.spectralregs)
        
class SpectraViewer():
    ''' Spectral plotter window for Auger survey or multiplex  spectra '''
    def __init__(self,root, parent):
        self.root = root
        self.parent = parent
        self.figure = mpl.figure.Figure(figsize=PLOT_SIZE, dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
        self.canvas = FigureCanvasTkAgg(self.figure,self.root)
        # just use standard toolbar
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,self.root)
        self.toolbar.update()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.show()
    
    def plot_multiplex(self, extracted, energy, elemdata, currxy, **pkwargs):
        ''' Add variable number of subplots
        called from GUIrois
        extracted is 1D numpy array (same len as full multiplex) -- either deriv
            or direct is passed 
        energy is full multiplex range ev vals for extracted spectrum (as list)
        elemdata has peak stop/start indices for plots -- only active elements
        currxy is X, Y of extracted spectrum (or avg x,y of lassoed ROI)
        '''
        # since # of subplots can change, need to destroy and recreate
        try:
            self.canvas.get_tk_widget().destroy() # destroy previous plot
            self.toolbar.destroy()
        except:
            pass
        plottype=pkwargs.get('type') # integ or deriv
        vals=pkwargs.get('vals') # list with scatter points/backfits/etc.
        # plot from elemdata[i][ holds indices
        numcols=min(len(elemdata),2) # 1 or 2 columns
        numrows=math.ceil(len(elemdata)/2)
        self.figure = mpl.figure.Figure(figsize=PLOT_SIZE, dpi=100)
        self.figure.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
        self.ax=[]
        for i, elemd in enumerate(elemdata):
            self.ax.append(self.figure.add_subplot(numrows,numcols,i+1))
            [lowind, junk]=elemd[3]
            [junk,highind]=elemd[4]
            idealev=elemd[8] # ideal peak eV 
            symbol=elemd[0] # name of element/peak
            self.ax[i].plot(energy[lowind:highind], extracted[lowind:highind])
            self.ax[i].axvline(x=idealev)
            energy=[int(i) for i in energy] # ensure ints

            # for deriv vals is list of dfs w/ scatter points
            if plottype=='deriv':
                # derxvals and deryvals passed np arrays to add pospeak/negpeak
                # as scatter plot 
                [elem, xvals, yvals, ampl]=vals[i]
                self.ax[i].scatter(xvals, yvals, color='r')
                # add elem and amplitude as text label
                tempstr=symbol+' Ampl:'+ "%.2f" % ampl
                self.ax[i].set_title(tempstr, fontsize=10)
            elif plottype=='integ':
                # elem, peak energy (integration center), integcnts, slope/ intercept of backfit
                [elem, peakev, integcnts, slope, intercept]=vals[i] 
                # Scatter point at integration center
                yvals= extracted[energy.index(peakev)]
                self.ax[i].scatter(peakev, yvals, color='r')
                # Plot background fit line
                x=np.linspace(min(energy[lowind:highind]), max(energy[lowind:highind]), 100)
                self.ax[i].plot(x, x*slope+intercept, color='r')
                # add elem symbol and integcounts as subplot title
                tempstr=symbol+' Integcnts:'+str(integcnts)
                self.ax[i].set_title(tempstr, fontsize=10)
                # label with integrat
        # add vertical lines at ideal position
        labelstr='X: '+str(currxy[0])+' Y: '+str(currxy[1])
        self.ax[0].text(0.05,0.95, labelstr, transform=self.ax[0].transAxes, fontsize=12)
        # recreate and pack
        self.canvas = FigureCanvasTkAgg(self.figure, self.root)        
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,self.root)
        self.toolbar.update()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.show()
    
    def label_quant(self, elems, vals):
        ''' Add quant text label with active elems and at. % values '''
        if self.EDXfile is None:return
        # Add vertical lines at known element energies
        fullstr=''
        for i, (elem,val) in enumerate(zip(elems, vals)):
            tempstr=r'$%s_{%.0f}$' %(elem, float(val))
            fullstr+=tempstr            
        # transform=ax.transAxes uses coords from 0 to 1 (instead of true x and y vals)
        self.ax.text(0.05,0.95, fullstr, fontsize=30, verticalalignment='top', transform=self.ax.transAxes)
        self.canvas.show()

