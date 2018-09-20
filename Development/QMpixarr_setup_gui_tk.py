# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:41:13 2017

@author: tkc
"""
import os
import tkinter as tk
import tkinter.messagebox as tkmess
from tkinter import filedialog
import matplotlib as mpl # using path, figure, rcParams
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.widgets import Lasso
from matplotlib import path
import pandas as pd
import numpy as np
from PIL import Image, ImageDraw
import shutil, sys, fileinput

PLOT_SIZE = (8,6)  # better on laptop
PLOT_SIZE = (7,6) # width, height
AESQUANTPARAMFILE='C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\AESquantparams.csv'
# Couldn't figure out proper encoding on some .phi file so just open/overwrite existing
AREAFILE='C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\spatial_areas_sample_min.phi'
BEAMDEFLECT='FJS-10kVdeflected' # Setting to deflect beam off sample after quant map

def pixgui():
    ''' Launcher function for Auger quantmap pixarray creator '''
    root = tk.Tk()
    root.wm_title("QM pixarray setup")
    GUIMain_pixarray(root)
    root.mainloop()
    return

class GUIMain_pixarray():
    ''' Main container for plotter, options (at right), and fileloader (bottom) 
    pass current working directory as default directory '''
    def __init__(self, root):
        self.root = root
        self.root.wm_title("QM pixarray setup")
        self.left_frame = tk.Frame(self.root)
        self.image_frame = tk.Frame(self.left_frame)
        self.image_frame.pack(side=tk.TOP)
        self.opts_frame = tk.Frame(self.left_frame) 
        self.opts_frame.pack(side=tk.BOTTOM)
        self.left_frame.pack(side=tk.LEFT)
        
        self.multiplex_frame = tk.Frame(self.root) 
        self.multiplex_frame.pack(side=tk.LEFT)

        self.pixviewer = Pixarrviewer(self.image_frame,self) # root,parent
        self.pixopts = Pixarropts(self.opts_frame,self)
        self.multparams = Multiplexopts(self.multiplex_frame, self)
        
        # Menubars
        self.menubar=tk.Menu(self.root)
        filemenu = tk.Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="Load image", command=self.load_image)
        filemenu.add_command(label="Save QM pixarray", 
                             command=self.pixopts.save_pixarray)
        filemenu.add_command(label="Save image", command=lambda: 
            self.args_popup_menu({'command':'save_image', 
                'entry':['filename',self.pixopts.pixarrname.get()],
                'radio':['imtype',['full field','crop']],
                'radio2':['overlay',['outline','grid','none']]}))
        filemenu.add_command(label="Exit", command=self.on_quitapp)
        self.menubar.add_cascade(label="File", menu=filemenu)
        self.root.config(menu=self.menubar)
    
    def on_quitapp(self):
        msg = "Quitting:\nUnsaved progress will be lost.\nDo you wish to Continue?"
        if tkmess.askokcancel("Quantmap",msg):
            self.root.destroy()
       
    def load_image(self):
        ''' Load standard QM (quantmap) file (paramlog with data returned to a DataManager '''
        fullpath= filedialog.askopenfilename(title='Select sem image', 
                filetypes=[("jpg file","*.jpg"), ("Auger sem file","*.sem")])        
        self.pixviewer.load_image(fullpath)
     
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
                    myargs.update({val[0]:tkvars[i].get()})
                else:
                    myargs.update({'command':kwargs.get('command')})
            self.pixviewer.runcmd(**myargs)
            t.destroy()
        t = tk.Toplevel(self.root) # open new toplevel window
        tkvars=[] # Display and alter params passed in kwargs

        # Key gives type of tkinter object
        for i, (key, val) in enumerate(kwargs.items()):
            if 'rad' in key: # Make radiobutton w/ choices list 
                prow=tk.Frame(t)
                [param, choices]=kwargs.get(key,[])
                tk.Label(prow, text=param).pack(side=tk.LEFT)
                tkvars.append(tk.StringVar()) # single common variable for chosen radiobutton
                for j, val in enumerate(choices): # list of opts for radiobutton
                    tk.Radiobutton(prow, text=val, value=val, variable=tkvars[i]).pack(side=tk.LEFT)
                prow.pack(side=tk.TOP)
            elif 'chk' in key: # each dict val has param name, default bool val as 2 item list
                prow=tk.Frame(t)
                [param, val]=kwargs.get(key,['',''])
                tkvars.append(tk.BooleanVar())
                tkvars[i].set(val)
                tk.Checkbutton(prow, text=param, variable=tkvars[i]).pack(side=tk.LEFT)
                prow.pack(side=tk.TOP)
            elif 'ent' in key:
                prow=tk.Frame(t)
                [param, val]=kwargs.get(key,[])
                tk.Label(prow, text=param).pack(side=tk.LEFT)
                tkvars.append(tk.StringVar())
                tk.Entry(prow, textvariable=tkvars[i]).pack(side=tk.LEFT)
                prow.pack(side=tk.TOP)
            elif key=='command': # put command name at top? 
                topframe=tk.Frame(t)
                tkvars.append(tk.StringVar()) # unused dummy
                tk.Label(topframe, text=key).pack(side=tk.LEFT)
                topframe.pack(side=tk.TOP)
        # Row for abort & run buttons
        prow=tk.Frame(t)
        tk.Button(prow, text='Abort', command=abort).pack(side=tk.LEFT)
        mystr='Run '+kwargs.get('command','')
        tk.Button(prow, text=mystr, command=runcmd).pack(side=tk.LEFT)
        prow.pack(side=tk.TOP)

class NavMapToolbar(NavigationToolbar2TkAgg):
    ''' Custom matplotlib toolbar w/ lasso pt remover and point picker
    parent is GUIplotter
    '''
    def __init__(self, canvas, root, parent):
        self.canvas = canvas
        self.root   = root
        self.parent = parent # plotter is Pixarrviewer
        self.ax= self.parent.ax # axes needed for interaction
        self.xys = None # for xy vals later associated with plot
        self.selected = None # Holding area for selected indices

        # Generic mpl toolbar using tkagg (with standard buttons)
        NavigationToolbar2TkAgg.__init__(self, canvas, root)
        # Create lasso and link to multi_select_callback
        self.lasso_button= tk.Button(master=self, text='Select region', padx=2, pady=2, command=self.startlasso)
        self.lasso_button.pack(side=tk.LEFT,fill="y")
        self.showgrid_button= tk.Button(master=self, text='Show/hide grid', padx=2, pady=2, 
            command=self.parent.toggle_grid) # show/hide toggle in pixviewer
        self.showgrid_button.pack(side=tk.LEFT,fill="y")

    def startlasso(self):
        ''' Activated by lasso menu bar button on click; disconnects prior IDs, prep for lasso button press
        '''
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpresslasso)

    def onpresslasso(self, event):
        ''' Create lasso when button pressed on active canvas/axes '''
        # ensure that current dataset is active
        self.xys = self.parent.xys # pixels list passed from plotter (parent)
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callbacklasso)
        #  self.canvas.widgetlock(self.lasso)  # skip... gives ValueError: already locked

    def callbacklasso(self, verts):
        print('Verts length is', len(verts))
        # Verts is a list of x,y coords describing drawn path
        p = path.Path(verts)
        # true/false array
        ind = p.contains_points(self.xys)
        self.selected=[i for i in range(0, len(self.xys)) if ind[i]==True]
        self.canvas.draw_idle()
        # self.canvas.widgetlock.release(self.lasso) # valueerror you don't own this lock
        del self.lasso
        self.canvas.mpl_disconnect(self.cid) # disconnect lasso tool
        
    def extract_spectrum(self):
        ''' Map only lassoed (or circular regTake single pixel or lassoed pixels, generate extracted spectrum  '''
        # Use method in GUIrois (which in turn calls QMfile method)
        self.parent.parent.opts.map_custom(self.selected)
        # TODO selected is probably sequential index # so convert to X,Y list

class Pixarrviewer():
    '''Window for display of image and overlaid pixarray
    lasso to set active/ inactive pixels 
    pixviewer can handle image --load, save w/ overlay, save/crop w/ overlay
    '''
    def __init__(self,root, parent):
        self.root = root
        self.parent = parent
        self.figure = mpl.figure.Figure(figsize=PLOT_SIZE, dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
        self.canvas = FigureCanvasTkAgg(self.figure, self.root)
        self.xys = None # used by lasso
        self.image = None # base 512x512 image for superimposition of array
        self.directory = None
        self.uniquename = None
        self.gridbool = False # show/hide for overlaid array/grid
        #  TODO is grid necessary or just pass and modify full pixarray?
        self.gridlist = [] # visual representation of current pixarray
        # xys used for lasso selection
        self.xys=[[i,j] for i in range(0,512) for j in range(0,512)]
        # Custom navselecttoolbar w/ interactive buttons
        self.toolbar = NavMapToolbar(self.canvas, self.root,self)
        self.toolbar.update()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.show()
    
    def pass_pixarray(self, pixarr):
        ''' Pass created pix array to viewer .. make grid for ImageDraw overlay '''
        self.gridlist=[]
        # pixarray already includes built in margin 
        for index, row in pixarr.iterrows():
            self.gridlist.append([row.X1,row.Y1,row.X2,row.Y2])
        
    def show_image(self):
        ''' Reload of base image (and optional grid overlay) into pixarrviewer window '''
        try:
            self.canvas.get_tk_widget().destroy() # destroy previous plot
            self.toolbar.destroy()
        except:
            pass
        print('Starting show_image.')
        self.figure = mpl.figure.Figure(figsize=PLOT_SIZE, dpi=100)
        self.figure.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
        self.ax = self.figure.add_subplot(111)
        self.ax.imshow(self.image)
        if self.gridbool:
            draw=ImageDraw.Draw(self.image)
            for i,[x1,y1,x2,y2] in enumerate(self.gridlist):
                draw.rectangle((x1,y1,x2,y2), outline='red')
        self.canvas = FigureCanvasTkAgg(self.figure, self.root)        
        self.toolbar = NavMapToolbar(self.canvas, self.root,self)
        self.toolbar.update()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.show()
        print('Show image completed.')

    def load_image(self, fullpath):
        ''' Load and resize base image'''
        image=Image.open(fullpath)
        if image.size[0]!=512:
            image=image.resize((512,512), resample=Image.NEAREST)
        (directory, filename)=os.path.split(fullpath)
        self.directory = directory
        self.uniquename = ".".join(filename.split('.')[0:-1])
        self.image= image
        self.show_image()
    
    def mod_grid(self):
        ''' Modify grid on/off pixel regions based on lasso, rect selector 
        passed back to pixarropts before save of pixarray '''
        # TODO finish lasso or rectangle selector mods of grid
        pass
        
    def toggle_grid(self):
        ''' Toggle visibility of overlaid QM pixarray grid '''
        if self.gridbool:
            self.gridbool=False
        else:
            self.gridbool=True
        self.show_image() # clear and reload main view

    def runcmd(self, **pkwargs):
        ''' Menu/main lauched with pop up arguments entered in toplevel
        passed via pkwargs... can define multiple '''
        command=pkwargs.get('command')
        if command=='save_image':
            self.save_image(pkwargs)
    
    def save_image(self, pkwargs):
        ''' Various types of image saves 
        imtype:  fullfield or crop (subtract margins)
        overlay:  outline, grid or none '''
        print('pkwargs are:', pkwargs)
        thisim=self.image
        imtype=pkwargs.get('imtype','full field')
        overlay=pkwargs.get('overlay','none')
        fname=pkwargs.get('filename','none')
        if overlay=='grid':
            draw=ImageDraw.Draw(thisim)
            for i,[x1,y1,x2,y2] in enumerate(self.gridlist):
                draw.rectangle((x1,y1,x2,y2), outline='red')
        elif overlay=='outline':
            # find min x1, max x2, min y1, max y2
            minx1=min([i[0] for i in self.gridlist])
            maxx2=max([i[2] for i in self.gridlist])
            miny1=min([i[1] for i in self.gridlist])
            maxy2=max([i[3] for i in self.gridlist])
            # single rectangular outline of quant mapped area
            draw=ImageDraw.Draw(thisim)
            draw.rectangle((minx1,miny1,maxx2,maxy2), outline='red')
        #  Perform crop at end 
        if imtype=='crop':
            minx1=min([i[0] for i in self.gridlist])
            maxx2=max([i[2] for i in self.gridlist])
            miny1=min([i[1] for i in self.gridlist])
            maxy2=max([i[3] for i in self.gridlist])
            # perform crop to quant mapped region
            thisim.crop([minx1, miny1, maxx2, maxy2])
        imname=self.directory+'/'+fname
        thisim.save(imname)
            
class Pixarropts():
    ''' Parent is GUImain, manages pixarray parameter displayed in pixarrviewer
    handles modifications (on/off) for 
    and save of pixarray file, autotool, spatial array files
    '''
    def __init__(self,root,parent):
        self.root = root
        self.parent = parent
        # Lots of tk variables 
        self.arrsize=tk.IntVar() # Number of pixels for QM scan (usually differs from 512 pix of imaging)
        self.arrsize.set(100)
        self.margin=tk.DoubleVar()
        self.margin.set(0.2)
        self.regint=tk.IntVar()
        self.regint.set(2)
        self.atbool=tk.BooleanVar() # make new Autotool file
        self.atbool.set(False)
        self.areabool=tk.BooleanVar() # make new spatial area files
        self.areabool.set(False)
        self.scramblebool=tk.BooleanVar() 
        self.scramblebool.set(False)
        self.basename=tk.StringVar()
        self.basename.set('100x100array20marg')
        self.pixarrname=tk.StringVar()
        self.pixarrname.set('samplename')
        # Display/ entry of quantmap pixarray parameters
        self.button_frame = tk.Frame(self.root, pady=10)
        self.button_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        
        self.pixarr_frame = tk.Frame(self.root, pady=10)
        self.pixarr_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # Display/entry of multiplex scan parameters
        self.multiplex_frame = tk.Frame(self.root, pady=10)
        self.multiplex_frame.pack(side=tk.TOP,fill=tk.X,expand=1)

        self.pixarr= None # created by makesquarearray
        self.show_pixarrvars()
                
        # permanent buttons in button_frame
        rowframe=tk.Frame(self.button_frame)
        tk.Button(rowframe, text='Make QM pixarray', 
                  command=self.makesquarearray).pack(side=tk.LEFT,fill=tk.X,expand=1)
        rowframe.pack(fill=tk.X, expand=1)
        
    def show_pixarrvars(self):
        ''' Set up all the tkvars into pixarr_frame '''
        rownum=0
        tk.Label(self.pixarr_frame, text='Array Size (~100 x 100 max)').grid(row=rownum, column=1)
        tk.Entry(self.pixarr_frame, textvariable=self.arrsize).grid(row=rownum, column=0)
        rownum+=1
        tk.Label(self.pixarr_frame, text='Unmapped margin (fraction)').grid(row=rownum, column=1)
        tk.Entry(self.pixarr_frame, textvariable=self.margin).grid(row=rownum, column=0)
        rownum+=1
        tk.Label(self.pixarr_frame, text='Image reg interval (Autotool)').grid(row=rownum, column=1)
        tk.Entry(self.pixarr_frame, textvariable=self.regint).grid(row=rownum, column=0)
        rownum+=1
        tk.Label(self.pixarr_frame, text='Pix array name').grid(row=rownum, column=1)
        tk.Entry(self.pixarr_frame, textvariable=self.pixarrname).grid(row=rownum, column=0)
        rownum+=1
        tk.Label(self.pixarr_frame, text='Basename Autotool & spatial areas').grid(row=rownum, column=1)
        tk.Entry(self.pixarr_frame, textvariable=self.basename).grid(row=rownum, column=0)
        rownum+=1
        tk.Checkbutton(self.pixarr_frame, variable=self.atbool, text='Make new Autotool file?').grid(row=rownum, column=0)
        tk.Checkbutton(self.pixarr_frame, variable=self.areabool, text='Create new spatial area files?').grid(row=rownum, column=1)
        rownum+=1
        tk.Checkbutton(self.pixarr_frame, variable=self.scramblebool, text='hop beam to minimize charging?').grid(row=rownum, column=0)
        rownum+=1
        tk.Label(self.pixarr_frame, text='If same parameters, existing Autotool and spatial area files can be reused.').grid(row=rownum, column=0, columnspan=3)
        rownum+=1    
    
    def makesquarearray(self):
        ''' Divide up 512x512 pixels in map into n areas and format correctly for spatial areas phi files
        (which are loaded using Autotool loops into PHI Smartsoft); 
        Mapping proceed horizontally (x horizontal, y vertical) 
        ARGS:
        arraysize -- # of pixels in one direction
        margin - % of image field that is unmapped (20% means 10% of 512 field at both edges (aka 51 pixels))
            is unmapped
        basename - basename for area definition files 
            e.g. "50x50array20m" basemname makes files 50x50array20m1, 50x50array20m2, etc.
            
        KWARGS: 'regint' - interval at which to build in image registration into autotool loop; val of 1 means
               register every 20 pixels (since each area file holds 20 defined spatial areas);  passed to makeautotool
               * this is best way to incorporate image registration in quantmap process... more flexible
               interval allowed; if instead one builds in image reg into multiplex itself, one has to 
               run image registration much more frequently which unnecessarily increases acquisition time
            'writeareas' (bool) - write new spatial area definition files (old ones can be reused if same arraysize 
                 and margin)
            'writeAutotool' --  write of new Autotool sequence (can be reused if same arraysize/regint)
            scrambled
        '''
        print('Making square array')
        pix=512  # default value for Auger instrument
        arraysize=self.arrsize.get()
        width=(pix*(1-self.margin.get())/arraysize) # width/height of scan pixel in terms of 512x512 field
        startxy=int(pix*self.margin.get()/2) # split margin between top/bottom, left/right
        mycols=['Xindex','Yindex','Areanumber','PHIname','Subnumber','X1','Y1','X2','Y2', 'Width', 'Height']
        dim=arraysize**2
        # square is the pixarray file correlating 101.spe with associated pixel in quantmap 
        self.pixarr=pd.DataFrame(index=np.arange(0,dim), columns=mycols)
        # x is horizontal axis and mapping proceeds by going across top row
        for index,row in self.pixarr.iterrows():
            xindex=index//arraysize # remainder is row (0th is top row)
            yindex=index%arraysize # mod is correct column (0th is left column)
            self.pixarr.loc[index]['Xindex']=xindex # remainder is row
            self.pixarr.loc[index]['Yindex']=yindex # mod is correct column
            left=int(width*yindex+startxy) # left-right position depends on column 
            self.pixarr.loc[index]['X1']=left
            right=int(width*yindex+startxy+width)
            self.pixarr.loc[index]['X2']=right
            top=int(width*xindex+startxy) 
            self.pixarr.loc[index]['Y1']=top # top-bottom position depends on row
            bottom=int(width*xindex+startxy+width)
            self.pixarr.loc[index]['Y2']=bottom
            self.pixarr.loc[index]['Width']=right-left # variations due to rounding error
            self.pixarr.loc[index]['Height']=bottom-top
            # true area number describing pix position after file combination
            self.pixarr.loc[index]['Areanumber']=index+1
            # max of 20 areas allowed per spatial area .phi file
            self.pixarr.loc[index]['Subnumber']=index%20 # Before combination w/ 20 areas per file
            filenum=index//20+1 # filenumber of multiplex 
            self.pixarr.loc[index]['PHIname']=self.basename.get()+str(filenum)
        filelist=self.pixarr.PHIname.unique()
        filelist=np.ndarray.tolist(filelist)
        # map/move beam non-sequentially to minimize local charging
        if self.scramblebool.get():
            areanums=np.ndarray.tolist(self.pixarr.Areanumber.unique())
            areanums=np.random.permutation(areanums)
            self.pixarr['Areanumber']=pd.Series(areanums)
            self.pixarr=self.pixarr.sort_values(['Areanumber'])
            # need to reassign subnumber and PHIname
            self.pixarr=self.pixarr.reset_index(drop=True)
            for index,row in self.pixarr.iterrows():
                self.pixarr=self.pixarr.set_value(index,'Subnumber', index%20)
                self.pixarr=self.pixarr.set_value(index,'PHIname', self.basename.get()+str(index//20+1))
        self.parent.pixviewer.pass_pixarray(self.pixarr)
        print('Square array created.')

    def save_pixarray(self):
        '''Menu/main lauched save of pixarray file (after linking with
        underlying data files '''
        if self.pixarr is None:
            return
        filelist=self.pixarr.PHIname.unique()
        filelist=np.ndarray.tolist(filelist)
        if self.areabool.get():
            for i, fname in enumerate(filelist):
                thisfile=self.pixarr[self.pixarr['PHIname']==fname]
                self.writeAESareas(thisfile, fname) # writes each list of 20 areas to separate .phi text file
            print('New spatial area files saved.. e.g.', fname)
        if self.atbool.get(): # option write to new Autotool file
            atframe=self.makeautotool(filelist)
            ATname='AT'+self.basename.get()+'.phi'
            self.writeautotool(atframe, ATname)
            print('Saving new autotool file', ATname)
        # Instead of C:\Temp copy multiplex and spatial areas files to Smartsoft settings folders
        mydir=self.parent.pixviewer.directory # get directory from pixviewer window
        if mydir==None: # in case image hasn't been loaded
            mydir=filedialog.askdirectory('Find data directory')
        fname=mydir+'/'+self.pixarrname.get()+'_pixarr.csv'
        self.pixarr.to_csv(fname, index=False)
        print('Pixarray file saved', fname)
    
    def writeAESareas(self, df, PHIname):
        ''' Optional write of stage positions to .phi spatial areas file in chunks of 25 positions
        Some weird encoding so just read and modify existing pff file '''    
        datastr=''
        datastr+='[SpatialArea]\nArea Count='
        datastr+=str(len(df)) # number of areas
        datastr+='\n'
        for i in range(0,len(df)):
            datastr+='Area Active '
            datastr+='%d' % i
            datastr+='=True\n'
        for i in range(0,len(df)):
            datastr+='Area Mode '
            datastr+='%d' % i
            datastr+='=Area\n'
        for i in range(0,len(df)):
            datastr+='Area Left '
            datastr+='%d' % i
            datastr+='='
            val=df.iloc[i]['X1']
            datastr+='%d' % val
            datastr+='\n'
        for i in range(0,len(df)):
            datastr+='Area Top '
            datastr+='%d' % i
            datastr+='='
            val=df.iloc[i]['Y1']
            datastr+='%d' % val
            datastr+='\n' 
        for i in range(0,len(df)):
            datastr+='Area Right '
            datastr+='%d' % i
            datastr+='='
            val=df.iloc[i]['X2']
            datastr+='%d' % val
            datastr+='\n'         
        for i in range(0,len(df)):
            datastr+='Area Bottom '
            datastr+='%d' % i
            datastr+='='
            val=df.iloc[i]['Y2']
            datastr+='%d' % val
            datastr+='\n'      
        for i in range(0,len(df)):
            datastr+='Area Width '
            datastr+='%d' % i
            datastr+='='
            val=df.iloc[i]['Width']
            datastr+='%d' % val
            datastr+='\n'      
        for i in range(0,len(df)):
            datastr+='Area Height '
            datastr+='%d' % i
            val=df.iloc[i]['Height']
            datastr+='='
            datastr+='%d' % val
            datastr+='\n' 
        # Write this chunk of files to .phi spatial areas file (done w/ file replace method since encoding is weird unknown type)
        filename=PHIname+'.phi'
        try:
            shutil.copyfile(AREAFILE, filename)
        except:
            mydir=filedialog.askdirectory(title='Select directory with spatial area file example.')
            shutil.copyfile(mydir+'/'+'spatial_areas_sample_min.phi', filename)
        for line in fileinput.input(filename, inplace=1):
            sys.stdout.write(datastr)
        
    def makeautotool(self, filelist, multacq='QM_multiplex.phi'):
        '''Generates df with Autotool commands and data values (for generating Autotool phi file)
        7/10 Modified with image reg insertion 
        
        spatial area files have to be in Auger's default directory for them..
        attempting load from elsewhere (i.e. from C:\Temp) doesn't seem to work
        
        kwarg: regint - interval for insertion of image registrations 
        multiplex breaks -- used in shifted situation (series of files at different shifts 
        multacq -- name of multiplex file to load... not currently implemented
        '''
        mycols=['Command','Data']
        atframe=pd.DataFrame(columns=mycols)
        atframe['Command']=['AES:Register Image','SEM:Photo']
        # drop load of first multiplex ... user can just load it before start

        ''' TODO multiplex break not yet implemented
        allows load of different multiplex file in middle of scans (in cases of
        spatially varying charging across mapped region)
        # CODE SECTION
        if 'multibreaks' in kwargs:
            multibreaks=kwargs.get('multibreaks',[])
            multinames=kwargs.get('multinames',[])
            multacq=multinames[0] # set to first shifted multiplex to load
            if 0 in multibreaks: # first multiplex file loaded w/ multacq
                multibreaks.remove(0)
        '''
        # Add first multiplex load (rest are below)
        ''' # just load multiplex manually (unless multibreaks is implemented)
        newrow=pd.DataFrame(index=np.arange(0,1), columns=mycols)
        newrow=newrow.set_value(0,'Command','AES:Load Multiplex Setting...')
        newrow=newrow.set_value(0,'Data', multacq)
        atframe=pd.concat([atframe,newrow], ignore_index=True)
        '''
        for i, file in enumerate(filelist):
            newrow=pd.DataFrame(index=np.arange(0,2), columns=mycols)
            newrow=newrow.set_value(0,'Command','AES:Load Area Define Setting...')
            newrow=newrow.set_value(0,'Data', file)
            newrow=newrow.set_value(1,'Command','AES:Multiplex Acquire')
            atframe=pd.concat([atframe,newrow], ignore_index=True)
        ''' TODO Multiplex breaks section 
        # Now add load of next shifted multiplex file 
        if 'multibreaks' in kwargs:
            if i in multibreaks:
                lindex=multibreaks.index(i)
                newrow=pd.DataFrame(index=np.arange(0,1), columns=mycols)
                newrow=newrow.set_value(0,'Command','AES:Load Multiplex Setting...')
                # multfile must be in settings/multiplex acquire (no file extensions in Autotool data cols)
                newrow=newrow.set_value(0,'Data', multinames[lindex+1].replace('.phi',''))
                atframe=pd.concat([atframe,newrow], ignore_index=True)
        '''
        # Tag on ending to autotool file
        newrow=pd.DataFrame(index=np.arange(0,1), columns=mycols)
        newrow['Command']=['SEM:Photo']
        atframe=pd.concat([atframe,newrow], ignore_index=True) # ending SEM photo
        newrow=pd.DataFrame(index=np.arange(0,1), columns=mycols)
        newrow['Command']=['SEM:SEM Load By Ref...']
        newrow['Data']=[BEAMDEFLECT]
        atframe=pd.concat([atframe,newrow], ignore_index=True) # add final beam deflection
        return atframe
    
    def writeautotool(self, atframe, atname):
        ''' Write of standard autotool loop for quantmap  
        atframe created by makeautotool... name is passed 
        weird encoding so just read and modify existing pff file '''

        datastr=''
        datastr+='[AutoTool]\nTask Count='
        datastr+=str(len(atframe)) # 
        datastr+='\n'
        for index, row in atframe.iterrows():
            datastr+='Task '
            datastr+='%d' % index
            command=atframe.loc[index]['Command']
            datastr+='='
            datastr+=command
            datastr+='\n'
        datastr+='Data Count='
        datastr+=str(len(atframe)) # 
        datastr+='\n'    
        for index, row in atframe.iterrows():
            datastr+='Data '
            datastr+='%d' % index
            datastr+='='
            val=atframe.loc[index]['Data']
            if str(val)!='nan':
                datastr+=str(val) # could be int in some cases
            datastr+='\n'
            # Write this chunk of files to .phi spatial areas file (done w/ file replace method since encoding is weird unknown type)
        shutil.copyfile(AREAFILE, atname)
        for line in fileinput.input(atname, inplace=1):
            sys.stdout.write(datastr)

    def plot_maps(self):
        ''' display 2D arrays of various types in mapviewer '''
        activeelems=[]
        plotmaps=[]
        title=''
        for i, tkbool in enumerate(self.tkelems):
            if tkbool.get():
                if self.plottype.get()=='Shiftmap':
                    if self.QMfile.shiftmaps[i] is not None:
                        activeelems.append(self.QMfile.elements[i])
                        plotmaps.append(self.QMfile.shiftmaps[i])                                
                        title='Peak shift'
                elif self.plottype.get()=='Amplmap':
                    if self.QMfile.amplmaps[i] is not None:
                        activeelems.append(self.QMfile.elements[i])
                        plotmaps.append(self.QMfile.amplmaps[i])                                
                        title='Peak amplitude'
                elif self.plottype.get()=='Elemmap':
                    if self.QMfile.amplmaps[i] is not None:
                        activeelems.append(self.QMfile.elements[i])
                        plotmaps.append(self.QMfile.elemmaps[i])                                
                        title='Element map'
        print("Plotting ", self.plottype.get(), "for elements", ".".join(activeelems))
        self.parent.mapviewer.replot_maps(activeelems, plotmaps, title)

class Multiplexopts():
    ''' Parent is GUImain, manages multiplex element and other parameter 
    for creation of QM multiplex file 
    also includes time estimate updates for display somewhere 
    based on existing QMmultiplex_setup (spyder launched tk version)
    '''
    def __init__(self, root, parent):
        self.root = root
        self.parent = parent
        # Load AESquantparams
        self.aesquantparams = pd.DataFrame()
        self.loadAESquantparams()
        self.elements=[]
        self.elemparams=[] # standard QM vals from aesquantparams
        # common multiplex scan parameters 
        self.dwell=tk.IntVar()
        self.dwell.set(20)
        self.numcycles=tk.IntVar()
        self.numcycles.set(3)
        self.peakshift=tk.IntVar() # globally appliable peak shift (and no local shift option)
        self.peakshift.set(0)
        self.peakwidth=tk.IntVar() # globally appliable peak shift
        self.peakwidth.set(7)
        # self.regint=tk.IntVar() # image registration interval (if done w/in multiplex)
        # self.regint.set(0)
        self.phiname=tk.StringVar()
        self.phiname.set('QMmultiplex.phi')
        self.timeest=tk.DoubleVar()

        # lists of tkvars with multiplex parameters 
        self.peaks=[]
        self.sweeps=[]
        self.widths=[]
        self.lowers=[] # lower ev of scan range
        self.uppers=[]

        # Shows all scan regions within multiplex
        self.top_frame= tk.Frame(self.root)
        # show/alter common multiplex parameters
        self.param_frame = tk.Frame(self.top_frame, pady=10)
        self.param_frame.pack(side=tk.LEFT,fill=tk.X,expand=1)
        # button frame below global params
        self.button_frame = tk.Frame(self.top_frame, pady=10)
        self.button_frame.pack(side=tk.LEFT,fill=tk.X,expand=1)
        self.top_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        self.mult_frame = tk.Frame(self.root, pady=10)
        self.mult_frame.pack(side=tk.BOTTOM,fill=tk.X,expand=1)
        
        tk.Button(self.button_frame, text='Choose elements', 
                  command=self.chooseelems).pack(side=tk.TOP,fill=tk.X,expand=1)
        tk.Button(self.button_frame, text='Apply global shift/width', 
                  command=self.applyglobals).pack(side=tk.TOP,fill=tk.X,expand=1)
        # Add a button label
        tk.Label(self.button_frame, text='Make local peak changes').pack(side=tk.TOP,fill=tk.X,expand=1)
        tk.Button(self.button_frame, text='Change ranges using widths', 
                  command=self.updatepeaks).pack(side=tk.TOP,fill=tk.X,expand=1)
                  command=self.updatewidths).pack(side=tk.TOP,fill=tk.X,expand=1)
        tk.Button(self.button_frame, text='Recalc/update', 
                  command=self.recalcupdate).pack(side=tk.TOP,fill=tk.X,expand=1)
        tk.Button(self.button_frame, text='Save multiplex', 
                  command=self.save_multiplex).pack(side=tk.TOP,fill=tk.X,expand=1)
        self.display_multparams()
    
    def save_multiplex(self):
        ''' Saves new multiplex to file '''
        multdf=self.makemultdf()
        # Save modified multiplex scan for QM
        self.writemultiplex(multdf)
    
    def makemultdf(self):
        ''' Reconstruct dataframe holding altered multiplex scan parameters 
        then feed to writemultiplex'''       
        mycols=['AtomNum', 'Elem', 'Active', 'Sweeps', 'EVstep', 'Lower', 'Upper',
           'Range', 'Lowpeak', 'Peak', 'Hipeak', 'Back1', 'Back2']
        multdf=pd.DataFrame(columns=mycols)
        atnos=[]
        regs=[]
        peaks=self.convert_tklist(self.peaks)
        sweeps=self.convert_tklist(self.sweeps)
        widths=self.convert_tklist(self.widths)
        lowers=self.convert_tklist(self.lowers)
        uppers=self.convert_tklist(self.uppers)
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(self.elemparams):
            atnos.append(atno)
            atnos.append(atno)
            atnos.append(atno)
            regs.append(elem+'L')
            regs.append(elem)
            regs.append(elem+'H')
        multdf['AtomNum']=atnos
        multdf['Elem']=regs
        multdf['Active']='Y'
        multdf['Sweeps']=sweeps
        multdf['Sweeps']=sweeps
        multdf['EVstep']=1.0
        multdf['Peak']=peaks
        multdf['Lower']=lowers
        multdf['Upper']=uppers
        multdf['Lowpeak']=multdf['Lower']+2
        multdf['Hipeak']=multdf['Upper']-2
        multdf['Back1']=multdf['Lower']
        multdf['Back2']=multdf['Upper']
        # convert half-widths to full widths
        widths=[2*i+1 for i in widths]
        multdf['Range']=widths
        # Eliminate any overlaps in scan range (remove from L or H not from main peak)
        # overlaps should already be gone b/c of tk gui methods 
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(self.elemparams):
            lowend=multdf.iloc[i]['Upper']
            mainstart=multdf.iloc[i+1]['Lower']
            mainend=multdf.iloc[i+1]['Upper']
            highstart=multdf.iloc[i+2]['Lower']
            # print(str(lowend), str(mainstart), str(mainend),str(highstart))
            if mainstart<lowend:
                multdf=multdf.set_value(multdf.index[i],'Upper', mainstart-1)
            if highstart<mainend:
                multdf=multdf.set_value(multdf.index[i+2],'Lower', mainend+1)     
        return multdf
    
    def writemultiplex(self, multdf):
        ''' Write of multiplex settings file (max 20 regions) after interactive param setting
        Some weird encoding so just read and modify existing pff file 
        image registration choices are cycles or areas
        kwargs:
        regmode - Areas or  (if not present, image registration done using autotool not in multiplex)
        reginterval - 2 (or whatever)'''
        phiname=self.phiname.get()
        dwell=self.dwell.get()
        numcycles=self.numcycles.get()
        
        datastr='' # long data string for multiplex file
        datastr+='[MultiplexAcq]\n'
        ''' image registration option within multiplex (better in autotool)
        # (also can be done in autotool at lower frequencies)
        if int(self.regint.get())>0:  # ensure this is integer
            regmode='Areas'
            reginterval=int(self.regint.get())
            datastr+='Register Image=True\nImage Registration Interval='
            datastr+=str(reginterval)
            datastr+='\n'
            datastr+='Image Registration Mode='
            datastr+=regmode
        else:
        '''
    
        datastr+='Register Image=False'
        # likely can skip other params if this is false
        datastr+='\nTime Per Step (ms)='
        datastr+='%0.1f' % dwell
        datastr+='\nNegative Values=Allow\nNumber of Regions='
        datastr+=str(len(multdf)) # number of multiplex regions
        datastr+='\nAtomic Number List Count='
        datastr+=str(len(multdf))
        datastr+='\n'
        for i in range(0,len(multdf)):
            datastr+='Atomic Number List '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['AtomNum']
            datastr+='%d' % val
            datastr+='\n'
        datastr+='Element List Count='
        datastr+=str(len(multdf))
        datastr+='\n'
        for i in range(0,len(multdf)):
            datastr+='Element List '
            datastr+='%d' % i
            datastr+='='
            strval=multdf.iloc[i]['Elem']
            datastr+=strval
            datastr+='\n'
        datastr+='Active Count='
        datastr+=str(len(multdf))
        datastr+='\n'
        for i in range(0,len(multdf)):
            datastr+='Active '
            datastr+='%d' % i
            datastr+='=True\n'  # won't be present in Df if not active
        datastr+='Sweeps Count='
        datastr+=str(len(multdf))
        datastr+='\n' 
        for i in range(0,len(multdf)):
            datastr+='Sweeps '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Sweeps']
            datastr+='%d' % val
            datastr+='\n'
        datastr+='Lower Acq Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Lower Acq '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Lower']
            datastr+='%0.1f' % val # float with tenths place precision
            datastr+='\n'
        datastr+='Upper Acq Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Upper Acq '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Upper']
            datastr+='%0.1f' % val
            datastr+='\n'        
        datastr+='Acq Range Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Acq Range '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Range']
            datastr+='%0.1f' % val
            datastr+='\n'    
        datastr+='Lower Analysis Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Lower Analysis '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Lowpeak']
            datastr+='%0.1f' % val
            datastr+='\n'     
        datastr+='Upper Analysis Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Upper Analysis '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Hipeak']
            datastr+='%0.1f' % val
            datastr+='\n'    
        datastr+='eV Per Step Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='eV Per Step '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['EVstep']
            datastr+='%0.1f' % val
            datastr+='\n'   
        datastr+='Peak Energy Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Peak Energy '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Peak']
            datastr+='%0.1f' % val
            datastr+='\n'  
        datastr+='Background 1 Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Background 1 '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Back1']
            datastr+='%0.1f' % val
            datastr+='\n'  
        datastr+='Background 2 Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Background 2 '
            datastr+='%d' % i
            datastr+='='
            val=multdf.iloc[i]['Back2']
            datastr+='%0.1f' % val
            datastr+='\n'  
        datastr+='Number Of Channels Count='
        datastr+=str(len(multdf))
        datastr+='\n'     
        for i in range(0,len(multdf)):
            datastr+='Number Of Channels '
            datastr+='%d' % i
            datastr+='=1 thru 8\n'
        datastr+='Number of Cycles='
        datastr+=str(numcycles)
        datastr+='\n'
        # Write this chunk of files to .phi spatial areas file (done w/ file replace method since encoding is weird unknown type)
        try:
            shutil.copyfile(AREAFILE,phiname)
        except:
            mydir=filedialog.askdirectory(title='Select directory with spatial area file example.')
            shutil.copyfile(mydir+'/'+'spatial_areas_sample_min.phi',phiname)
        for line in fileinput.input(phiname, inplace=1):
            sys.stdout.write(datastr)
        print(phiname,' saved to current working directory')
                   
    def convert_tklist(self, tklist):
        ''' convert list of tk inter variables to normal python list '''
        newlist=[]
        for i, val in enumerate(tklist):
            newlist.append(val.get())
        return newlist
    
    def display_multparams(self):
        ''' Display of global params for multiplex scan '''
            # Initialize w/ shift of zero
        # elem name, lower region center, peak center, upper reg center, # sweeps
        shift=0
        self.elemparams=self.getQMparams(shift)
        # Display of array size and other params
        tk.Label(self.param_frame, text='Dwell time (ms)').grid(row=0, column=6)
        tk.OptionMenu(self.param_frame, self.dwell, 5, 10, 20) # drop-down 
        tk.Entry(self.param_frame, textvariable=self.dwell).grid(row=0, column=7)
        tk.Label(self.param_frame, text='# of cycles').grid(row=1, column=6)
        tk.Entry(self.param_frame, textvariable=self.numcycles).grid(row=1, column=7)
        tk.Label(self.param_frame, text='Peak shift').grid(row=2, column=6)
        tk.Entry(self.param_frame, textvariable=self.peakshift).grid(row=2, column=7)
        tk.Label(self.param_frame, text='Half width').grid(row=3, column=6)
        tk.Entry(self.param_frame, textvariable=self.peakwidth).grid(row=3, column=7)
        tk.Label(self.param_frame, text='Multiplex name').grid(row=4, column=6)
        tk.Entry(self.param_frame, textvariable=self.phiname).grid(row=4, column=7)
        tk.Label(self.param_frame, text='Time (hrs)').grid(row=5, column=6)
        tk.Entry(self.param_frame, textvariable=self.timeest).grid(row=5, column=7)
    
        # tk var not needed for peak/ element regions ... just display

    def init_peaks(self):
        ''' all tkvars for multiplex scan params '''
        print('initializing peak parameters for elements',",".join(self.elements))
        # clear all existing
        self.peaks=[]
        self.sweeps=[]
        self.widths=[]
        self.lowers=[] # lower ev of scan range
        self.uppers=[]
        # Create lists of tk string variables
        for i in range(0,3*len(self.elements)):
            self.peaks.append(tk.DoubleVar()) # center of scan window
            self.lowers.append(tk.DoubleVar()) # lower limit of scan window
            self.uppers.append(tk.DoubleVar())
            self.sweeps.append(tk.IntVar())
            self.widths.append(tk.IntVar())
        # initialize based on elemparams
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(self.elemparams):
            self.peaks[3*i].set(int(self.elemparams[i][1])) # lower window for element
            self.lowers[3*i].set(int(self.elemparams[i][1]-7))
            self.uppers[3*i].set(int(self.elemparams[i][1]+7))
            self.widths[3*i].set(7)
            self.sweeps[3*i].set(1)
            self.peaks[3*i+1].set(int(self.elemparams[i][2])) # peak itself
            self.lowers[3*i+1].set(int(self.elemparams[i][2]-7))
            self.uppers[3*i+1].set(int(self.elemparams[i][2]+7))
            self.widths[3*i+1].set(7)
            self.sweeps[3*i+1].set(int(self.elemparams[i][4]))
            self.peaks[3*i+2].set(int(self.elemparams[i][3])) # upper window for element
            self.lowers[3*i+2].set(int(self.elemparams[i][3]-7))
            self.uppers[3*i+2].set(int(self.elemparams[i][3]+7))
            self.widths[3*i+2].set(7)
            self.sweeps[3*i+2].set(1)
    
    def display_peaks(self):
        ''' Display all multiplex peak scan params in mult_frame '''
        print('display peaks called')
        for child in self.mult_frame.winfo_children():
            child.destroy()
        # First col are energy region/element labels (not alterable)
        tk.Label(self.mult_frame, text='Region').grid(row=0, column=0)
        tk.Label(self.mult_frame, text='Peak').grid(row=0, column=1)
        tk.Label(self.mult_frame, text='Lower').grid(row=0, column=2)
        tk.Label(self.mult_frame, text='Upper').grid(row=0, column=3)
        tk.Label(self.mult_frame, text='Sweeps').grid(row=0, column=4)
        tk.Label(self.mult_frame, text='Half-width').grid(row=0, column=5)
        # Display peak params (centers/ lower ev/ upper ev/ width/
        rownum=1 # place peaks after column header row
        print('Elements are',",".join(self.elements))
        for i, elem in enumerate(self.elements):
            # first list element regions (i.e. SiL, Si, SiH)
            tempstr=elem+'L'
            tk.Label(self.mult_frame, text=tempstr).grid(row=rownum, column=0)
            tk.Entry(self.mult_frame, textvariable=self.peaks[3*i], width=7).grid(row=rownum, column=1)
            tk.Entry(self.mult_frame, textvariable=self.lowers[3*i], width=7).grid(row=rownum, column=2)
            tk.Entry(self.mult_frame, textvariable=self.uppers[3*i], width=7).grid(row=rownum, column=3)
            tk.Entry(self.mult_frame, textvariable=self.sweeps[3*i], width=7).grid(row=rownum, column=4)
            tk.Entry(self.mult_frame, textvariable=self.widths[3*i], width=7).grid(row=rownum, column=5)
            rownum+=1
            # now handle peak line itself
            tk.Label(self.mult_frame, text=elem).grid(row=rownum, column=0)
            tk.Entry(self.mult_frame, textvariable=self.peaks[3*i+1], width=7).grid(row=rownum, column=1)
            tk.Entry(self.mult_frame, textvariable=self.lowers[3*i+1], width=7).grid(row=rownum, column=2)
            tk.Entry(self.mult_frame, textvariable=self.uppers[3*i+1], width=7).grid(row=rownum, column=3)
            tk.Entry(self.mult_frame, textvariable=self.sweeps[3*i+1], width=7).grid(row=rownum, column=4)
            tk.Entry(self.mult_frame, textvariable=self.widths[3*i+1], width=7).grid(row=rownum, column=5)
            rownum+=1
            tempstr=elem+'H'
            tk.Label(self.mult_frame, text=tempstr).grid(row=rownum, column=0)
            tk.Entry(self.mult_frame, textvariable=self.peaks[3*i+2], width=7).grid(row=rownum, column=1)
            tk.Entry(self.mult_frame, textvariable=self.lowers[3*i+2], width=7).grid(row=rownum, column=2)
            tk.Entry(self.mult_frame, textvariable=self.uppers[3*i+2], width=7).grid(row=rownum, column=3)
            tk.Entry(self.mult_frame, textvariable=self.sweeps[3*i+2], width=7).grid(row=rownum, column=4)
            tk.Entry(self.mult_frame, textvariable=self.widths[3*i+2], width=7).grid(row=rownum, column=5)
            rownum+=1
    
    def getQMparams(self, shift):
        ''' retrieve direct peak location, low & high energy regions (background below and 
        above peak), and default # of sweeps for passed elements '''
        self.elemparams=[]
        for i, elem in enumerate(self.elements):
            match=self.aesquantparams[self.aesquantparams['element']==elem]
            if len(match)==1:
                lowback=round(match.iloc[0]['QMlow']+shift,1)
                hiback=round(match.iloc[0]['QMhigh']+shift,1)
                # Direct peak location val is relative to negpeak
                peak=round(match.iloc[0]['negpeak']+match.iloc[0]['integpeak']+shift,1)
                # Default # of sweeps for actual peak
                sweeps=int(match.iloc[0]['QMsweeps'])
                atno=int(match.iloc[0]['atno']) # atomic number
                self.elemparams.append([elem, lowback, peak, hiback, sweeps, atno])
            else:
                print("Couldn't find element", elem)

    def recalc(self):
        ''' Update time estimates after any value changes (not tied to button) but called from 
        within every method
        '''
        # Get array size from pixarropts
        print(self.parent.pixopts.arrsize.get())
        # Acquisition time est. in hours
        numchannels=0
        for i in range(0,len(self.widths)):
            # 2 * half-width * # sweeps * num cycles is # of ev channels scanned
            numchannels+=int(self.numcycles.get()) * int(self.sweeps[i].get()) * (int(self.widths[i].get())*2+1)
        # time for single multiplex (per area)
        singletime=numchannels*(self.dwell.get())/1000
        numareas=self.parent.pixopts.arrsize.get()**2
        fulltime= singletime*numareas/3600
        tempstr = "%.2f" % fulltime
        self.timeest.set(tempstr)

    def updatewidth(self):
        ''' Use current lowers/uppers values to calculate each local width '''
        # first reset centers/ranges to initial values 
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(self.elemparams):
            if int(self.uppers[3*i].get() - self.lowers[3*i].get())/2 > 0:
                self.widths[3*i].set(int(self.uppers[3*i].get() - self.lowers[3*i].get())/2 )
            else:
                self.widths[3*i].set(0)
            if int(self.uppers[3*i+1].get() - self.lowers[3*i+1].get())/2 >0:   
                self.widths[3*i+1].set(int(self.uppers[3*i+1].get() - self.lowers[3*i+1].get())/2)
            else:
                self.widths[3*i+1].set(0)
            if int(self.uppers[3*i+2].get() - self.lowers[3*i+2].get())/2 >0:
                self.widths[3*i+2].set(int(self.uppers[3*i+2].get() - self.lowers[3*i+2].get())/2)
            else:
                self.widths[3*i+2].set(0)
    
    def adjustpeakcenter(self):
        ''' Squeezing of background regions changes effective center  '''
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(self.elemparams):
            self.peaks[3*i].set((self.uppers[3*i].get() + self.lowers[3*i].get())/2) # lower window for element
            self.peaks[3*i+2].set((self.uppers[3*i+2].get() + self.lowers[3*i+2].get())/2)
 
    def checkoverlaps(self):
        ''' After any changes adjust lower and upper backgrounds if overlapping w/ main region '''
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(self.elemparams):
            if self.uppers[3*i].get() >= self.lowers[3*i+1].get():
                self.uppers[3*i].set(self.lowers[3*i+1].get() - 1)
            if self.lowers[3*i+2].get() <= self.uppers[3*i+1].get():
                self.lowers[3*i+2].set(self.uppers[3*i+1].get() + 1)
        # TODO width of OL or OH cannot be negative (constraint in updatewidth?)
        # Changing lowers/ uppers change resulting width
        self.updatewidth()
        # changing lowers uppers for lowback and highback regions can change peakcenter
        self.adjustpeakcenter()
    
    def applyglobals(self):
        ''' Update each peak positions and widths based on global peak shift and 
        global width 
        button launched'''
        for i, elem in enumerate(self.elemparams):
            self.peaks[3*i].set(int(self.elemparams[i][1]+self.peakshift.get()))
            self.lowers[3*i].set(int(self.elemparams[i][1]-self.peakwidth.get() +self.peakshift.get()))
            self.uppers[3*i].set(int(self.elemparams[i][1] + self.peakwidth.get() +self.peakshift.get()))
            self.widths[3*i].set(self.peakwidth.get())
            self.peaks[3*i+1].set(int(self.elemparams[i][2]+self.peakshift.get()))
            self.lowers[3*i+1].set(int(self.elemparams[i][2]-self.peakwidth.get() +self.peakshift.get()))
            self.uppers[3*i+1].set(int(self.elemparams[i][2] + self.peakwidth.get() +self.peakshift.get()))
            self.widths[3*i+1].set(self.peakwidth.get())
            self.peaks[3*i+2].set(int(self.elemparams[i][3]+self.peakshift.get()))
            self.lowers[3*i+2].set(int(self.elemparams[i][3]-self.peakwidth.get() +self.peakshift.get()))
            self.uppers[3*i+2].set(int(self.elemparams[i][3] + self.peakwidth.get() +self.peakshift.get()))
            self.widths[3*i+2].set(self.peakwidth.get())
        self.checkoverlaps()
        self.recalc() # update time estimate
    
    def updatepeaks(self):
        ''' Use local widths (and global shift) to update peak centers and ranges 
        no vals for local energy shift but could change peaks/lowers/uppers '''
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(self.elemparams):
            self.peaks[3*i].set(int(self.elemparams[i][1] + self.peakshift.get() )) # lower window for element
            self.lowers[3*i].set(int(self.elemparams[i][1]+ self.peakshift.get() - self.widths[3*i].get()))
            self.uppers[3*i].set(int(self.elemparams[i][1]+ self.peakshift.get() + self.widths[3*i].get()))
            self.peaks[3*i+1].set(int(self.elemparams[i][2] + self.peakshift.get() )) # peak itself
            self.lowers[3*i+1].set(int(self.elemparams[i][2] + self.peakshift.get() - self.widths[3*i+1].get()))
            self.uppers[3*i+1].set(int(self.elemparams[i][2] + self.peakshift.get() + self.widths[3*i+1].get()))
            self.peaks[3*i+2].set(int(self.elemparams[i][3]  + self.peakshift.get())) # upper window for element
            self.lowers[3*i+2].set(int(self.elemparams[i][3]+ self.peakshift.get() - self.widths[3*i+2].get()))
            self.uppers[3*i+2].set(int(self.elemparams[i][3]+ self.peakshift.get() + self.widths[3*i+2].get()))
        
        self.checkoverlaps()
        # acquisition time est. in hours
        self.recalc()

    def updatewidths(self):
        ''' Use current lowers/uppers values to calculate each local width '''
        # first reset centers/ranges to initial values 
        for i, [elem, lowreg, peak, hireg, sweep, atno] in enumerate(self.elemparams):
            self.widths[3*i].set(int(self.uppers[3*i].get() - self.lowers[3*i].get())/2 )
            self.widths[3*i+1].set(int(self.uppers[3*i+1].get() - self.lowers[3*i+1].get())/2)
            self.widths[3*i+2].set(int(self.uppers[3*i+2].get() - self.lowers[3*i+2].get())/2)
        self.checkoverlaps()
        # acquisition time est. in hours
        self.recalc()

    def recalcupdate(self):
        ''' Update time estimates, widths, adjust overlap boundaries
        '''
        self.checkoverlaps()
        self.recalc()
        
    def chooseelems(self):
        ''' Select elements using pop-up toplevel;  all available peaks from AESquantparams.csv '''
        tframe = tk.Toplevel(self.root)
        # Subset of elements selected (on) by default
        elemdict={'S':1,'C':1,'Ti':1,'O':1,'Fe1':1,'Fe2':1,'Na':1,'Mg':1,'Al':1,'Si':1,'Fe':1,'Ca':1}
        preset1={'C':1,'O':1,'Si':1,'Fe':1,'Mg':1,'Ca':1}
        preset2={'O':1,'Mg':1,'Si':1,'Fe':1}
        
        # All available elements/peaks are those with entries in Aesquantparams.csv
        elems=np.ndarray.tolist(self.aesquantparams.element.unique()) 
        varlist=[] # list of tkinter IntVars
        for i, col in enumerate(elems): # set up string variables
            varlist.append(tk.IntVar())
            val=elemdict.get(col,0) # set to 1 or 0 based on above default dictionary
            varlist[i].set(val) # set default value based on elemdict
            
        tk.Label(tframe, text='Select elements for plotting or quant').grid(row=0,column=0)
        
        def clearall():
            ''' Set all tkinter vars to zero ''' 
            for i, col in enumerate(elems): # set up string variables
                varlist[i].set(0) # set default value based on elemdict
        def choose1():
            ''' Have available preset defaults and adjust checkbox values '''
            # preset1={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1,'Fe2':1}
            # Still have to pass these through as tkinter ints
            for i, col in enumerate(elems): # set up string variables
                val=preset1.get(col,0) # set to 1 or 0 based on above default dictionary
                varlist[i].set(val) # set default value based on elemdict
        def choose2():
            ''' Have available preset defaults and adjust checkbox values '''
            # preset2={'S':1,'Mg':1,'Si':1,'Fe':1,'Ca':1,'Fe2':1}
            # Still have to pass these through as tkinter ints
            for i, col in enumerate(elems): # set up string variables
                val=preset2.get(col,0) # set to 1 or 0 based on above default dictionary
                varlist[i].set(val) # set default value based on elemdict
        
        def selectelems():
            '''Choose checked elems and close popup '''
            self.elements=[] # list of strings with plot number and x or y
            for i, val in enumerate(varlist): # result in normal string, not tkinter StringVar
                if val.get()==1: # this element is active
                    self.elements.append(elems[i]) # add element if box is checked
            print('Selected elements are',",".join(self.elements))
            self.getQMparams(int(self.peakshift.get())) # get QM params with current global peak shift 
            self.init_peaks()
            self.display_peaks() # now show peaks with chosen elements
            tframe.destroy()
                
        for i, col in enumerate(elems):
            # choose row, col grid position (starting row 1)
            thisrow=i%3+1 # three column setup
            thiscol=i//3
            ent=tk.Checkbutton(tframe, text=elems[i], variable=varlist[i])
            ent.grid(row=thisrow, column=thiscol)
        # Add preset 1 button (defined above)
        els=list(preset1)
        mystr=', '.join(els)
        c=tk.Button(tframe, text=mystr, command=choose1)
        lastrow=len(elems)%3+2
        c.grid(row=lastrow, column=0)
        # Add preset 2 button
        els=list(preset2)
        mystr=', '.join(els)
        d=tk.Button(tframe, text=mystr, command=choose2)
        lastrow=len(elems)%3+3
        d.grid(row=lastrow, column=0)
        # clear all 
        e=tk.Button(tframe, text='Clear all', command=clearall)
        lastrow=len(elems)%3+4
        e.grid(row=lastrow, column=0)
        g=tk.Button(tframe, text='Done', command=selectelems)
        lastrow=len(elems)%3+7
        g.grid(row=lastrow, column=0)
        # Do further initialization once toplevel is destroyed
        
    def loadAESquantparams(self):
        ''' Loads standard values of Auger quant parameters 
        TODO what about dealing with local shifts '''
        # Checkbutton option for local (or standard) AESquantparams in file loader?
        try:
            self.aesquantparams=pd.read_csv(AESQUANTPARAMFILE, encoding='utf-8')
        except:
            dir=filedialog.askdirectory()
            self.aesquantparams=pd.read_csv(dir+'/'+'aesquantparams.csv', encoding='utf-8')
        print('AESquantparams loaded')
