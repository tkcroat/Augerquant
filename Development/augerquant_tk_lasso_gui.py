# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 15:27:05 2017

@author: tkc
"""
import os
import pandas as pd
import sys
import glob
import numpy as np
import tkinter as tk
import tkinter.messagebox as tkmess
from tkinter import filedialog
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.widgets import Lasso
from PIL import ImageTk,Image

PLOT_SIZE = (8,5)
MPL_STYLE = {
    "text.color":"lightblue",
    "axes.labelcolor":"lightblue",
    "axes.edgecolor":"black",
    "axes.facecolor":"0.4",
    "xtick.color": "lightblue",
    "ytick.color": "lightblue",
    "figure.facecolor":"black",
    "figure.edgecolor":"black",
    "text.usetex":False
}
mpl.rcParams.update(MPL_STYLE)

#-----------Misc.-------------#
DEFAULT_DIRECTORY = os.getcwd()

class GUIMain():
    ''' Main container for plotter, options (at right), and fileloader (bottom)  '''
    def __init__(self,root):
        self.root = root
        self.root.wm_title("Auger spectra viewer ")
        self.top_frame = tk.Frame(self.root)
        self.top_frame.pack(side=tk.TOP) 
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM)
        self.plot_frame = tk.Frame(self.top_frame)
        self.plot_frame.pack(side=tk.LEFT)
        self.options_frame = tk.Frame(self.top_frame) 
        self.options_frame.pack(side=tk.LEFT)
        self.loader_frame = tk.Frame(self.bottom_frame)
        self.loader_frame.pack(side=tk.LEFT,fill=tk.BOTH)
        self.plotter    = GUIPlotter(self.plot_frame,self)
        self.options    = GUIOptions(self.options_frame,self)
        self.loader  = GUIFileloader(self.loader_frame,self)
        
class NavSelectToolbar(NavigationToolbar2TkAgg): 
    ''' Custom matplotlib toolbar w/ lasso and picker called by GUIplotter '''
    def __init__(self, canvas,root,parent):
        self.canvas = canvas
        self.root   = root
        self.parent = parent
        NavigationToolbar2TkAgg.__init__(self, canvas,root)
        self.lasso_button = self._custom_button(text="lasso",command=lambda: self.lasso(
                lambda inds: self.parent.multi_select_callback(inds),"lasso"))
        self.pick_button = self._custom_button(text="select",command=lambda: self.picker(
                lambda ind: self.parent.single_select_callback(ind),"select"))

    def _custom_button(self, text, command, **kwargs):
        button = tk.Button(master=self, text=text, padx=2, pady=2, command=command, **kwargs)
        button.pack(side=tk.LEFT,fill="y")
        return button

    def contains_points(self,verts,callback):
        ''' '''
        xys = self.parent.xy
        if xys is not None:
            p = mpl.path.Path(verts)
            ind = [ii for ii,xy in enumerate(xys) if p.contains_point(xy)]
            callback(ind)
        self.canvas.draw_idle()
        del self.lasso_obj
            
    def press_lasso(self,event,callback):
        ''' button press with lasso active returning lasso object '''
        if event.inaxes is None: return
        self.lasso_obj = Lasso(event.inaxes, (event.xdata, event.ydata),
                           lambda verts: self.contains_points(verts,callback))

    def press_picker(self,event,callback):
        ind = event.ind
        if ind: callback(event.ind[0])
        
    def _disconnect_all_ids(self):
        self.canvas.widgetlock.release(self)
        if self._idPress is not None:
            self._idPress = self.canvas.mpl_disconnect(self._idPress)
            self.mode = ''
        if self._idRelease is not None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
            self.mode = ''
        
    def lasso(self,callback,msg):
        ''' Activated by lasso menu bar button '''
        self._active = "LASSO"
        self._disconnect_all_ids()
        self._idPress = self.canvas.mpl_connect(
            'button_press_event', lambda event: self.press_lasso(event,callback))
        self.mode = msg
        self.canvas.widgetlock(self)
        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self._active)
        self.set_message(self.mode)
        
    def picker(self,callback,msg):
        ''' Activated by picker menu bar button  '''
        self._active = "PICKER"
        self._disconnect_all_ids()
        self._idPress = self.canvas.mpl_connect(
            'pick_event', lambda event: self.press_picker(event,callback))
        self.mode = msg
        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self._active)
        self.set_message(self.mode)

class GUIPlotter():
    def __init__(self,root,parent):
        self.root = root
        self.parent = parent
        self.figure = mpl.figure.Figure(figsize=PLOT_SIZE, dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(bottom=0.15,right=0.95,top=0.95)
        self.canvas = FigureCanvasTkAgg(self.figure,self.root)
        self.toolbar = NavSelectToolbar(self.canvas,self.root,self)
        self.toolbar.update()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        # placeholder for important data attributes
        self.data_manager = None
        self.spelist= None
        self.Augerparamlog
        self.active_viewers = {}
        self.toolbar.pick_button.invoke()
        self.canvas.show()

    def set_data_manager(self,data_manager):
        ''' method called by guifileloader ''' 
        # Data_manager ... list of pandas dataframes 
        self.data_manager = data_manager

    def plot(self,**kwargs):
        if self.data_manager is None:return
        self.data_manager.toggle_types(self.parent.options.toggle_buttons)
        self.ax.cla()
        x_field  = self.parent.options.xaxis_key_var.get()
        y_field  = self.parent.options.yaxis_key_var.get()
        xscale = self.parent.options.xaxis_scale_var.get()
        yscale = self.parent.options.yaxis_scale_var.get()
        x = self.data_manager.cdata[x_field]
        y = self.data_manager.cdata[y_field]
        if x.min() <= 0.0 and xscale == "log":
            self.parent.options.xaxis_scale_check.invoke()
            xscale = self.parent.options.xaxis_scale_var.get()
        if y.min() <= 0.0 and yscale == "log":
            self.parent.options.yaxis_scale_check.invoke()
            yscale = self.parent.options.yaxis_scale_var.get()
        self.xy = np.vstack((x,y)).transpose()
        self.ax.set_xscale(xscale)
        self.ax.set_yscale(yscale)
        self.ax.set_xlabel(x_field)
        self.ax.set_ylabel(y_field)
        self.current_plot = self.ax.scatter(x,y,picker=True,**kwargs)
        self.update_colors()
        self.canvas.show()

    def multi_select_callback(self,inds):
        ''' Linked to lasso_button in NavSelectToolbar.. can launch GUI multiviewer 
        inds are lasso selected points '''
        if self.data_manager is None:return
        if not inds:return
        print('Lasso launcher called.')
        ''' TODO finish processing
        self.launch_multi_viewer(inds)
        '''
            
    def single_select_callback(self,ind):
        ''' Linked to pick_button in NavSelectToolbar.. can launch GUIviewersingle '''
        if self.data_manager is None:return
        print('point selector called')
        print(ind) # indices returned by mpl point selector
        ''' TODO finish processing of data selector
        self.launch_viewer(ind)
        '''
        
    def launch_viewer(self,ind):
        ''' Placeholder for some process launched by single point data selector '''
        if self.data_manager is None:return
        if ind in self.active_viewers.keys():
            self.active_viewers[ind].root.destroy()
        self.active_viewers[ind] = GUIViewerSingle(self,ind)
    
class GUIFileloader():
    ''' Picks directory and loads main Auger param files '''
    def __init__(self,root,parent):
        self.root = root
        self.parent = parent
        self.top_frame = tk.Frame(self.root)
        self.top_frame.pack(side=tk.TOP,anchor=tk.W)
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM,fill=tk.BOTH,expand=1)
        tk.Label(self.top_frame,text="Directory:",padx=8,pady=2,height=1).pack(side=tk.LEFT,anchor=tk.W)
        self.directory_entry = tk.Entry(self.top_frame,width=90,bg="lightblue",
                                    fg="black",highlightcolor="lightblue",insertbackground="black",
                                    highlightthickness=2)
        self.directory_entry.pack(side=tk.LEFT,fill=tk.BOTH,expand=1,anchor=tk.W)
        self.directory_entry.insert(0,DEFAULT_DIRECTORY)
        tk.Button(self.top_frame,text="Browse",command=self.launch_dir_finder
                  ).pack(side=tk.LEFT,fill=tk.BOTH,expand=1,anchor=tk.W)
        self.load_button = tk.Button(self.bottom_frame,text="Load Auger spectra",width=60,
                                       command=self.load_AESdataset)
        self.load_button.pack(side=tk.BOTTOM,expand=1,anchor=tk.CENTER)
        
    def load_AESdataset(self):
        ''' Load standard AES files (paramlogwith data returned to a DataManager '''
        directory = self.directory_entry.get()
        AESdata = AESdataset(directory)
        # populate spectra selector checkboxes (in GUIoptions)
        self.parent.options.populate_specselector(AESdata.spelist)
        # data_manager = DataManager(output)
        # self.parent.plotter.set_data_manager(AESdata)
        
    
    def launch_dir_finder(self):
        directory = filedialog.askdirectory()
        self.directory_entry.delete(0,tk.END)
        self.directory_entry.insert(0,directory)

class AESspectrafile():
    ''' Single instance of AES spectra file created from row of spelist (child of AESdataset) '''
    def __init__(self, parent, row):
        # open files from directory arg
        self.parent=parent
        self.path=self.parent.path # same path as AESdataset parent
        self.numareas=row.Numareas
        self.evbreaks=row.Evbreaks
        
    def savecsv():
        ''' Save any changes to underlying csv file '''

    
class AESdataset():
    ''' Loads all dataframes with Auger parameters from current project folder 
    holds param files, backfit files, smdifpeaklogs, etc.  '''
    def __init__(self, directory):
        # open files from directory arg
        self.path=directory
        os.chdir(directory)
        self.Augerparamlog, self.spelist, self.Smdifpeakslog, self.Backfitlog, self.Integquantlog, self.AESquantparams=self.open_main_files()
        print(str(len(self.spelist)),' spectral files loaded from Auger dataset.')
    
    def open_main_files(self):
        ''' Auto loads Auger param files from working directory including 
        EDXparalog- assorted params associated w/ each SEM-EDX or TEM-EDX emsa file 
        Backfitparamslog - ranges and parameters for EDX background fits
        Integquantlog - subtracted and corrected counts for chosen elements
        Peakfitlog - params of gaussian fits to each element (xc, width, peakarea, Y0, rsquared)'''
        if os.path.exists('Augerparamlog.csv'):
            Augerparamlog=pd.read_csv('Augerparamlog.csv', encoding='cp437')
            spelist=Augerparamlog[pd.notnull(Augerparamlog['Areas'])]
        else:
            Augerparamlog=pd.DataFrame()
            spelist=pd.DataFrame()
        if os.path.exists('Smdifpeakslog.csv'):
            Smdifpeakslog=pd.read_csv('Smdifpeakslog.csv', encoding='cp437')
        else:
            Smdifpeakslog=pd.DataFrame()
        if os.path.exists('Backfitlog.csv'):
            Backfitlog=pd.read_csv('Backfitlog.csv', encoding='cp437')
        else:
            Backfitlog=pd.DataFrame()
        if os.path.exists('Integquantlog.csv'):
            Integquantlog=pd.read_csv('Integquantlog.csv', encoding='cp437')
        else:
            Integquantlog=pd.DataFrame()
        # Print TEM or SEM to console based on beam kV
        try:
            AESquantparams=pd.read_csv('C:\\Users\\tkc\\Documents\\Python_Scripts\\Augerquant\\Params\\AESquantparams.csv', encoding='utf-8')
        except:
            AESquantparams=pd.DataFrame()
        return Augerparamlog, spelist, Smdifpeakslog, Backfitlog, Integquantlog, AESquantparams
    
    def savechanges(self):
        ''' TODO need method of saving alterations made to these log files '''
        print('Save changes method of AESdataset called')
        pass
    
class GUIOptions():
    ''' parent is GUImain '''
    def __init__(self,root,parent):
        self.root = root
        self.parent = parent
        self.specselect_frame = tk.Frame(self.root,pady=10)
        self.specselect_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        self.misc_opts_frame = tk.Frame(self.root,pady=10)
        self.misc_opts_frame.pack(side=tk.TOP,fill=tk.X,expand=1)
        # list of checkbuttons one per spectrum in specselect frame
        self.spec_tklist=[]  # list of tk bools (active or inactive)

        self.replot_button = tk.Button(
            self.misc_opts_frame,text="Replot",command=self.filter_files,
            padx=2, pady=6)
        self.replot_button.pack(side=tk.BOTTOM,fill=tk.X,expand=1)

        self.save_button = self._custom_button(
            self.misc_opts_frame,"Save AESdata changes", self.save)
        
        self.quit_button = self._custom_button(
            self.misc_opts_frame,"Quit", self.quitapp)

    def populate_specselector(self, spelist):
        ''' On project load, regenerate list of tk bools from spelist, update specselect frame view '''
        self.spec_tklist=[]
        for i, name in self.spelist:
            self.spec_tklist.append(tk.BooleanVar())
            self.spec_tklist[i].set(0) # Default unselected
            # Fill spectra selector frame w/ associated checkbuttons
            tk.Checkbutton(self.specselect_frame, text=name, variable=self.spec_tklist[i]).pack(side=tk.TOP)

    def save(self):
        ''' TODO set this to call save method of data object '''
        if self.parent.plotter.data_manager is None:return
        dump_name = filedialog.asksaveasfilename()
        if dump_name:
            np.save(dump_name,self.parent.plotter.data_manager)
        
    def filter_files(self):
        ''' Evaluate tkbools for each spe, then call plot method in GUIplotter'''
        if self.parent.plotter.data_manager is None:return
        # Filter and pass selected subset of files to plot method in GUIplotter)
        newlist=[]
        for i, val in enumerate(self.spec_tklist):
            print(val)
            if val==1:
                newlist.append(self.spec_tklist[i].get())

    def quitapp(self):
        msg = "Quitting:\nUnsaved progress will be lost.\nDo you wish to Continue?"
        if tkmess.askokcancel("Combustible Lemon",msg):
            self.parent.root.destroy()
