# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 17:49:52 2017

@author: tkc
"""

class AESshift_app(QtWidgets.QMainWindow):
    ''' Main app/GUI for SEM-EDX background fit corrections 
    old version with EDXcanvas split from EDX_app '''
    def __init__(self, parent=None):
        super(AESshift_app, self).__init__(parent)

        self.setGeometry(50,50,1050,650) # main window sizing
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Auger shift determination")

        self.main_widget = QtWidgets.QWidget(self) # define main matplotlib widget                
        self.vbl = QtWidgets.QVBoxLayout(self.main_widget) # simple layout w widgets stacked

        self.main_widget.setFocus() # keyboard input to this widgets
        self.setCentralWidget(self.main_widget)

        self.EDXdata=EDXdataset() # instance w/ 6 pandas dataframes with params

        # sets EDX canvas (mpl plot and spinbox) as main widget
        sc = EDXCanvas(self, self.main_widget, self.EDXdata)
        self.vbl.addWidget(sc)
        self.setLayout(self.vbl) # finalizes layout 
        self.show()
        
    def file_Quit(self):
        ''' Done button within child (EDXCanvas) calls quit for main window 
        '''
        self.close()

