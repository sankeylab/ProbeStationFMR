# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 18:27:51 2017

@author: admin
"""
import spinmob.egg as egg
import cv2 as _cv2

############################################
### Camera interface
############################################

def image_to_rgb_data(image): 
    """
    Converts the array from cv2 to an array of data for our plotter.
    """
    return image.transpose((1,0,2)).astype(float)


def data_to_image(image, rescale=False):
    """
    Converts our plotty data to an image array for cv2.
    """    
    if rescale: 
        imax = max(image)
        imin = min(image)        
        return ((image.transpose((1,0)) - imin) * 256.0/(imax-imin)).astype(int) 
    else: 
        return image.transpose((1,0)).astype(int)

class ImageWithButtons(egg.gui.GridLayout):
    
    def __init__(self, window):
        """
        This object is a grid layout containing an image with save / load 
        buttons.
        
        You must supply a window object so that it can connect buttons to 
        actions, etc.
        """

        # initialize the grid layout
        egg.gui.GridLayout.__init__(self)    
        
        # no need for more margins
        self._layout.setContentsMargins(0,0,0,0)

        # store the window object
        self._window = window

        # add the save, load, and image
        #self.button_save = self.place_object(egg.gui.Button("Save"))
        #self.button_load = self.place_object(egg.gui.Button("Load"))
        imgobj = egg.pyqtgraph.ImageView()
        imgobj.ui.histogram.hide()
        imgobj.ui.roiBtn.hide()
        imgobj.ui.menuBtn.hide()
        self.image       = self.place_object(imgobj, 0,1, column_span=3, alignment=0)
       
        
        # sha-clacky the buttons together
        self.set_column_stretch(2,10)        
        
        # data
        self.data = 0.0
        

        # connect the buttons to the functionality
        #self._window.connect(self.button_save.signal_clicked, self.button_save_clicked)
        #self._window.connect(self.button_load.signal_clicked, self.button_load_clicked)
        
    def set_data(self, data, **kwargs):
        """
        Sets the image view data to the supplied array.
        """
        self.image.setImage(data, **kwargs)
        self.data = data
        
    def set_levels(self, minvalue, maxvalue, minlevel, maxlevel):
        """
        Sets the minimum and maximum values of the histogram as well as the levelbars. 
        """
        self.image.setLevels(minlevel, maxlevel)
        self.image.ui.histogram.setHistogramRange(minvalue, maxvalue)
    
    def save_image(self, path="ask"):
        """
        Saves the image.
        """
        # get a valid path
        if path=="ask": path = egg.dialogs.save("*.png")
        if not path: return

        # save the image
        _cv2.imwrite(path, data_to_image(self.image.image))
    
    def button_save_clicked(self, *a): self.save_image()

    def load_image(self, path="ask"):
        """
        Loads an image.
        """
        # get a path
        if path=="ask": path = egg.dialogs.open_single("*.png")
        if not path: return
        
        # load the image
        rgb = image_to_rgb_data(_cv2.imread(path))
    
        # assume r+g+b / 3 by default.
        self.set_data((rgb[:,:,0]+rgb[:,:,1]+rgb[:,:,2]) / 3.0)

    def button_load_clicked(self, *a): 
        self.load_image()