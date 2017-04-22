# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:11:20 2017

@author: admin
"""
import numpy    as _n
import spinmob as _spinmob
_d = _spinmob.data
_g = _spinmob.egg._gui
import PyQt4.QtGui as _qt

# import pyqtgraph and create the App.
import pyqtgraph.widgets.MatplotlibWidget as _mplw
#_a = _g.mkQApp()

class ColorWFPlot(_g.GridLayout):
    """
    Creates tabbed color and waterfall plots for provided 3D data. Creates a 2D 
    databox with user-specified dimensions while storing the two independent 
    variables in the header.
    """
    def __init__(self, size, bnds = [0,1,0,1], file_type="*.dat", **kwargs):

        # Initialize widget and databox
        _g.GridLayout.__init__(self, margins=False)
        #_d.databox.__init__(self, **kwargs)

        # Initialize data array and headers
        self.dbs = []
        self.set_size_and_bounds(size, bnds)
        
        # Top row is main controls
        self.place_object(_g.Label("Raw Data:"), alignment=1)
        self.button_load     = self.place_object(_g.Button("Load")                .set_width(50), alignment=1)
        self.button_save     = self.place_object(_g.Button("Save")                .set_width(50), alignment=1)
        self.place_object(_g.Label("Plotted Databox:"), alignment=1)
        self.dbs_combobox = self.place_object(_g.ComboBox())

        # keep the buttons shaclackied together
        self.set_column_stretch(5)
        
        # Second row adds the plots
        self.tabs = self.place_object(_g.TabArea(False), 0,1, row_span=4, column_span=5, alignment=0)
        
        # Configure the color plot
        self.colortab = self.tabs.add_tab("Color Plot")
        self.colorfig = _mplw.MatplotlibWidget()
        self.colorax  = self.colorfig.getFigure().add_subplot(111)
        
        # Initialize color plot
        imgplotpad = _n.zeros((size[0], size[1]))
        self.colorplot = self.colorax.imshow(imgplotpad, aspect='auto', extent = self.data_bounds, \
            interpolation = 'none', cmap = 'bwr')
        self.colorplot_cb = self.colorfig.getFigure().colorbar(self.colorplot)

        # Draw.        
        self.colorfig.getFigure().tight_layout()        
        self.colortab.place_object(self.colorfig)
        
        # Configure the waterfall plot
        self.wftab = self.tabs.add_tab("Waterfall Plot")
        self.wffig  = _mplw.MatplotlibWidget()
        self.wfax = self.wffig.getFigure().add_subplot(111)
        self.wfax.autoscale_view(True,True,True)
        self.wfplots = [] # Object containing plots corresponding to each data point        
        
        # Draw.
        self.wffig.getFigure().tight_layout()        
        self.wftab.place_object(self.wffig, column_span = 4)
        
        # Spinner for waterfall plots spacing
        self.wftab.new_autorow()
        self.wfautoscale = self.wftab.place_object(_g.Button('Autoscale').set_checkable(True).set_checked(True))
        self.wftab.place_object(_g.Label("Relative Spacing:").set_style(text_align='right'))
        self.wfspinner = self.wftab.place_object(_g.NumberBox(value = 1, step = 0.1))
        self.wfspacing = self.wftab.place_object(_g.CheckBox().set_checked(True))

        # file type (e.g. *.dat) for the file dialogs
        self._file_type = file_type

        # Button functionality
        self.button_save      .signal_clicked.connect(self._button_save_clicked)
        self.button_load      .signal_clicked.connect(self._button_load_clicked)
        self.wfspinner        .signal_changed.connect(self.update_wfplot)
        self.wfspacing        .signal_changed.connect(self.update_wfplot)
        self.wfautoscale      .signal_toggled.connect(self.update_wfplot)
        self.dbs_combobox     .signal_changed.connect(self.update_plots)
        
    def _button_save_clicked(self, *a):
        """
        Called whenever the button is clicked.
        """
        self.save_file()

    def _button_load_clicked(self, *a):
        """
        Called whenever the button is clicked.
        """
        self.load_file()        

    def save_file(self, path="ask", force_overwrite=False, just_settings=False):
        """
        Saves the data in the different databoxes to files.

        just_settings=True means only save the configuration of the controls
        """

        # if it's just the settings file, make a new databox
        if just_settings: d = _d.databox()

        # otherwise use the internal databoxes
        else: d = self.dbs

        # add all the controls settings
        for x in self._autosettings_controls: self._store_gui_setting(d, x)

        # save the file
        if path == "ask": 
            path = str(_qt.QFileDialog.getSaveFileName(None,'Select where to save the databoxes'))
            
        for k in range(len(self.dbs)):
            print (path, self.dbs_combobox.get_text(k))
            _d.databox.save_file(d[k], path + "." + self.dbs_combobox.get_text(k), self._file_type, force_overwrite)

    def load_file(self, path="ask", just_settings=False):
        """
        Loads a data file. After the file is loaded, calls self.after_load_file(self),
        which you can overwrite if you like!

        just_settings=True will only load the configuration of the controls,
        and will not plot anything or run after_load_file
        """
        # if it's just the settings file, make a new databox
        #if just_settings:
        #    d = _d.databox()
        #    header_only = True

        # otherwise use the internal databox
        #else:
        #    d = self.dbs[self.dbs_combobox.get_current_index()]
        #    header_only = False

        # Currently, loading a databox means that you flush out all the previous ones.
        self.dbs = []
        cbox_len = len(self.dbs_combobox.get_all_items())  
 
        # Prompt the user to fetch the databoxes       
        paths = _qt.QFileDialog.getOpenFileNames(None, 'Select databoxes to import')
        
        # Figure out the data size and the bounds from a single dummy databox   
        d = _d.databox()
        _d.databox.load_file(d, paths[0], filters=self._file_type, quiet=just_settings)
            
        self.set_size_and_bounds([d.headers['xsize'], d.headers['ysize']], \
                                 [d.headers['xbounds'][0], d.headers['xbounds'][1], \
                                  d.headers['ybounds'][0], d.headers['ybounds'][1]])        
        
        # Load the data
        for k in range(len(paths)):
            self.dbs.append(_d.databox())        
            self.dbs_combobox.add_item(paths[k].split('.')[-1])
            _d.databox.load_file(self.dbs[k], paths[k], filters=self._file_type, quiet=just_settings)
            
        # Wipe out the previous combobox entries    
        for k in range(cbox_len): self.dbs_combobox.remove_item(0)
            
        # Update the plots
        self.update_plots()


        # import the settings if they exist in the header
       # if not None == _d.databox.load_file(d, path, filters=self._file_type, header_only=header_only, quiet=just_settings):
        #    # loop over the autosettings and update the gui
         #   for x in self._autosettings_controls: self._load_gui_setting(d, x)


    def after_load_file(self,*args):
        """
        Called after a file is loaded. Does nothing. Feel free to overwrite!

        The first argument is the DataboxPlotInstance so your function can
        tell which instance loaded a file.
        """
        return


    def add_databox(self, name = None):
        """
        Adds a new databox containing another variable of data.
        """
        self.dbs.append(_d.databox())
        self.dbs_combobox.add_item(str(name))
        #self.dbs[-1].insert_header('xvals', self.xvals)
        #self.dbs[-1].insert_header('yvals', self.yvals)
        
        self.update_headers()
       
    def set_databox(self, index = 1):
        """
        Sets the currently plotted databox.
        """
        self.dbs_combobox.set_current_index(index)
        self.update_plots()
        
    def clear_databoxes(self):
        """
        Clears all databoxes.
        """
        for k in self.dbs: 
            k.clear()
       
    def set_size_and_bounds(self, size, bnds = [0,1,0,1]):
        """
        Sets the size and bounds of the final data array.
        Format: size = [xsteps, ysteps], bnds = [xmin, xmax, ymin, ymax]
        """
        self.clear_databoxes()        
        
        self.array_size  = size
        self.data_bounds = bnds
        
        # Build databox
        #self.clear()
        #for n in range(size[0]): self["f"+str(n)] = []
            
        # Initialize arrays for independent variables and store in header
        self.xvals = _n.linspace(bnds[0], bnds[1], size[0])
        self.yvals = _n.linspace(bnds[2], bnds[3], size[1])
#        for k in range(len(self.dbs)):        
#            self.dbs[k].insert_header("xvals", self.xvals)
#            self.dbs[k].insert_header("yvals", self.yvals)
            
        # Update headers
        self.update_headers()

            
    def insert_header(self, hkey, value):
        """
        Inserts a header in each databox.
        """
        for k in self.dbs:
            k.insert_header(hkey, value)
                
        
    def update_headers(self):
        """
        Updates the size and bounds headers in all databoxes.
        """
        for k in self.dbs:
            k.update_headers({'xbounds': self.data_bounds[0:2], 'ybounds': self.data_bounds[2:], \
                              'xsize': self.array_size[0], 'ysize': self.array_size[1]})
    
    
    def set_labels(self, labels):
        """
        Sets the labels on the plots. Four-string list of the format
        [xlabel, colorax_ylabel, wfax_ylabel, title = '']
        """
        self.colorax.set_xlabel(labels[0])
        self.colorax.set_ylabel(labels[1])
        if len(labels) == 4: self.colorax.set_title(labels[3])        
        
        self.wfax.set_xlabel(labels[0])
        self.wfax.set_ylabel(labels[2])
        if len(labels) == 4: self.wfax.set_title(labels[3])


    def append_data_row(self, db, newdata):
        """
        Appends a new row of data to the databox.
        """
        self.dbs[db].append_data_point(newdata)
        self.update_plots()  
    
    
    def update_data_row(self, db, newdata, index):
        """
        Updates the nth data point with new data.
        """
        self.dbs[db].pop_data_point(index)
        self.dbs[db].insert_data_point(newdata,index)

        
    def fetch_data(self, db = None):
        """
        Fetches all the data in a certain databox.
        """
        if db == None: db = self.dbs_combobox.get_current_index()
            
        # Fetch the data
        try:
            data = _n.array(self.dbs[db][0::]).transpose()
            datarows = data.size/self.array_size[0]
        except:
            data = []
            datarows = 0

        return data, datarows
        
 
    def update_colorplot(self, data = None):
        """
        Updates the color plot
        """
        # fetch data in current databox
        data, datarows = self.fetch_data()
        
        # pad with zeros for the image plot and update colorbar
        imgplotpad = _n.zeros( (self.array_size[1] - datarows, self.array_size[0]) )
        if data != []:
#            self.colorplot.remove()
            self.colorplot = self.colorax.imshow(_n.flipud(_n.vstack((data,imgplotpad))), aspect='auto', extent = self.data_bounds, \
                interpolation = 'none', cmap = 'bwr')
            self.colorplot_cb.set_clim(vmin = data[_n.isfinite(data)].min(), vmax = data[_n.isfinite(data)].max())
            self.colorplot_cb.draw_all()        

        # Redraw
        self.colorfig.getFigure().tight_layout() 
        self.colorfig.canvas.draw()  
  
  
    def update_wfplot(self):
        """
        Updates the waterfall plot
        """        
        # fetch data in current databox
        data, datarows = self.fetch_data()
        data = _n.array(data)

        # Remove previous plots
        for k in range(len(self.wfplots)):
            self.wfplots[k].remove()
        self.wfplots = []
        self.wfax.set_color_cycle(None) # So the plots don't change color on every update

        # Figure out if we offset the plots
        hasbuffer = self.wfspacing.is_checked()/2 # True = 2 for some reason...

        # Redraw waterfall plot
        wf_buffer = 0     # Adjusts the relative spacing between the plots
        for k in range(0, datarows):
            # Fetch all (valid) values of the dataset for plotting
            xvals = self.xvals[_n.isfinite(data[k,:])]                      
            yvals = data[k,_n.isfinite(data[k,:])]
            try:
                ymin, ymax = _n.min(yvals), _n.max(yvals)
            except ValueError:
                ymin, ymax = 0, 0
            
            self.wfplots.append(self.wfax.plot(xvals, yvals + hasbuffer*(-ymin + wf_buffer) )[0])
            wf_buffer += (ymax - ymin)*float(self.wfspinner.get_value())
                
        # Redraw
        if self.wfautoscale.is_checked() and datarows != 0: # Autoscale if checked
            plot_ydata = []            
            for k in range(0, datarows):
                plot_ydata.append(self.wfplots[k].get_data()[1])
            self.wfax.set_ylim(_n.min(plot_ydata), _n.max(plot_ydata))
        self.wffig.getFigure().tight_layout() 
        self.wffig.canvas.draw()


    def update_plots(self):
        """
        Updates the plots.
        """
        self.update_colorplot()
        self.update_wfplot()