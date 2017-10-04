from optparse import OptionParser
import sys, os
import bag2hdf5, h5py

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np

import pickle
    
class CalibrationWidget(object):
    def __init__(self, regions, t, force, calibration_filename):
        self.regions = regions
        self.t = t
        self.force = force
        self.calibration_filename = calibration_filename
        
        ## create GUI
        self.app = QtGui.QApplication([])
        self.w = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.w.setLayout(self.layout)
        
        btn = QtGui.QPushButton('save calibration')
        btn.pressed.connect(self.save_calibration)
        self.layout.addWidget(btn, 0, 0)
        
        self.p1 = pg.PlotWidget(title="Basic array plotting", x=t, y=force)
        self.p1.enableAutoRange('xy', False)
        self.layout.addWidget(self.p1, 0, 1)
        
        for r, region in self.regions.items():
            lr = pg.LinearRegionItem(values=region['values'])
            f = 'update_linear_region_' + str(r)
            lr.sigRegionChanged.connect(self.__getattribute__(f))
            self.p1.addItem(lr)
            
    def update_linear_region(self, r, linear_region):
        self.regions[r]['values'] = linear_region.getRegion()
        indices = [0, 0]
        indices[0] = np.argmin( np.abs( self.t - self.regions[r]['values'][0] ) )
        indices[1] = np.argmin( np.abs( self.t - self.regions[r]['values'][1] ) )
        self.regions[r]['indices'] = indices
        self.regions[r]['mean'] = np.mean(self.force[self.regions[r]['indices'][0]: self.regions[r]['indices'][-1]])
        self.regions[r]['std'] = np.std(self.force[self.regions[r]['indices'][0]: self.regions[r]['indices'][-1]])
        print 'region: ', r, ' :', self.regions[r]['mean']
        
    def update_linear_region_1(self, linear_region):
        self.update_linear_region('1', linear_region)
    
    def update_linear_region_2(self, linear_region):
        self.update_linear_region('2', linear_region)
        
    def update_linear_region_3(self, linear_region):
        self.update_linear_region('3', linear_region)


    def run(self):

        ## Display the widget as a new window
        self.w.show()

        ## Start the Qt event loop
        self.app.exec_()
        
    def save_calibration(self):
        f = open(self.calibration_filename, 'w')
        pickle.dump(self.regions, f)
        f.close()
        print 'calibration saved'
        

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
        
    ## Read data #############################################################
    parser = OptionParser()
    parser.add_option('--filename', type=str, help="the .bag file")
    parser.add_option('--max_strlen', type=int, default=255,
                        help="maximum length of encoded strings")
    parser.add_option('--out', type=str, default=None,
                        help="name of output file")
    parser.add_option('--topic', type=str, nargs='*',
                        help="topic name to convert. defaults to all. "
                        "multiple may be specified.")
    (options, args) = parser.parse_args()
           
    fname = os.path.splitext(options.filename)[0]
            
    if options.out is not None:
        output_fname = options.out
    else:
        output_fname = fname + '.hdf5'
            
    print 'Output name: ', output_fname
            
    if not os.path.exists(output_fname):
        if not os.path.exists(options.filename):
            print >> sys.stderr, 'No file %s' % options.filename
            sys.exit(1)
        bag2hdf5.bag2hdf5(options.filename,
                 output_fname,
                 max_strlen=options.max_strlen,
                 topics='/phidgets_daq/force')            
    
    f = h5py.File(output_fname, 'r')  
            
    t = f['phidgets_daq']['force']['t']
    force = f['phidgets_daq']['force']['value']
    
    print 'loaded time and force'


    calibration_filename = fname + '_calibration.pickle'
    if os.path.exists(calibration_filename):
        cal = open(calibration_filename)
        regions = pickle.load(cal)
        cal.close()
    else:
        # guess time ranges:
        region_length = 0.5
        
        region1_start = t[0]
        region1_end = region1_start + region_length
        
        region2_start = t[int(len(t)/2.)]
        region2_end = region2_start + region_length
        
        region3_start = t[-1]
        region3_end = region3_start - region_length

        #############################################################################
        
        regions = {'1':  {  'values': [region1_start, region1_end],
                            'indices': [0, 0],
                            'mean': 0,
                            'std': 0},
                   '2':  {  'values': [region2_start, region2_end],
                            'indices': [0, 0],
                            'mean': 0,
                            'std': 0},
                   '3':  {  'values': [region3_start, region3_end],
                            'indices': [0, 0],
                            'mean': 0,
                            'std': 0}}
                            
    #############################################################################
    
    calibrationwidget = CalibrationWidget(regions, t, force, calibration_filename)
    calibrationwidget.run()
