# alkali_flies_of_mono_lake
Data and software associated with the paper "Super-hydrophobic diving flies (Ephydra hians) and the hypersaline waters of Mono Lake"

## What you need to run our analysis
* Ubuntu (we used Ubuntu 12-16)
* Python (2.7)
* ROS (Robot Operating System, Kinetic): http://wiki.ros.org/kinetic/Installation/Ubuntu
* apt-get repositories: git python-pip python-scipy python-h5py python-progressbar python-sympy python-networkx
* Manual downloads: 
..* https://github.com/taynaud/python-louvain/
..* http://www.pyqtgraph.org/
* pip installs: pandas (0.19), statsmodels
* My packages:
..* FigureFirst: https://github.com/FlyRanch/figurefirst
..* FlyPlotLib: https://github.com/florisvb/FlyPlotLib
..* DataFit: https://github.com/florisvb/DataFit
..* FlyStat: https://github.com/florisvb/FlyStat
* Inkscape

You may wish to do all of this in a virtual environment.

## Downloading the data

## Downloading the code
All of the code is contained in this repository, but we depend on the repos listed above in "what you need to run our analysis". Git clone this repo. 

## Making the data automatically accessible to the analysis
In mono_paper_locations/mono_paper_locations, create a duplicate of data_locations_analysiscavetech_organized.py. Edit the file so that the paths correspond to the data 

## Processing raw data
Raw data is saved as .bag files, which contain raw force measurements, lvdt, and movie images, all time-synced. See http://wiki.ros.org/ROS/Tutorials/Recording%20and%20playing%20back%20data

Our first step was to manually segment the data into the portions where the fly is entering the water, stable and submerged, or exiting the water. We did this using a pyqtgraph gui. You can view our segments for data (and change them!) as follows:

Directory: raw_analysis
Command: [code]python align_force_data.py --file=FILENAME.bag[/code]

(to change the selection, drag the vertical lines, and click save)

These time segments are saved (along with mean and std values) in a pickle file associated with each filename.

Calibration data was preprocessed in a similar way; we selected stable segments of before, during, and after calibration weight placement and used the mean values of these segments to determine the calibration for each fly, which is also saved as a separate pickle file. 

Directory: raw_analysis
Command: [code]python extract_calibration_data.py --file=FILENAME.bag[/code]

## 
