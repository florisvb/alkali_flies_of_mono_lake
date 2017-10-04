# Alkali Flies of Mono Lake
Data and software associated with the paper "Super-hydrophobic diving flies (Ephydra hians) and the hypersaline waters of Mono Lake"

This readme assumes working knowledge of Ubuntu and python.

## What you need to run our analysis
* Ubuntu (we used Ubuntu 12-16)
* Python (2.7)
* ROS (Robot Operating System, Kinetic): http://wiki.ros.org/kinetic/Installation/Ubuntu
* apt-get repositories: git python-pip python-scipy python-h5py python-progressbar python-sympy python-networkx
* Manual downloads: 
  * https://github.com/taynaud/python-louvain/
  * http://www.pyqtgraph.org/
* pip installs: pandas (0.19), statsmodels
* My packages:
  * FigureFirst: https://github.com/FlyRanch/figurefirst
  * FlyPlotLib: https://github.com/florisvb/FlyPlotLib
  * DataFit: https://github.com/florisvb/DataFit
  * FlyStat: https://github.com/florisvb/FlyStat
* Inkscape

You may wish to do all of this in a virtual environment.

## Downloading the data

## Downloading the code
All of the code is contained in this repository, but we depend on the repos listed above in "what you need to run our analysis". Git clone this repo. 

## Making the data automatically accessible to the analysis
We ran our analysis on several different computers, so to keep track of everything, we created a python package that points to the data and figure template locations. In order to run our analysis, you will need to add your machine and local paths to this repository. 

In `mono_paper_locations/mono_paper_locations`, create a duplicate of `data_locations_analysiscavetech_organized.py`, e.g. `data_locations_yourname.py`. Edit the file so that the paths correspond to the data locations on your machine. Next, you will need to create an environmental variable called `mono_paper_locations` (e.g. type `export mono_paper_locations mono_paper_yourname` in any terminal window in which you plan to run our code, or add that to your .bashrc). Add an `elif` statement to the `__init__.py` file in `mono_paper_locations` that matches, for example, `mono_paper_yourname` to  `data_locations_yourname.py` as in the other if and elif statements.

In `mono_paper_locations/mono_paper_locations`, edit the file `figure_template_locations.py`, so that the paths match your system.

Install the package (from `mono_paper_locations` type `python setup.py install`). 

## Processing raw data
Raw data is saved as .bag files, which contain raw force measurements, lvdt, and movie images, all time-synced. See http://wiki.ros.org/ROS/Tutorials/Recording%20and%20playing%20back%20data

Our first step was to manually segment the data into the portions where the fly is entering the water, stable and submerged, or exiting the water. We did this using a pyqtgraph gui. You can view our segments for data (and change them!) as follows:

* Directory: raw_analysis
* Command: `python align_force_data.py --file=FILENAME.bag`

To change the selection, drag the vertical lines, and click save. These time segments are saved (along with mean and std values) in a pickle file associated with each filename.

Calibration data was preprocessed in a similar way; we selected stable segments of before, during, and after calibration weight placement and used the mean values of these segments to determine the calibration for each fly, which is also saved as a separate pickle file. 

* Directory: raw_analysis
* Command: `python extract_calibration_data.py --file=FILENAME.bag`

## Installing the analysis

Install the following analysis packages within this repository (`python setup.py install`)
* gcms_analysis
* mono_analysis

## Running the analysis

In each "figure" folder (except figure1) there is a make_figureX.py file. Run this file (`python ./make_figureX.py`) to rerun the analysis and update the associated svg figure files in that directory. You can use this to trace backwards our analysis, most of which can be found in mono_analysis/plot_raw_traces.py

## Species and colloquial names

In our analysis and data we refer to the species by the following alternative names:

Fucelia rufitibia: blue kelp fly
Coelopa vanduzeei: black kelp fly
Ephydra hians: alkali fly / mono lake fly
Ephydra sp: santa ana fly
Helaeomyia petrolei: oil fly / petroleum fly
Drosophila melanogaster: melanogaster
Drosophila virilis: virilis
