# forestlas
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![LiDAR derived vertical profiles](http://www2.geog.ucl.ac.uk/~ucfaptv/3_PLOTS_geo_trees_bar_spectra_v5.png)
Python code for generating metrics of forest vertical structure from airborne LiDAR data. This code was developed as 
part of my PhD (completed in 2016, can be viewed 
<a href=https://www.researchgate.net/publication/290436021_Assessment_of_forest_canopy_vertical_structure_with_multi-scale_remote_sensing_from_the_plot_to_the_large_area>here</a>) 
and was developed over the forests of Victoria, Australia.
The aim was to develop a suite of metrics that are robust to forest type i.e. can be applied without prior information of 
forest structure.

There are a number of methods available, check this <a href=https://github.com/philwilkes/forestlas/blob/master/forestlas_intro.ipynb>
Jupyter notebook</a> for an introduction.
Functions include reading `.las` files to numpy array, writing to `.las` as well as a number of methods to dice, slice and tile 
LiDAR data.
The main set of functions found in `forestlas.canopyComplexity`.
These allow you to derive metrics of vertical canopy structure such as <i>Pgap</i> and also estimate number of canopy layers.
More information can be found in this paper <a href=https://doi.org/10.1111/2041-210X.12510>Wilkes, P. et al. (2016). Using discrete-return airborne laser scanning to 
quantify number of canopy strata across diverse forest types. Methods in Ecology and Evolution, 7(6), 700–712</a>. 


#### Funding
This research was funded by the Australian Postgraduate Award, Cooperative Research Centre for Spatial Information 
under Project 2.07, TERN/AusCover and Commonwealth Scientific and IndustrialResearch Organisation (CSIRO) Postgraduate 
Scholarship.
