'''Compute 3d THR on IVT. Process multiple data files too large to fit into
RAM.

Difference from compute_thr_singlefile.py:

    this is for computing THR process on multiple data files that are too
    large to fit into RAM. Like a moving average, the THR process is
    less accurate at the ends of the data domain. For data saved into multiple
    files (usually separated by time), this may create some discontinuity
    at the ends of each file. To overcome this, 2 files are read in at a time,
    data are concatenated in time dimension, therefore the transition
    between these 2 data files is "smooth". Then the 3rd file is read in to
    form a smooth transition between the 2nd and the 3rd files. Then the same
    process rotates on. The outputs are saved one for a year.

Input:
    IVT (integrated Vapor Transport) in netCDF format, one file for a year
    (or a month, depending on your specific organization of data files).
    The file name is assumed to have a format that there is a field that
    specifies the year, e.g. "ivt_m1-60_6_%d_cln.nc". The %d field is replaced
    by a year, e.g. 2000.

    Data are assumed to be in the format: (time, level, latitude, longitude)
    dimensions (level dimension is optional, if present, should be a
    singleton axis of length 1).

    NOTE: data should have proper time, latitude and longitude axes.

Optional input:

    Orographic data providing the surface terrain elevations, that correspond
    to the IVT data. This is used to perform some extra computations over
    high terrain regions to enhance the inland penetration of ARs. The mostly
    affected area is the western coast of North America. Other areas are mostly
    not affected.

Usage:

    Change global parameters in the Globals section to point to the storage
    location of IVT data, and specificy an output folder to save results.

    Specify the latitudinal domain in LAT1, LAT2.

    The KERNEL parameter specifies the $t$ and $s$ parameters of the
    structuring element size.
    $t$ is in number time steps, $s$ is number of grid cells.
    See paper for more details, but basically the choices of $t$ and $s$
    should correspond to the synoptic temporal and spatial scales.

    SHIFT_LON shifts the longitude by a given degree of longitudes, so
    that the Pacific and Atlantic basins can be centered.

    Run the script as:
        ```
        python compute_thr_multifile.py
        ```

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-07-22 10:23:35.
'''

from __future__ import print_function


#--------------Globals------------------------------------------
#--------------------Year range--------------------
YEARS=range(2003,2006)

#-----------IVT data----------------------
SOURCEDIR1='/home/guangzhi/datasets/erai/erai_qflux'
FILE1_BASE='ivt_m1-60_6_%d_cln.nc'
VARIN='ivt'          # data id in nc file

LAT1=20; LAT2=60      # degree, latitude domain

#-------Structuring element for erosion (E)-------
KERNEL=[8,6,6]   # half length of time (time steps), and half length of spatial (number of grids)

SHIFT_LON=80  # shift longitudinally to center Pacific and Altantic

# Orographic file, providing surface terrain elevation info.
# This is optional, can be used to enhance the continent-penetration
# of landfalling ARs.
ORO_FILE='/home/guangzhi/datasets/oro_s_a_1900_erai.nc'
HIGH_TERRAIN=600 # surface height (in m) above which land surface is defined
                 # as high terrain. Extra computations are performed over
                 # high terrain areas to enhance continent-penetration of
                 # landfalling ARs.

#------------------Output folder------------------
OUTPUTDIR='/home/guangzhi/datasets/quicksave2/THR/'







#--------Import modules-------------------------
import os
from ipart.utils import funcs
from ipart import thr




#-------------Main---------------------------------
if __name__=='__main__':


    filelist=[]
    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    #--------------------Read in orographic data--------------------
    oroNV=funcs.readNC(ORO_FILE, 'z')
    oroNV.data=oroNV.data/9.8
    oroNV=oroNV.sliceData(LAT1, LAT2, axis=1)
    oroNV=oroNV.shiftLon(SHIFT_LON).squeeze()

    for year in YEARS:
        #-----------Read in data----------------------
        file_in_name=FILE1_BASE %(year)
        abpath_in=os.path.join(SOURCEDIR1,file_in_name)
        filelist.append(abpath_in)

    if len(filelist)<2:
        raise Exception("Need to give at least 2 files. For single file, use compute_thr_singlefile.py")

    selector=funcs.Selector(LAT1, LAT2, axis=2)
    thr.rotatingTHR(filelist, VARIN, KERNEL, OUTPUTDIR,
            oroNV=oroNV, selector=selector, high_terrain=HIGH_TERRAIN)

