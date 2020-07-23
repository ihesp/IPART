'''Compute 3d THR on IVT. Process a single data file.

Input:
    IVT (integrated Vapor Transport) in netCDF format.

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
        python compute_thr_singlefile.py
        ```

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-07-22 10:22:08.
'''

from __future__ import print_function


#--------------Globals------------------------------------------

#-----------IVT data----------------------
IVT_FILE='/home/guangzhi/datasets/erai_qflux/ivt_m1-60_6_2007_cln-cea-proj.nc'
VARIN='ivt'          # data id in nc file

LAT1=0; LAT2=90      # degree, latitude domain

#-------Structuring element for erosion (E)-------
KERNEL=[10,6,6]   # half length of time (time steps), and half length of spatial (number of grids)

SHIFT_LON=80  # shift longitudinally to center Pacific and Altantic

# Orographic file, providing surface terrain elevation info.
# This is optional, can be used to enhance the continent-penetration
# of landfalling ARs.
ORO_FILE='/home/guangzhi/datasets/oro_s_a_1900_erai-cea-proj.nc'
HIGH_TERRAIN=600 # surface height (in m) above which land surface is defined
                 # as high terrain. Extra computations are performed over
                 # high terrain areas to enhance continent-penetration of
                 # landfalling ARs.

#------------------Output folder------------------
OUTPUTDIR='/home/guangzhi/datasets/quicksave2/THR'







#--------Import modules-------------------------
import os
from ipart.utils import funcs
from ipart import thr



#-------------Main---------------------------------
if __name__=='__main__':


    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    #-----------Read in data----------------------
    print('\n### <compute_thr_singlefile>: Read in file:\n',IVT_FILE)
    var=funcs.readNC(IVT_FILE, VARIN)
    var=var.sliceIndex(0, 200, axis=0)

    #--------------------Read in orographic data--------------------
    oroNV=None

    #-----------------Shift longitude-----------------
    var=var.sliceData(LAT1, LAT2, 2)
    var=var.shiftLon(SHIFT_LON)

    #----------------------Do THR----------------------
    ivt, ivtrec, ivtano=thr.THR(var, KERNEL, oroNV=oroNV,
            high_terrain=HIGH_TERRAIN)

    #--------Save------------------------------------
    fname=os.path.split(IVT_FILE)[1]
    file_out_name='%s-THR-kernel-t%d-s%d.nc'\
            %(os.path.splitext(fname)[0], KERNEL[0], KERNEL[1])

    abpath_out=os.path.join(OUTPUTDIR,file_out_name)
    print('\n### <testrotatingfilter>: Saving output to:\n',abpath_out)
    funcs.saveNC(abpath_out, ivt, 'w')
    funcs.saveNC(abpath_out, ivtrec, 'a')
    funcs.saveNC(abpath_out, ivtano, 'a')


