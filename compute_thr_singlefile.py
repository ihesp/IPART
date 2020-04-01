'''Compute 3d THR on IVT. Process a single data file.

Input:
    IVT (integrated Vapor Transport) in netCDF format.

    Data are assumed to be in the format: (time, level, latitude, longitude)
    dimensions (level dimension is optional, if present, should be a
    singleton axis of length 1).

    NOTE: data should have proper time, latitude and longitude axes.

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

Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
Update time: 2020-04-01 12:14:48.
'''

from __future__ import print_function


#--------------Globals------------------------------------------

#-----------IVT data----------------------
IVT_FILE='/home/guangzhi/datasets/erai/ERAI_AR_THR/ivt_m1-60_6_1984_crop.nc'
VARIN='ivt'          # data id in nc file

LAT1=0; LAT2=90      # degree, latitude domain

#-------Structuring element for erosion (E)-------
KERNEL=[16,6,6]   # half length of time (time steps), and half length of spatial (number of grids)

SHIFT_LON=80  # shift longitudinally to center Pacific and Altantic

#------------------Output folder------------------
OUTPUTDIR='/home/guangzhi/datasets/erai/ERAI_AR_THR/'







#--------Import modules-------------------------
import os
import cdms2 as cdms
import MV2 as MV
import numpy as np
from skimage import morphology
from utils import funcs




def filterData(ivt,kernel,verbose=True):
    """Perform THR filtering process on 3d data

    Args:
        ivt (TransientVariable): 3D input IVT data, with dimensions (time, lat, lon) or
                                 (time, level, lat, lon).
        kernel (list or tuple): list/tuple of integers specifying the shape of the kernel/structuring
                                element used in the gray erosion process.
    Returns:
        ivt (TransientVariable): 3D array, input <ivt> squeezed.
        ivtrec (TransientVariable): 3D array, the reconstruction component from the THR process.
        ivtano (TransientVariable): 3D array, the difference between input <ivt> and <ivtrec>.
    """

    ndim=np.ndim(ivt)
    ivt=ivt(squeeze=1)

    #-------------------3d ellipsoid-------------------
    ele=funcs.get3DEllipse(*kernel)
    #dt=kernel[0] # half length in time dimesion

    #################### use a cube to speed up ##############
    if kernel[0]>=10 or kernel[1]>=6:
        ele=np.ones(ele.shape)
    ##########################################################

    if verbose:
        print('\n# <filterData>: Computing erosion ...')

    lm=morphology.erosion(ivt.data,selem=ele)

    if verbose:
        print('\n# <filterData>: Computing reconstruction ...')

    ivtrec=morphology.reconstruction(lm,ivt,method='dilation')
    ivtrec=MV.array(ivtrec)
    ivtano=ivt-ivtrec
    ivtano=MV.array(ivtano)
    ivtrec=MV.array(ivtrec)

    if ndim==4:
        levax=cdms.createAxis([0,])
        levax.designateLevel()
        levax.id='z'
        levax.name='level'
        levax.units=''

        ivt=funcs.addExtraAxis(ivt,levax,1)
        ivtano=funcs.addExtraAxis(ivtano,levax,1)
        ivtrec=funcs.addExtraAxis(ivtrec,levax,1)

    axislist=ivt.getAxisList()
    ivtano.setAxisList(axislist)
    ivtrec.setAxisList(axislist)

    ivtrec.id='ivt_rec'
    ivtrec.long_name='Integrated moisture transport, minimal reconstruction'
    ivtrec.standard_name=ivtrec.long_name
    ivtrec.title=ivtrec.long_name
    ivtrec.units=ivt.units

    ivtano.id='ivt_ano'
    ivtano.long_name='Integreated moisture transport, anomaly wrt minimal reconstruction'
    ivtano.standard_name=ivtano.long_name
    ivtano.title=ivtano.long_name
    ivtano.units=ivt.units

    return ivt, ivtrec, ivtano





#-------------Main---------------------------------
if __name__=='__main__':


    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    #-----------Read in data----------------------
    var=funcs.readVar(IVT_FILE, 'ivt')

    #-----------------Shift longitude-----------------
    var=var(latitude=(LAT1, LAT2))
    var=var(longitude=(SHIFT_LON,SHIFT_LON+360))

    #----------------------Do THR----------------------
    ivt, ivtrec, ivtano=filterData(var, KERNEL)

    #--------Save------------------------------------
    fname=os.path.split(IVT_FILE)[1]
    file_out_name='%s-minimal-rec-ano-kernel-t%d-s%d.nc'\
            %(os.path.splitext(fname)[0], KERNEL[0], KERNEL[1])

    abpath_out=os.path.join(OUTPUTDIR,file_out_name)
    print('\n### <testrotatingfilter>: Saving output to:\n',abpath_out)
    fout=cdms.open(abpath_out,'w')
    fout.write(ivt,typecode='f')
    fout.write(ivtrec,typecode='f')
    fout.write(ivtano,typecode='f')
    fout.close()


