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

Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
Update time: 2020-04-01 12:14:48.
'''

from __future__ import print_function


#--------------Globals------------------------------------------

#-----------IVT data----------------------
IVT_FILE='/home/guangzhi/datasets/erai_qflux/ivt_m1-60_6_2007_cln-cea-proj.nc'
VARIN='ivt'          # data id in nc file

LAT1=0; LAT2=90      # degree, latitude domain

#-------Structuring element for erosion (E)-------
KERNEL=[16,6,6]   # half length of time (time steps), and half length of spatial (number of grids)

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
OUTPUTDIR='/home/guangzhi/datasets/erai/ERAI_AR_THR/'







#--------Import modules-------------------------
import os
import cdms2 as cdms
import MV2 as MV
import numpy as np
from skimage import morphology
from utils import funcs



def filterData(ivt, kernel, oro=None, high_terrain=600, verbose=True):
    """Perform THR filtering process on 3d data

    Args:
        ivt (TransientVariable): 3D input IVT data, with dimensions
            (time, lat, lon) or (time, level, lat, lon).
        kernel (list or tuple): list/tuple of integers specifying the shape of
            the kernel/structuring element used in the gray erosion process.
    Keyword Args:
        oro (TransientVariable): 2D array, surface orographic data in meters.
            This additional surface height info is used to perform a separate
            reconstruction computation for areas with high elevations, and
            the results can be used to enhance the continent-penetration
            ability of landfalling ARs. Sensitivity in landfalling ARs is
            enhanced, other areas are not affected. Needs to have compatible
            shape as <ivt>.
        high_terrain (float): minimum orographic height to define as high
            terrain area, within which a separate reconstruction is performed.
            Only used if <oro> is not None.
    Returns:
        ivt (TransientVariable): 3D array, input <ivt> squeezed.
        ivtrec (TransientVariable): 3D array, the reconstruction component from the THR process.
        ivtano (TransientVariable): 3D array, the difference between input <ivt> and <ivtrec>.
    """

    ndim=np.ndim(ivt)
    ivt=ivt(squeeze=1)

    #-------------------3d ellipsoid-------------------
    ele=funcs.get3DEllipse(*kernel)

    #################### use a cube to speed up ##############
    # empirical
    if kernel[0]>=16 or kernel[1]>=6:
        ele=np.ones(ele.shape)
    ##########################################################

    # reconstruction element: a 6-connectivity element
    rec_ele=np.zeros([3,3,3])
    rec_ele[0,1,1]=1
    rec_ele[1,:,1]=1
    rec_ele[1,1,:]=1
    rec_ele[2,1,1]=1

    if verbose:
        print('\n# <filterData>: Computing erosion ...')

    lm=morphology.erosion(ivt.data,selem=ele)

    if verbose:
        print('\n# <filterData>: Computing reconstruction ...')

    ivtrec=morphology.reconstruction(lm, ivt, method='dilation', selem=rec_ele)

    # perform an extra reconstruction over land
    if oro is not None:
        oro_rs=MV.where(oro>=high_terrain, 1, 0)
        oro_rs=funcs.addExtraAxis(oro_rs,axis=0)
        oro_rs=np.repeat(oro_rs, len(ivt), axis=0)

        ivtrec_oro=morphology.reconstruction(lm*oro_rs, ivt, method='dilation',
                selem=rec_ele)
        ivtano=MV.maximum(ivt-ivtrec, (ivt-ivtrec_oro)*oro_rs)
    else:
        ivtano=ivt-ivtrec

    ivtrec=ivt-ivtano

    if ndim==4:
        levax=cdms.createAxis([0,])
        levax.designateLevel()
        levax.id='z'
        levax.name='level'
        levax.units=''

        ivt=funcs.addExtraAxis(ivt,levax,1)
        ivtrec=funcs.addExtraAxis(ivtrec,levax,1)
        ivtano=funcs.addExtraAxis(ivtano,levax,1)

    axislist=ivt.getAxisList()
    ivtrec.setAxisList(axislist)
    ivtano.setAxisList(axislist)

    ivtrec.id='ivt_rec'
    ivtrec.long_name='Integrated moisture transport, THR reconstruction'
    ivtrec.standard_name=ivtrec.long_name
    ivtrec.title=ivtrec.long_name
    ivtrec.units=ivt.units

    ivtano.id='ivt_ano'
    ivtano.long_name='Integreated moisture transport, anomaly wrt THR reconstruction'
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

    #--------------------Read in orographic data--------------------
    oro=funcs.readVar(ORO_FILE, 'oro')

    #-----------------Shift longitude-----------------
    var=var(latitude=(LAT1, LAT2))
    var=var(longitude=(SHIFT_LON,SHIFT_LON+360))
    oro=oro(latitude=(LAT1, LAT2))
    oro=oro(longitude=(SHIFT_LON,SHIFT_LON+360))(squeeze=1)

    #----------------------Do THR----------------------
    ivt, ivtrec, ivtano=filterData(var, KERNEL, oro=oro,
            high_terrain=HIGH_TERRAIN)

    #--------Save------------------------------------
    fname=os.path.split(IVT_FILE)[1]
    file_out_name='%s-THR-kernel-t%d-s%d.nc'\
            %(os.path.splitext(fname)[0], KERNEL[0], KERNEL[1])

    abpath_out=os.path.join(OUTPUTDIR,file_out_name)
    print('\n### <testrotatingfilter>: Saving output to:\n',abpath_out)
    fout=cdms.open(abpath_out,'w')
    fout.write(ivt,typecode='f')
    fout.write(ivtrec,typecode='f')
    fout.write(ivtano,typecode='f')
    fout.close()

