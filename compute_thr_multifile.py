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

Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
Update time: 2019-12-06 22:55:40.
'''

from __future__ import print_function


#--------------Globals------------------------------------------
#--------------------Year range--------------------
YEARS=range(2003,2006)

#-----------IVT data----------------------
SOURCEDIR1='/home/guangzhi/datasets/erai/erai_qflux'
FILE1_BASE='ivt_m1-60_6_%d_cln.nc'
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
OUTPUTDIR='/home/guangzhi/datasets/erai/ivt_thr/'







#--------Import modules-------------------------
import os
import cdms2 as cdms
import MV2 as MV
import numpy as np
from cdms2.selectors import Selector
from utils import funcs
from compute_thr_singlefile import filterData



def rotatingFiltering(filelist, varin, selector, kernel, outputdir, oro=None,
        high_terrain=600, verbose=True):
    '''Compute time filtering on data in different files.

    Args:
        filelist (list): list of abs paths to data files. User is responsible to
                    make sure files in list have correct chronological order.
                    Note that time axis in data files should be at the 1st axis.
        varin (str): variable id in files.
        selector: selector obj to select subset of data.
        outputdir (str): path to folder to save outputs.

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

    Designed to perform temporal filtering on data that are too large to fit
    into memory, e.g. high-resolution data across multiple decades.

    Function will read in 2 files at once, call the filtering function on the
    concatenated data, and shift 1 step to the next 2 files. If at the begining,
    pad 0s to the left end. If in the mid, pad filtered data in the mid of
    the concatenated data in the previous step. If at the end, pad 0s to
    the right end.

    The filtering function <func> is assumed to apply a filtering window with
    odd length n, and truncates (n-1)/2 points from both ends. If the function
    doesn't truncate data, will raise an exception.
    '''

    #----------------Check input files----------------
    funcs.checkFiles(filelist)

    for ii, fii in enumerate(filelist[:-1]):

        if ii==0:
            var1=funcs.readVar(fii, varin)
            var1=var1(selector)
            var1=var1(longitude=(SHIFT_LON,SHIFT_LON+360))
        else:
            var1=var2
            del var2

        fii2=filelist[ii+1]
        var2=funcs.readVar(fii2, varin)
        var2=var2(selector)
        var2=var2(longitude=(SHIFT_LON,SHIFT_LON+360))

        timeidx=funcs.interpretAxis('time',var1)
        if timeidx!=0:
            raise Exception("Time axis in variable is not at axis 0.")
        timeidx=funcs.interpretAxis('time',var2)
        if timeidx!=0:
            raise Exception("Time axis in variable is not at axis 0.")

        n1=var1.shape[0]

        vartmp=funcs.cat(var1,var2,axis=0)
        vartmp, vartmp_rec, vartmp_ano=filterData(vartmp, kernel, oro=oro,
            high_terrain=high_terrain)

        # crop end points
        dt=kernel[0]
        vartmp_rec=vartmp_rec[dt:-dt]
        vartmp_ano=vartmp_ano[dt:-dt]

        if dt<=0:
            raise Exception("dt<=0 not supported yet")

        if verbose:
            print('\n# <rotatingFiltering>: Concatenated var shape:',vartmp.shape)
            print('# <rotatingFiltering>: Filtered var shape:',vartmp_rec.shape)
            print('# <rotatingFiltering>: Length difference:',dt)

        if ii==0:
            #----------------------Pad 0s----------------------
            left_rec=MV.zeros((dt,)+var1.shape[1:])
            left_rec.mask=True
            left_ano=MV.zeros((dt,)+var1.shape[1:])
            left_ano.mask=True
        else:
            left_rec=mid_left_rec
            left_ano=mid_left_ano

        rec_pad=funcs.cat(left_rec,vartmp_rec,axis=0)
        rec1=rec_pad[:n1]

        ano_pad=funcs.cat(left_ano,vartmp_ano,axis=0)
        ano1=ano_pad[:n1]

        var1=vartmp[:n1]

        if verbose:
            print('\n# <rotatingFiltering>: Shape of left section after padding:', rec1.shape)

        rec1.setAxisList(var1.getAxisList())
        rec1.id=vartmp_rec.id
        rec1.long_name=var1.long_name+' rotating filtered'
        rec1.standard_name=rec1.long_name
        rec1.title=rec1.long_name
        rec1.units=var1.units

        ano1.setAxisList(var1.getAxisList())
        ano1.id=vartmp_ano.id
        ano1.long_name=var1.long_name+' rotating filtered'
        ano1.standard_name=ano1.long_name
        ano1.title=ano1.long_name
        ano1.units=var1.units

        # left to pad in next iteration
        mid_left_rec=vartmp_rec[n1-dt:n1]
        mid_left_ano=vartmp_ano[n1-dt:n1]

        #-----------------------Save-----------------------
        fname=os.path.split(fii)[1]
        file_out_name='%s-THR-kernel-t%d-s%d.nc'\
                %(os.path.splitext(fname)[0], kernel[0], kernel[1])

        abpath_out=os.path.join(outputdir,file_out_name)
        print('\n### <testrotatingfilter>: Saving output to:\n',abpath_out)
        fout=cdms.open(abpath_out,'w')
        fout.write(var1,typecode='f')
        fout.write(rec1,typecode='f')
        fout.write(ano1,typecode='f')
        fout.close()

        # save the right section for the last file
        if ii==len(filelist)-2:
            right=MV.zeros((dt,)+var2.shape[1:])
            right.mask=True
            rec2=rec_pad[n1:]
            rec2=funcs.cat(rec2,right,axis=0)

            ano2=ano_pad[n1:]
            ano2=funcs.cat(ano2,right,axis=0)

            var2=vartmp[n1:]

            if verbose:
                print('\n# <rotatingFiltering>: Shape of last section after padding:', ano2.shape)

            rec2.setAxisList(var2.getAxisList())
            rec2.id=vartmp_rec.id
            rec2.long_name=var2.long_name+' rotating filtered'
            rec2.standard_name=rec2.long_name
            rec2.title=rec2.long_name
            rec2.units=var2.units

            ano2.setAxisList(var2.getAxisList())
            ano2.id=vartmp_ano.id
            ano2.long_name=var2.long_name+' rotating filtered'
            ano2.standard_name=ano2.long_name
            ano2.title=ano2.long_name
            ano2.units=var2.units

            #-----------------------Save-----------------------
            fname=os.path.split(fii2)[1]
            file_out_name='%s-THR-kernel-t%d-s%d.nc'\
                    %(os.path.splitext(fname)[0], kernel[0], kernel[1])
            abpath_out=os.path.join(outputdir,file_out_name)
            print('\n### <testrotatingfilter>: Saving output to:\n',abpath_out)
            fout=cdms.open(abpath_out,'w')
            fout.write(var2,typecode='f')
            fout.write(rec2,typecode='f')
            fout.write(ano2,typecode='f')
            fout.close()


    return




#-------------Main---------------------------------
if __name__=='__main__':


    filelist=[]
    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    #--------------------Read in orographic data--------------------
    oro=funcs.readVar(ORO_FILE, 'oro')
    oro=oro(latitude=(LAT1, LAT2))
    oro=oro(longitude=(SHIFT_LON,SHIFT_LON+360))(squeeze=1)

    for year in YEARS:
        #-----------Read in data----------------------
        file_in_name=FILE1_BASE %(year)
        abpath_in=os.path.join(SOURCEDIR1,file_in_name)
        filelist.append(abpath_in)

    if len(filelist)<2:
        raise Exception("Need to give at least 2 files. For single file, use compute_thr_singlefile.py")

    selector=Selector(latitude=(LAT1,LAT2))
    rotatingFiltering(filelist, VARIN, selector, KERNEL, OUTPUTDIR,
            oro=oro, high_terrain=HIGH_TERRAIN)

