"""Detect atmospheric rivers (ARs) from integrated water vapor transports
(IVTs) using Top-hat by reconstruction (THR) algorithm.

In this script, read in some input data, detect ARs, collect all results
and save them to disk in one go.

See also detect_ARs_generator_version.py, detection results are yielded
once for each time step, and the results are saved to disk as long as
they are computed.

# Input data

1. uflux, vflux:

    Instantaneous vertically integrated moisture flux, in kg/m/s.
    Data should be formatted into 4D (time, singleton_z, latitude, longitude),
    or 3D (time, latitude, longitude).

2. IVT:

    This is the magnitude of vertically integrated moisture flux, i.e.
    IVT^2 = uflux^2 + vflux^2.

3. THR results:

    Instantaneous THR reconstruction and anomalies of IVT, in kg/m/s.
    This is the outcome of the THR process. See compute_thr_singlefile.py,
    compute_thr_multifile.py for scripts to perform this process.
    The THR algorithm is implemented in ipart.thr.THR.

uflux, vflux, IVT, and THR result data should be formatted into 4D
(time, singleton_z, latitude, longitude) or 3D (time, latitude, longitude),
and have compatible shapes.

The user also needs to provide time, latitude and longitude axes metadata.
In this script, these data are read in from the netCDF files using the
netCDF4 package. If you are using some other package, e.g.
CDAT, xarray, iris or something else, please adjust the relevant code
accordingly.

# Domain

Take only northern hemisphere, shift longitude to 80 E.

# Output data

1. labels:

    Labels of detected ARs. This is a 3D ndarray with dimension
    (time, lat, lon). At each time point, a unique integer label is assign
    to each detected AR, and the AR region is filled with the label value in
    the (lat, lon) map.

2. angles:

    Orientation differences between AR axes and horizontal moisture fluxes,
    measured in degrees.

3. crossfluxes:

    Cross-sectional moisture fluxes in all ARs, in kg/m/s, computed as
    the projection of the total moisture flux onto the local AR axis.

labels, angles and crossfluxes have the same dimension and are saved into a
netCDF file.

4. result_df:

    Table of detected AR records. The columns of the table includes:

        id, time, centroid_x, centroid_y, axis_x, axis_y, ... etc.

    This table is saved to a .csv file.

5. AR detection result plots (optional):

    If set `PLOT=True`, will also plot out the IVT, THR reconstruction and
    THR anomaly distributions at each time point when any AR is detected.
    The boundary of all detected ARs are also marked out.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-07-22 10:17:58.
"""

from __future__ import print_function


#######################################################################
#                               Globals                               #
#######################################################################

#--------------------Time range--------------------
YEAR=2007
TIME_START='%d-01-01 00:00:00' %YEAR
TIME_END='%d-01-10 18:00:00' %YEAR

#-----------u-qflux----------------------
UQ_FILE_NAME='/home/guangzhi/datasets/erai_qflux/uflux_m1-60_6_%d_cln-cea-proj.nc' %YEAR
UQ_VAR='uflux'

#-----------v-qflux----------------------
VQ_FILE_NAME='/home/guangzhi/datasets/erai_qflux/vflux_m1-60_6_%d_cln-cea-proj.nc' %YEAR
VQ_VAR='vflux'

#-----------------ivt reconstruction and anomalies-----------------
IVT_FILE_NAME='/home/guangzhi/datasets/quicksave2/THR/ivt_m1-60_6_%d_cln-cea-proj-THR-kernel-t10-s6.nc' %YEAR

#------------------Output folder------------------
OUTPUTDIR='/home/guangzhi/datasets/quicksave2/THR'
LABEL_FILE_OUT_NAME='ar_label-angle-flux_%d.nc' %YEAR
RECORD_FILE_OUT_NAME='ar_records_%d.csv' %YEAR



PLOT=True          # create maps of found ARs or not

LAT1=0; LAT2=90      # degree, latitude domain
SHIFT_LON=80          # degree, shift left bound to longitude.

PARAM_DICT={
    # kg/m/s, define AR candidates as regions >= than this anomalous ivt.
    # If None is given, compute a threshold based on anomalous ivt data. See
    # the docstring of ipart.AR_detector.determineThresLow() for details.
    'thres_low' : 1,
    # km^2, drop AR candidates smaller than this area.
    'min_area': 50*1e4,
    # km^2, drop AR candidates larger than this area.
    'max_area': 1800*1e4,
    # float, min length/width ratio.
    'min_LW': 2,
    # degree, exclude systems whose centroids are lower than this latitude.
    # NOTE this is the absolute latitude for both NH and SH. For SH, systems
    # with centroid latitude north of -20 will be excluded.
    'min_lat': 20,
    # degree, exclude systems whose centroids are higher than this latitude.
    # NOTE this is the absolute latitude for both NH and SH. For SH, systems
    # with centroid latitude south of -80 will be excluded.
    'max_lat': 80,
    # km, ARs shorter than this length is treated as relaxed.
    'min_length': 2000,
    # km, ARs shorter than this length is discarded.
    'min_length_hard': 1500,
    # degree lat/lon, error when simplifying axis using rdp algorithm.
    'rdp_thres': 2,
    # grids. Remove small holes in AR contour.
    'fill_radius': None,
    # do peak partition or not, used to separate systems that are merged
    # together with an outer contour.
    'single_dome': False,
    # max prominence/height ratio of a local peak. Only used when single_dome=True
    'max_ph_ratio': 0.7,
    # minimal proportion of flux component in a direction to total flux to
    # allow edge building in that direction
    'edge_eps': 0.4,
    # bool, if True, treat the data as zonally cyclic (e.g. entire hemisphere
    # or global). ARs covering regions across the longitude bounds will be
    # correctly treated as one. If your data is not zonally cyclic, or a zonal
    # shift of the data can put the domain of interest to the center, consider
    # doing the shift and setting this to False, as it will save computations.
    'zonal_cyclic': True,
    }



#--------Import modules-------------------------
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import num2date

from ipart.utils import funcs
from ipart.utils import plot
from ipart.AR_detector import plotAR, findARs



#-------------Main---------------------------------
if __name__=='__main__':

    #-----------Read in flux data----------------------
    quNV=funcs.readNC(UQ_FILE_NAME, UQ_VAR)
    qvNV=funcs.readNC(VQ_FILE_NAME, VQ_VAR)

    #-----------------Shift longitude-----------------
    quNV=quNV.shiftLon(SHIFT_LON)
    qvNV=qvNV.shiftLon(SHIFT_LON)

    #-------------------Read in ivt and THR results-------------------
    ivtNV=funcs.readNC(IVT_FILE_NAME, 'ivt')
    ivtrecNV=funcs.readNC(IVT_FILE_NAME, 'ivt_rec')
    ivtanoNV=funcs.readNC(IVT_FILE_NAME, 'ivt_ano')

    #--------------------Slice data--------------------
    quNV=quNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=2).squeeze()
    qvNV=qvNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=2).squeeze()
    ivtNV=ivtNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=2).squeeze()
    ivtrecNV=ivtrecNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=2).squeeze()
    ivtanoNV=ivtanoNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=2).squeeze()

    #--------------------Data shape check--------------------
    if np.ndim(quNV.data)!=3 or np.ndim(qvNV.data)!=3:
        raise Exception("<qu> and <qv> should be 3D data.")
    if quNV.shape!=qvNV.shape or ivtNV.shape!=quNV.shape:
        raise Exception("Data shape dismatch: qu.shape=%s; qv.shape=%s; ivt.shape=%s"\
                %(quNV.shape, qvNV.shape, ivtNV.shape))

    #-----------------Get coordinates-----------------
    latax=quNV.getLatitude()
    lonax=quNV.getLongitude()
    timeax=ivtNV.getTime()

    #######################################################################
    #                           Detect ARs                                #
    #######################################################################
    time_idx, labels, angles, crossfluxes, result_df = findARs(ivtNV.data,
            ivtrecNV.data, ivtanoNV.data, quNV.data, qvNV.data, latax, lonax,
            times=timeax, **PARAM_DICT)


    #######################################################################
    #                          Save all records                           #
    #######################################################################
    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    if PLOT:
        plot_dir=os.path.join(OUTPUTDIR, 'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

    #-----------------Save label file-----------------
    abpath_out=os.path.join(OUTPUTDIR, LABEL_FILE_OUT_NAME)
    print('\n### <detect_ARs>: Saving output to:\n',abpath_out)
    funcs.saveNC(abpath_out, labels, 'w')
    funcs.saveNC(abpath_out, angles, 'a')
    funcs.saveNC(abpath_out, crossfluxes, 'a')

    #-----------------Save record file-----------------
    abpath_out=os.path.join(OUTPUTDIR, RECORD_FILE_OUT_NAME)
    print('\n### <detect_ARs>: Saving output to:\n',abpath_out)
    # Necessary: to remove ... in csv file
    if sys.version_info.major==2:
        np.set_printoptions(threshold=np.inf)
    elif sys.version_info.major==3:
        np.set_printoptions(threshold=sys.maxsize)
    result_df.to_csv(abpath_out,index=False)


    #-------------------Create plots and save------------------------
    if PLOT:
        print('\n# <detect_ARs>: Plotting ...')
        label_timeax=num2date(labels.getTime(), labels.axislist[0].units,
                only_use_cftime_datetimes=False)

        for (ii, timett) in zip(time_idx, label_timeax):

            timett_str=str(timett)

            slab=ivtNV.data[ii]
            slabrec=ivtrecNV.data[ii]
            slabano=ivtanoNV.data[ii]
            ardf=result_df[result_df.time==timett]

            plot_vars=[slab,slabrec,slabano]
            titles=['IVT', 'THR_recon', 'THR_ano']
            iso=plot.Isofill(plot_vars,12,1,1,min_level=0,max_level=800)

            figure=plt.figure(figsize=(12,10),dpi=100)

            for jj in range(len(plot_vars)):
                ax=figure.add_subplot(3,1,jj+1)
                pobj=plot.plot2(plot_vars[jj],iso,ax,projection='cyl',
                        xarray=lonax, yarray=latax,
                        title='%s %s' %(timett_str, titles[jj]),
                        fix_aspect=False)

                bmap=pobj.bmap
                plotAR(ardf,ax,bmap)

            #----------------- Save plot------------
            plot_save_name='ar_%s' %(timett_str)
            plot_save_name=os.path.join(plot_dir,plot_save_name)
            print('\n# <detect_ARs>: Save figure to', plot_save_name)
            figure.savefig(plot_save_name+'.png',dpi=100,bbox_inches='tight')

            plt.close('all')


