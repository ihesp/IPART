"""Locate atmospheric rivers (ARs) from integrated water vapor transports
(IVTs) using Top-hat by reconstruction (THR) algorithm.

# Input data

1. uflux, vflux:

    6-hourly vertically integrated moisture flux, in kg/m/s.
    Data should be formatted into 4D (time, singleton_z, latitude, longitue),
    or 3D (time, latitude, longitude).

2. IVT:

    6-hourly THR reconstruction and anomalies of IVT, in kg/m/s.
    File names:

        ivt_m1-60_6_<year>_cln-minimal-rec-ano-kernel-t<T>-s<S>.nc

    Where:
        <year>: year in yyyy.
        <T>: kernel size in time dimension.
        <S>: kernel size in space (horizontal) dimension.

    E.g.:

        ivt_m1-60_6_2000_cln-minimal-rec-ano-kernel-t12-s8.nc

    These are generated in compute_thr_singlefile.py or compute_thr_multifile.py

    Data should be formatted into 4D (time, singleton_z, latitude, longitue),
    or 3D (time, latitude, longitude).

NOTE: data should have proper time, latitude and longitude axes!

# Domain

Take only northern hemisphere, shift longitude to 80 E.

Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
Update time: 2020-04-01 12:15:12.
"""

from __future__ import print_function



#######################################################################
#                               Globals                               #
#######################################################################

#--------------------Time range--------------------
YEAR=2007
TIME_START='%d-01-01 00:00:00' %YEAR
TIME_END='%d-01-05 18:00:00' %YEAR

#-----------u-qflux----------------------
SOURCEDIR1='/home/guangzhi/datasets/erai_qflux/'
UQ_FILE_NAME='uflux_m1-60_6_%d_cln-cea-proj.nc' %YEAR
UQ_VAR='uflux'

#-----------v-qflux----------------------
SOURCEDIR2='/home/guangzhi/datasets/erai_qflux'
VQ_FILE_NAME='vflux_m1-60_6_%d_cln-cea-proj.nc' %YEAR
VQ_VAR='vflux'

#-----------------ivt reconstruction and anomalies-----------------
SOURCEDIR3='/home/guangzhi/datasets/erai/ERAI_AR_THR/'
IVT_FILE_NAME='ivt_m1-60_6_%d_cln-cea-proj-THR-kernel-t16-s6.nc' %YEAR


#------------------Output folder------------------
OUTPUTDIR='/home/guangzhi/datasets/erai/ERAI_AR_THR/%d/' %YEAR
LABEL_FILE_OUT_NAME='ar_s_6_%d_label-angle-flux.nc' %YEAR
RECORD_FILE_OUT_NAME='ar_records_%d.csv' %YEAR



PLOT=True          # create maps of found ARs or not

LAT1=0; LAT2=90      # degree, latitude domain
# NOTE: this has to match the domain seletion in compute_thr_singlefile.py

RESO=0.75             # degree, (approximate) horizontal resolution of input data.
SHIFT_LON=80          # degree, shift left bound to longitude. Should match
                      # that used in compute_thr_singlefile.py

PARAM_DICT={
    # kg/m/s, define AR candidates as regions >= than this anomalous ivt.
    'thres_low' : 1,
    # km^2, drop AR candidates smaller than this area.
    'min_area': 50*1e4,
    # km^2, drop AR candidates larger than this area.
    'max_area': 1800*1e4,
    # float, isoperimetric quotient. ARs larger than this (more circular in shape) is treated as relaxed.
    'max_isoq': 0.6,
    # float, isoperimetric quotient. ARs larger than this is discarded.
    'max_isoq_hard': 0.7,
    # degree, exclude systems whose centroids are lower than this latitude.
    'min_lat': 20,
    # degree, exclude systems whose centroids are higher than this latitude.
    'max_lat': 80,
    # km, ARs shorter than this length is treated as relaxed.
    'min_length': 2000,
    # km, ARs shorter than this length is discarded.
    'min_length_hard': 1500,
    # degree lat/lon, error when simplifying axis using rdp algorithm.
    'rdp_thres': 2,
    # grids. Remove small holes in AR contour.
    'fill_radius': max(1,int(4*0.75/RESO)),
    # do peak partition or not, used to separate systems that are merged
    # together with an outer contour.
    'single_dome': False,
    # max prominence/height ratio of a local peak. Only used when single_dome=True
    'max_ph_ratio': 0.6,
    # minimal proportion of flux component in a direction to total flux to
    # allow edge building in that direction
    'edge_eps': 0.4
    }





#--------Import modules-------------------------
import os
import sys
import cdms2 as cdms
import numpy as np
import matplotlib.pyplot as plt

from AR_tracker.utils import funcs,plot
from AR_tracker.AR_detector import plotAR, findARsGen




#-------------Main---------------------------------
if __name__=='__main__':


    #-----------Read in flux data----------------------
    file_in_name=UQ_FILE_NAME
    abpath_in=os.path.join(SOURCEDIR1,file_in_name)
    qu=funcs.readVar(abpath_in, UQ_VAR)

    file_in_name=VQ_FILE_NAME
    abpath_in=os.path.join(SOURCEDIR2,file_in_name)
    qv=funcs.readVar(abpath_in, VQ_VAR)

    #-----------------Shift longitude-----------------
    qu=qu(longitude=(SHIFT_LON,SHIFT_LON+360))
    qv=qv(longitude=(SHIFT_LON,SHIFT_LON+360))

    #-------------------Read in ivt and THR results-------------------
    file_in_name=IVT_FILE_NAME
    abpath_in=os.path.join(SOURCEDIR3,file_in_name)
    print('\n### <detect_ARs2>: Read in file:\n',abpath_in)
    fin=cdms.open(abpath_in,'r')
    ivt=fin('ivt')
    ivtrec=fin('ivt_rec')
    ivtano=fin('ivt_ano')
    fin.close()

    #--------------------Slice data--------------------
    qu=qu(time=(TIME_START,TIME_END), latitude=(LAT1, LAT2))(squeeze=1)
    qv=qv(time=(TIME_START,TIME_END), latitude=(LAT1, LAT2))(squeeze=1)
    ivt=ivt(time=(TIME_START,TIME_END))(squeeze=1)
    ivtrec=ivtrec(time=(TIME_START,TIME_END))(squeeze=1)
    ivtano=ivtano(time=(TIME_START,TIME_END))(squeeze=1)

    #--------------------Data shape check--------------------
    if np.ndim(qu)!=3 or np.ndim(qv)!=3:
        raise Exception("<qu> and <qv> should be 3D data.")
    if qu.shape!=qv.shape or ivt.shape!=qu.shape:
        raise Exception("Data shape dismatch: qu.shape=%s; qv.shape=%s; ivt.shape=%s"\
                %(qu.shape, qv.shape, ivt.shape))

    #-----------------Get coordinates-----------------
    latax=qu.getLatitude()
    lonax=qu.getLongitude()
    timeax=ivt.getTime().asComponentTime()
    timeax=['%d-%02d-%02d %02d:00' %(timett.year,timett.month,\
                timett.day,timett.hour) for timett in timeax]
    timeax=None

    #-----------------Prepare outputs-----------------
    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    if PLOT:
        plot_dir=os.path.join(OUTPUTDIR, 'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

    # nc file to save AR location labels
    abpath_out=os.path.join(OUTPUTDIR, LABEL_FILE_OUT_NAME)
    print('\n### <detect_ARs2>: Saving output to:\n',abpath_out)
    ncfout=cdms.open(abpath_out,'w')

    # csv file to save AR record table
    abpath_out=os.path.join(OUTPUTDIR, RECORD_FILE_OUT_NAME)
    print('\n### <detect_ARs2>: Saving output to:\n',abpath_out)
    # Necessary: to remove ... in csv file
    if sys.version_info.major==2:
        np.set_printoptions(threshold=np.inf)
    elif sys.version_info.major==3:
        np.set_printoptions(threshold=sys.maxsize)

    with open(abpath_out, 'a') as dfout:

        #############################################################
        #                     Start processing                      #
        #############################################################
        finder_gen = findARsGen(ivt, ivtrec, ivtano, qu, qv, latax, lonax,
                PARAM_DICT, times=timeax)
        next(finder_gen)  # prime the generator to prepare metadata

        for (tidx, timett, label, angle, cross, result_df) in finder_gen:

            #------------------Save record to csv file------------------
            result_df.to_csv(dfout, header=dfout.tell()==0, index=False)

            #-------------------Save labels to nc file-------------------
            ncfout.write(label)
            ncfout.write(angle,typecode='f')
            ncfout.write(cross,typecode='f')

            #-------------------Plot------------------------
            if PLOT:

                timett_str=str(timett)

                slab=ivt[tidx]
                slabrec=ivtrec[tidx]
                slabano=ivtano[tidx]

                plot_vars=[slab,slabrec,slabano]
                titles=['IVT', 'Reconstruction', 'THR']
                iso=plot.Isofill(plot_vars,12,1,1,min_level=0,max_level=800)

                figure=plt.figure(figsize=(12,10),dpi=100)

                for jj in range(len(plot_vars)):
                    ax=figure.add_subplot(3,1,jj+1)
                    pobj=plot.plot2(plot_vars[jj],iso,ax,projection='cyl',
                            title='%s %s' %(timett_str, titles[jj]),
                            fix_aspect=False)

                    bmap=pobj.bmap
                    plotAR(result_df,ax,bmap)

                #----------------- Save plot------------
                plot_save_name='ar_%s' %(timett_str)
                plot_save_name=os.path.join(plot_dir,plot_save_name)
                print('\n# <detect_ARs2>: Save figure to', plot_save_name)
                figure.savefig(plot_save_name+'.png',dpi=100,bbox_inches='tight')

                plt.close('all')


    #######################################################################
    #                          Save all records                           #
    #######################################################################
    ncfout.close()



