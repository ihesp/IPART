'''Locate atmospheric rivers (ARs) from integrated water vapor transports
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
Update time: 2019-12-02 11:31:01.
'''

from __future__ import print_function



#######################################################################
#                               Globals                               #
#######################################################################

#--------------------Time range--------------------
YEAR=2004
TIME_START='%d-01-01 00:00:00' %YEAR
TIME_END='%d-03-31 18:00:00' %YEAR

#-----------u-qflux----------------------
SOURCEDIR1='/home/guangzhi/datasets/erai/erai_qflux'
UQ_FILE_NAME='uflux_m1-60_6_%d_cln.nc'
UQ_VAR='uflux'

#-----------v-qflux----------------------
SOURCEDIR2='/home/guangzhi/datasets/erai/erai_qflux'
VQ_FILE_NAME='vflux_m1-60_6_%d_cln.nc'
VQ_VAR='vflux'

#-----------------ivt reconstruction and anomalies-----------------
SOURCEDIR3='/home/guangzhi/datasets/erai/ivt_thr/'
IVT_FILE_NAME='ivt_m1-60_6_%d_cln-minimal-rec-ano-kernel-t16-s6.nc'

#------------------Output folder------------------
OUTPUTDIR='/home/guangzhi/datasets/erai/ERAI_AR_THR/%d/' %YEAR



PLOT=True          # create maps of found ARs or not
SINGLE_DOME=False  # do peak partition or not

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
    'min_length_hard': 800,
    # degree lat/lon, error when simplifying axis using rdp algorithm.
    'rdp_thres': 2,
    # grids. Remove small holes in AR contour.
    'fill_radius': max(1,int(4*0.75/RESO)),
    # max prominence/height ratio of a local peak.
    'max_ph_ratio': 0.4,
    # minimal proportion of flux component in a direction to total flux to
    # allow edge building in that direction
    'edge_eps': 0.4
    }





#--------Import modules-------------------------
import os
import cdms2 as cdms
import MV2 as MV
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from skimage import morphology

from utils import funcs,plot
from river_tracker1_funcs import areaFilt, maskToGraph, getARAxis, cropMask,\
        partPeaks, getARData, uvDecomp, save2DF, plotAR
from compute_thr_singlefile import readVar






def findARs(anoslab, quslab, qvslab, areas, costhetas, sinthetas, param_dict):
    '''Find ARs from a time snap

    Args:
        anoslab (cdms.TransientVariable): (n * m) 2D anomalous IVT slab, in kg/m/s.
        quslab (cdms.TransientVariable): (n * m) 2D u-flux slab, in kg/m/s.
        qvslab (cdms.TransientVariable): (n * m) 2D v-flux slab, in kg/m/s.
        areas (cdms.TransientVariable): (n * m) 2D grid cell area slab, in km^2.
        costhetas (cdms.TransientVariable): (n * m) 2D slab of grid cell shape:
                                          cos=dx/sqrt(dx^2+dy^2).
        sinthetas (cdms.TransientVariable): (n * m) 2D slab of grid cell shape:
                                          sin=dy/sqrt(dx^2+dy^2).
        param_dict (dict): parameter dict defined in Global preamble.

    Returns:
        masks (list): list of 2D binary masks, each with the same shape as
                      <anoslab> etc., and with 1s denoting the location of a
                      found AR.
        axes (list): list of AR axis coordinates. Each coordinate is defined
                     as a Nx2 ndarray storing (y, x) indices of the axis
                     (indices defined in the matrix of corresponding mask
                     in <masks>.)
        mask2 (ndarray): 2D binary mask showing all ARs in <masks> merged into
                         one map.
        axismask (ndarray): 2D binary mask showing all axes in <axes> merged
                            into one map. Overlaying <axismask> over <mask2>
                            would show all the ARs at this time step, with
                            their axes.

        If no ARs are found, <masks> and <axes> will be []. <mask2> and
        <axismask> will be zeros.
    '''

    # fetch parameters
    thres_low=param_dict['thres_low']
    min_area=param_dict['min_area']
    max_area=param_dict['max_area']
    min_lat=param_dict['min_lat']
    max_lat=param_dict['max_lat']
    fill_radius=param_dict['fill_radius']
    max_ph_ratio=param_dict['max_ph_ratio']
    edge_eps=param_dict['edge_eps']
    max_isoq_hard=param_dict['max_isoq_hard']

    mask0=np.where(anoslab>thres_low,1,0)
    mask0=areaFilt(mask0,areas,min_area,max_area)

    # prepare outputs
    masks=[]
    axes=[]
    mask2=np.zeros(mask0.shape)
    axismask=np.zeros(mask0.shape)

    if mask0.max()==0:
        return masks, axes, mask2, axismask


    #---------------Separate local peaks---------------
    labels=measure.label(mask0,connectivity=2)
    mask1=np.zeros(mask0.shape)
    latax=anoslab.getLatitude()

    for ii in range(labels.max()):
        ii+=1
        maskii=np.where(labels==ii,1,0)

        #-------------Skip if latitude too low or too high---------
        rpii=measure.regionprops(maskii, intensity_image=np.array(anoslab))[0]
        centroidy,centroidx=rpii.weighted_centroid
        centroidy=latax[int(centroidy)]
        min_lat_idx=np.argmin(np.abs(latax[:]-min_lat))

        if (centroidy<=min_lat and maskii[:min_lat_idx].sum()/float(maskii.sum())>=0.5)\
                or centroidy>=max_lat:
            continue

        if SINGLE_DOME:
            cropmask,cropidx=cropMask(maskii)
            maskii2=partPeaks(cropmask,cropidx,anoslab,max_ph_ratio)
            mask1=mask1+maskii2
        else:
            mask1=mask1+maskii

    mask1=areaFilt(mask1,areas,min_area,max_area)

    if mask1.max()==0:
        return masks, axes, mask2, axismask

    #--------------------Find axes--------------------
    labels=measure.label(mask1,connectivity=2)
    filldisk=morphology.disk(fill_radius)

    for ii in range(labels.max()):
        maskii=np.where(labels==ii+1,1,0)

        #-------------Skip if latitude too low or too high---------
        rpii=measure.regionprops(maskii, intensity_image=np.array(anoslab))[0]
        centroidy,centroidx=rpii.weighted_centroid
        centroidy=latax[int(centroidy)]
        min_lat_idx=np.argmin(np.abs(latax[:]-min_lat))

        if (centroidy<=min_lat and maskii[:min_lat_idx].sum()/float(maskii.sum())>=0.5)\
                or centroidy>=max_lat:
            continue

        # filter by isoperimetric quotient
        isoquoii=4*np.pi*rpii.area/rpii.perimeter**2

        if isoquoii>=max_isoq_hard:
            continue

        #-----------------Fill small holes-----------------
        maskii=morphology.closing(maskii,selem=filldisk)

        #----------Convert mask to directed graph----------
        gii=maskToGraph(maskii,quslab,qvslab,costhetas,sinthetas,edge_eps)

        #--------------Get AR axis from graph--------------
        axisarrii,axismaskii=getARAxis(gii,quslab,qvslab,maskii)
        axes.append(axisarrii)
        masks.append(maskii)
        axismask=axismask+axismaskii
        mask2=mask2+maskii

    return masks, axes, mask2, axismask



#-------------Main---------------------------------
if __name__=='__main__':


    #-----------Read in flux data----------------------
    file_in_name=UQ_FILE_NAME %YEAR
    abpath_in=os.path.join(SOURCEDIR1,file_in_name)
    qu=readVar(abpath_in, UQ_VAR)

    file_in_name=VQ_FILE_NAME %YEAR
    abpath_in=os.path.join(SOURCEDIR2,file_in_name)
    qv=readVar(abpath_in, VQ_VAR)

    #-----------------Shift longitude-----------------
    qu=qu(longitude=(SHIFT_LON,SHIFT_LON+360))
    qv=qv(longitude=(SHIFT_LON,SHIFT_LON+360))

    #-------------------Read in ivt-------------------
    file_in_name=IVT_FILE_NAME %YEAR
    abpath_in=os.path.join(SOURCEDIR3,file_in_name)
    print('\n### <river_tracker1>: Read in file:\n',abpath_in)
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

    dxs=funcs.dLongitude(qu,R=6371)
    dys=funcs.dLatitude(qu,R=6371)
    areamap=dxs*dys # km^2
    costhetas=dxs/MV.sqrt(dxs**2+dys**2)
    sinthetas=dys/MV.sqrt(dxs**2+dys**2)

    timeax=ivt.getTime().asComponentTime()

    #-----------------Prepare outputs-----------------
    # save found AR records:
    #     key: time str in 'yyyy-mm-dd hh:00'
    #     value: pandas dataframe. See getARData().
    result_dict={}

    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    if PLOT:
        plot_dir=os.path.join(OUTPUTDIR, 'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

    file_out_name='ar_s_6_%d_label-angle-flux.nc' %YEAR
    abpath_out=os.path.join(OUTPUTDIR,file_out_name)
    print('\n### <river_tracker1>: Saving output to:\n',abpath_out)
    ncfout=cdms.open(abpath_out,'w')



    #######################################################################
    #                          Start processing                           #
    #######################################################################

    #----------------Loop through time----------------
    for ii, timett in enumerate(timeax):

        timett_str='%d-%02d-%02d %02d:00' %(timett.year,timett.month,\
                timett.day,timett.hour)

        print('\n# <river_tracker1>: Processing time: %s' %timett_str)

        slab=ivt[ii]
        slabano=ivtano[ii]
        slabrec=ivtrec[ii]
        quslab=qu[ii]
        qvslab=qv[ii]

        # decompose background-transient
        qurec,quano,qvrec,qvano=uvDecomp(quslab,qvslab,slabrec,slabano)

        # find ARs
        mask_list,axis_list,armask,axismask=findARs(slabano,quano,qvano,
                areamap,costhetas,sinthetas,PARAM_DICT)

        # skip if none
        if armask.sum()==0:
            continue

        # fetch AR related data
        labels,angles,crossfluxes,ardf=getARData(
                slab,quslab,qvslab,
                slabano,quano,qvano,
                areamap,
                mask_list,axis_list,timett_str,PARAM_DICT,SHIFT_LON,
                False,OUTPUTDIR)

        # prepare nc output
        timeaxii=cdms.createAxis([timett.torel('days since 1900-1-1').value])
        timeaxii.designateTime()
        timeaxii.id='time'
        timeaxii.units='days since 1900-1-1'

        labels=funcs.addExtraAxis(labels,timeaxii)
        angles=funcs.addExtraAxis(angles,timeaxii)
        crossfluxes=funcs.addExtraAxis(crossfluxes,timeaxii)

        # save to disk
        ncfout.write(labels,typecode='f')
        ncfout.write(angles,typecode='f')
        ncfout.write(crossfluxes,typecode='f')

        result_dict[timett_str]=ardf

        #-------------------Plot------------------------
        if PLOT:
            plot_vars=[slab,slabrec,slabano]
            iso=plot.Isofill(plot_vars,12,1,1,min_level=0,max_level=800)

            figure=plt.figure(figsize=(12,10),dpi=100)

            for jj in range(len(plot_vars)):
                ax=figure.add_subplot(3,1,jj+1)
                pobj=plot.plot2(plot_vars[jj],iso,ax,projection='cyl',
                        title=timett_str,
                        fix_aspect=False)

                bmap=pobj.bmap
                plotAR(ardf,ax,bmap)

            #----------------- Save plot------------
            plot_save_name='ar_%s' %(timett_str)
            plot_save_name=os.path.join(plot_dir,plot_save_name)
            print('\n# <river_tracker1>: Save figure to', plot_save_name)
            figure.savefig(plot_save_name+'.png',dpi=100,bbox_inches='tight')

            plt.close('all')


    #######################################################################
    #                          Save all records                           #
    #######################################################################

    ncfout.close()

    result_df=save2DF(result_dict)
    file_out_name='ar_records_%s-%s.csv'\
            %(TIME_START.replace(' ','_').replace(':','-'),
            TIME_END.replace(' ','_').replace(':','-'))

    abpath_out=os.path.join(OUTPUTDIR,file_out_name)
    print('\n### <river_tracker1>: Saving output to:\n',abpath_out)
    # Necessary: to remove ... in csv file
    np.set_printoptions(threshold='nan')
    result_df.to_csv(abpath_out,index=False)


