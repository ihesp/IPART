"""Detect atmospheric rivers (ARs) from integrated water vapor transports
(IVTs) using Top-hat by reconstruction (THR) algorithm.

In this script, read in some input data, detect ARs, and the detection results
are yielded once for each time step, and the results are saved to disk as long
as they are computed.

See also detect_ARs.py, where the detection results are collected and saved
to disk in one go.

# Input data

1. uflux, vflux:

    Instantaneous vertically integrated moisture flux, in kg/(m*s).
    u-moisture flux component's standard_name:
        "eastward_atmosphere_water_transport_across_unit_distance".
    v-moisture flux component's standard_name:
        "northward_atmosphere_water_transport_across_unit_distance".

    Data should be formatted into 4D (time, singleton_z, latitude, longitude),
    or 3D (time, latitude, longitude).

2. IVT:

    This is the magnitude of vertically integrated moisture flux, i.e.
    IVT^2 = uflux^2 + vflux^2.

3. THR results:

    Instantaneous THR reconstruction and anomalies of IVT, in kg/(m*s).
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

    Cross-sectional moisture fluxes in all ARs, in kg/(m*s), computed as
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
Update time: 2020-07-22 10:13:31.
"""

from __future__ import print_function


#######################################################################
#                               Globals                               #
#######################################################################
ADD_DUMMY_COMP=False
#--------------------Time range--------------------
YEAR=1979
TIME_START='%d-01-01 00:00:00' %YEAR
TIME_END='%d-12-31 23:59:00' %YEAR

PROJECTION = 'nplaea'
BLAT0 = -20 if PROJECTION == 'splaea' else 20

THR_FILE_NAME_BASE = 'ivt_s_6_%s_cln-%s-proj[b1c889-ebf31ba]-THR-kernel-t0-s19-latweight.nc'

PARAM_SUFFIX = THR_FILE_NAME_BASE.index('THR')
PARAM_SUFFIX = THR_FILE_NAME_BASE[PARAM_SUFFIX:-3]
OUTPUTDIR='/run/media/guangzhi/ERA5/ERA5_AR/ARsno_oro_polar/' + PARAM_SUFFIX
OUTPUTDIR='/home/guangzhi/datasets/era5/ARsno_oro_polar/' + PARAM_SUFFIX

LAT1 = 23
PLOT=True          # create maps of found ARs or not

#LAT1=0; LAT2=90      # degree, latitude domain
#SHIFT_LON=80          # degree, shift left bound to longitude.

PARAM_DICT={
    # kg/m/s, define AR candidates as regions >= than this anomalous ivt.
    # If None is given, compute a threshold based on anomalous ivt data. See
    # the docstring of ipart.AR_detector.determineThresLow() for details.
    'thres_low' : 8.9,
    # km^2, drop AR candidates smaller than this area.
    #'min_area': 50*1e4,
    'min_area': 40*1e4,
    # km^2, drop AR candidates larger than this area.
    'max_area': 1800*1e4,
    # float, min length/width ratio.
    'min_LW': 2,
    # degree, exclude systems whose centroids are lower than this latitude.
    # NOTE this is the absolute latitude for both NH and SH. For SH, systems
    # with centroid latitude north of -20 will be excluded.
    'min_lat': 23,
    # degree, exclude systems whose centroids are higher than this latitude.
    # NOTE this is the absolute latitude for both NH and SH. For SH, systems
    # with centroid latitude south of -80 will be excluded.
    'max_lat': 90,
    # km, ARs shorter than this length is treated as relaxed.
    'min_length': 1500,
    # km, ARs shorter than this length is discarded.
    'min_length_hard': 1000,
    # degree lat/lon, error when simplifying axis using rdp algorithm.
    'rdp_thres': 2,
    # grids. Remove small holes in AR contour.
    'fill_radius': 2,
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
    'zonal_cyclic': False,
    }



#--------Import modules-------------------------
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import cartopy.crs as ccrs
from netCDF4 import Dataset

from ipart.utils import funcs
#from ipart.utils import plot
from tools import plot
from ipart.AR_detector import findARsGenPolar


def dummy():
    x=np.random.random([500,500])
    np.linalg.inv(x)*x.T

def plotVar(var, xx, yy, method, ax, bmap, min_lat, projection, title=None, labels=None,
        label_color='yellow'):

    ax.patch.set_color('0.5')

    if method.ext_1 is False and method.ext_2 is False:
        extend='neither'
    elif method.ext_1 is True and method.ext_2 is False:
        extend='min'
    elif method.ext_1 is False and method.ext_2 is True:
        extend='max'
    else:
        extend='both'

    if projection == 'nplaea':
        lats=np.arange(20,90,30)
    else:
        lats=np.arange(-80,0,30)

    lons=np.arange(-180,181,30)

    if labels is None:
        meridians=[1,1,1,1]
    else:
        meridians=labels

    bmap.drawcoastlines(linewidth=0.5)
    bmap.drawparallels(lats,linewidth=0.5)
    #bmap.drawparallels([min_lat,], linewidth=2.0, color='yellow', dashes=[1,0])
    min_latx = np.arange(-180, 181, 1)
    min_laty = np.ones(len(min_latx)) * min_lat
    min_latx, min_laty = bmap(min_latx, min_laty)
    ax.plot(min_latx, min_laty, color='w', linestyle='--', lw=1.4)
    bmap.drawmeridians(lons,labels=meridians,linewidth=0.5,\
            labelstyle='+/-',fontsize=10)

    if method.method=='boxfill':
        cs=bmap.imshow(var,cmap=method.cmap,ax=ax,
                vmin=method.vmin,vmax=method.vmax,interpolation='nearest')
    else:
        cs=bmap.contourf(xx,yy,np.array(var),method.levels,\
                cmap=method.cmap,extend=extend,latlon=True)

    for ll in lats:
        xll,yll=bmap(180., ll)
        ax.text(xll, yll, u'%+d\N{DEGREE SIGN}' %ll,
                color=label_color,
                fontsize=9,
                horizontalalignment='center',
                verticalalignment='center')

    if title is not None:
        ax.set_title(title,loc='left', pad=20)

    return cs,extend


def breakCurveAtEdge(xs, ys, dx=10):
    '''Segment curve coordinates at the left, right edges
    '''

    idx=[]  # break point indices
    new_xs=[] # result list for x coordinates segments
    new_ys=[]
    diff = abs(np.diff(xs))
    idx = np.where(diff>dx)[0]

    if len(idx)==0:
        new_xs.append(xs)
        new_ys.append(ys)
    else:
        idx += 1
        idx = list(idx)
        idx.insert(0,0)
        idx.append(len(xs))

        for i1, i2 in zip(idx[:-1], idx[1:]):
            new_xs.append(xs[i1:i2])
            new_ys.append(ys[i1:i2])

    return new_xs, new_ys

def plotAR(ardf, ax, bmap):
    '''Helper function to plot the regions and axes of ARs

    Args:
        ardf (pandas.DataFrame): table containing AR records.
        ax (matplotlib axis): axis to plot onto.
        lonax (ndarray): 1d array of the longitude axis the plot is using.
    '''

    for ii in range(len(ardf)):

        vv=ardf.iloc[ii]
        isrelaxkk=vv['is_relaxed']

        # plot contour
        px=vv['contour_x']
        py=vv['contour_y']

        px_segs, py_segs=breakCurveAtEdge(px, py)
        # note that the GeoAxes (of cartopy) doesn't seem to carry these info.

        for xjj, yjj in zip(px_segs, py_segs):
            if len(xjj) < 2:
                continue

            linewidth=1.0 if isrelaxkk else 1.0
            linestyle=':' if isrelaxkk else '-'
            xjj,yjj = bmap(xjj,yjj)
            ax.plot(xjj,yjj,color='k',linestyle=linestyle,linewidth=linewidth)

        # plot axis
        px=vv['axis_x']
        py=vv['axis_y']

        px, py = bmap(px, py)
        ax.plot(px,py,'g-',linewidth=1.3)

    return


def sliceLat(var, lats, min_lat):

    lat_idx = np.where(abs(lats.data) >= abs(min_lat))
    ymin = np.min(lat_idx[0])
    ymax = np.max(lat_idx[0])
    xmin = np.min(lat_idx[1])
    xmax = np.max(lat_idx[1])
    var.data = var.data[..., ymin:ymax+1, xmin:xmax+1]

    return var




#-------------Main---------------------------------
if __name__=='__main__':

    for MONTH in range(1, 13):

        YEARMONTH='%d%s' %(YEAR, str(MONTH).rjust(2, '0'))

        #-----------------Get coordinates-----------------
        COORD_FILE_NAME='/run/media/guangzhi/ERA5/ERA5_AR/IVT_polar/area_s_a_1900_%s-proj[b1c889-ebf31ba].nc'\
                %(PROJECTION)
        latax=funcs.readNC(COORD_FILE_NAME, 'lat')
        lonax=funcs.readNC(COORD_FILE_NAME, 'lon')
        dlats=funcs.readNC(COORD_FILE_NAME, 'dlat')
        dlons=funcs.readNC(COORD_FILE_NAME, 'dlon')
        areas=funcs.readNC(COORD_FILE_NAME, 'area')

        #-----------u-qflux----------------------
        UQ_FILE_NAME='/run/media/guangzhi/ERA5/ERA5_AR/IVT_polar/uq_s_6_%s_cln-%s-proj[b1c889-ebf31ba].nc'\
                %(YEARMONTH, PROJECTION)
        UQ_VAR='uq'

        #-----------v-qflux----------------------
        VQ_FILE_NAME='/run/media/guangzhi/ERA5/ERA5_AR/IVT_polar/vq_s_6_%s_cln-%s-proj[b1c889-ebf31ba].nc'\
                %(YEARMONTH, PROJECTION)
        VQ_VAR='vq'

        #-----------------ivt reconstruction and anomalies-----------------
        IVT_FILE_NAME='/run/media/guangzhi/ERA5/ERA5_AR/THR_polar/' + \
                THR_FILE_NAME_BASE %(YEARMONTH, PROJECTION)

        #------------------Output folder------------------
        LABEL_FILE_OUT_NAME='ar_label-angle-flux_%s.nc' %YEARMONTH
        RECORD_FILE_OUT_NAME='ar_records_%s.csv' %YEARMONTH


        #-----------Read in flux data----------------------
        quNV=funcs.readNC(UQ_FILE_NAME, UQ_VAR)
        qvNV=funcs.readNC(VQ_FILE_NAME, VQ_VAR)

        #-------------------Read in ivt and THR results-------------------
        ivtNV=funcs.readNC(IVT_FILE_NAME, 'ivt')
        #ivtrecNV=funcs.readNC(IVT_FILE_NAME, 'ivt_rec')
        ivtanoNV=funcs.readNC(IVT_FILE_NAME, 'ivt_ano')
        ivtrec = ivtNV.data - ivtanoNV.data
        ivtrecNV = funcs.NCVAR(ivtrec, 'ivt_rec', ivtNV.axislist, None)

        #--------------------Slice data--------------------
        quNV=quNV.sliceData(TIME_START, TIME_END).squeeze()
        qvNV=qvNV.sliceData(TIME_START, TIME_END).squeeze()
        ivtNV=ivtNV.sliceData(TIME_START, TIME_END).squeeze()
        ivtrecNV=ivtrecNV.sliceData(TIME_START, TIME_END).squeeze()
        ivtanoNV=ivtanoNV.sliceData(TIME_START, TIME_END).squeeze()

        quNV = sliceLat(quNV, latax, LAT1)
        qvNV = sliceLat(qvNV, latax, LAT1)
        lonax = sliceLat(lonax, latax, LAT1)
        dlats = sliceLat(dlats, latax, LAT1)
        dlons = sliceLat(dlons, latax, LAT1)
        areas = sliceLat(areas, latax, LAT1)
        latax = sliceLat(latax, latax, LAT1)

        #--------------------Data shape check--------------------
        if np.ndim(quNV.data)!=3 or np.ndim(qvNV.data)!=3:
            raise Exception("<qu> and <qv> should be 3D data.")
        if quNV.shape!=qvNV.shape or ivtNV.shape!=quNV.shape:
            raise Exception("Data shape dismatch: qu.shape=%s; qv.shape=%s; ivt.shape=%s"\
                    %(quNV.shape, qvNV.shape, ivtNV.shape))

        timeax=ivtNV.getTime()

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
        ncfout=Dataset(abpath_out, 'w')

        # csv file to save AR record table
        abpath_out=os.path.join(OUTPUTDIR, RECORD_FILE_OUT_NAME)
        print('\n### <detect_ARs2>: Saving output to:\n',abpath_out)
        # Necessary: to remove ... in csv file
        if sys.version_info.major==2:
            np.set_printoptions(threshold=np.inf)
        elif sys.version_info.major==3:
            np.set_printoptions(threshold=sys.maxsize)

        if os.path.exists(abpath_out):
            print('\n# <cppl_detect_ARs>: Removing existing file:', abpath_out)
            os.remove(abpath_out)

        with open(abpath_out, 'a') as dfout:

            #############################################################
            #                     Start processing                      #
            #############################################################
            finder_gen = findARsGenPolar(ivtNV.data, ivtrecNV.data, ivtanoNV.data,
                    quNV.data, qvNV.data, latax.data, lonax.data,
                    dlats.data, dlons.data, areas.data, times=timeax, **PARAM_DICT)
            next(finder_gen)  # prime the generator to prepare metadata

            for (tidx, timett, label, result_df) in finder_gen:

                #if ADD_DUMMY_COMP:
                    #dummy()

                #------------------Save record to csv file------------------
                result_df.to_csv(dfout, header=dfout.tell()==0, index=False)

                #-------------------Save labels to nc file-------------------
                funcs.saveNCDims(ncfout, label.axislist)
                funcs._saveNCVAR(ncfout, label, 'int')
                #funcs._saveNCVAR(ncfout, angle)
                #funcs._saveNCVAR(ncfout, cross)

                if len(result_df) == 0:
                    continue

                #-------------------Plot------------------------
                if PLOT:

                    bmap=Basemap(projection=PROJECTION,
                            boundinglat=BLAT0,lon_0=180,
                            resolution='l',
                            fix_aspect=True)

                    timett_str=str(timett)

                    slab=ivtNV.data[tidx]
                    slabrec=ivtrecNV.data[tidx]
                    slabano=ivtanoNV.data[tidx]

                    plot_vars=[slab,slabrec,slabano]
                    titles=['(a) IVT', '(b) Reconstruction', '(c) THR']
                    labels=[[1,0,1,1], [0,0,1,1], [0,1,1,1]]
                    iso=plot.Isofill(plot_vars,16,1,1,min_level=0,max_level=800)

                    figure=plt.figure(figsize=(15, 5),dpi=100)
                    min_lat = PARAM_DICT['min_lat']
                    min_lat *= 1 if PROJECTION=='nplaea' else -1

                    for jj in range(len(plot_vars)):
                        ax=figure.add_subplot(1, len(plot_vars), jj+1)
                        cs, extend=plotVar(plot_vars[jj], lonax.data, latax.data,
                                iso, ax, bmap, min_lat, PROJECTION, titles[jj], labels[jj])
                        plotAR(result_df,ax,bmap)

                    # colorbar
                    cax=figure.add_axes([0.15,0.05,0.7,0.02])

                    cbar=plt.colorbar(cs,cax=cax,orientation='horizontal',\
                            ticks=None,drawedges=False,extend=extend,
                            aspect=35)

                    cbar.set_label('kg/(m s)', fontsize=9)

                    #----------------- Save plot------------
                    plot_save_name='ar_%s' %(timett_str)
                    plot_save_name=os.path.join(plot_dir,plot_save_name)
                    print('\n# <detect_ARs2>: Save figure to', plot_save_name)
                    figure.savefig(plot_save_name+'.png',dpi=100,bbox_inches='tight')

                    plt.close('all')


        #----------------Close the nc file----------------
        ncfout.close()



