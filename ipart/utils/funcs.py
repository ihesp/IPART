'''Utility functions

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-06-05 23:23:22.
'''
from __future__ import print_function


class DataError(Exception):
    def __init__(self, string=None):
        if string != None:
            self.message = string
    def __str__(self):
        return self.message


def get3DEllipse(t,y,x):
    '''Get a binary 3D ellipse structuring element

    Args:
        t (int): ellipse axis length in the t (1st) dimension.
        y (int): ellipse axis length in the y (2nd) dimension.
        x (int): ellipse axis length in the x (3rd) dimension.
        Note that the axis length is half the size of the ellipse
        in that dimension.

    Returns:
        result (ndarray): 3D binary array, with 1s side the ellipse
            defined as (T/t)^2 + (Y/y)^2 + (X/x)^2 <= 1.
    '''
    import numpy as np

    at=np.arange(-t,t+1)
    ax=np.arange(-x,x+1)
    ay=np.arange(-y,y+1)
    T,Y,X=np.meshgrid(at,ay,ax,indexing='ij')
    dd=(X/float(x))**2+(Y/float(y))**2+(T/float(t))**2
    result=np.where(dd<=1,1,0)

    return result


def getQuantiles(slab,percents=None,verbose=False):
    '''Find quantiles of a slab

    <slab>: ndarray, whose quantiles will be found.
    <percents>: float or a list of floats, left percentage(s). Right quantiles
                will be found by (1-percentage).

    Return <quantiles>: nested list of left and right quantiles for corresponding
                       percentages.
    '''

    import numpy

    if percents is None:
        percents=numpy.array([0.001,0.005,0.01,0.025,0.05,0.1])
    percents=numpy.array(percents)
    if percents.ndim!=1:
        raise Exception("<percents> needs to be a 1D array.")

    #-------Remove nans and masked values--------
    mask=getMissingMask(slab)
    slab=numpy.array(slab)
    slab=slab[numpy.where(mask==False)]

    flatten=slab.flatten()
    flatten.sort()
    n=len(flatten)

    qlidx=(n*percents).astype('int')
    qridx=(n*(1-percents)).astype('int')
    ql=flatten[qlidx]
    qr=flatten[qridx]

    quantiles=list(zip(ql,qr))

    if verbose:
        for ii,pii in enumerate(percents):
            print('# <getQuantiles>: %0.3f left quantile: %f.  %0.3f right quantile: %f.'\
                    %(pii,ql[ii],1-pii,qr[ii]))

    return quantiles


#-------Copies selected attributes from source object to dict--
def attribute_obj2dict(source_object,dictionary=None,verbose=False):
    '''Copies selected attributes from source object to dict
    to <dictionary>.

    <source_object>: object from which attributes are copied.
    <dictionary>: None or dict. If None, create a new dict to store
                  the result. If a dict, use attributes from <source_object>
                  to overwrite or fill the dict.

    Update time: 2016-01-18 11:00:55.
    '''


    if dictionary is None:
        dictionary={}

    #------------------Attribute list------------------
    att_list=['name','id','dataset','source','title','long_name','standard_name',\
            'units','syno','end','harms','filename','comments','description']

    #-----------------Copy attributes-----------------
    for att in att_list:
        if hasattr(source_object,att):
            dictionary[att]=getattr(source_object,att).strip()
            if verbose:
                print('\n# <attribute_obj2dict>: %s: %s' %(att, dictionary[att]))

    return dictionary


#-------------Copy attributes from dict to target object----------
def attribute_dict2obj(dictionary,target_object,verbose=False):
    '''Copies attributes from dictionary to target object.

    <dictionary>: dict, contains attributes to copy.
    <target_object>: obj, attributes are copied to.

    Return <target_object>: target object with new attributes.

    Update time: 2016-01-18 11:31:25.
    '''

    for att in dictionary.keys():
        setattr(target_object,att,dictionary[att])
        if verbose:
            print('\n# <attribute_dict2obj>: Copy attribute: %s = %s' %(att,dictionary[att]))

    return target_object

#-------------------Add an extra axis to a data slab -------------
def addExtraAxis(slab,newaxis=None,axis=0,verbose=False):
    """Adds an extra axis to a data slab.

    <slab>: variable to which the axis is to insert.
    <newaxis>: axis object, could be of any length. If None, create a dummy
               singleton axis.
    <axis>: index of axis to be inserted, e.g. 0 if <newaxis> is inserted
            as the 1st dimension.

    Return: <slab2>.
    Update time: 2013-10-09 12:34:32.
    """

    import cdms2 as cdms
    import MV2 as MV

    if newaxis is None:
        newaxis=cdms.createAxis([1,])
        newaxis.units=''

    # add new axis to axis list of input <slab>
    axislist=slab.getAxisList()
    axislist.insert(axis,newaxis)

    #----------------Reshape----------------
    shape=list(slab.shape)
    shape.insert(axis,len(newaxis))
    slab2=MV.reshape(slab,shape)

    #------------Create variable------------
    att_dict=attribute_obj2dict(slab)
    slab2=cdms.createVariable(slab2,axes=axislist,attributes=att_dict,\
            typecode='f')
    slab2.id=slab.id

    if verbose:
        print('\n# <addExtraAxis>: Originial variable shape:',slab.shape)
        print('# <addExtraAxis>: New variable shape:',slab2.shape)

    return slab2


#-------------Concatenate transient variables---------------------
def cat(var1,var2,axis=0,verbose=False):
    '''Concatenate 2 variables along axis.

    <var1>,<var2>: Variables to be concatenated, in the order of \
            <var1>, <var2>;
    <axis>: int, index of axis to be concatenated along.

    Return <result>
    '''

    import MV2 as MV
    import numpy

    try:
        order=var1.getAxisListIndex()
    except:
        order=numpy.arange(var1.ndim)  # if var1 is np.ndarray

    var1=MV.array(var1)
    var2=MV.array(var2)

    try:
        attdict=attribute_obj2dict(var1)
        hasatt=True
    except:
        hasatt=False

    if not hasattr(var1.getAxis(axis),'units'):
        ax=var1.getAxis(axis)
        ax.units=''
        var1.setAxis(axis,ax)
    if not hasattr(var2.getAxis(axis),'units'):
        ax=var2.getAxis(axis)
        ax.units=''
        var2.setAxis(axis,ax)

    if verbose:
        print('# <cat>: Original order:',order)

    if axis!=0:
        #----Switch order------
        order[axis]=0
        order[0]=axis
        if verbose:
            print('# <cat>: New order:',order)

        var1=var1(order=order)
        var2=var2(order=order)

        result=MV.concatenate((var1,var2))
        #result=numpy.concatenate((var1,var2),axis=0)

        #NOTE: There seems to be some problems with MV.concatenate() when axis
        # is not 0, but can not remember what the problem is. That is why this function
        # is written.
        # And also some issues regards to the re-ordering and MV.concatenate()
        # method defined here. When I concatenated something along the 2nd
        # axis and do a MV.std(var,axis=2) (and numpy.std(), an attributeError was raised.
        # But other times it works ok. Maybe because of some attributes of my
        # variable is gone when putting into MV.std(). No idea why.
        # That problem was solved by replacing MV.concatenate() with numpy.concatenate().
        # But this will cause the output to be numpy.ndarray rather than MV.transientVariable.
        # So be aware that this function may cause some errors if inputs <var1>,<var2>
        # are numpy.ndarray.

        #-------Switch back----------
        result=result(order=order)
    else:
        result=MV.concatenate((var1,var2))

    if hasatt:
        result=attribute_dict2obj(attdict,result)

    return result


#------Interpret and convert an axis id to index----------
def interpretAxis(axis,ref_var,verbose=True):
    '''Interpret and convert an axis id to index

    <axis>: axis option, integer or string.
    <ref_var>: reference variable.

    Return <axis_index>: the index of required axis in <ref_var>.

    E.g. index=interpretAxis('time',ref_var)
         index=0

         index=interpretAxis(1,ref_var)
         index=1

    Update time: 2013-09-23 13:36:53.
    '''

    import sys
    import numpy

    if isinstance(axis,(int,numpy.integer)):
        return axis

    # interpret string dimension
    #elif type(axis)==type('t'):
    elif isinstance(axis,str if sys.version_info[0]>=3 else basestring):
        axis=axis.lower()

        if axis in ['time', 'tim', 't']:
            dim_id = 'time'
        elif axis in ['level', 'lev','z']:
            dim_id = 'level'
        elif axis in ['latitude', 'lat','y']:
            dim_id = 'latitude'
        elif axis in ['longitude', 'long', 'lon','x']:
            dim_id = 'longitude'
        else:
            dim_id = axis

        dim_index = ref_var.getAxisIndex(dim_id)

        if dim_index==-1:
            raise Exception("Required dimension not in <var>.")

        return dim_index

    else:
        raise Exception("<axis> type not recognized.")


#----------Check exsitance of files in file list-----------
def checkFiles(file_list,verbose=True):
    '''Check existance of files in a list.

    <file_list>: a list of ABSOLUTE paths to be checked;

    Usefull before a long list of iteration to make sure every data
    file are ready on the disk.

    Function prompts enquiry if any file is missing in the list.
    '''

    import os
    import sys
    if sys.version_info.major==3:
        from builtins import input as input # py2 py3 compatible
    else:
        input=raw_input

    for fileii in file_list:
        if os.path.exists(fileii)==False:
            print('# <checkFiles>: File not found.',fileii)
            input("Press Enter to continue...")

    return


#----Get mask for missing data (masked or nan)----
def getMissingMask(slab):
    '''Get a bindary denoting missing (masked or nan).

    <slab>: nd array, possibly contains masked values or nans.

    Return <mask>: nd bindary, 1s for missing, 0s otherwise.
    '''
    import numpy

    nan_mask=numpy.where(numpy.isnan(slab),1,0)

    if not hasattr(slab,'mask'):
        mask_mask=numpy.zeros(slab.shape)
    else:
        if slab.mask.size==1 and slab.mask==False:
            mask_mask=numpy.zeros(slab.shape)
        else:
            mask_mask=numpy.where(slab.mask,1,0)

    mask=numpy.where(mask_mask+nan_mask>0,1,0)

    return mask


#-------Retrieve required axis from variable-------
def getAxis(axis,ref_var,verbose=True):
    dim_idx=interpretAxis(axis,ref_var)
    try:
        ax=ref_var.getAxis(dim_idx)
    except:
        raise Exception("<axis> %s not found in variable." %str(axis))

    if ax is None:
        raise Exception("<axis> %s not found in variable." %str(axis))

    return ax


def greatCircle(lat1,lon1,lat2,lon2,r=None,verbose=False):
    '''Compute the great circle distance on a sphere

    <lat1>, <lat2>: scalar float or nd-array, latitudes in degree for
                    location 1 and 2.
    <lon1>, <lon2>: scalar float or nd-array, longitudes in degree for
                    location 1 and 2.

    <r>: scalar float, spherical radius.

    Return <arc>: great circle distance on sphere.

    <arc> is computed by:

        arc = r * dsigma
        dsigma = arctan( sqrt(A) / B)
        A = (cos(<lat2>) * sin(<dlon>))^2 +
            (cos(<lat1>) * sin(<lat2>) - sin(<lat1>) * cos(<lat2>) * cos(<don>))^2
        B = sin(<lat1>) * sin(<lat2>) + cos(<lat1>) * cos(<lat2>) * cos(<dlon>)
        dlon = abs(lon1 - lon2)

    For details see wiki page:
    http://en.wikipedia.org/wiki/Great-circle_distance

    Update time: 2014-08-11 20:02:05.
    '''

    import numpy as np
    from numpy import sin, cos

    if r is None:
        r=6371000. #m

    d2r=lambda x:x*np.pi/180
    lat1,lon1,lat2,lon2=map(d2r,[lat1,lon1,lat2,lon2])
    dlon=abs(lon1-lon2)

    numerator=(cos(lat2)*sin(dlon))**2 + \
            (cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon))**2
    numerator=np.sqrt(numerator)
    denominator=sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(dlon)

    dsigma=np.arctan2(numerator,denominator)
    arc=r*dsigma

    if verbose:
        print('\n# <greatCircle>: <dsigma>:',dsigma)
        print('# <greatCircle>: <arc>:', arc)

    return arc


#----------------------Get a slab from a variable----------------------
def getSlab(var,index1=-1,index2=-2,verbose=True):
    '''Get a slab from a variable

    <var>: nd array with dimension >=2.
    <index1>,<index2>: str, indices denoting the dimensions from which a slab is to slice.

    Return <slab>: the (1st) slab from <var>.
                   E.g. <var> has dimension (12,1,241,480), getSlab(var) will
                   return the 1st time point with singleton dimension squeezed.

    Update time: 2015-07-14 19:23:42.
    '''
    import numpy

    ndim=numpy.ndim(var)
    if ndim<2:
        raise DataError('Dimension in <var> is smaller than 2.')
    if ndim==2:
        return var

    slab='dummy'
    slicer=['0',]*ndim
    slicer[index1]=':'
    slicer[index2]=':'
    string='slab=var[%s]' %','.join(slicer)
    exec(string)
    return slab


#-------------Change latitude axis to south-to-north-----------------
def increasingLatitude(slab,verbose=False):
    '''Changes a slab so that is always has latitude running from
    south to north.

    <slab>: input transientvariable. Need to have a proper latitude axis.

    Return: <slab2>, if latitude axis is reversed, or <slab> otherwise.

    If <slab> has a latitude axis, and the latitudes run from north to south, a
    copy <slab2> is made with the latitudes reversed, i.e., running from south
    to north.

    Update time: 2016-01-18 11:58:11.
    '''

    latax=getAxis('lat',slab)
    '''
    try:
        latax=slab.getLatitude()
    except:
        raise DataError('Failed to obtain latitude axis from <slab>.')

    if latax is None:
        raise DataError('Failed to obtain latitude axis from <slab>.')
    '''

    #-----Reverse latitudes if necessary------------------
    if latax[0]>latax[-1]:
        if verbose:
            print('\n# <increasingLatitude>: Reversing latitude axis.')
        slab2=slab(latitude=(latax[-1],latax[0]))
        return slab2
    else:
        if verbose:
            print('\n# <increasingLatitude>: Latitude axis correct. Not changing.')
        return slab


#----------Delta_Latitude----------------------------
def dLongitude(var,side='c',R=6371000):
    '''Return a slab of longitudinal increment (meter) delta_x.

    Args:
        var (TransientVariable): variable from which latitude axis is obtained;
        side (str): 'n': northern boundary of each latitudinal band;
                    's': southern boundary of each latitudinal band;
                    'c': central line of latitudinal band;

                 -----     'n'
                /-----\     'c'
               /_______\     's'


        R (float): radius of Earth;

    Returns:
        delta_x (TransientVariable): a 2-D slab with grid information copied from <var>.

    UPDATE: 2014-08-05 11:12:27:
        In computing <delta_x>, the longitudinal increment should be taken
        from the actual longitude axis (bounds).
        Fortunately this is not affecting any previous computations which are all
        globally.

    '''

    import numpy
    import MV2 as MV

    latax=getAxis('lat',var)
    lonax=getAxis('lon',var)

    #----------Get axes---------------------
    var=increasingLatitude(var)
    lonax=var.getLongitude()

    latax_bounds=latax.getBounds()
    lonax_bounds=lonax.getBounds()
    lon_increment=numpy.ptp(lonax_bounds,axis=1)*numpy.pi/180.

    if side=='n':
        lats=latax_bounds.max(axis=1)
    elif side=='c':
        lats=latax[:]
    elif side=='s':
        lats=latax_bounds.min(axis=1)

    lats=abs(lats)*numpy.pi/180.
    delta_x=R*numpy.cos(lats)[:,None]*lon_increment[None,:]
    delta_x=MV.where(delta_x<=1e-8,1,delta_x)

    delta_x.setAxisList((latax,lonax))


    return delta_x


def dLongitude2(lats, lons, side='c', R=6371000):
    '''Return a slab of longitudinal increment (meter) delta_x.

    Args:
        lats (ndarray): 1d array, latitude coordinates in degrees.
        lons (ndarray): 1d array, longitude coordinates in degrees.
        side (str): 'n': northern boundary of each latitudinal band;
                    's': southern boundary of each latitudinal band;
                    'c': central line of latitudinal band;

                 -----     'n'
                /-----\     'c'
               /_______\     's'


    Keyword Args:
        R (float): radius of Earth.

    Returns:
        delta_x (TransientVariable): a 2-D slab with grid information copied from <var>.

    New in v2.0.
    '''

    import numpy
    import cdms2 as cdms

    latax=cdms.createAxis(lats)
    latax.designateLatitude()
    latax.id='y'
    latax.name='latitude'
    latax.units='degree north'

    lons=numpy.sort(lons)
    lonax=cdms.createAxis(lons)
    lonax.designateLongitude()
    lonax.id='x'
    lonax.name='longitude'
    lonax.units='degree east'

    #----------Get bounds---------------------
    latax_bounds=latax.getBounds()
    lonax_bounds=lonax.getBounds()
    lon_increment=numpy.ptp(lonax_bounds,axis=1)*numpy.pi/180.

    if side=='n':
        lats=latax_bounds.max(axis=1)
    elif side=='c':
        lats=latax[:]
    elif side=='s':
        lats=latax_bounds.min(axis=1)

    lats=abs(lats)*numpy.pi/180.
    delta_x=R*numpy.cos(lats)[:,None]*lon_increment[None,:]
    delta_x=numpy.where(delta_x<=1e-8,1,delta_x)

    #delta_x.setAxisList((latax,lonax))

    return delta_x


#----------Delta_Longitude----------------------------
def dLatitude(var,R=6371000,verbose=True):
    '''Return a slab of latitudinal increment (meter) delta_y.

    Args:
        var (TransientVariable): variable from which latitude axis is abtained;

    Keyword Args:
        R (float): radius of Earth;

    Returns:
        delta_y (TransientVariable): a 2-D slab with grid information copied from\
            <var>.
    '''

    import numpy
    import MV2 as MV

    latax=getAxis('lat',var)
    lonax=getAxis('lon',var)

    #---------Get axes and bounds-------------------
    latax_bounds=latax.getBounds()
    delta_y=latax_bounds.ptp(axis=1)*numpy.pi/180.*R

    #-------Repeat array to get slab---------------
    delta_y=MV.repeat(delta_y[:,None],len(lonax),axis=1)
    delta_y.setAxisList((latax,lonax))

    return delta_y

#----------Delta_Longitude----------------------------
def dLatitude2(lats, lons, R=6371000, verbose=True):
    '''Return a slab of latitudinal increment (meter) delta_y.

    Args:
        lats (ndarray): 1d array, latitude coordinates in degrees.
        lons (ndarray): 1d array, longitude coordinates in degrees.

    Keyword Args:
        R (float): radius of Earth;

    Returns:
        delta_y (TransientVariable): a 2-D slab with grid information copied from\
            <var>.

    New in v2.0.
    '''

    import numpy
    import cdms2 as cdms

    latax=cdms.createAxis(lats)
    latax.designateLatitude()
    latax.id='y'
    latax.name='latitude'
    latax.units='degree north'

    lons=numpy.sort(lons)
    lonax=cdms.createAxis(lons)
    lonax.designateLongitude()
    lonax.id='x'
    lonax.name='longitude'
    lonax.units='degree east'

    #---------Get axes and bounds-------------------
    latax_bounds=latax.getBounds()
    delta_y=latax_bounds.ptp(axis=1)*numpy.pi/180.*R

    #-------Repeat array to get slab---------------
    delta_y=numpy.repeat(delta_y[:,None],len(lonax),axis=1)
    #delta_y.setAxisList((latax,lonax))

    return delta_y


#--------------------Get contour from a binary mask--------------------
def getBinContour(mask,lons=None,lats=None,return_largest=True):
    '''Get contour from a binary mask

    <mask>: 2d array, binary mask.
    <lons>,<lats>: 1d array, x, y coordinates for <mask>.

    Return <cont>: Nx2 array, coordinates of the contour of the largest
                   continuous region represented by 1s in <mask>.
    '''

    import matplotlib.pyplot as plt
    import numpy

    assert numpy.ndim(mask)==2, "<mask> needs to be 2D."
    if lons is not None:
        assert numpy.ndim(lons)==1, "<lons> needs to be 1D."
        assert len(lons)==mask.shape[1], "<lons> doesn't match <mask> shape."
    if lats is not None:
        assert numpy.ndim(lats)==1, "<lats> needs to be 1D."
        assert len(lats)==mask.shape[0], "<lats> doesn't match <mask> shape."

    fig,ax=plt.subplots()
    if lons is None:
        lons=numpy.arange(mask.shape[1])
    if lats is None:
        lats=numpy.arange(mask.shape[0])

    cs=ax.contourf(lons,lats,mask,[0.9,1.1]).collections
    conts=cs[0].get_paths()
    if return_largest:
        conts.sort(key=lambda x:len(x.vertices))
        #cont=conts[-1].vertices
        cont=conts[-1]
    else:
        cont=conts
    ax.cla()
    plt.close(fig)

    return cont

#-----------Find index of value in array-----------
def findIndex(x,a):
    '''Find index of value in array
    <x>: scalar, value to search.
    <a>: 1d array.

    Return <idx>: int, index in <a> that a[idx] is closest to <x>.
                  If <idx> is 0 or len(a)-1, and <x> is too far from the
                  closest value, return None.
    '''

    import numpy
    if not numpy.isscalar(x):
        raise Exception("<x> needs to be scalar.")
    if numpy.ndim(a)!=1:
        raise Exception("<a> needs to be 1d array.")

    idx=numpy.argmin(abs(x-a))
    if idx==0 and abs(a[0]-x) > abs(a[1]-a[0]):
        idx=None
        #raise Exception("<x> not in range of <a>.")
    if idx==len(a)-1 and abs(x-a[-1]) > abs(a[-1]-a[-2]):
        idx=None
        #raise Exception("<x> not in range of <a>.")
    return idx


def getBearing(lat1,lon1,lat2,lon2):
    '''Compute bearing from point 1 to point2

    <lat1>, <lat2>: scalar float or nd-array, latitudes in degree for
                    location 1 and 2.
    <lon1>, <lon2>: scalar float or nd-array, longitudes in degree for
                    location 1 and 2.

    Return <theta>: (forward) bearing in degree.
    NOTE that the bearing from P1 to P2 is in general not the same as that
    from P2 to P1.
    '''
    import numpy as np
    from numpy import sin, cos

    d2r=lambda x:x*np.pi/180
    lat1,lon1,lat2,lon2=map(d2r,[lat1,lon1,lat2,lon2])
    dlon=lon2-lon1

    theta=np.arctan2(sin(dlon)*cos(lat2),
            cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon))
    theta=theta/np.pi*180
    theta=(theta+360)%360

    return theta


def getCrossTrackDistance(lat1,lon1,lat2,lon2,lat3,lon3,r=None):
    '''Compute cross-track distance

    Args:
        lat1, lon1 (float): scalar float or nd-array, latitudes and longitudes in
                        degree, start point of the great circle.
        lat2, lon2 (float): scalar float or nd-array, latitudes and longitudes in
                        degree, end point of the great circle.
        lat3, lon3 (float): scalar float or nd-array, latitudes and longitudes in
                        degree, a point away from the great circle.

    Returns:
        dxt (float): great cicle distance between point P3 to the closest point
                  on great circle that connects P1 and P2.

                  NOTE that the sign of dxt tells which side of the 3rd point
                  P3 is on.

    See also getCrossTrackPoint(), getAlongTrackDistance().
    '''
    import numpy as np
    from numpy import sin

    if r is None:
        r=6371000.  #m

    # get angular distance between P1 and P3
    delta13=greatCircle(lat1,lon1,lat3,lon3,r=1.)
    # bearing between P1, P3
    theta13=getBearing(lat1,lon1,lat3,lon3)*np.pi/180
    # bearing between P1, P2
    theta12=getBearing(lat1,lon1,lat2,lon2)*np.pi/180

    dtheta=np.arcsin(sin(delta13)*sin(theta13-theta12))
    dxt=r*dtheta

    return dxt


def readVar(abpath_in, varid):
    '''Read in netcdf variable

    Args:
        abpath_in (str): absolute file path to nc file.
        varid (str): id of variable to read.
    Returns:
        var (TransientVariable): 4d TransientVariable.
    '''

    import cdms2 as cdms
    import numpy as np
    import MV2 as MV

    print('\n# <readVar>: Read in file:\n',abpath_in)
    fin=cdms.open(abpath_in,'r')
    var=fin(varid)

    if np.ndim(var) not in [3, 4]:
        raise Exception("Rank of <var> needs to be 3 or 4. Either (time, lat, lon) or (time, lev, lat, lon)")

    #-----------------Try get latitude-----------------
    try:
        latax=var.getLatitude()
        if latax is None:
            raise Exception("latax is None")
    except:
        raise Exception("Failed to fetch latitude axis in data.")

    #-----------------Try get longitude-----------------
    try:
        lonax=var.getLongitude()
        if lonax is None:
            raise Exception("lonax is None")
    except:
        raise Exception("Failed to fetch longitude axis in data.")

    #----------------Try get time axis----------------
    try:
        timeax=var.getTime()
        if timeax is None:
            raise Exception("timeax is None")
    except:
        raise Exception("Failed to fetch time axis in data.")

    levax=cdms.createAxis([0,])
    levax.designateLevel()
    levax.id='z'
    levax.name='level'
    levax.units=''

    ndim=np.ndim(var)
    axislist=[timeax, levax, latax, lonax]

    if ndim==3:
        var=MV.reshape(var, (len(timeax), 1, len(latax), len(lonax)))

    var.id=varid
    var.setAxisList(axislist)
    fin.close()

    return var


def getTimeAxis(times, ntime, ref_time='days since 1900-01-01'):
    '''Create a time axis

    Args:
        times (list or tuple or array or None): array of strings giving the
            time stamps. It is assumed to be in chronological order.
            If None, default to create a dummy time axis, with 6-hourly time
            step, starting from <ref_time>, with a length of <ntime>.
        ntime (int): length of the time axis. If <times> is not None, it is
            checked to insure that length of <times> equals <ntime>. If <times>
            is None, <ntime> is used to create a dummy time axis with a length
            of <ntime>.

    Keyword Args:
        ref_time (str): reference time point. If <times> is not None, used to
            create the numerical values of the resulant time axis with this
            reference time. If <times> is None, used to create a dummy time
            axis with this reference time.

    Returns:
        result (cdms.axis.TransientAxis): a cdms.axis.TransientAxis obj. If
            <times> is not None, the axis is using the provided time stamps.
            Otherwise, it is a dummy time axis, with 6-hourly time
            step, starting from <ref_time>, with a length of <ntime>.

    New in v2.0.
    '''

    import sys
    import pandas as pd
    import numpy as np
    import cdtime
    import cdms2 as cdms
    import warnings

    if times is None:
        times=np.arange(ntime)/4.
        warnings.warn("No time stamps are given. Default to 6-houly time steps since 1900-01-01.")

    if not isinstance(times, (tuple, list, pd.core.series.Series,\
            np.ndarray)):
        times=[times,]

    if len(times) != ntime:
        raise Exception("Length of <times> doesn't equal <ntime>.")

    isstr=lambda x: isinstance(x, str if sys.version_info[0]>=3 else basestring)
    datetime2Comp=lambda x: cdtime.s2c('%4d-%2d-%2d %2d:%2d:%2d'\
            %(x.year, x.month, x.day, x.hour, x.minute, x.second))

    if all([isstr(ii) for ii in times]):
        # times is all strs
        # try convert to cdtime
        try:
            dt=[pd.to_datetime(ii) for ii in times]
            result=[datetime2Comp(ii) for ii in dt]
        except:
            raise Exception("Failed to convert times to time axis. Pls check your <times> input.")
        else:
            result=[ii.torel(ref_time).value for ii in result]
            result=cdms.createAxis(result)
            result.designateTime()
            result.id='time'
            result.units=ref_time
    else:
        try:
            result=cdms.createAxis(times)
            result.designateTime()
            result.id='time'
            result.units=ref_time
        except:
            raise Exception("Failed to convert times to time axis. Pls check your <times> input.")

    return result


def breakCurveAtEdge(xs, ys, left_bound, right_bound):
    '''Segment curve coordinates at the left, right edges

    Args:
        xs (ndarray): 1d array of x- coordinates.
        ys (ndarray): 1d array of y- coordinates.
        left_bound (float): left most bound of the map domain.
        right_bound (float): right most bound of the map domain.
    Returns:
        new_xs (list): list of 1d arrays, each being a segment of the
            original input <xs>.
        new_ys (list): list of 1d arrays, each being a segment of the
            original input <ys>.

    This function segment a curve's coordinates into a number of segments
    so that when plotted in a basemap plot, a zonally cyclic curve won't
    be plotted as jumping straight lines linking the left and right bounds.
    '''

    idx=[]  # break point indices
    new_xs=[] # result list for x coordinates segments
    new_ys=[]
    for ii, xii in enumerate(xs[:-1]):
        xii2=xs[ii+1]
        dx=abs(xii2-xii)  # direct x-length from p_i to p_i+1
        dx2=min(abs(xii-left_bound) + abs(right_bound-xii2),
                abs(xii-right_bound) + abs(xii2-left_bound))
        # dx2 is the x-length if going from p_i to an edge then to p_i+1
        # if dx > dx2, going through the map edge is shorter, then need to
        # break at here
        if dx>dx2:
            idx.append(ii+1)

    if len(idx)==0:
        new_xs.append(xs)
        new_ys.append(ys)
    else:
        idx.insert(0,0)
        idx.append(len(xs))

        for i1, i2 in zip(idx[:-1], idx[1:]):
            new_xs.append(xs[i1:i2])
            new_ys.append(ys[i1:i2])

    return new_xs, new_ys

