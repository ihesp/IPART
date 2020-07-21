'''Utility functions

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-06-05 23:23:22.
'''
from __future__ import print_function
import copy
from datetime import datetime, timedelta
import numpy as np
#import pandas as pd
from netCDF4 import Dataset, date2num, num2date


class DataError(Exception):
    def __init__(self, string=None):
        if string != None:
            self.message = string
    def __str__(self):
        return self.message

def isInteger(x):
    """Check an input is integer

    Args:
        x (unknow type): input

    Returns: True if <x> is integer type, False otherwise.
    """

    import sys
    import numpy

    if sys.version_info.major==2:
        return isinstance(x, (int, long, numpy.integer))
    else:
        return isinstance(x, (int, numpy.integer))

class Selector(object):
    def __init__(self, v1, v2, axis=0):
        #if v1>v2:
            #v1, v2= v2, v1
        self.v1=v1
        self.v2=v2
        self.axis=axis

class NCVAR(object):
    def __init__(self, data, id=None, axislist=None, attributes=None):

        if id is None:
            id='unnamed'

        self.id=str(id)
        self.data=data

        if axislist is None:
            axislist=[]
            for ii in range(data.ndim):
                axisii=createAxis(ii, np.arange(data.shape[ii]).astype('f'))
                axislist.append(axisii)

        if attributes is None:
            attributes={'id': str(id),
                    'name': str(id),
                    'long_name': '',
                    'standard_name': '',
                    'title': '',
                    'units': ''}

        self._axislist=axislist
        self._attributes=attributes
        for kk, vv in attributes.items():
            setattr(self, kk, vv)

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.data.__repr__()

    def __getitem__(self, idx):
        if isInteger(idx):
            return self.sliceIndex(idx, idx+1)
        elif isinstance(idx, (list, tuple)):
            if all([isinstance(ii, slice) for ii in idx]):
                result=self
                for ii, sii in enumerate(idx):
                    if sii.start is None and sii.stop is None:
                        continue
                    result=result.sliceIndex(sii.start, sii.stop, axis=ii)
                return result

    def __call__(self, selectors):
        if selectors is None:
            return self
        if isinstance(selectors, Selector):
            selectors=(selectors,)
        if isinstance(selectors, (list, tuple)) and\
                all([isinstance(ii, Selector) for ii in selectors]):

            result=self
            for ii, sii in enumerate(selectors):
                if sii.v1 is None and sii.v2 is None:
                    continue
                result=result.sliceData(sii.v1, sii.v2, sii.axis)

            return result
        else:
            raise Exception("Input needs to be a tuple of Selectors.")



    @property
    def axislist(self):
        return copy.deepcopy(self._axislist)

    @axislist.setter
    def axislist(self, newlist):
        self._axislist=newlist

    @property
    def attributes(self):
        return copy.deepcopy(self._attributes)

    @attributes.setter
    def attributes(self, att):
        self._attributes=att

    @property
    def dims(self):
        return [aa.id for aa in self._axislist]

    @property
    def dtype(self):
        return str(self.data.dtype)

    @property
    def ndim(self):
        return self.data.ndim

    @property
    def shape(self):
        return self.data.shape

    def __len__(self):
        return self.data.shape[0]

    def info(self):
        result=['### Description of slab ###',
                '  id: %s' %self.id,
                '  shape: %s' %str(self.shape),
                '  filename: %s' %str(self.attributes.get('filename')),
                '  missing_value: %s' %str(self.attributes.get('missing_value')),
                '  comments: %s' %str(self.attributes.get('comments')),
                '  grid_name: %s' %str(self.attributes.get('grid_name')),
                '  grid_type: %s' %str(self.attributes.get('grid_type')),
                '  long_name: %s' %str(self.attributes.get('long_name')),
                '  units: %s' %str(self.attributes.get('units')),
                '  standard_name: %s' %str(self.attributes.get('standard_name')),
                '  Order: %s' %str([aa.id for aa in self.axislist]),]
        axes=[' ']
        for ii, axisii in enumerate(self.axislist):
            attrii=axisii.attributes
            keys=list(attrii.keys())
            axes.append('** Dimension %d **' %(ii+1))
            keys.sort()
            for kk in keys:
                vv=str(attrii[kk])
                lkk='   %s: %s' %(kk, vv)
                axes.append(lkk)
            axes.append('   length: %d' %len(axisii.data))
            axes.append('   1st: %s'  %str(axisii.data[0]))
            axes.append('   lst: %s'  %str(axisii.data[-1]))

        result.extend(axes)
        result.append('### End of description ###')
        result=[ll+'\n' for ll in result]
        result=''.join(result)
        print(result)

        return

    def sliceIndex(self,idx1,idx2,axis=0,squeeze=True):

        #if idx1>idx2:
            #idx1,idx2=idx2,idx1

        axisobj=self._axislist[axis]
        axisdata=axisobj.data
        if idx1 is None:
            idx1=0
        if idx2 is None:
            idx2=len(axisdata)

        newaxisdata=axisdata[idx1:idx2]

        if len(newaxisdata)==0:
            raise Exception("No data found in [idx1,idx2).")

        slicer=[slice(None),]*self.ndim
        slicer[axis]=slice(idx1,idx2)

        newdata=self.data[slicer]

        if idx2-idx1==1 and squeeze:
            newdata=np.squeeze(newdata, axis=axis)

        #--------------Create a new data obj--------------
        axislist=self.axislist
        if idx2-idx1==1 and squeeze:
            axislist.pop(axis)
        else:
            newaxis=NCVAR(newaxisdata, axisobj.id, [], axisobj.attributes)
            axislist[axis]=newaxis
        result=NCVAR(newdata, self.id, axislist, self.attributes)
        for ii in axislist:
            setattr(result, ii.id, ii)

        return result


    def sliceData(self, v1, v2, axis=0):
        if v1>v2:
            v1,v2=v2,v1

        axisobj=self._axislist[axis]
        axisdata=axisobj.data
        idx=np.ma.where((axisdata>=v1) & (axisdata<v2))[0]
        slicer=[slice(None),]*self.ndim

        if len(idx)>0:
            slicer[axis]=idx
        else:
            raise Exception("No data found in [v1,v2].")

        newdata=self.data[slicer]   # if var is np.ndarray

        #--------------Create a new data obj--------------
        newaxisdata=axisdata[idx]
        newaxis=NCVAR(newaxisdata, axisobj.id, [], axisobj.attributes)
        axislist=self.axislist
        axislist[axis]=newaxis
        result=NCVAR(newdata, self.id, axislist, self.attributes)
        for ii in axislist:
            #exec('result.%s=ii' %str(ii.id))
            setattr(result, ii.id, ii)
        return result

    def squeeze(self):
        axislist=[]
        for ii in self.axislist:
            if len(ii.data)!=1:
                axislist.append(ii)
        self.data=np.squeeze(self.data)
        self.axislist=axislist
        result=NCVAR(self.data, self.id, self.axislist, self.attributes)
        return result

    def getAxis(self, idx):
        if idx not in range(self.ndim):
            raise Exception("<idx> not in data shape.")
        return self.axislist[idx]

    def getLatitude(self):

        for axisii in self.axislist:
            if axisii.id.lower() in ['y', 'lat', 'latitude']:
                return axisii.data
        return

    def getLongitude(self):

        for axisii in self.axislist:
            if axisii.id.lower() in ['x', 'lon', 'longitude']:
                return axisii.data
        return

    def getTime(self):

        for axisii in self.axislist:
            if axisii.id.lower() in ['t', 'time']:
                return axisii.data
        return

    def getLevel(self):

        for axisii in self.axislist:
            if axisii.id.lower() in ['z', 'level']:
                return axisii.data
        return

    def shiftLon(self, dx):

        lonidx=interpretAxis('longitude', self)
        if lonidx<0:
            raise Exception("Longitude axis not found in var.")

        lonax=self.axislist[lonidx]
        lons=lonax.data
        dxidx=np.argmin(abs(lons-dx))
        lons=np.roll(lons, -dxidx)
        lonax.data=lons
        self.data=np.roll(self.data, -dxidx, axis=lonidx)

        return self

def squeezeTo3D(vv):

    if np.ndim(vv) not in [3, 4]:
        raise Exception("Input <ivt>, <ivtrec>, <ivtano>, <qu> and <qv> should be 3D or 4D.")
    if np.ndim(vv)==4 and vv.shape[0]!=1:
        vv=np.squeeze(vv)
    elif np.ndim(vv)==4 and vv.shape[0]==1:
        vv=vv[:,0,:,:]
    elif np.ndim(vv)==3 and vv.shape[0]==1:
        pass

    return vv

def createAxis(id, data, attributes=None):

    if attributes is None:
        attributes={'id': str(id),
                'name': str(id),
                'long_name': 'axis_%s' %str(id),
                'units': ''}

    result=NCVAR(data, id, [], attributes)
    return result

def getBounds(axisdata, width=1.):
    # only for axis
    if np.ndim(axisdata)>1:
        raise Exception("Only work for axis data.")
    if len(axisdata)==1:
        delta=width/2.
        return np.array([[axisdata[0]-delta, axisdata[0]+delta]])

    bound1=0.5*(axisdata[1:]+axisdata[:-1])
    b_lower=np.hstack([axisdata[0]*2.-bound1[0], bound1])
    b_upper=np.hstack([bound1, axisdata[-1]*2.-bound1[-1]])
    return np.vstack([b_lower,b_upper]).T

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

    if percents is None:
        percents=np.array([0.001,0.005,0.01,0.025,0.05,0.1])
    percents=np.array(percents)
    if percents.ndim!=1:
        raise Exception("<percents> needs to be a 1D array.")

    #-------Remove nans and masked values--------
    mask=getMissingMask(slab)
    slab=np.array(slab)
    slab=slab[np.where(mask==False)]

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

    for kk, vv in source_object.attributes.items():
        dictionary[kk]=vv
        if verbose:
            print('\n# <attribute_obj2dict>: %s: %s' %(kk, vv))

    return dictionary


#-------------Copy attributes from dict to target object----------
def attribute_dict2obj(dictionary,target_object,verbose=False):
    '''Copies attributes from dictionary to target object.

    <dictionary>: dict, contains attributes to copy.
    <target_object>: obj, attributes are copied to.

    Return <target_object>: target object with new attributes.

    Update time: 2016-01-18 11:31:25.
    '''

    for kk, vv in dictionary.items():
        setattr(target_object, kk, vv)
        target_object.attributes[kk]=vv
        if verbose:
            print('\n# <attribute_dict2obj>: Copy attribute: %s = %s' %(kk,vv))

    return target_object


def addExtraAxis(slab, newaxis, axis=0, verbose=False):
    """Adds an extra axis to a data slab.

    <slab>: variable to which the axis is to insert.
    <newaxis>: axis object, could be of any length. If None, create a dummy
               singleton axis.
    <axis>: index of axis to be inserted, e.g. 0 if <newaxis> is inserted
            as the 1st dimension.

    Return: <slab2>.
    Update time: 2013-10-09 12:34:32.
    """

    if newaxis is None:
        newaxis=NCVAR(np.array([0.]), 'newaxis', [], {'id': 'newaxis',
            'name': 'newaxis', 'units': ''})

    # add new axis to axis list of input <slab>
    axislist=slab.axislist
    axislist.insert(axis, newaxis)

    #----------------Reshape----------------
    shape=list(slab.shape)
    shape.insert(axis,len(newaxis))
    slab2=np.reshape(slab.data, shape)

    #------------Create variable------------
    att_dict=slab.attributes
    slab2=NCVAR(slab2, id=slab.id, axislist=axislist, attributes=att_dict)

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

    attdict=getattr(var1, 'attributes', None)
    axis1=var1.getAxis(axis)
    axis2=var2.getAxis(axis)
    newaxis=np.r_[axis1.data, axis2.data]
    newaxis=NCVAR(newaxis, axis1.id, [], axis1.attributes)

    result=np.ma.concatenate((var1.data, var2.data), axis=axis)
    axislist=var1.axislist
    axislist[axis]=newaxis
    result=NCVAR(result, id=var1.id, axislist=axislist, attributes=attdict)

    return result


def concatenate(var_list, axis=0, verbose=False):

    if not isinstance(var_list, (list, tuple)):
        raise Exception("<var_list> needs to be a tuple of NCVAR variables.")

    if not all([isinstance(vii, NCVAR) for vii in var_list]):
        raise Exception("<var_list> needs to be a tuple of NCVAR variables.")

    if len(var_list)==0:
        raise Exception("Needs at least 1 variable.")
    if len(var_list)==1:
        return var_list[0]

    result=var_list[0]
    resultdata=result.data
    newaxis=result.getAxis(axis)
    newaxisdata=newaxis.data
    axislist=result.axislist
    attdict=getattr(result, 'attributes', None)

    for vii in var_list[1:]:
        axis2=vii.getAxis(axis)
        newaxisdata=np.r_[newaxisdata, axis2.data]
        resultdata=np.ma.concatenate([resultdata, vii.data], axis=axis)

    newaxis=NCVAR(newaxisdata, newaxis.id, [], newaxis.attributes)
    axislist[axis]=newaxis
    result=NCVAR(resultdata, id=result.id, axislist=axislist, attributes=attdict)

    return result

#------Interpret and convert an axis id to index----------
def interpretAxis(axis, ref_var, verbose=True):
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

    if isinstance(axis, (int, np.integer)):
        return axis

    # interpret string dimension
    elif isinstance(axis, str if sys.version_info[0]>=3 else basestring):
        axis=axis.lower()
        dim_idx = -1

        for ii, axisii in enumerate(ref_var.axislist):
            idii=axisii.id

            if axis in ['time', 'tim', 't']:
                if idii in ['time', 'tim', 't']:
                    dim_idx = ii
            elif axis in ['level', 'lev','z']:
                if idii in ['level', 'lev','z']:
                    dim_idx = ii
            elif axis in ['latitude', 'lat','y']:
                if idii in ['latitude', 'lat','y']:
                    dim_idx = ii
            elif axis in ['longitude', 'long', 'lon','x']:
                if idii in ['longitude', 'long', 'lon','x']:
                    dim_idx = ii

        if dim_idx==-1:
            raise Exception("Required dimension not in <var>.")
        return dim_idx
    else:
        raise Exception("<axis> type not recognized.")

def getAttributes(var):
    attr={}
    for kk in var.ncattrs():
        attr[kk]=var.getncattr(kk)
    return attr

def readNC(abpath_in, varid):

    fin=Dataset(abpath_in, 'r')
    var=fin.variables[varid]

    dims=var.dimensions
    axislist=[]
    for dd in dims:
        ncaxis=fin.variables[dd]
        if dd=='time':
            try:
                axisii=num2date(ncaxis[:], ncaxis.units,
                        only_use_cftime_datetimes=False,
                        only_use_python_datetimes=True)
            except:
                axisii=num2date(ncaxis[:], ncaxis.units,
                        only_use_cftime_datetimes=False)
            #axisii=pd.to_datetime(axisii)
        else:
            axisii=ncaxis[:]
        axisattr=getAttributes(ncaxis)
        axisattr['isunlimited']=fin.dimensions[dd].isunlimited()
        axisii=NCVAR(axisii, dd, [], axisattr)
        axislist.append(axisii)

    attributes=getAttributes(var)
    data=var[:]

    ncvar=NCVAR(data, varid, axislist, attributes)
    for ii in axislist:
        #exec('ncvar.%s=ii' %str(ii.id))
        setattr(ncvar, ii.id, ii)

    return ncvar


def saveNC(abpath_out, var, mode='w', dtype=None):

    with Dataset(abpath_out, mode) as fout:

        saveNCDims(fout, var.axislist)
        _saveNCVAR(fout, var, dtype)


def _saveNCVAR(fout, var, dtype=None):

    #-----------------Create variable-----------------
    if dtype is None:
        dtype=np.float32

    if var.id not in fout.variables.keys():
        varout=fout.createVariable(var.id, dtype, var.dims, zlib=True)
        #varout.set_collective(True)
        #----------------Create attributes----------------
        for kk, vv in var.attributes.items():
            try:
                varout.setncattr(kk, vv)
            except:
                pass
        varout[:]=var.data
    else:
        varout=fout.variables[var.id]
        '''
        if np.all(varout[-1].mask):
            varout[-1]=var.data
        else:
            newdata=np.concatenate([varout[:], var.data])
            varout[:]=newdata
        timeax=fout.variables['time']
        if timeax[-1].mask:
            timeax[-1]=var.getTime()
        '''
        timeax=fout.variables['time']
        tlen=len(timeax)
        tnew=var.getTime()
        t0=tnew[0]
        tidx=np.where(timeax==t0)[0]
        if len(tidx)>0:
            tidx=tidx[0]
        else:
            tidx=tlen
        timeax[tidx:tidx+len(tnew)]=tnew
        varout[tidx:tidx+len(tnew)]=var.data


    return


def saveNCDims(fout, axislist):

    #----------------Create dimensions----------------
    for aa in axislist:

        if aa.id in fout.dimensions.keys():
            continue

        if aa.id in ['t', 'time', 'Time', 'T']:
            lenii=None

            if all([isinstance(ii, datetime) for ii in aa.data]):
                axisdata=date2num(aa.data, aa.units)
            else:
                axisdata=aa.data
        else:
            lenii=len(aa.data)
            axisdata=aa.data

        if hasattr(aa, 'isunlimited') and aa.isunlimited:
            lenii=None

        fout.createDimension(aa.id, lenii)
        axisvar=fout.createVariable(aa.id, np.float32, (aa.id,), zlib=True)
        axisvar[:]=axisdata

        for kk,vv in aa.attributes.items():
            if kk!='isunlimited':
                axisvar.setncattr(kk, vv)



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

    if isinstance(slab, NCVAR):
        slab=slab.data

    nan_mask=np.where(np.isnan(slab),1,0)

    if not hasattr(slab,'mask'):
        mask_mask=np.zeros(slab.shape)
    else:
        if slab.mask.size==1 and slab.mask==False:
            mask_mask=np.zeros(slab.shape)
        else:
            mask_mask=np.where(slab.mask,1,0)

    mask=np.where(mask_mask+nan_mask>0,1,0)

    return mask


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

    ndim=np.ndim(var)
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
def increasingLatitude(slab, axis, lats=None, verbose=False):
    '''Changes a slab so that is always has latitude running from
    south to north.

    <slab>: input transientvariable. Need to have a proper latitude axis.

    Return: <slab2>, if latitude axis is reversed, or <slab> otherwise.

    If <slab> has a latitude axis, and the latitudes run from north to south, a
    copy <slab2> is made with the latitudes reversed, i.e., running from south
    to north.

    Update time: 2016-01-18 11:58:11.
    '''

    if isinstance(slab, NCVAR):
        lats_data=slab.getLatitude()
        if lats_data is None:
            if lats is None:
                raise Exception("No latitude axis info found.")
            lats_data=lats
        data=slab.data
        isnc=True
    elif isinstance(slab, np.ndarray):
        if lats is None:
            raise Exception("Need to provide latitude axis data.")
        lats_data=lats
        data=slab
        isnc=False


    #-----Reverse latitudes if necessary------------------
    if lats_data[0]>lats_data[-1]:
        if verbose:
            print('\n# <increasingLatitude>: Reversing latitude axis.')

        data=np.flip(data, axis=axis)

        if not isnc:
            return data
        else:
            lats_data=lats_data[::-1]
            axislist=slab.axislist
            latax=createAxis(axislist[axis].id, lats_data,
                    axislist[axis].attributes)
            axislist[axis]=latax
            slab2=NCVAR(data, slab.id, axislist, slab.attributes)
        return slab2
    else:
        if verbose:
            print('\n# <increasingLatitude>: Latitude axis correct. Not changing.')
        return slab

def increasingLatitude2(slab, axis, lats, verbose=False):
    '''Changes a slab so that is always has latitude running from
    south to north.

    <slab>: input transientvariable. Need to have a proper latitude axis.

    Return: <slab2>, if latitude axis is reversed, or <slab> otherwise.

    If <slab> has a latitude axis, and the latitudes run from north to south, a
    copy <slab2> is made with the latitudes reversed, i.e., running from south
    to north.

    Update time: 2016-01-18 11:58:11.
    '''

    #-----Reverse latitudes if necessary------------------
    if lats[0]>lats[-1]:
        if verbose:
            print('\n# <increasingLatitude>: Reversing latitude axis.')

        slab=np.flip(slab, axis=axis)
        lats=lats[::-1]
        return slab, lats
    else:
        if verbose:
            print('\n# <increasingLatitude>: Latitude axis correct. Not changing.')
        return slab, lats




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

    lons=np.sort(lons)

    #----------Get bounds---------------------
    latax_bounds=getBounds(lats)
    lonax_bounds=getBounds(lons)
    lon_increment=np.ptp(lonax_bounds,axis=1)*np.pi/180.

    if side=='n':
        lats=latax_bounds.max(axis=1)
    elif side=='c':
        pass
    elif side=='s':
        lats=latax_bounds.min(axis=1)

    lats=abs(lats)*np.pi/180.
    delta_x=R*np.cos(lats)[:,None]*lon_increment[None,:]
    delta_x=np.where(delta_x<=1e-8,1,delta_x)

    return delta_x


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

    lons=np.sort(lons)

    #---------Get axes and bounds-------------------
    latax_bounds=getBounds(lats)
    delta_y=latax_bounds.ptp(axis=1)*np.pi/180.*R

    #-------Repeat array to get slab---------------
    delta_y=np.repeat(delta_y[:,None],len(lons),axis=1)

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

    assert np.ndim(mask)==2, "<mask> needs to be 2D."
    if lons is not None:
        assert np.ndim(lons)==1, "<lons> needs to be 1D."
        assert len(lons)==mask.shape[1], "<lons> doesn't match <mask> shape."
    if lats is not None:
        assert np.ndim(lats)==1, "<lats> needs to be 1D."
        assert len(lats)==mask.shape[0], "<lats> doesn't match <mask> shape."

    fig,ax=plt.subplots()
    if lons is None:
        lons=np.arange(mask.shape[1])
    if lats is None:
        lats=np.arange(mask.shape[0])

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

    if not np.isscalar(x):
        raise Exception("<x> needs to be scalar.")
    if np.ndim(a)!=1:
        raise Exception("<a> needs to be 1d array.")

    idx=np.argmin(abs(x-a))
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
        times (list or tuple or array): array of datetime objs, or strings,
            giving the time stamps. It is assumed to be in chronological order.
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
        result (list): a list of datetime objs. If <times> is not None,
            the datetime objs are using the provided time stamps.
            Otherwise, it is a new time series, with 6-hourly time
            step, starting from <ref_time>, with a length of <ntime>.

    New in v2.0.
    '''

    import sys
    import numpy as np
    import warnings

    #---If no time stamps given, create a dummy one---
    if times is None:
        warnings.warn("No time stamps are given. Default to 6-houly time steps since 1900-01-01.")
        t1=datetime.strptime(ref_time.split('since')[1].strip(),
                '%Y-%m-%d')
        dt=timedelta(hours=6)
        result=[t1+dt*ii for ii in range(ntime)]

        return result

    if not isinstance(times, (tuple, list, np.ndarray)):
        times=[times,]

    if len(times) != ntime:
        raise Exception("Length of <times> doesn't equal <ntime>.")

    isstr=lambda x: isinstance(x, str if sys.version_info[0]>=3 else basestring)

    #------If times are already datetime, return------
    if all([isinstance(ii, datetime) for ii in times]):
        result=times

    #---If times are in strings, convert to datetime---
    if all([isstr(ii) for ii in times]):
        # times is all strs # try convert to datetime
        try:
            result=[datetime.strptime(ii, '%Y-%m-%d %H:%M') for ii in times]
        except:
            raise Exception("Failed to convert times to time axis. Pls check your <times> input is in the '1989-06-04 12:00' format.")

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


#------ Get binary grids inside a contour ------------
def getGridsInContour(contour,x,y):
    '''Get binary grids inside a contour
    <contour>: Nx2 ndarray, (x,y) coordinates of a contour.
               Or a matplotlib Path obj. If the latter, function can
               remove holes if exists in the contour.
    <x>,<y>: 1d array, x and y coordinates of the grid.

    Return <validcoords>: 2d array of shape (len(y), len(x)), binary slab
                          with 1s inside of <contour> and 0s elsewhere.
    '''
    from matplotlib.path import Path

    #---------------Create a dummy path---------------
    dummy=Path([[0,0],[1,1]])

    #-----------Check type of input contour-----------
    if type(contour)==type(dummy):
        ispath=True
    else:
        ispath=False
        assert np.ndim(contour)==2, "<contour> needs to be 2D."
        assert contour.shape[1]==2, "<contour> needs to be Nx2."

    assert np.ndim(x)==1, "<x> needs to be 1D."
    assert np.ndim(y)==1, "<y> needs to be 1D."

    if ispath:
        cxy=contour.vertices
    else:
        cxy=contour

    #--------Get a bounding box around contour--------
    xmin,ymin=np.min(cxy,axis=0)
    xmax,ymax=np.max(cxy,axis=0)

    xmin_idx=findIndex(xmin,x) or 0
    ymin_idx=findIndex(ymin,y) or 0
    xmax_idx=findIndex(xmax,x) or len(x)-1
    ymax_idx=findIndex(ymax,y) or len(y)-1

    xmin_idx=max(0,xmin_idx-3) # enlarge a bit
    ymin_idx=max(0,ymin_idx-3) # enlarge a bit
    xmax_idx=min(len(x)-1,xmax_idx+3) # enlarge a bit
    ymax_idx=min(len(y)-1,ymax_idx+3) # enlarge a bit

    X2,Y2=np.meshgrid(x[xmin_idx:xmax_idx],y[ymin_idx:ymax_idx])
    coords=np.array(zip(X2.flat, Y2.flat))

    #------Create a path from vertices------
    if ispath:
        path=contour
    else:
        path=Path(contour)
    validcoords=path.contains_points(coords)
    validcoords=np.reshape(validcoords,X2.shape).astype('int')

    #-------------------Subtract holes-------------------
    if ispath:
        segs=contour.to_polygons()
        if len(segs)>1:
            areas=[polygonArea(sii[:,0],sii[:,1]) for sii in segs]
            winner_id=np.argmax(areas)

            for ii in range(len(segs)):
                if ii!=winner_id:
                    contii=Path(segs[ii])
                    holeii=contii.contains_points(coords)
                    holeii=np.reshape(holeii,X2.shape).astype('int')
                    validcoords=validcoords-holeii

    #------------------Paste box back------------------
    mask=np.zeros([len(y),len(x)])
    mask[ymin_idx:ymax_idx,xmin_idx:xmax_idx]=validcoords

    return mask

def polygonArea(x,y):

    '''
    def isClosed(xs,ys):
        if np.alltrue([np.allclose(xs[0],xs[-1]),\
            np.allclose(ys[0],ys[-1]),xs.ptp(),ys.ptp()]):
            return True
        else:
            return False

    if not isClosed(x,y):
        # here is a minor issue: isclosed() on lat/lon can be closed,
        # but after projection, unclosed. Happens to spurious small
        # contours usually a triangle. just return 0.
        return 0
    area=np.sum(y[:-1]*np.diff(x)-x[:-1]*np.diff(y))
    return np.abs(0.5*area)
    '''
    return np.abs(signedArea(x, y))

def signedArea(x, y):

    x=np.asarray(x)
    y=np.asarray(y)
    def isClosed(xs,ys):
        if np.alltrue([np.allclose(xs[0],xs[-1]),\
            np.allclose(ys[0],ys[-1]),xs.ptp(),ys.ptp()]):
            return True
        else:
            return False

    if not isClosed(x, y):
        x=np.r_[x, x[0]]
        y=np.r_[y, y[0]]

    area = np.sum(-y[:-1] * np.diff(x) + x[:-1] * np.diff(y)) * 0.5

    return area

def averager(slab, axis=-1, weights='equal', coords=None, keepdims=True,
        verbose=True):

    if not isinstance(axis, (tuple, list)):
        axis=(axis,)

    if any([ii > np.ndim(slab)-1 for ii in axis]):
        raise Exception("<axis> exceeds slab's dimension.")

    if weights=='equal':
        result=np.ma.mean(slab, axis=axis, keepdims=True)

    elif weights=='generate':
        if coords is None:
            result=np.ma.mean(slab, axis=axis, keepdims=True)
        else:
            if not isinstance(coords, (tuple, list)):
                coords=[coords,]

            if len(axis) != len(coords):
                raise Exception("Length of <axis> doesn't equal that of <coords>.")

            result=slab
            for idxii, cii in zip(axis, coords):
                boundii=getBounds(cii)
                wii=np.abs(np.diff(boundii, axis=1))
                wii=arrayGrow(wii, result.shape, axis=idxii)
                wii=np.ma.array(wii)
                wii.mask=result.mask

                result=np.ma.sum(result*wii, axis=idxii, keepdims=True)/\
                        np.ma.sum(wii, axis=idxii, keepdims=True)

    if not keepdims:
        result=np.squeeze(result, axis=axis)

    return result


def areaAverage(slab, axis=(-2, -1), weights='equal', lats=None, lons=None,
        keepdims=True, verbose=True):

    if not isinstance(axis, (tuple, list)):
        axis=(axis,)

    if any([ii > np.ndim(slab)-1 for ii in axis]):
        raise Exception("<axis> exceeds slab's dimension.")

    if weights=='equal':
        result=np.ma.mean(slab, axis=axis, keepdims=True)

    elif weights=='generate':
        if lats is None and lons is None:
            result=np.ma.mean(slab, axis=axis, keepdims=True)
        else:
            if lats is None or lons is None:
                raise Exception("Need to provide lat and lon coordinates.")

            dlat=dLatitude2(lats, lons)
            dlon=dLongitude2(lats, lons)
            area=dlat*dlon

            area=slabGrow(area, slab.shape)

            if hasattr(slab, 'mask'):
                area.mask=slab.mask

            num=np.ma.sum(slab*area, axis=axis, keepdims=True)
            den=np.ma.sum(area, axis=axis, keepdims=True)
            result=num/den

            '''
            result=slab
            for idxii, cii in zip(axis, [dlat, dlon]):
                wii=arrayGrow(cii, result.shape, axis=idxii)
                wii=np.ma.array(wii)
                if hasattr(slab, 'mask'):
                    wii.mask=slab.mask

                result=np.ma.sum(result*wii, axis=idxii, keepdims=True)/\
                        np.ma.sum(wii, axis=idxii, keepdims=True)
            '''

    if not keepdims:
        result=np.squeeze(result, axis=axis)

    return result


def areaStd(slab, axis=(-2, -1), weights='equal', lats=None, lons=None,
        keepdims=True, verbose=True):

    if not isinstance(axis, (tuple, list)):
        axis=(axis,)

    if any([ii > np.ndim(slab)-1 for ii in axis]):
        raise Exception("<axis> exceeds slab's dimension.")

    if weights=='equal':
        result=np.ma.std(slab, axis=axis, keepdims=True)

    elif weights=='generate':
        if lats is None and lons is None:
            result=np.ma.std(slab, axis=axis, keepdims=True)
        else:
            if lats is None or lons is None:
                raise Exception("Need to provide lat and lon coordinates.")

            dlat=dLatitude2(lats, lons)
            dlon=dLongitude2(lats, lons)
            area=dlat*dlon

            area=slabGrow(area, slab.shape)

            if hasattr(slab, 'mask'):
                area.mask=slab.mask

            num=np.ma.sum(slab*area, axis=axis, keepdims=True)
            den=np.ma.sum(area, axis=axis, keepdims=True)
            mean=num/den
            ano=slab-mean
            result=np.ma.sum(ano*ano*area, axis=axis, keepdims=True)/den
            result=np.ma.sqrt(result)

    if not keepdims:
        result=np.squeeze(result, axis=axis)

    return result

#---------------Grow a 2-D slab to 3-D or 4-D--------------------
def slabGrow(slab, shape, verbose=False):
    '''Grow an 2-D slab to 3-D or 4-D slab, according to a reference variable.

    <slab> An MV 2-D slab to be repeated;
    <ref_var> Multi-dimentional (3-D or 4-D) reference variable;

    Return <slab>.
    Usually used to extend one coeffient to match the size of a variable.
    e.g. <slab> is a 2-D land-sea-mask, the function then grows the 2-D slab\
            to match dimention of reference variable <ref_var>, either \
            (t,z,y,x) or (t,y,x).
    '''

    slabshape=slab.shape

    match=False
    for ii, jj in zip(shape[:-1], shape[1:]):
        if ii==slabshape[0] and jj==slabshape[1]:
            match=True
            break

    if not match:
        raise Exception("Doesn't find matching slab in shape.")

    order=range(len(shape))

    order_remain=order[0:-2]
    order_remain.reverse()

    for indexii in order_remain:
        slab=np.ma.reshape(slab, (1,)+slab.shape)
        slab=np.ma.repeat(slab, shape[indexii], axis=0)

    slab=np.ma.array(slab)

    if verbose:
        print('# <slabGrow>: New slab shape:', slab.shape)

    return slab


#------------Grow an 1-D array to 2-D, 3-D or 4-D-------------------
def arrayGrow(array,shape,axis=None,verbose=True):
    '''Grow an 1-D array to a 2-D, 3-D or 4-D slab, according to a\
            reference variable.

    <array>: 1-D array to be repeated;
    <ref_var>: nd (2-D, 3-D or 4-D) reference variable;
    <axis>: an optional argument to specify the axis of the array using\
            indices (0,1,2,or 3).
            usefull when 2 axes in ref_var have same length.\

    Return <newarray>.

    Usually used to extend one coefficient to match the size of a variable.
    e.g. <array> is longitudinal increment calculated using latitude values,\
            the function then grows the 1-D array to match dimentions of reference
            variable <ref_var>, which could be (t,z,y,x) or (t,y,x).

    NOTE: difference to arrayGrowOld(): <ref_var> can be an numpy.ndarray rather
          than transientVariable.

    Update time: 2016-01-18 16:30:43.
    '''

    import numpy

    #-----------------Check dimension-----------------
    array=np.atleast_1d(np.squeeze(array))
    if numpy.ndim(array)>1:
        raise Exception("<array> needs to be 1d array.")

    length=len(array)

    if axis is None:
        try:
            index=shape.index(length)
        except:
            raise Exception('Shape unmatch along the given axis')
    else:
        index=axis
        if shape[index]!=length:
            raise Exception('Shape unmatch along the given axis')

    #-----------Switch this axis to the end of shape list------
    ndim=len(shape)
    order=range(ndim)

    if index!=len(shape)-1:
        order[index]=order[-1]
        order[-1]=index

    tmpshape=numpy.array(shape)[order]

    #---------------Expand by repeating---------------
    array=numpy.array(array)
    newarray=array[numpy.newaxis,:]
    newarray=numpy.repeat(newarray, numpy.prod(shape)/length, axis=0)
    newarray=numpy.reshape(newarray,tmpshape)
    newarray=numpy.transpose(newarray,order)

    return newarray


