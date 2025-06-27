'''Utility functions

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-07-22 09:27:36.
'''
from __future__ import print_function
import sys, os
import copy
from datetime import datetime, timedelta
import cftime
import numpy as np
import pandas as pd
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
    Returns:
        True if <x> is integer type, False otherwise.
    """
    if sys.version_info.major==2:
        return isinstance(x, (int, long, np.integer))
    else:
        return isinstance(x, (int, np.integer))

class Selector(object):
    def __init__(self, v1, v2, axis=0):
        '''Similar to slice to slice NCVAR variable'''
        self.v1=v1
        self.v2=v2
        self.axis=axis


class NCVAR(object):
    def __init__(self, data, id=None, axislist=None, attributes=None):
        '''An object to store netcdf data and metadata

        Args:
            data (ndarray): data array.
        Keyword Args:
            id (str or None): id of data. If None, use a default "unnamed".
            axislist (tuple or list or None): a tuple/list of NCVAR objs
                acting as axes of the data. The length should equal the rank
                of <data>. If None, create one axis for each dimension using
                integer indices.
            attributes (dict or None): variable attribute dictionary. If None,
                create a new dict filled with some basic attributes.

        This is a rudimentary attempt to mimic the API of CDAT package in
        handling netcdf data. The netCDF4 API, imo, goes against the philosiphy
        of of netcdf format. The data and metadata, after read in using
        netCDF4, are saved to different variables, rather than bandled
        together. One has to maintain these metadata himself, and it gets
        bloody annoying quickly when things get complicated.

        This class will try to save data and metadata together, and provide
        some utility methods for interacting via metadata.

        Examples:
            * to create a variable:
                ```
                var = NCVAR(data_array, 'sst', (times, lats, lons),
                     {'name': 'sst,
                      'long_name': 'sea surface temperature',
                      'standard_name': 'sea_surface_temperature',
                      'units': 'K'})
                ```
            * to get the np data array:
                ```var.data```
            * to query the shape:
                ```
                var.shape, var.ndim, len(var)
                ```
            * to get the time axis:
                ```var.getTime()```
            * to get the latitude axis:
                ```var.getLatitude()```
            * to slice by latitude:
                ```slab = var.sliceData(10, 60, axis=1)```
            * to get the 1st time slice:
                ```slab = var[0]```
            * to slice by indices:
                ```slab = var.sliceIndex(1, 10, axis=2)```
                Alternatively, can slice by indices using a Selector:
            * to create a selector:
                ```sel = Selector(1, 10, axis=2)```
            * to slice using the selector:
                ```slab = var(sel)```
                This will slice the data as ```data[:,:,1:10]```, the 2nd
                axis will be sliced accordingly.

        The native CDAT API is more intuitive.

        Unlike CDAT, the obj can not be used directly as ndarray. To do
        computations using data, one has to access the "data" attribute:
            ```
            res = np.sin(var.data)
            ```

        It can be confusing when something is an NCVAR or an ndarray. To help
        distinguish, will use a naming convention that appends a "NV" after
        a variable name to indicate that this is an NCVAR: ivtNV, sstNV, etc..
        Their correponding ndarray are named ivt, sst, etc..
        '''

        if id is None:
            id='unnamed'

        self.id=str(id)
        self.data=data

        # create axes using indices
        if axislist is None:
            axislist=[]
            for ii in range(data.ndim):
                axisii=createAxis(ii, np.arange(data.shape[ii]).astype('f'))
                axislist.append(axisii)

        # fill some basic attributes
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
        '''Indexing method

        Args:
            idx (int or tuple/list): if int, indexing the 1st dimension by
                the given index.
                If tuple/list, a tuple/list of slice objs to slice the data.
        Returns:
            result (NCVAR): sliced result.
        '''
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
        '''Indexing using selectors

        Args:
            selectors (Selector or tuple/list or None): if None, return
                self.
                If Selector or tuple/list of Selectors, slice the data using
                ranges defined in Selectors.
        Returns:
            result (NCVAR): sliced result.
        '''
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
        '''Print some summary info of the variable'''

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
        '''Slice the data by start, end indices along a given dimension

        Args:
            idx1, idx2 (int or None): start, end indices. If idx1 is None,
                changed to 0. If idx2 is None, changed to the length of the
                dimension specified by <axis>.
            axis (int): axis to do slicing.
        Keyword Args:
            squeeze (bool): if True, will squeeze singleton dimensions.
        Returns:
            result (NCVAR): sliced data.
        '''

        #--------------------Slice axis--------------------
        axisobj=self._axislist[axis]
        axisdata=axisobj.data
        if idx1 is None:
            idx1=0
        if idx2 is None:
            idx2=len(axisdata)

        newaxisdata=axisdata[idx1:idx2]
        if len(newaxisdata)==0:
            raise Exception("No data found in [idx1,idx2).")

        #--------------------Slice data--------------------
        slicer=[slice(None),]*self.ndim
        slicer[axis]=slice(idx1,idx2)
        newdata=self.data[tuple(slicer)]
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
        '''Slice the data by start, end values along a given dimension

        Args:
            v1, v2 (float): start, end dimension values within which to slice
                the data.
            axis (int): axis to do slicing.
        Returns:
            result (NCVAR): sliced data.
        '''
        if v1>v2:
            v1,v2=v2,v1

        try:
            if axis==interpretAxis('time', self):
                v1=pd.to_datetime(v1).to_pydatetime()
                v2=pd.to_datetime(v2).to_pydatetime()
        except:
            pass

        #--------------------Slice axis--------------------
        axisobj=self._axislist[axis]
        axisdata=axisobj.data
        idx=np.ma.where((axisdata>=v1) & (axisdata<v2))[0]
        slicer=[slice(None),]*self.ndim

        if len(idx)>0:
            slicer[axis]=idx
        else:
            raise Exception("No data found in [v1,v2].")

        newdata=self.data[tuple(slicer)]   # if var is np.ndarray

        #--------------Create a new data obj--------------
        newaxisdata=axisdata[idx]
        newaxis=NCVAR(newaxisdata, axisobj.id, [], axisobj.attributes)
        axislist=self.axislist
        axislist[axis]=newaxis
        result=NCVAR(newdata, self.id, axislist, self.attributes)
        for ii in axislist:
            setattr(result, ii.id, ii)
        return result

    def squeeze(self):
        '''Squeeze singleton dimensions
        '''
        axislist=[]
        for ii in self.axislist:
            if len(ii.data)!=1:
                axislist.append(ii)
        self.data=np.squeeze(self.data)
        self.axislist=axislist
        result=NCVAR(self.data, self.id, self.axislist, self.attributes)

        return result

    def getAxis(self, idx):
        '''Get the dimension obj given by axis index

        Args:
            idx (int): axis number.
        Returns:
            result (NCVAR): NCVAR obj acting as the dimension.
        '''
        if idx not in range(self.ndim):
            raise Exception("<idx> not in data shape.")
        return self.axislist[idx]

    def getLatitude(self):
        '''Get latitude dimension values
        '''
        for axisii in self.axislist:
            if axisii.id.lower() in ['y', 'lat', 'latitude']:
                return axisii.data
        return

    def getLongitude(self):
        '''Get longitude dimension values
        '''
        for axisii in self.axislist:
            if axisii.id.lower() in ['x', 'lon', 'longitude']:
                return axisii.data
        return

    def getTime(self):
        '''Get time dimension values
        '''
        for axisii in self.axislist:
            if axisii.id.lower() in ['t', 'time']:
                return axisii.data
        return

    def getLevel(self):
        '''Get level dimension values
        '''
        for axisii in self.axislist:
            if axisii.id.lower() in ['z', 'level']:
                return axisii.data
        return

    def shiftLon(self, dx):
        '''Shift data along longitude dimension

        Args:
            dx (float): target longitude to shift to.
        Returns:
            result (NCVAR): new variable with longitude shifted so that it
                starts from the lonigutde of <dx>.
        '''
        lonidx=interpretAxis('longitude', self)
        if lonidx<0:
            raise Exception("Longitude axis not found in var.")

        lonax=self._axislist[lonidx]
        lons=lonax.data
        dxidx=np.argmin(abs(lons-dx))
        left=lons[:dxidx]
        right=lons[dxidx:]
        #lons=np.roll(lons, -dxidx)
        lons=np.r_[right, left+360]
        lonax.data=lons
        self.data=np.roll(self.data, -dxidx, axis=lonidx)

        return self

def squeezeTo3D(vv):
    '''Squeeze ndarray to 3D

    Args:
        vv (ndarray): 3d or 4d ndarray. If 3d, return as it is. If 4d and 1st
            dimension length is not 1, squeeze all singleton dimensions.
            If 4d and 1st dimension length 1, squeeze the 2nd dimension.
    Returns:
        vv (ndarray): 3d ndarray.
    '''

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
    '''Create a NCVAR for axis

    Args:
        id (str): axis id, e.g. 'time', 'latitude'.
        data (ndarray): 1d array, axis data.
    Keyword Args:
        attributes (dict or None): attribute dictionary.
    Returns:
        resultNV (NCVAR): an NCVAR obj acting as the axis for a variable.
    '''
    if attributes is None:
        attributes={'id': str(id),
                'name': str(id),
                'long_name': 'axis_%s' %str(id),
                'units': ''}
    resultNV=NCVAR(data, id, [], attributes)
    return resultNV

def getBounds(axisdata, width=1.):
    '''Compute axis bounds

    Args:
        axisdata (ndarray): 1d array, axis data.
    Keyword Args:
        width (float): default data spacing, only used when len(axisdata)==1.
    Returns:
        result (ndarray): Nx2 ndarray. 1st column is the lower bounds for
            each axis value, 2nd column the upper bounds.

    Notes: bounds are computed as the mid points between axis points. The
    2 outer bounds are computed using linear extrapolation using the 2 end
    points.
    '''
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


def addExtraAxis(slabNV, newaxis, axis=0, verbose=False):
    """Adds an extra axis to a data slab.

    Args:
        slabNV (NCVAR): variable to which the axis is to insert.
        newaxis (NCVAR): axis object, could be of any length. If None,
            create a dummy singleton axis.
    Keyword Args:
        axis (int): index of axis to be inserted, e.g. 0 if <newaxis> is
            inserted as the 1st dimension.
    Returns:
        slab2NV (NCVAR): variable with an extra axis inserted.
    """

    if newaxis is None:
        newaxis=NCVAR(np.array([0.]), 'newaxis', [], {'id': 'newaxis',
            'name': 'newaxis', 'units': ''})

    # add new axis to axis list of input <slab>
    axislist=slabNV.axislist
    axislist.insert(axis, newaxis)

    #----------------Reshape----------------
    shape=list(slabNV.shape)
    shape.insert(axis,len(newaxis))
    slab2NV=np.reshape(slabNV.data, shape)

    #------------Create variable------------
    att_dict=slabNV.attributes
    slab2NV=NCVAR(slab2NV, id=slabNV.id, axislist=axislist, attributes=att_dict)

    if verbose:
        print('\n# <addExtraAxis>: Originial variable shape:',slabNV.shape)
        print('# <addExtraAxis>: New variable shape:',slab2NV.shape)

    return slab2NV


#-------------Concatenate transient variables---------------------
def cat(var1NV, var2NV, axis=0, verbose=False):
    '''Concatenate 2 NCVAR variables along axis.

    Args:
        var1NV, var2NV (NCVAR): variables to be concatenated, in the order of
            <var1NV>, <var2NV>.
    Keyword Args:
        axis (int): index of axis to be concatenated along.
    Returns:
        result (NCVAR): concatenated variable.

    Result will use the metadata of the 1st variable.
    '''

    attdict=getattr(var1NV, 'attributes', None)
    axis1=var1NV.getAxis(axis)
    axis2=var2NV.getAxis(axis)
    newaxis=np.r_[axis1.data, axis2.data]
    newaxis=NCVAR(newaxis, axis1.id, [], axis1.attributes)

    result=np.ma.concatenate((var1NV.data, var2NV.data), axis=axis)
    axislist=var1NV.axislist
    axislist[axis]=newaxis
    result=NCVAR(result, id=var1NV.id, axislist=axislist, attributes=attdict)

    return result


def concatenate(var_list, axis=0, verbose=False):
    '''Concatenate multiple NCVAR variables along axis.

    Args:
        var_list (list/tuple): variables to be concatenated.
    Keyword Args:
        axis (int): index of axis to be concatenated along.
    Returns:
        resultNV (NCVAR): concatenated variable.

    Result will use the metadata of the 1st variable.
    '''

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
    resultNV=NCVAR(resultdata, id=result.id, axislist=axislist,
            attributes=attdict)

    return resultNV

#------Interpret and convert an axis id to index----------
def interpretAxis(axis, ref_varNV, verbose=True):
    '''Interpret and convert an axis id to index

    Args:
        axis (int or str): axis option, integer or string.
        ref_varNV (NCVAR): NCVAR variable.
    Returns:
        axis_index (int): the index of required axis in <ref_varNV>.

    E.g. index=interpretAxis('time', ref_varNV)
         index=0

         index=interpretAxis(1, ref_varNV)
         index=1
    '''

    if isinstance(axis, (int, np.integer)):
        return axis

    # interpret string dimension
    elif isinstance(axis, str if sys.version_info[0]>=3 else basestring):
        axis=axis.lower()
        dim_idx = -1

        for ii, axisii in enumerate(ref_varNV.axislist):
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
    '''Get attribute dict from an NETCDF obj'''
    attr={}
    for kk in var.ncattrs():
        attr[kk]=var.getncattr(kk)
    return attr

def num2dateWrapper(values, units, calendar='standard'):
    '''A wrapper of netCDF4.num2date to work for different versions

    Args:
        values (ndarray): date values.
        units (str): date units.
    Keyword Args:
        calendar (str): calendar used. Default to 'standard'.
    Returns:
        axis (ndarray): date values converted to python datetime objs.

    NOTE this is to work for different netcdf4 version. Versions prior to
    1.4 do not seem to have the kwarg "only_use_cftime_datetimes" and
    "only_use_python_datetimes".
    '''

    try:
        axis=num2date(values, units, calendar=calendar,
                only_use_cftime_datetimes=False,
                only_use_python_datetimes=True)
    except:
        try:
            axis=num2date(values, units, calendar=calendar,
                    only_use_cftime_datetimes=False)
        except:
            axis=num2date(values, units, calendar=calendar)

    return axis

def readNC(abpath_in, varid, calendar='standard'):
    '''Read in a variable from an netcdf file

    Args:
        abpath_in (str): absolute file path to the netcdf file.
        varid (str): id of variable to read.
    Keyword Args:
        calendar (str): calendar used. Default to 'standard'.
    Returns:
        ncvarNV (NCVAR): variable stored as an NCVAR obj.
    '''

    fin=Dataset(abpath_in, 'r')
    var=fin.variables[varid]

    #---------------------Get axes---------------------
    dims=var.dimensions
    axislist=[]
    for dd in dims:
        ncaxis=fin.variables[dd]
        if dd.lower() in ['time', 't']:
            # convert time to datetime objs
            axisii=num2dateWrapper(ncaxis[:], ncaxis.units, calendar)
        else:
            axisii=ncaxis[:]
        axisattr=getAttributes(ncaxis)
        axisattr['isunlimited']=fin.dimensions[dd].isunlimited()
        axisii=NCVAR(axisii, dd.lower(), [], axisattr)
        axislist.append(axisii)

    attributes=getAttributes(var)
    data=var[:]

    # create NCVAR var
    ncvarNV=NCVAR(data, varid, axislist, attributes)
    for ii in axislist:
        setattr(ncvarNV, ii.id, ii)

    ncvarNV=increasingLatitude(ncvarNV, interpretAxis('latitude', ncvarNV))

    return ncvarNV


def saveNC(abpath_out, varNV, mode='w', dtype=None):
    '''Save data to netcdf file

    Args:
        abpath_out (str): absolute file path.
        varNV (NCVAR): NCVAR obj.
    Keyword Args:
        mode (str): writing mode.
        dtype (None or str): date type. Default to np.float32.
    '''

    with Dataset(abpath_out, mode) as fout:
        saveNCDims(fout, varNV.axislist)
        _saveNCVAR(fout, varNV, dtype)


def _saveNCVAR(fout, varNV, dtype=None):
    '''Write variable data to netcdf file

    Args:
        fout (netcdf file handle): opened netcdf file handle.
        varNV (NCVAR): NCVAR obj.
    Keyword Args:
        dtype (None or str): date type. Default to np.float.
    '''

    #-----------------Create variable-----------------
    if dtype is None:
        dtype=np.float32

    if varNV.id not in fout.variables.keys():
        varout=fout.createVariable(varNV.id, dtype, varNV.dims, zlib=True)
        #varout.set_collective(True)
        #----------------Create attributes----------------
        for kk, vv in varNV.attributes.items():
            try:
                varout.setncattr(kk, vv)
            except:
                pass
        # assign data
        varout[:]=varNV.data
    else:
        varout=fout.variables[varNV.id]
        timeax=fout.variables['time']
        tlen=len(timeax)
        tnew=varNV.getTime()
        t0=tnew[0]
        tidx=np.where(timeax==t0)[0]
        if len(tidx)>0:
            # if time already exists
            tidx=tidx[0]
        else:
            # new time point
            tidx=tlen
        timeax[tidx:tidx+len(tnew)]=tnew
        varout[tidx:tidx+len(tnew)]=varNV.data

    return


def saveNCDims(fout, axislist):
    '''Write dimension info to netcdf file

    Args:
        fout (netcdf file handle): opened netcdf file handle.
        axislist (tuple or list): tuple/list of NCVAR objs acting as axes.
    '''

    #----------------Create dimensions----------------
    for aa in axislist:

        if aa.id in fout.dimensions.keys():
            continue
        if aa.id in ['t', 'time', 'Time', 'T']:
            lenii=None
            if all([isinstance(ii, (datetime, cftime.datetime)) for ii in aa.data]):
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
            if kk!='isunlimited' and kk!='_FillValue':
                axisvar.setncattr(kk, vv)

    return


#----------Check exsitance of files in file list-----------
def checkFiles(file_list, verbose=True):
    '''Check existance of files in a list.

    Args:
        file_list (list): a list of ABSOLUTE paths to be checked.

    Usefull before a long list of iteration to make sure every data
    file are ready on the disk.

    Function prompts enquiry if any file is missing in the list.
    '''

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

    Args:
        slab (ndarray): ndarray, possibly contains masked values or nans.
    Returns:
        mask (ndarray): bindary, 1s for missing, 0s otherwise.
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
def getSlab(var, index1=-1, index2=-2, verbose=True):
    '''Get a 2d slab from a variable

    Args:
        var (ndarray): ndarray with dimension >=2.
        index1,index2 (int): indices denoting the dimensions from which a 2d
            slab is to slice.
    Returns:
        slab (ndarray): the (1st) 2d slab from <var>.
               E.g. <var> has dimension (12,1,241,480), getSlab(var) will
               return the 1st time point with singleton dimension squeezed.
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
    string='slab=var[tuple(%s)]' %','.join(slicer)
    exec(string)
    return slab

#-------------Change latitude axis to south-to-north-----------------
def increasingLatitude(slab, axis, lats=None, verbose=False):
    '''Changes a slab so that is always has latitude running from
    south to north.

    Args:
        slab (NCVAR or ndarray): input data.
        axis (int): axis index for the latitude dimension.
    Keyword Args:
        lats (None or ndarray): if 1darray, the latitude values. Use None
            when <slab> is an NCVAR which comes with its latitude axis.
    Returns:
        slab2 (NCVAR or ndarray): data with latitude orientated from south
            to north. It is an NCVAR if <slab> is NCVAR, ndarray if <slab>
            is ndarray.

    See also increasingLatitude2()
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

    Args:
        slab (ndarray): input data.
        axis (int): axis index for the latitude dimension.
        lats (ndarray): 1darray, the latitude values.
    Returns:
        slab (ndarray): data with latitude orientated from south
            to north.
        lats (ndarray): 1darray, oriented latitude values.

    See also increasingLatitude()
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


def dLongitude(lats, lons, side='c', R=6371000):
    '''Return a slab of longitudinal increment (meter) delta_x.

    Args:
        lats (ndarray): 1d array, latitude coordinates in degrees.
        lons (ndarray): 1d array, longitude coordinates in degrees.
    Keyword Args:
        side (str): 'n': northern boundary of each latitudinal band;
                    's': southern boundary of each latitudinal band;
                    'c': central line of latitudinal band;

                 -----     'n'
                /-----\     'c'
               /_______\     's'

        R (float): radius of Earth.
    Returns:
        delta_x (ndarray): 2d array, longitudinal increments.
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
def dLatitude(lats, lons, R=6371000, verbose=True):
    '''Return a slab of latitudinal increment (meter) delta_y.

    Args:
        lats (ndarray): 1d array, latitude coordinates in degrees.
        lons (ndarray): 1d array, longitude coordinates in degrees.
    Keyword Args:
        R (float): radius of Earth;
    Returns:
        delta_x (ndarray): 2d array, latitudinal increments.
            <var>.
    '''

    lons=np.sort(lons)

    #---------Get axes and bounds-------------------
    latax_bounds=getBounds(lats)
    delta_y=np.ptp(latax_bounds,axis=1)*np.pi/180.*R

    #-------Repeat array to get slab---------------
    delta_y=np.repeat(delta_y[:,None],len(lons),axis=1)

    return delta_y

#--------------------Get contour from a binary mask--------------------
def getBinContour(mask,lons=None,lats=None,return_largest=True):
    '''Get contour from a binary mask

    Args:
        mask (ndarray): 2d array, binary mask.
    Keyword Args:
        lons,lats (1darray or None): if 1d array, x, y coordinates for <mask>.
            if None, create a coordinate array with indices.
    Returns:
        cont (ndarray): Nx2 array, coordinates of the contour of the largest
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

    if lons is None:
        lons=np.arange(mask.shape[1])
    if lats is None:
        lats=np.arange(mask.shape[0])

    cs=plt.contourf(lons,lats,mask,[0.9,1.1]).collections
    conts=cs[0].get_paths()
    if return_largest:
        conts.sort(key=lambda x:len(x.vertices))
        #cont=conts[-1].vertices
        cont=conts[-1]
    else:
        cont=conts
    plt.clf()
    plt.close()
    return cont

#-----------Find index of value in array-----------
def findIndex(x,a):
    '''Find index of value in array

    Args:
        x (float or int): scalar, value to search.
        a (1d array): array to search from.
    Returns:
        idx (int): index in <a> that a[idx] is closest to <x>.
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
    if idx==len(a)-1 and abs(x-a[-1]) > abs(a[-1]-a[-2]):
        idx=None
    return idx


def getBearing(lat1,lon1,lat2,lon2):
    '''Compute bearing from point 1 to point2

    Args:
        lat1,lat2 (float or ndarray): scalar float or nd-array, latitudes in
            degree for location 1 and 2.
        lon1,lon2 (float or ndarray): scalar float or nd-array, longitudes in
            degree for location 1 and 2.
    Returns:
        theta (float or ndarray): (forward) bearing in degree.

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

    NOTE: deprecated, use netCDF4 instead of CDAT.
    '''

    """
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
    """
    pass


def getTimeAxis(times, ntime, ref_time='days since 1900-01-01'):
    '''Create a time axis

    Args:
        times (list or tuple or array): array of datetime objs, or strings,
            giving the time stamps in the format of 'yyyy-mm-dd HH:MM'. It is
            assumed to be in chronological order.
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
    if all([isinstance(ii, (datetime, cftime.datetime)) for ii in times]):
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
    so that when plotted in a basemap/cartopy plot, a zonally cyclic curve
    won't be plotted as jumping straight lines linking the left and right
    bounds.
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
def getGridsInContour(contour, x, y):
    '''Get binary grids inside a contour

    Args:
        contour (ndarray): Nx2 ndarray, (x,y) coordinates of a contour.
               Or a matplotlib Path obj. If the latter, function can
               remove holes if exists in the contour.
        x,y (ndarray): 1d array, x and y coordinates of the grid.
    Returns:
        validcoords (ndarray): 2d array of shape (len(y), len(x)), binary slab
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
    coords=np.array(list(zip(X2.flat, Y2.flat)))

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
    return np.abs(signedArea(x, y))

def signedArea(x, y):
    '''Get signed polygon area

    Args:
        x, y (ndarray): 1d array, x and y coordinates of the polygon.
    Returns:
        area (float): signed area of the polygon.
    '''

    x=np.asarray(x)
    y=np.asarray(y)
    def isClosed(xs,ys):
        if np.alltrue([np.allclose(xs[0],xs[-1]),\
            np.allclose(ys[0],ys[-1]),np.ptp(xs),np.ptp(ys)]):
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
    '''Compute weighted average along a given axis/axes.

    Args:
        slab (ndarray): data to compute average.
    Keyword Args:
        axis (int or tuple/list): if int, axis to compute average along.
            if tuple/list, multiple axes to compute averages along.
        weights (str): 'equal': compute normal average.
                       'generate': compute weighted average using axis
                       coordinates increments as weights.
        coords (ndarray or tuple/list or None): if ndarray, an 1d array
            giving the axis coordinates of the axis/axes given in <axis>.
            If None, the same as weights='equal', i.e. equal weights.
        keepdims (bool): if True, remove the singleton axis after averaging.
    Returns:
        result (ndarray): averaged result.
    '''

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
    '''Compute area-weighted average.

    Args:
        slab (ndarray): data to compute average.
    Keyword Args:
        axis (tuple/list): tuple/list of the axes to compute areal averages.
        weights (str): 'equal': compute normal average.
                       'generate': compute weighted average using grid cell
                       areas as weights.
        lats, lons (ndarray or None): if ndarray, an 1d array
            giving the y- and x- coordinates.
            If both None, the same as weights='equal', i.e. equal weights.
        keepdims (bool): if True, remove the singleton axis after averaging.
    Returns:
        result (ndarray): averaged result.
    '''

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

            dlat=dLatitude(lats, lons)
            dlon=dLongitude(lats, lons)
            area=dlat*dlon
            area=slabGrow(area, slab.shape)

            if hasattr(slab, 'mask'):
                area.mask=slab.mask

            num=np.ma.sum(slab*area, axis=axis, keepdims=True)
            den=np.ma.sum(area, axis=axis, keepdims=True)
            result=num/den

    if not keepdims:
        result=np.squeeze(result, axis=axis)

    return result


def areaStd(slab, axis=(-2, -1), weights='equal', lats=None, lons=None,
        keepdims=True, verbose=True):
    '''Compute area-weighted standard deviations.

    Args:
        slab (ndarray): data to compute std.
    Keyword Args:
        axis (tuple/list): tuple/list of the axes to compute std.
        weights (str): 'equal': compute normal std.
                       'generate': compute weighted std using grid cell
                       areas as weights.
        lats, lons (ndarray or None): if ndarray, an 1d array
            giving the y- and x- coordinates.
            If both None, the same as weights='equal', i.e. equal weights.
        keepdims (bool): if True, remove the singleton axis after averaging.
    Returns:
        result (ndarray): area-weighted standard deviations.
    '''

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

            dlat=dLatitude(lats, lons)
            dlon=dLongitude(lats, lons)
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

#---------------Grow a 2D slab to 3D or 4D--------------------
def slabGrow(slab, shape, verbose=False):
    '''Grow an 2D slab to 3D or 4D slab given target shape

    Args;
        slab (ndarray): 2d array to be repeated.
        shape (tuple/list): target shape to broadcast to.
    Returns:
        slab (ndarray): slab repeatd to target shape.

    Usually used to extend one coeffient to match the size of a variable.
    e.g. <slab> is a 2D land-sea-mask, the function then grows the 2D slab
            to match dimention of a variable, either (t,z,y,x) or (t,y,x).
    '''

    slabshape=slab.shape

    match=False
    for ii, jj in zip(shape[:-1], shape[1:]):
        if ii==slabshape[0] and jj==slabshape[1]:
            match=True
            break
    if not match:
        raise Exception("Doesn't find matching slab in shape.")

    order=list(range(len(shape)))
    order_remain=order[0:-2]
    order_remain.reverse()
    for indexii in order_remain:
        slab=np.ma.reshape(slab, (1,)+slab.shape)
        slab=np.ma.repeat(slab, shape[indexii], axis=0)
    slab=np.ma.array(slab)

    if verbose:
        print('# <slabGrow>: New slab shape:', slab.shape)

    return slab


#------------Grow an 1D array to 2D, 3D or 4D-------------------
def arrayGrow(array, shape, axis=None, verbose=True):
    '''Grow an 1D array to a 2D, 3D or 4D slab given target shape

    Args:
        array (ndarray): 1d array to be repeated.
        shape (tuple/list): target shape to broadcast to.
    Keyword Args:
        axis (int or None): an optional argument to specify the axis index of
            the array after broadcast. Usefull when 2 axes in <shape> have
            same length.
            If None, use the length of <array> to find the axis index.
            E.g. len(array)=100, shape=(4, 100, 200). Then axis is 1.
            If len(array)=100, shape=(3, 100, 100), can use axis=1 or axis=2
            to specify which axis should match <array>.
    Returns:
        newarray (ndarray): array repeatd to target shape.

    Usually used to extend one coefficient to match the size of a variable.
    e.g. <array> is longitudinal increment calculated using latitude values,
    the function then grows the 1d array to match dimentions of reference
    variable <ref_var>, which could be (t,z,y,x) or (t,y,x).
    '''

    #-----------------Check dimension-----------------
    array=np.atleast_1d(np.squeeze(array))
    if np.ndim(array)>1:
        raise Exception("<array> needs to be 1d array.")

    #--------------------Find axis--------------------
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
    tmpshape=np.array(shape)[order]

    #---------------Expand by repeating---------------
    array=np.array(array)
    newarray=array[np.newaxis,:]
    newarray=np.repeat(newarray, np.prod(shape)/length, axis=0)
    newarray=np.reshape(newarray,tmpshape)
    newarray=np.transpose(newarray,order)

    return newarray
