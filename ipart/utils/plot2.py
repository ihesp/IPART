'''Plotting Functions.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2017-10-09 14:32:05.
'''

#--------Import modules--------------
import numpy
import matplotlib.pyplot as plt
from . import funcs2 as functions

DEFAULT_CMAP=plt.cm.PRGn
DEFAULT_CMAP=plt.cm.RdBu_r
#DEFAULT_CMAP=plt.cm.bwr




def mkscale(n1,n2,nc=12,zero=1):
    '''Copied from vcs/util.py

    Function: mkscale

    Description of function:
    This function return a nice scale given a min and a max

    option:
    nc # Maximum number of intervals (default=12)
    zero # Not all implemented yet so set to 1 but values will be:
           -1: zero MUST NOT be a contour
            0: let the function decide # NOT IMPLEMENTED
            1: zero CAN be a contour  (default)
            2: zero MUST be a contour
    Examples of Use:
    >>> vcs.mkscale(0,100)
    [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    >>> vcs.mkscale(0,100,nc=5)
    [0.0, 20.0, 40.0, 60.0, 80.0, 100.0]
    >>> vcs.mkscale(-10,100,nc=5)
    [-25.0, 0.0, 25.0, 50.0, 75.0, 100.0]
    >>> vcs.mkscale(-10,100,nc=5,zero=-1)
    [-20.0, 20.0, 60.0, 100.0]
    >>> vcs.mkscale(2,20)
    [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
    >>> vcs.mkscale(2,20,zero=2)
    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]

    '''
    if n1==n2: return [n1]
    import numpy
    nc=int(nc)
    cscale=0  # ???? May be later
    #import pdb; pdb.set_trace()
    min=numpy.min([n1,n2])
    max=numpy.max([n1,n2])
    if zero>1.:
        if min>0. : min=0.
        if max<0. : max=0.
    rg=float(max-min)  # range
    delta=rg/nc # basic delta
    # scale delta to be >10 and <= 100
    lg=-numpy.log10(delta)+2.
    il=numpy.floor(lg)
    delta=delta*(10.**il)
    max=max*(10.**il)
    min=min*(10.**il)
    if zero>-0.5:
        if delta<=20.:
            delta=20
        elif delta<=25. :
            delta=25
        elif delta<=40. :
            delta=40
        elif delta<=50. :
            delta=50
        elif delta<=60.:
            delta=60
        elif delta<=70.:
            delta=70
        elif delta<=80.:
            delta=80
        elif delta<=90.:
            delta=90
        elif delta<=101. :
            delta=100
        first = numpy.floor(min/delta)-1.
    else:
        if delta<=20.:
            delta=20
        elif delta<=25. :
            delta=25
        elif delta<=40. :
            delta=40
        elif delta<=50. :
            delta=50
        elif delta<=60.:
            delta=60
        elif delta<=70.:
            delta=70
        elif delta<=80.:
            delta=80
        elif delta<=90.:
            delta=90
        elif delta<=101. :
            delta=100
        first=numpy.floor(min/delta)-1.5

    scvals=delta*(numpy.arange(2*nc)+first)
    a=0
    for j in range(len(scvals)):
        if scvals[j]>min :
            a=j-1
            break
    b=0
    for j in range(len(scvals)):
        if scvals[j]>=max :
            b=j+1
            break
    if cscale==0:
        cnt=scvals[a:b]/10.**il
    else:
        #not done yet...
        raise Exception('ERROR scale not implemented in this function')
    return list(cnt)


#---Translate an integer index to letter index---
def index2Letter(index,verbose=True):
    '''Translate an integer index to letter index

    <index>: integer index for a subplot.

    Return <letter>: corresponding letter index for <index>.

    <index> to letter indexing is defined as following:
    ----------------------------------------------------
    <index>                     letter index
    ----------------------------------------------------
       1                            (a)
       2                            (b)
       3                            (c)
       ...                          ...
    ----------------------------------------------------
    '''

    #---------------Create <index> dict---------------
    index_dict={\
            1:'(a)',\
            2:'(b)',\
            3:'(c)',\
            4:'(d)',\
            5:'(e)',\
            6:'(f)',\
            7:'(g)',\
            8:'(h)',\
            9:'(i)',\
            10:'(j)',\
            11:'(k)',\
            12:'(l)',\
            13:'(m)',\
            14:'(n)',\
            15:'(o)',\
            16:'(p)',\
            17:'(q)',\
            18:'(r)',\
            19:'(s)',\
            20:'(t)',\
            21:'(u)',\
            22:'(v)',\
            23:'(w)',\
            24:'(x)',\
            25:'(y)',\
            26:'(z)'
            }

    #-------------------Check inputs-------------------
    if index<=0:
        raise Exception("<index> needs to be positive.")
    if index>26:
        return str(index)
    else:
        return index_dict[index]


#--------------Get a matplotlib blue-to-red colormap obj--------------
def colormapBR(reverse=False,verbose=True):
    '''Get a matplotlib blue-to-red colormap obj

    <reverse>: bool, if True, colormap goes from red-white-blue.

    Return <cmap>: matplotlib colormap.
    Update time: 2015-10-25 16:09:53.
    '''

    import matplotlib.pyplot as plt

    if reverse:
        return plt.cm.bwr_r
    else:
        return plt.cm.bwr


#-------------Map color indices----------------------
def mapColor(levels,cmap,split=0,verbose=True):
    '''Map a list of levels linearly onto a colormap.

    <levels>: a list of level values;
    <cmap>: matplotlib colormap obj;
    <split>: int, control behavior of negative and positive values\
            0: Do not split negatives and positives, map onto entire\
               range of [0,1];
            1: split only when vmin<0 and vmax>0, otherwise\
               map onto entire range of [0,1];
               If split is to be applied, negative values are mapped onto first half [0,0.5],
               and postives onto second half (0.5,1].
            2: force split, if vmin<0 and vmax>0, same as <split>==1;
                If vmin<0 and vmax<0, map onto 1st half [0,0.5];
                If vmin>0 and vmax>0, map onto 2nd half (0.5,1].

    Return: <colors>: a list of matplotlib colors sampled from <cmap>,
                      each corresponds to a level in <levels>.

    Usage:
        cmap=plt.cm.jet
        xx=numpy.arange(12)
        colors=plot.mapColor(range(10),cmap,split=0)
        for ii in range(10):
            yy=numpy.sin(xx)+ii
            plt.plot(xx,yy,color=colors[ii])
        plt.show()

    Update time: 2015-10-26 14:03:44.
    '''

    import numpy

    levels=numpy.sort(levels)
    minlevel=levels[0]
    maxlevel=levels[-1]

    #-------------Shift colormap if needed-------------
    cmap=remappedColorMap(cmap,minlevel,maxlevel,split)
    colors=[cmap(1.*ii/len(levels)) for ii in range(len(levels))]

    return colors


def remappedColorMap(cmap,vmin,vmax,split=2,name='shiftedcmap'):
    '''Re-map the colormap to split positives and negatives.

    <cmap> : The matplotlib colormap to be altered
    <vmin> : float, minimal level in data.
    <vmax> : float, maximal level in data.
    <split>: int, control behavior of negative and positive values\
            0: Do not split negatives and positives, map onto entire\
               range of [0,1];
            1: split only when vmin<0 and vmax>0, otherwise\
               map onto entire range of [0,1];
               If split is to be applied, negative values are mapped onto first half [0,0.5],
               and postives onto second half (0.5,1].
            2: force split, if vmin<0 and vmax>0, same as <split>==1;
                If vmin<0 and vmax<0, map onto 1st half [0,0.5];
                If vmin>0 and vmax>0, map onto 2nd half (0.5,1].
    '''
    import numpy
    import matplotlib
    import matplotlib.pyplot as plt

    #-----------Return cmap if not splitting-----------
    if split==0:
        return cmap
    if split==1 and vmin*vmax>=0:
        return cmap

    #------------Resample cmap if splitting------------
    cdict = {\
    'red': [],\
    'green': [],\
    'blue': [],\
    'alpha': []\
    }

    vmin,vmax=numpy.sort([vmin,vmax]).astype('float')

    shift_index=numpy.linspace(0,1,256)

    #-------------------Force split-------------------
    if vmin<0 and vmax<=0 and split==2:
        if vmax<0:
            idx=numpy.linspace(0,0.5,256,endpoint=False)
        else:
            idx=numpy.linspace(0,0.5,256,endpoint=True)

    elif vmin>=0 and vmax>0 and split==2:
        '''
        if vmin>0:
            idx=numpy.linspace(0.5,1,256,endpoint=False)
        else:
            idx=numpy.linspace(0.5,1,256,endpoint=True)
        '''
        idx=numpy.linspace(0.5,1,256,endpoint=True)

    #--------------------Split -/+--------------------
    if vmin*vmax<0:
        mididx=int(abs(vmin)/(vmax+abs(vmin))*256)
        idx=numpy.hstack([numpy.linspace(0,0.5,mididx),\
                numpy.linspace(0.5,1,256-mididx)])

    #if vmin*vmax==0:
        #return cmap

    #-------------------Map indices-------------------
    for ri,si in zip(idx,shift_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


#--------------Get levels for isofill--------------
def isoLevels(vars,num=15,zero=1,max_level=None,min_level=None,\
        ql=None,qr=None):
    '''Get levels for isofill plot.

    <vars>: one or a list of variables, from which the <minlevel>
            and <maxlevel> is obtained.
            If <vars> has more than 1 variables, then function calculates
            the minimum and maximum of all variables, thus the color
            legend would be unified for all subplots.
            Note that if <max_level> and/or <min_level> are given,
            they will override the maximum/minimal levels derived
            from <vars>.
    <Num> is the (maximum) number of isoline levels;
    <zero>: controls 0 in created <levels>:\
            -1: 0 is NOT allowed to be a level;
            1: 0 is permitted to be a level;
            2: 0 forced to be a level.
    <max_level>,<min_level>: the max/min limits to be plotted out, values
                            outside the range will be grouped into the
                            last level intervals on both ends.
    <ql>,<qr>: extreme percentiles for the lower and upper boundaries.
               Could be one of the values in the list:

                    percents=[0.001,0.005,0.01,0.025,0.05,0.1]

               E.g. ql=0.001; qr=0.005
               means that the 0.1% percentile will be set to the minimum level;
               0.005% (from the right most, or 99.95% from the left most) percentil will
               be set to the maximum level.

               If both <ql> and <min_level> are given, use <min_level>.
               If both <qr> and <max_level> are given, use <max_level>.

    ###NOTE:###
    <minlevel> and <maxlevel> are better derived using numpy.min() and\
            numpy.max(). MV.min() and MV.max() may have problems.

    *Change: Mar-05-2013:
             Now function uses vcs.mkscale, instead of levelByNum() to
             generate isofill levels.
    *Change: Mar-07-2013:
             Now function receives a list of variables (<*var>) to obtain
             <minlevel> and <maxlevel>, for the purpose of unifying legend
             for multiple plots.
    Update time: 2013-08-15 14:42:48: added ext_1, ext_2 control.
    Update time: 2014-06-16 20:21:17: added <ql>, <qr> percentiles.
    '''


    if type(vars) is list or type(vars) is tuple:
        pass
    else:
        vars=[vars,]

    #-------------------Get max/min-------------------
    vmin,vmax=getRange(vars,min_level,max_level,ql,qr)
    data_min,data_max=getRange(vars,verbose=False)

    #--------------------Make scale--------------------
    levels=mkscale(vmin,vmax,num,zero)

    #---------Use actual data range if overflow---------
    '''
    levels=[ii if ii>=data_min else data_min for ii in levels]
    levels=[ii if ii<=data_max else data_max for ii in levels]
    levels=list(set(levels))
    levels.sort()
    '''

    ext1=True if data_min<numpy.min(levels) else False
    ext2=True if data_max>numpy.max(levels) else False


    return levels,ext1,ext2


#----------------Get min/max from 1 or more variables----------------
def getRange(vars,min_level=None,max_level=None,\
        ql=None,qr=None,verbose=True):
    '''Get min/max value

    <vars>: a list of variables.
    '''

    import numpy

    #-----------------Print some stats-----------------
    percents=[0.001,0.005,0.01,0.025,0.05,0.1]
    if ql is not None or qr is not None:
        for ii,vii in enumerate(vars):
            vii=vii.flatten()
            #vii=numpy.reshape(vii,(1,)+vii.shape)
            if ii==0:
                var_all=vii
            else:
                var_all=numpy.concatenate((var_all,vii),axis=0)

        if verbose:
            print('\n# <getRange>: Get quantiles for vars')
        quantiles=functions.getQuantiles(var_all,percents,verbose)

    #-------Get min/max from all vars-----------------
    min_list=[numpy.nanmin(i) for i in vars]
    max_list=[numpy.nanmax(i) for i in vars]
    minlevel=min(min_list)
    maxlevel=max(max_list)

    #----------------Set lower boundary----------------
    if min_level is not None and ql is None:
        vmin=max(minlevel,min_level)
    elif min_level is None and ql is not None:
        if ql not in percents:
            raise Exception('<ql> should be one of'+repr(percents)+'.')
        vmin=quantiles[percents.index(ql)][0]
    elif min_level is not None and ql is not None:
        vmin=max(minlevel,min_level)
    else:
        vmin=minlevel

    #----------------Set upper boundary----------------
    if max_level is not None and qr is None:
        vmax=min(maxlevel,max_level)
    elif max_level is None and qr is not None:
        if qr not in percents:
            raise Exception('<qr> should be one of'+repr(percents)+'.')
        vmax=quantiles[percents.index(qr)][1]
    elif max_level is not None and qr is not None:
        vmax=min(maxlevel,max_level)
    else:
        vmax=maxlevel

    if vmax<vmin:
        vmin, vmax=vmax, vmin

    return vmin,vmax


class Isofill(object):
    def __init__(self,vars,num=15,zero=1,split=2,max_level=None,min_level=None,\
        ql=None,qr=None,cmap=None,verbose=True):
        '''Return an isofill object with specified color scheme.

        Args:
            vars: one or a list of variables, from which the <minlevel>
                    and <maxlevel> is obtained.
                    If <vars> has more than 1 variables, then function calculates
                    the minimum and maximum of all variables, thus the color
                    legend would be unified for all subplots.
                    Note that if <max_level> and/or <min_level> are given,
                    they will override the maximum/minimal levels derived
                    from <vars>.
        Keyword Args:
            Num (int): is the (maximum) number of isoline levels;
            zero (int): controls 0 in created <levels>:\
                    -1: 0 is NOT allowed to be a level;
                    1: 0 is permitted to be a level;
                    2: 0 forced to be a level.
            split (int): int, control behavior of negative and positive values\
                    0: Do not split negatives and positives, map onto entire\
                       range of [0,1];
                    1: split only when vmin<0 and vmax>0, otherwise\
                       map onto entire range of [0,1];
                       If split is to be applied, negative values are mapped onto first half [0,0.5],
                       and postives onto second half (0.5,1].
                    2: force split, if vmin<0 and vmax>0, same as <split>==1;
                        If vmin<0 and vmax<0, map onto 1st half [0,0.5];
                        If vmin>0 and vmax>0, map onto 2nd half (0.5,1].
            max_level,min_level (float): the max/min limits to be plotted out, values
                                    outside the range will be grouped into the
                                    last level intervals on both ends.
            ql,qr (float): extreme percentiles for the lower and upper boundaries.
                       Could be one of the values in the list:

                            percents=[0.001,0.005,0.01,0.025,0.05,0.1]

                       E.g. ql=0.001; qr=0.005
                       means that the 0.1% percentile will be set to the minimum level;
                       0.005% (from the right most, or 99.95% from the left most) percentil will
                       be set to the maximum level.

                       If both <ql> and <min_level> are given, use <min_level>.
                       If both <qr> and <max_level> are given, use <max_level>.
            cmap: specify a color map. Could be:
                    1) None, a default blue-white-red (bwr) color map will be created.
                    2) a matplotlib cmap obj.
                    3) a string name of a matplotlib cmap obj, to list a few:

                        'bwr':      blue-white-red
                        'RdBu':     darkred-white-blue
                        'RdYlBu':   red-yellow-white-blue
                        'RdYlGn':   red-yellow-white-green
                        'spectral': purple-yellow-cyan-blue
                        'seismic':  darkblue-white-darkred
                        'jet':      rainbow darkblue-darkred
                        'rainbow':  rainbow purple-red

                    Append '_r' to get the reversed colormap.

        Note:
            <minlevel> and <maxlevel> are better derived using numpy.min() and\
                    numpy.max(). MV.min() and MV.max() may have problems.
            Iso levels are computed using vcs function (mkscale()), and a matplotlib
            colormap is created (if not given), and the colormap will be changed so
            positive/negative splits (if required) is achieved.

        Update time: 2015-04-27 14:55:33
        '''

        self.functions=functions
        self.plt=__import__('matplotlib.pyplot')

        self.vars=vars
        self.num=num
        self.zero=zero
        self.split=split
        self.max_level=max_level
        self.min_level=min_level
        self.ql=ql
        self.qr=qr
        self.cmap=cmap
        self.ext_1=False
        self.ext_2=False
        self.method='isofill'

        #---------------Get 1st time points---------------
        if type(vars) is list or type(vars) is tuple:
            pass
        else:
            vars=[vars,]

        #self.vars=[get1stTime(ii) for ii in vars]
        if self.split not in [0,1,2]:
            raise Exception("<split> not in [0,1,2].")

        #--------------------Get levels--------------------
        self.levels=self.getLevels()

        #-------------------Get colormap-------------------
        self.cmap=self.getColormap()


        return

    #--------------Create contour levels--------------
    def getLevels(self):

        levels,ext1,ext2=isoLevels(self.vars,self.num,self.zero,\
                        self.max_level,self.min_level,self.ql,self.qr)

        #----------------Set extends----------------
        self.ext_1=ext1
        self.ext_2=ext2

        return levels




    #-----------------Create colormap-----------------
    def getColormap(self):
        if self.cmap is None:
            #cmap=self.plt.cm.bwr
            cmap=DEFAULT_CMAP
        elif self.cmap is not None and type(self.cmap) is str:
            cmpstr='cmap=self.plt.cm.'+self.cmap
            try:
                exec(cmpstr)
            except:
                raise Exception("Color map name wrong.")
        else:
            cmap=self.cmap

        #-------------Shift colormap if needed-------------
        cmap=remappedColorMap(cmap,self.levels[0],self.levels[-1],self.split)

        return cmap


class Boxfill(object):
    def __init__(self,vars,zero=1,split=2,max_level=None,min_level=None,\
        ql=None,qr=None,cmap=None,verbose=True):
        '''Return an isofill object with specified color scheme.

        Args:
            vars: one or a list of variables, from which the <minlevel>
                    and <maxlevel> is obtained.
                    If <vars> has more than 1 variables, then function calculates
                    the minimum and maximum of all variables, thus the color
                    legend would be unified for all subplots.
                    Note that if <max_level> and/or <min_level> are given,
                    they will override the maximum/minimal levels derived
                    from <vars>.
        Keyword Args:
            zero (int): controls 0 in created <levels>:\
                    -1: 0 is NOT allowed to be a level;
                    1: 0 is permitted to be a level;
                    2: 0 forced to be a level.
            split (int): int, control behavior of negative and positive values\
                    0: Do not split negatives and positives, map onto entire\
                       range of [0,1];
                    1: split only when vmin<0 and vmax>0, otherwise\
                       map onto entire range of [0,1];
                       If split is to be applied, negative values are mapped onto first half [0,0.5],
                       and postives onto second half (0.5,1].
                    2: force split, if vmin<0 and vmax>0, same as <split>==1;
                        If vmin<0 and vmax<0, map onto 1st half [0,0.5];
                        If vmin>0 and vmax>0, map onto 2nd half (0.5,1].
        Keyword Args:
            max_level,min_level (float): the max/min limits to be plotted out, values
                                    outside the range will be grouped into the
                                    last level intervals on both ends.
            ql,qr (float): extreme percentiles for the lower and upper boundaries.
                       Could be one of the values in the list:

                            percents=[0.001,0.005,0.01,0.025,0.05,0.1]

                       E.g. ql=0.001; qr=0.005
                       means that the 0.1% percentile will be set to the minimum level;
                       0.005% (from the right most, or 99.95% from the left most) percentil will
                       be set to the maximum level.

                       If both <ql> and <min_level> are given, use <min_level>.
                       If both <qr> and <max_level> are given, use <max_level>.
            cmap: specify a color map. Could be:
                    1) None, a default blue-white-red (bwr) color map will be created.
                    2) a matplotlib cmap obj.
                    3) a string name of a matplotlib cmap obj, to list of few:

                        'bwr':      blue-white-red
                        'RdBu':     darkred-white-blue
                        'RdYlBu':   red-yellow-white-blue
                        'RdYlGn':   red-yellow-white-green
                        'spectral': purple-yellow-cyan-blue
                        'seismic':  darkblue-white-darkred
                        'jet':      rainbow darkblue-darkred
                        'rainbow':  rainbow purple-red

                    Append '_r' to get the reversed colormap.

        NOTE:
            <minlevel> and <maxlevel> are better derived using numpy.min() and\
                    numpy.max(). MV.min() and MV.max() may have problems.
            Iso levels are computed using vcs function (mkscale()), and a matplotlib
            colormap is created (if not given), and the colormap will be changed so
            positive/negative splits (if required) is achieved.

        Update time: 2015-04-27 14:55:33
        '''
        self.functions=functions
        self.plt=__import__('matplotlib.pyplot')

        self.vars=vars
        self.zero=zero
        self.split=split
        self.max_level=max_level
        self.min_level=min_level
        self.ql=ql
        self.qr=qr
        self.cmap=cmap
        self.ext_1=False
        self.ext_2=False

        self.method='boxfill'

        #---------------Get 1st time points---------------
        if type(vars) is list or type(vars) is tuple:
            pass
        else:
            vars=[vars,]

        self.vars=[get1stTime(ii) for ii in vars]

        #-------------------Get max/min-------------------
        self.vmin,self.vmax=getRange(self.vars,min_level,max_level,ql,qr)
        data_min,data_max=getRange(self.vars,verbose=False)

        self.ext_1=True if data_min<self.vmin else False
        self.ext_2=True if data_max>self.vmax else False

        #-------------------Get colormap-------------------
        self.cmap=self.getColormap()


        return




    #-----------------Create colormap-----------------
    def getColormap(self):
        if self.cmap is None:
            #cmap=self.plt.cm.bwr
            cmap=DEFAULT_CMAP
        elif self.cmap is not None and type(self.cmap) is str:
            cmpstr='cmap=self.plt.cm.'+self.cmap
            try:
                exec(cmpstr)
            except:
                raise Exception("Color map name wrong.")
        else:
            cmap=self.cmap

        #-------------Shift colormap if needed-------------
        cmap=remappedColorMap(cmap,self.vmin,self.vmax,self.split)

        return cmap


#-------------------Get first time point in variable-------------------
def get1stTimeold(var,verbose=True):
    '''Get first time point in variable

    <var>: nd variable. If is transient variable, slice
           its 1st time point.
           If numpy.ndarray, take its 1st slab along
           the 1st dimension.

    NOTE: if transient variable, var.squeeze() will squeeze
          but with axes info gone. Therefore treat ndarray and
          transient variable differently.
    Update time: 2016-02-05 12:59:55.
    '''
    import numpy
    import MV2 as MV
    import cdms2

    #if type(var)!=cdms2.tvariable.TransientVariable:
    if not cdms2.isVariable(var):
        var=numpy.squeeze(var)
        if numpy.ndim(var)>2:
            time_idx=0
            var=numpy.take(var,[0,],axis=time_idx)
            var=numpy.squeeze(var)

    else:
        var=var(squeeze=1)
        if numpy.ndim(var)>2:
            try:
                time_idx=functions.interpretAxis('time',var)
            except:
                time_idx=0
            var=MV.take(var,[0,],axis=time_idx)
            var=var(squeeze=1)

    return var

#-------------------Get first time point in variable-------------------
def get1stTime(var,verbose=True):
    '''Get first time point in variable

    <var>: nd variable. If is transient variable, slice
           its 1st time point.
           If numpy.ndarray, take its 1st slab along
           the 1st dimension.

    NOTE: if transient variable, var.squeeze() will squeeze
          but with axes info gone. Therefore treat ndarray and
          transient variable differently.
    Update time: 2016-02-05 12:59:55.
    '''

    from .funcs import NCVAR, interpretAxis
    import numpy as np

    if isinstance(var, NCVAR):
        vardata=np.squeeze(var.data)
        if np.ndim(vardata)>2:
            try:
                time_idx=interpretAxis('time', var)
            except:
                time_idx=0
            result=var.sliceData(var.getTime()[0], var.getTime()[1], time_idx)
            result=result.squeeze()

    else:
        result=np.squeeze(var)
        if np.ndim(result)>2:
            time_idx=0
            result=np.take(result,[0,],axis=time_idx)
            result=np.squeeze(result)

    return result


class Plot2D(object):

    ncvar=functions.NCVAR

    def __init__(self, var, method, ax=None, xarray=None, yarray=None,\
            title=None, latlon=True, latlongrid=False, legend='global',
            legend_ori='horizontal', clean=False):
        '''Utility class for 2D plots

        Args:
            var (ndarray): variable to plot. At least 2D.
            method: plotting method, could be an instance of Boxfill, Isofill.
        Keyword Args:
            ax: matplotlib axis obj. If None, create a new.
            xarray (ndarray): 1d array, the array values for the x-axis. If None, use
                              the int indices for the x-dimension.
            yarray (ndarray): 1d array, the array values for the y-axis. If None, use
                              the int indices for the y-dimension.
            title (str): title to plot at subtitle. If None, plot only an alphabetical index.
            latlon (bool): plot lat/lon axis labels or not.
            latlongrid (bool): plot lat/lon grid lines or not.
            legend (str): location of colorbar. Could be: 'global': all subplots share
                          the colorbar of the 1st subplot in figure. or
                          'local': each subplot in figure uses its own colorbar.
            legend_ori (str): 'horizontal' or 'vertical', colorbar orientation.
            clean (bool): if True, omit axis labels, colorbar, subtitle, continents, boundaries etc..
                          Useful to overlay plots.
        '''
        import matplotlib.pyplot as plt

        self.var=var
        self.method=method
        self.ax=ax or plt.subplot(111)
        self.xarray=xarray
        self.yarray=yarray

        self.title=title
        self.latlon=latlon
        self.latlongrid=latlongrid
        self.legend=legend
        self.legend_ori=legend_ori
        self.clean=clean

        #---------------------Get slab---------------------
        self.var=Plot2D.getSlab(self.var)

        #---------------------Get grid---------------------
        self.lonax,self.latax,self.lons,self.lats=self.getGrid()

        #---------------Get geo and fontsize---------------
        self.geo, self.subidx, self.font_size=self.getGeo()




    @classmethod
    def checkBasemap(cls,var,xarray,yarray):
        '''Check variable should be plotted using basemap or not'''
        import numpy

        if isinstance(var, cls.ncvar) and var.getLatitude() is not None\
                and var.getLongitude() is not None:
            return True
        elif isinstance(var, numpy.ndarray) and\
                numpy.ndim(numpy.squeeze(var))>1 and\
                xarray is not None and yarray is not None:
            return True
        else:
            return False





    #-------------------Get first time point in variable-------------------
    @classmethod
    def getSlab(cls, var):
        '''Get a 2D slab from variable

        Args:
            var (ndarray): nd variable. If is transient variable, try to slice
                   its 1st time point.
                   If numpy.ndarray, try to take a slab from its last 2 dimensions.
        Returns:
            ndarray: 2D slab from input <var> to plot.

        Note:
            if transient variable, var.squeeze() will squeeze but with axes
            info gone. Therefore treat ndarray and transient variable
            differently.

        Update time: 2016-02-05 12:59:55.
        '''
        import numpy

        #------------------Numpy ndarray------------------
        if not isinstance(var, functions.NCVAR):
            var=numpy.squeeze(var)
            if numpy.ndim(var)==2:
                return var

            if numpy.ndim(var)>2:
                return functions.getSlab(var)

        #----------------NCVAR variable----------------
        else:
            var=var.squeeze()
            if numpy.ndim(var.data)==2:
                return var.data

            if numpy.ndim(var.data)>2:
                yidx=-1; xidx=-2
                return functions.getSlab(var.data,yidx,xidx)


    #---------------------Get grid---------------------
    def getGrid(self):
        '''Get lat/lon grid info from data'''
        import numpy

        if self.yarray is None:
            latax=numpy.arange(self.var.shape[0])
        else:
            latax=numpy.array(self.yarray)

        if self.xarray is None:
            lonax=numpy.arange(self.var.shape[1])
        else:
            lonax=numpy.array(self.xarray)

        if len(latax)!=self.var.shape[0]:
            raise Exception("X-axis dimention does not match")
        if len(lonax)!=self.var.shape[1]:
            raise Exception("Y-axis dimention does not match")

        lons,lats=numpy.meshgrid(lonax,latax)

        return lonax,latax,lons,lats


    def getGeo(self):

        geo=self.ax.get_geometry()[:2]
        subidx=self.ax.get_geometry()[-1]
        scale=1./max(geo)
        font_size=7*scale+8

        return geo, subidx, font_size


    #-------------Plot according to method-------------
    def plot(self):

        self.cs=self._plot()
        self.plotAxes()
        self.cbar=self.plotColorbar()
        if self.cbar is not None:
            var_units=getattr(self.var,'units','')
            if var_units:
                # option1: plot colorbar units below
                #self.cbar.set_label(var_units,fontsize=self.font_size)

                # option2: plot colorbar units to the right or below, depending on
                # orientation
                if self.legend_ori=='horizontal':
                    cbar_ax=self.cbar.ax
                    cbar_ax.text(1.07, 0.5, var_units, fontsize=self.font_size,
                            transform=cbar_ax.transAxes,
                            horizontalalignment='left',
                            verticalalignment='center')
                elif self.legend_ori=='vertical':
                    cbar_ax=self.cbar.ax
                    cbar_ax.text(0.5, -0.15, var_units, fontsize=self.font_size,
                            transform=cbar_ax.transAxes,
                            horizontalalignment='center',
                            verticalalignment='top')

        self.plotTitle()

        return self.cs


    def _plot(self):

        #-------------------Contour fill/line-------------------
        if self.method.method in ['isofill','isoline']:

            if self.method.ext_1 is False and self.method.ext_2 is False:
                extend='neither'
            elif self.method.ext_1 is True and self.method.ext_2 is False:
                extend='min'
            elif self.method.ext_1 is False and self.method.ext_2 is True:
                extend='max'
            else:
                extend='both'

            if self.method.method=='isofill':
                # make masked value grey, otherwise they will be white
                self.ax.patch.set_color('0.5')
                cs=self.ax.contourf(self.lons, self.lats, self.var,
                        self.method.levels, cmap=self.method.cmap,
                        extend=extend,hatch='/')

            elif self.method.method=='isoline':
                if self.method.black:
                    colors=['k']*len(self.method.levels)
                    cs=self.ax.contour(self.lons, self.lats, self.var,
                            self.method.levels, colors=colors,extend=extend,
                            linewidth=self.method.linewidth,
                            alpha=self.method.alpha)
                else:
                    # Plot as colored lines
                    cs=self.ax.contour(self.lons, self.lats, self.var,
                            self.method.levels, cmap=self.method.cmap,
                            extend=extend, linewidth=self.method.linewidth,
                            alpha=self.method.alpha)

                #-----------------Set line styles-----------------
                # for some reason it's not giving me dashed line for negatives.
                # have to set myself
                for ii in range(len(cs.collections)):
                    cii=cs.collections[ii]
                    lii=cs.levels[ii]
                    if lii<0:
                        if self.method.dash_negative:
                            cii.set_linestyle('dashed')
                        else:
                            cii.set_linestyle('solid')

                # For some reason the linewidth keyword in contour doesn't
                # work, has to set again.
                for ii in range(len(cs.collections)):
                    cs.collections[ii].set_linewidth(self.method.linewidth)

                #----------------Thicken some lines----------------
                if self.method.bold_lines is not None:
                    multi=2.0
                    idx_bold=[]
                    for bii in self.method.bold_lines:
                        idxii=numpy.where(numpy.array(self.method.levels)==bii)[0]
                        if len(idxii)>0:
                            idx_bold.append(int(idxii))

                    for bii in idx_bold:
                        cs.collections[bii].set_linewidth(self.method.linewidth*multi)

        #---------------------Boxfill---------------------
        elif self.method.method == 'boxfill':

            self.ax.patch.set_color('0.5')
            cs=self.ax.imshow(self.var,cmap=self.method.cmap,origin='lower',\
                    vmin=self.method.vmin,vmax=self.method.vmax,
                    interpolation='nearest',
                    extent=[self.lonax.min(),self.lonax.max(),\
                            self.latax.min(),self.latax.max()],
                    aspect='auto')

        #-------------------Pcolor fill-------------------
        elif self.method.method == 'pcolor':

            self.ax.patch.set_color('0.5')
            cs=self.ax.pcolormesh(self.lons,self.lats,self.var,\
                    cmap=self.method.cmap,vmin=self.method.levels[0],\
                    vmax=self.method.levels[-1])

        #------------------Hatch contourf------------------
        elif self.method.method=='hatch':
            # Skip if none == 1
            if numpy.all(self.var==0):
                nlevel=1
            else:
                nlevel=3
            cs=self.ax.contourf(self.lons[0,:],self.lats[:,0],self.var,nlevel,\
                    colors='none',hatches=[None,self.method.hatch],
                    alpha=0.)

        return cs



    #---------------Draw lat, lon grids---------------
    def plotAxes(self):

        import numpy

        n_lat=int(6*1./self.geo[0]+3)
        n_lon=int(6*1./self.geo[1]+4)

        lat_labels=numpy.array(mkscale(self.latax[0],self.latax[-1],n_lat,1))
        idx=numpy.where((lat_labels>=numpy.min(self.latax)) & (lat_labels<=numpy.max(self.latax)))
        lat_labels=numpy.array(lat_labels)[idx]

        lon_labels=numpy.array(mkscale(self.lonax[0],self.lonax[-1],n_lon,1))
        idx=numpy.where((lon_labels>=numpy.min(self.lonax)) & (lon_labels<=numpy.max(self.lonax)))
        lon_labels=numpy.array(lon_labels)[idx]

        self.ax.set_xticks(lon_labels)
        self.ax.set_yticks(lat_labels)
        self.ax.tick_params(axis='both',which='major',labelsize=self.font_size)

        if self.latlongrid:
            self.ax.grid(True)

        #--------Turn off lat/lon labels if required--------
        if self.latlon is False or self.clean:
            self.ax.xaxis.set_ticklabels([])
            self.ax.yaxis.set_ticklabels([])

        return


    #------------------Draw color bar------------------
    def plotColorbar(self):
        import matplotlib.pyplot as plt
        import matplotlib.colorbar as mcbar

        if self.method.method=='hatch':
            return

        if self.legend is None:
            return

        if self.method.method in ['isofill','isoline']:
            isdrawedges=True
        else:
            isdrawedges=False

        if self.method.method in ['boxfill','pcolor']:
            if self.method.ext_1 is False and self.method.ext_2 is False:
                extend='neither'
            elif self.method.ext_1 is True and self.method.ext_2 is False:
                extend='min'
            elif self.method.ext_1 is False and self.method.ext_2 is True:
                extend='max'
            else:
                extend='both'

        #---------------Draw local colorbar---------------
        if self.legend=='local':
            ticks=getattr(self.method,'levels',None)
            if self.legend_ori=='horizontal':
                pad=0.18
            elif self.legend_ori=='vertical':
                pad=0.06

            if ticks is None:
                if self.method.method in ['boxfill','pcolor']:
                    #--------------Create a colorbar axis--------------
                    cax,kw=mcbar.make_axes_gridspec(self.ax,orientation=self.legend_ori,
                            shrink=0.75,pad=pad,fraction=0.07)
                    cbar=plt.colorbar(self.cs,cax=cax,orientation=self.legend_ori,
                           drawedges=isdrawedges,extend=extend)
                    #cbar=plt.colorbar(self.cs,ax=self.ax,orientation=self.legend_ori,\
                            #shrink=0.7,ticks=ticks,drawedges=isdrawedges,pad=0.15,\
                            #extend=extend)
                else:
                    cax,kw=mcbar.make_axes_gridspec(self.ax,orientation=self.legend_ori,
                            shrink=0.75,pad=pad,fraction=0.07)
                    cbar=plt.colorbar(self.cs,cax=cax,orientation=self.legend_ori,\
                                    ticks=ticks,drawedges=isdrawedges)
                    #cbar=plt.colorbar(self.cs,ax=self.ax,orientation=self.legend_ori,\
                            #shrink=0.7,ticks=ticks,drawedges=isdrawedges,pad=0.15)

                cbar.ax.tick_params(labelsize=self.font_size)
            else:
                if all([ll==int(ll) for ll in ticks]):
                    tick_format='%d'
                    ticks=[int(ll) for ll in ticks]
                else:
                    tick_format=None

                ltop=ticks[::2]       #labels on top
                lbot=ticks[1:][::2]   #labels at bottom

                if self.method.method in ['boxfill','pcolor']:
                    cax,kw=mcbar.make_axes_gridspec(self.ax,orientation=self.legend_ori,
                            shrink=0.75,pad=pad,fraction=0.07)
                    cbar=plt.colorbar(self.cs,cax=cax,orientation=self.legend_ori,
                           ticks=lbot,drawedges=isdrawedges,extend=extend,
                           format=tick_format)
                    #cbar=plt.colorbar(self.cs,ax=self.ax,orientation=self.legend_ori,\
                            #shrink=0.8,ticks=lbot,drawedges=isdrawedges,pad=0.15,\
                            #extend=extend)
                else:
                    cax,kw=mcbar.make_axes_gridspec(self.ax,orientation=self.legend_ori,
                            shrink=0.75,pad=pad,fraction=0.07)
                    cbar=plt.colorbar(self.cs,cax=cax,ax=self.ax,orientation=self.legend_ori,
                           ticks=lbot,drawedges=isdrawedges,
                           format=tick_format)
                    #cbar=plt.colorbar(self.cs,ax=self.ax,orientation=self.legend_ori,\
                            #shrink=0.8,ticks=lbot,drawedges=isdrawedges,pad=0.15)

                cbar.ax.tick_params(labelsize=self.font_size)

                #-------------Print bottom tick labels-------------
                cbar.ax.set_xticklabels(lbot)

                #--------------Print top tick labels--------------
                vmin=cbar.norm.vmin
                vmax=cbar.norm.vmax

                #---------Compute top line tick locations---------
                # need to take into account the overflow triangle
                # by default the triangle is 5% of axis length
                if cbar.extend=='min':
                    shift_l=0.05
                    scaling=0.95
                elif cbar.extend=='max':
                    shift_l=0.
                    scaling=0.95
                elif cbar.extend=='both':
                    shift_l=0.05
                    scaling=0.90
                else:
                    shift_l=0.
                    scaling=1.

                if self.legend_ori=='horizontal':
                    for ii,tii in enumerate(ltop):
                        cbar.ax.text(shift_l+scaling*(tii-vmin)/(vmax-vmin),
                                1.3, str(tii),\
                                transform=cbar.ax.transAxes,va='bottom',
                                ha='center',fontsize=self.font_size)
                elif self.legend_ori=='vertical':
                    for ii in ltop:
                        cbar.ax.text(0.02,shift_l+scaling*(ii-vmin)/(vmax-vmin),
                                str(ii),transform=cbar.ax.transAxes,va='center',
                                ha='right',fontsize=self.font_size)
            return cbar


        #-----Use the 1st subplot as global color bar-----
        elif self.legend=='global' and self.subidx==1:
            ticks=getattr(self.method,'levels',None)
            if self.legend_ori=='horizontal':
                cax=self.ax.get_figure().add_axes([0.18,0.02,0.65,0.015])
            else:
                cax=self.ax.get_figure().add_axes([0.93,0.10,0.02,0.7])

            if ticks is None:
                if self.method.method in ['boxfill','pcolor']:
                    cbar=plt.colorbar(self.cs,cax=cax,orientation=self.legend_ori,\
                            ticks=ticks,drawedges=isdrawedges,extend=extend)
                else:
                    cbar=plt.colorbar(self.cs,cax=cax,orientation=self.legend_ori,\
                            ticks=ticks,drawedges=isdrawedges)
                cbar.ax.tick_params(labelsize=self.font_size)
            else:

                if all([ll==int(ll) for ll in ticks]):
                    tick_format='%d'
                    ticks=[int(ll) for ll in ticks]
                else:
                    tick_format=None

                ltop=ticks[::2]       #labels on top
                lbot=ticks[1:][::2]   #labels at bottom

                if self.method.method in ['boxfill','pcolor']:
                    cbar=plt.colorbar(self.cs,cax=cax,orientation=self.legend_ori,\
                        ticks=lbot,drawedges=isdrawedges,extend=extend,
                        format=tick_format)
                else:
                    cbar=plt.colorbar(self.cs,cax=cax,orientation=self.legend_ori,\
                        ticks=lbot,drawedges=isdrawedges,
                        format=tick_format)
                cbar.ax.tick_params(labelsize=self.font_size)

                #-------------Print bottom tick labels-------------
                cbar.ax.set_xticklabels(lbot)

                #--------------Print top tick labels--------------
                vmin=cbar.norm.vmin
                vmax=cbar.norm.vmax
                '''
                for ii in ltop:
                    cbar.ax.text((ii-vmin)/(vmax-vmin), 1.3, str(ii),\
                        transform=cbar.ax.transAxes,va='bottom',
                        ha='center',fontsize=self.font_size)
                '''

                #---------Compute top line tick locations---------
                # need to take into account the overflow triangle
                # by default the triangle is 5% of axis length
                if cbar.extend=='min':
                    shift_l=0.05
                    scaling=0.95
                elif cbar.extend=='max':
                    shift_l=0.
                    scaling=0.95
                elif cbar.extend=='both':
                    shift_l=0.05
                    scaling=0.90
                else:
                    shift_l=0.
                    scaling=1.

                if self.legend_ori=='horizontal':
                    for ii,tii in enumerate(ltop):
                        cbar.ax.text(shift_l+scaling*(tii-vmin)/(vmax-vmin),
                                1.3, str(tii),\
                                transform=cbar.ax.transAxes,va='bottom',
                                ha='center',fontsize=self.font_size)
                elif self.legend_ori=='vertical':
                    for ii in ltop:
                        cbar.ax.text(0.02,shift_l+scaling*(ii-vmin)/(vmax-vmin),
                                str(ii),transform=cbar.ax.transAxes,va='center',
                                ha='right',fontsize=self.font_size)

            #var_units=getattr(self.var,'units','')
            #cbar.set_label(var_units,fontsize=self.font_size)

            return cbar
        else:
            return




    #-------------Plot subplot title/index-------------
    def plotTitle(self):

        if self.clean:
            return

        if self.title is not None and type(self.title) is str:

            if self.geo[0]*self.geo[1]>1:
                # force indexing
                import re
                rep=re.compile('^\\((.*?)\\)(.*)')
                match=rep.match(self.title)
                # overwrite indexing if title starts with (a), (b) or so
                if match is not None:
                    tidx_str='(%s)' %match.group(1)
                    tstr=match.group(2).strip()
                else:
                    tidx_str=index2Letter(self.subidx)
                    tstr=self.title
                self.ax.set_title('%s %s' %(tidx_str,tstr),loc='left',\
                        fontsize=self.font_size)
            else:
                self.ax.set_title('%s' %self.title,loc='left',\
                        fontsize=self.font_size)
        elif self.title is None and self.geo[0]*self.geo[1]>1:
            self.ax.set_title('%s' %index2Letter(self.subidx),loc='left',\
                    fontsize=self.font_size)

        return


class Plot2Basemap(Plot2D):
    def __init__(self,var,method,ax=None,legend='global',title=None,\
            xarray=None,yarray=None,\
            latlon=True,latlongrid=False,projection='merc',fill_color='0.8',
            legend_ori='horizontal',clean=False,fix_aspect=False):
        '''Utility class for 2D geographical plots using basemap

        Args:
            var (ndarray): variable to plot. At least 2D.
            method: plotting method, could be an instance of Boxfill, Isofill.
        Keyword Args:
            ax: matplotlib axis obj. If None, create a new.
            legend (str): location of colorbar. Could be: 'global': all subplots share
                          the colorbar of the 1st subplot in figure. or
                          'local': each subplot in figure uses its own colorbar.
            title (str): title to plot at subtitle. If None, plot only an alphabetical index.
            xarray (ndarray): 1d array, the array values for the x-axis. If None, use
                              the int indices for the x-dimension.
            yarray (ndarray): 1d array, the array values for the y-axis. If None, use
                              the int indices for the y-dimension.
            latlon (bool): plot lat/lon axis labels or not.
            latlongrid (bool): plot lat/lon grid lines or not.
            projection (str): map projection, used when plotting with basemap.
            fill_color: color to fill continent or masked regions.
            legend_ori (str): 'horizontal' or 'vertical', colorbar orientation.
            clean (bool): if True, omit axis labels, colorbar, subtitle, continents, boundaries etc..
                          Useful to overlay plots.
            fix_aspect (bool): passed to the basemap plotting function (e.g. contourf())
                               for control of aspect ratio.
        '''

        Plot2D.__init__(self,var,method,ax=ax,\
                xarray=xarray,yarray=yarray,\
                title=title,latlon=latlon,\
                latlongrid=latlongrid,legend=legend,
                legend_ori=legend_ori,clean=clean)

        self.projection=projection
        self.fill_color=fill_color
        self.fix_aspect=fix_aspect

        #------------------Check basemap------------------
        if not Plot2D.checkBasemap(var,self.xarray,self.yarray):
            raise Exception("<var> not suitable for basemap plot.")


    def getGrid(self):

        try:
            latax=self.var.getLatitude()[:]
            lonax=self.var.getLongitude()[:]
            lons,lats=numpy.meshgrid(lonax,latax)
        except:
            lonax,latax,lons,lats=super(Plot2Basemap,self).getGrid()

        return lonax,latax,lons,lats


    def _plot(self):

        from mpl_toolkits.basemap import Basemap
        from mpl_toolkits.basemap import addcyclic

        #------------------Create basemap------------------
        if self.method.method=='gis':
            self.projection='cyl'

        if self.projection in ['cyl','merc','cea']:
            bmap=Basemap(projection=self.projection,\
                    llcrnrlat=self.latax[0],llcrnrlon=self.lonax[0],\
                    urcrnrlat=self.latax[-1],urcrnrlon=self.lonax[-1],\
                    ax=self.ax,fix_aspect=self.fix_aspect)

        elif self.projection in ['npaeqd', 'nplaea', 'npstere']:

            self.var,self.lonax=addcyclic(self.var,self.lonax)
            self.lons,self.lats=numpy.meshgrid(self.lonax,self.latax)
            lat_0=numpy.min(self.latax)-5
            lon_0=180.

            bmap=Basemap(projection=self.projection,\
                    boundinglat=lat_0,lon_0=lon_0,
                    ax=self.ax,fix_aspect=self.fix_aspect)

        elif self.projection in ['spaeqd', 'splaea', 'spstere']:

            self.var,self.lonax=addcyclic(self.var,self.lonax)
            self.lons,self.lats=numpy.meshgrid(self.lonax,self.latax)
            lat_0=numpy.max(self.latax)+5
            lon_0=180.
            bmap=Basemap(projection=self.projection,\
                    boundinglat=lat_0,lon_0=lon_0,
                    ax=self.ax,fix_aspect=self.fix_aspect)

        self.bmap=bmap

        #-------------Plot according to method-------------
        #-------------------Contour fill/line-------------------
        if self.method.method in ['isofill','isoline']:

            if self.method.ext_1 is False and self.method.ext_2 is False:
                extend='neither'
            elif self.method.ext_1 is True and self.method.ext_2 is False:
                extend='min'
            elif self.method.ext_1 is False and self.method.ext_2 is True:
                extend='max'
            else:
                extend='both'

            if self.method.method=='isofill':
                cs=bmap.contourf(self.lons,self.lats,self.var,self.method.levels,latlon=True,\
                        cmap=self.method.cmap,ax=self.ax,extend=extend)

            elif self.method.method=='isoline':
                if self.method.black:
                    colors=['k']*len(self.method.levels)
                    cs=bmap.contour(self.lons,self.lats,self.var,self.method.levels,latlon=True,\
                            colors=colors,ax=self.ax,extend=extend,
                            linewidth=self.method.linewidth,
                            alpha=self.method.alpha)
                else:
                    cs=bmap.contour(self.lons,self.lats,self.var,self.method.levels,latlon=True,\
                            cmap=self.method.cmap,ax=self.ax,extend=extend,
                            alpha=self.method.alpha)

                #-----------------Set line styles-----------------
                # for some reason it's not giving me dashed line for negatives.
                # have to set myself
                for ii in range(len(cs.collections)):
                    cii=cs.collections[ii]
                    lii=cs.levels[ii]
                    if lii<0:
                        if self.method.dash_negative:
                            cii.set_linestyle('dashed')
                        else:
                            cii.set_linestyle('solid')

                # For some reason the linewidth keyword in contour doesn't
                # work, has to set again.
                for ii in range(len(cs.collections)):
                    cs.collections[ii].set_linewidth(self.method.linewidth)

                if self.method.bold_lines is not None:
                    multi=2.0
                    idx_bold=[]
                    for bii in self.method.bold_lines:
                        idxii=numpy.where(numpy.array(self.method.levels)==bii)[0]
                        if len(idxii)>0:
                            idx_bold.append(int(idxii))

                    for bii in idx_bold:
                        cs.collections[bii].set_linewidth(self.method.linewidth*multi)


        #---------------------Boxfill---------------------
        elif self.method.method == 'boxfill':

            cs=bmap.imshow(self.var,cmap=self.method.cmap,ax=self.ax,\
                    vmin=self.method.vmin,vmax=self.method.vmax,interpolation='nearest')

        #-------------------Pcolor fill-------------------
        elif self.method.method == 'pcolor':

            cs=bmap.pcolormesh(self.lons,self.lats,self.var,latlon=True,\
                    cmap=self.method.cmap,ax=self.ax,vmin=self.method.levels[0],\
                    vmax=self.method.levels[-1])

        #------------------Hatch contourf------------------
        elif self.method.method=='hatch':
            # Skip if none == 1
            if numpy.all(self.var==0):
                nlevel=1
            else:
                nlevel=3
            cs=bmap.contourf(self.lons,self.lats,self.var,nlevel,latlon=True,\
                    colors='none',ax=self.ax,hatches=[None,self.method.hatch],
                    alpha=0.)

        elif self.method.method=='gis':
            cs=bmap.arcgisimage(service='ESRI_Imagery_World_2D',xpixels=self.method.xpixels,
                    dpi=self.method.dpi,verbose=self.method.verbose)

        self.plotOthers()

        return cs



    def plotOthers(self):

        if self.clean:
            return

        #-------------------Draw others-------------------
        if not self.method.method=='gis':
            self.bmap.drawcoastlines(linewidth=0.5,linestyle='solid',color='k',\
                antialiased=True)
            self.bmap.drawcountries(linewidth=0.5,linestyle='solid',color='k',\
                    antialiased=True)
            #self.bmap.fillcontinents(color='w',lake_color=None,alpha=0.2)
            self.bmap.drawmapboundary(color='k',linewidth=1.0,fill_color=self.fill_color)
            #self.bmap.drawrivers(linewidth=0.5,linestyle='solid',color='b',\
                    #antialiased=True)

        return


    #---------------Draw lat, lon grids---------------
    def plotAxes(self):

        if self.clean:
            return

        def getLabelBool(geo,idx):
            if geo[0]*geo[1]==1:
                parallels=[1,1,0,0]
                meridians=[0,0,0,1]
            else:
                ridx,cidx=numpy.unravel_index(idx-1,geo)
                if cidx==0:
                    parallels=[1,0,0,0]
                elif cidx==geo[1]-1:
                    parallels=[0,1,0,0]
                else:
                    parallels=[0,0,0,0]
                if ridx==0:
                    meridians=[0,0,0,0]
                elif ridx==geo[0]-1:
                    meridians=[0,0,0,1]
                else:
                    meridians=[0,0,0,0]
            return parallels, meridians

        n_lat=int(6*1./self.geo[0]+3)
        n_lon=int(6*1./self.geo[1]+4)

        lat_labels=numpy.array(mkscale(self.latax[0],self.latax[-1],n_lat,1))
        idx=numpy.where((lat_labels>=self.latax[0]) & (lat_labels<=self.latax[-1]))
        lat_labels=numpy.array(lat_labels)[idx]

        lon_labels=numpy.array(mkscale(self.lonax[0],self.lonax[-1],n_lon,1))
        idx=numpy.where((lon_labels>=self.lonax[0]) & (lon_labels<=self.lonax[-1]))
        lon_labels=numpy.array(lon_labels)[idx]

        #-----------------Set axes labels-----------------
        if self.latlon is False or self.clean:
            parallels=[0,0,0,0]
            meridians=[0,0,0,0]
            #self.ax.xaxis.set_ticklabels([]) Doesn't work
        elif self.latlon=='all':
            parallels=[1,1,0,0]
            meridians=[0,0,0,1]
        else:
            parallels,meridians=getLabelBool(self.geo,self.subidx)

        #------------------Set grid lines------------------
        if self.latlongrid is False or self.clean:
            #linewidth=0
            zorder=-2
        else:
            #linewidth=0.5
            zorder=None

        #-----------------Draw axes/lables/ticks-----------------
        self.ax.set_yticks(lat_labels)
        self.ax.set_xticks(lon_labels)
        self.ax.tick_params(axis='both',which='major',labelsize=self.font_size)
        self.ax.xaxis.set_ticklabels([])
        self.ax.yaxis.set_ticklabels([])
        self.bmap.drawparallels(lat_labels,labels=parallels,
                zorder=zorder,
                #linewidth=linewidth,\
                labelstyle='+/-',fontsize=self.font_size)
        labels=self.bmap.drawmeridians(lon_labels,labels=meridians,
                zorder=zorder,
                #linewidth=linewidth,\
                labelstyle='+/-',fontsize=self.font_size)

        #--------------Fix label offset issue--------------
        y0, y1=self.ax.get_ylim()
        h=y1-y0
        yoffset=0.03
        for key, (lines,texts) in labels.items():
            for text in texts:
                x,y = text.get_position()
                text.set_position((x, y0-yoffset*h))

        #-------Label parallels in polar projections-------
        if self.projection in ['npaeqd', 'nplaea', 'npstere', 'spaeqd',
                'splaea', 'spstere']:
            for ll in lat_labels:
                xll,yll=self.bmap(180.,ll)
                self.ax.text(xll,yll,u'%+d\N{DEGREE SIGN}' %ll,
                        fontsize=max(4,int(self.font_size*0.7)),
                        horizontalalignment='center',
                        verticalalignment='center')

        return



    def plotColorbar(self):

        cbar=super(Plot2Basemap,self).plotColorbar()
        #if cbar is not None:
            #var_units=getattr(self.var,'units','')
            #cbar.set_label(var_units,fontsize=self.font_size)

        return cbar


def plot2(var,method,ax=None,legend='global',\
        xarray=None,yarray=None,\
        title=None,latlon=True,latlongrid=False,fill_color='0.8',
        projection='merc',legend_ori='horizontal',clean=False,
        isbasemap=True,
        fix_aspect=True,verbose=True):
    '''A helper function for quickly create 2D plots

    Args:
        var (ndarray): variable to plot. At least 2D.
        method: plotting method, could be an instance of Boxfill, Isofill.
    Keyword Args:
        ax: matplotlib axis obj. If None, create a new.
        xarray (ndarray): 1d array, the array values for the x-axis. If None, use
                          the int indices for the x-dimension.
        yarray (ndarray): 1d array, the array values for the y-axis. If None, use
                          the int indices for the y-dimension.
        title (str): title to plot at subtitle. If None, plot only an alphabetical index.
        latlon (bool): plot lat/lon axis labels or not.
        latlongrid (bool): plot lat/lon grid lines or not.
        fill_color: color to fill continent or masked regions.
        projection (str): map projection, used when plotting with basemap.
        legend (str): location of colorbar. Could be: 'global': all subplots share
                      the colorbar of the 1st subplot in figure. or
                      'local': each subplot in figure uses its own colorbar.
        legend_ori (str): 'horizontal' or 'vertical', colorbar orientation.
        clean (bool): if True, omit axis labels, colorbar, subtitle, continents, boundaries etc..
                      Useful to overlay plots.
        isbasemap (bool): plot using basemap or not. Usually used to force plot as a normal
                          2d plot instead of geographical plot using basemap.
        fix_aspect (bool): passed to the basemap plotting function (e.g. contourf())
                           for control of aspect ratio.
    '''

    if numpy.ndim(var)==1:
        raise Exception("<var> is 1D")

    if isbasemap and Plot2D.checkBasemap(var,xarray,yarray):
        try:
            var=functions.increasingLatitude(var)
        except:
            pass
        plotobj=Plot2Basemap(var,method,ax=ax,legend=legend,\
                xarray=xarray,yarray=yarray,\
                title=title,latlon=latlon,latlongrid=latlongrid,\
                fill_color=fill_color,projection=projection,
                legend_ori=legend_ori,clean=clean,fix_aspect=fix_aspect)
    else:
        plotobj=Plot2D(var,method,ax=ax,legend=legend,\
                xarray=xarray,yarray=yarray,\
                title=title,latlon=latlon,latlongrid=latlongrid,
                legend_ori=legend_ori,clean=clean)
    cs=plotobj.plot()

    return plotobj




