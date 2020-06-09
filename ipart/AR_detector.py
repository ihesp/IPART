'''AR detection functions.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-06-05 23:21:24.
'''

from __future__ import print_function
import numpy as np
import pandas as pd
import networkx as nx
from skimage import measure
from skimage import morphology
from scipy import ndimage
import matplotlib.pyplot as plt

import cdms2 as cdms
import MV2 as MV
from genutil import statistics as stats
import cdutil

from ipart.utils import rdp
from ipart.utils import funcs
from ipart.utils import peak_prominence2d as pp2d

NX_VERSION=nx.__version__[0]



def plotGraph(graph, ax=None, show=True):
    '''Helper func to plot the graph of an AR coordinates
    '''

    if ax is None:
        fig=plt.figure()
        ax=fig.add_subplot(111)

    pos=[(ii[1],ii[0]) for ii in graph.nodes()]  # x,y
    pos_dict=dict(zip(graph.nodes(),pos))
    nx.draw(graph,ax=ax,pos=pos_dict,node_size=15,node_color='darkgray',
            edge_color='dimgray')
    if show:
        plt.show(block=False)
    return


def areaFilt(mask, area, min_area=None, max_area=None):
    '''Filter AR binary masks by region areas

    Args:
        mask (ndarray): 2D binary mask with detected objects shown as 1s.
        area (ndarray): 2D map showing grid cell areas in km^2.
    Keyword Args:
        min_area (float or None): if not None, minimum area to filter objects
                                  in <mask>.
        max_area (float or None): if not None, maximum area to filter objects
                                  in <mask>.

    Returns:
        result (ndarray): 2D binary mask with objects area-filtered.
    '''

    if min_area is None and max_area is None:
        return mask

    labels=measure.label(mask,connectivity=1)
    n=labels.max()+1
    areas=ndimage.sum(area,labels,np.arange(n))
    sel=np.ones(n,bool)
    if min_area is not None:
        sel=np.where(areas<min_area,0,sel)
    if max_area is not None:
        sel=np.where(areas>max_area,0,sel)
    # remove background area
    sel[0]=0
    result=sel[labels]

    return result


def spherical2Cart(lat, lon):
    clat=(90-lat)*np.pi/180.
    lon=lon*np.pi/180.
    x=np.cos(lon)*np.sin(clat)
    y=np.sin(lon)*np.sin(clat)
    z=np.cos(clat)

    return np.array([x,y,z])


def cart2Spherical(x, y, z, shift_lon):

    r=np.sqrt(x**2+y**2+z**2)
    clat=np.arccos(z/r)/np.pi*180
    lat=90.-clat
    lon=np.arctan2(y,x)/np.pi*180
    lon=(lon+360)%360

    lon=np.where((lon>=0) & (lon<shift_lon),lon+360,lon)

    return np.array([lat,lon,np.ones(lat.shape)])


def computeTheta(p1, p2):
    r'''Tangent line to the arc \|p1-p2\|

    Args:
        p1,p2 (float): (lat,lon) coordinates
    '''
    p1=spherical2Cart(p1[0],p1[1])
    p2=spherical2Cart(p2[0],p2[1])
    theta=p2-np.dot(p1,p2)*p1
    norm=np.linalg.norm(theta)
    if norm>0:
        theta=theta/norm
    return theta


def wind2Cart(u, v, lats, lons):
    '''Convert u,v winds to Cartesian, consistent with spherical2Cart.
    '''

    latsr=lats*np.pi/180
    lonsr=lons*np.pi/180
    vh=v*np.sin(latsr)

    ux=-u*np.sin(lonsr) - vh*np.cos(lonsr)
    uy=u*np.cos(lonsr) - vh*np.sin(lonsr)
    uz=v*np.cos(latsr)
    vs=np.array([ux,uy,uz])

    return vs


def cart2Wind(vs, lats, lons):
    '''Convert winds in Cartesian to u,v, inverse to wind2Cart.
    '''

    latsr=lats*np.pi/180
    lonsr=lons*np.pi/180

    #ux=vs[0]
    uy=vs[1]
    uz=vs[2]

    u=uy/np.cos(lonsr) + uz*np.tan(latsr)*np.tan(lonsr)
    v=uz/np.cos(latsr)

    return u,v


def maskToGraph(mask, quslab, qvslab, costhetas, sinthetas, edge_eps,
        connectivity=2):
    '''Create graph from AR mask

    Args:
        mask (ndarray): 2D binary map showing the location of an AR with 1s.
        quslab (cdms.TransientVariable): 2D map of u-flux.
        qvslab (cdms.TransientVariable): 2D map of v-flux.
        costhetas (cdms.TransientVariable): (n * m) 2D slab of grid cell shape:
                                          cos=dx/sqrt(dx^2+dy^2).
        sinthetas (cdms.TransientVariable): (n * m) 2D slab of grid cell shape:
                                          sin=dy/sqrt(dx^2+dy^2).
        edge_eps (float): float in (0,1), minimal proportion of flux component
                          in a direction to total flux to allow edge building
                          in that direction. Defined in Global preamble.
        connectivity (int): 1 or 2. 4- or 8- connectivity in defining neighbor-
                            hood relationship in a 2D square grid.
    Returns:
        g (networkx.DiGraph): directed planar graph constructed from AR mask
                              and flows.
    '''


    quslab=np.array(quslab)
    qvslab=np.array(qvslab)
    wsslab=np.sqrt(quslab**2+qvslab**2)

    g=nx.DiGraph()

    # 1 connectivity edges
    # the list approach
    '''
    y,x=np.where(mask)
    zipcor=zip(y,x)
    right=[(yi,xi) for yi,xi in zipcor if (yi,xi+1) in zipcor]
    left=[(yi,xi) for yi,xi in zipcor if (yi,xi-1) in zipcor]
    up=[(yi,xi) for yi,xi in zipcor if (yi+1,xi) in zipcor]
    down=[(yi,xi) for yi,xi in zipcor if (yi-1,xi) in zipcor]

    # nodes to the right/left/up/down
    right0=[(yi,xi+1) for yi,xi in right]
    left0=[(yi,xi-1) for yi,xi in left]
    up0=[(yi+1,xi) for yi,xi in up]
    down0=[(yi-1,xi) for yi,xi in down]
    '''

    # the shifting approach
    right=np.roll(mask, -1, axis=1)*mask
    left=np.roll(mask, 1, axis=1)*mask
    up=np.roll(mask, -1, axis=0)*mask
    down=np.roll(mask, 1, axis=0)*mask

    def addWeightedEdges2(nodes1,speedslab,d):
        '''Add directed edges to graph. For shifting approach
        '''

        ratio=np.where(np.isclose(wsslab, 0., atol=1e-5), 0, speedslab/wsslab)
        idx=np.where(ratio>=edge_eps, 1, 0)*nodes1
        idx=zip(*np.where(idx>0))

        for ii, (yii,xii) in enumerate(idx):
            # nii: start, nii2: end
            yii=int(yii)
            xii=int(xii)
            nii=(yii,xii)
            if d=='r':
                nii2=(yii,xii+1)
            elif d=='l':
                nii2=(yii,xii-1)
            elif d=='u':
                nii2=(yii+1,xii)
            elif d=='d':
                nii2=(yii-1,xii)
            elif d=='tr':
                nii2=(yii+1,xii+1)
            elif d=='br':
                nii2=(yii-1,xii+1)
            elif d=='tl':
                nii2=(yii+1,xii-1)
            elif d=='bl':
                nii2=(yii-1,xii-1)

            meanivt=speedslab[yii,xii]

            g.add_edge(nii,nii2,
                    weight=np.exp(-meanivt/1e2),
                    ivt=meanivt)


    def addWeightedEdges(nodes1,nodes2,speedslab):
        '''Add directed edges to graph. For the list approach
        '''

        # nii: start, nii2: end
        for nii,nii2 in zip(nodes1,nodes2):
            if speedslab[nii]/wsslab[nii]>=edge_eps:
                meanivt=speedslab[nii]

                g.add_edge(nii,nii2,
                        weight=np.exp(-meanivt/1e2),
                        ivt=meanivt)

    # add 1 connectivity edges
    #addWeightedEdges(right,right0,quslab)
    #addWeightedEdges(left,left0,-quslab)
    #addWeightedEdges(up,up0,qvslab)
    #addWeightedEdges(down,down0,-qvslab)
    addWeightedEdges2(right,quslab,'r')
    addWeightedEdges2(left,-quslab,'l')
    addWeightedEdges2(up,qvslab,'u')
    addWeightedEdges2(down,-qvslab,'d')

    # 2 connectivity edges
    if connectivity==2:
        # the list approach
        '''
        tr=[(yi,xi) for yi,xi in zipcor if (yi+1,xi+1) in zipcor]
        br=[(yi,xi) for yi,xi in zipcor if (yi-1,xi+1) in zipcor]
        tl=[(yi,xi) for yi,xi in zipcor if (yi+1,xi-1) in zipcor]
        bl=[(yi,xi) for yi,xi in zipcor if (yi-1,xi-1) in zipcor]

        tr0=[(yi+1,xi+1) for yi,xi in tr]
        br0=[(yi-1,xi+1) for yi,xi in br]
        tl0=[(yi+1,xi-1) for yi,xi in tl]
        bl0=[(yi-1,xi-1) for yi,xi in bl]
        '''

        # the shifting approach
        tr=np.roll(np.roll(mask, -1, axis=0), -1, axis=1)*mask
        br=np.roll(np.roll(mask, 1, axis=0), -1, axis=1)*mask
        tl=np.roll(np.roll(mask, -1, axis=0), 1, axis=1)*mask
        bl=np.roll(np.roll(mask, 1, axis=0), 1, axis=1)*mask

        # add 2 connectivity edges
        #addWeightedEdges(tr,tr0,quslab*costhetas+qvslab*sinthetas)
        #addWeightedEdges(br,br0,quslab*costhetas-qvslab*sinthetas)
        #addWeightedEdges(tl,tl0,-quslab*costhetas+qvslab*sinthetas)
        #addWeightedEdges(bl,bl0,-quslab*costhetas-qvslab*sinthetas)
        addWeightedEdges2(tr,quslab*costhetas+qvslab*sinthetas,'tr')
        addWeightedEdges2(br,quslab*costhetas-qvslab*sinthetas,'br')
        addWeightedEdges2(tl,-quslab*costhetas+qvslab*sinthetas,'tl')
        addWeightedEdges2(bl,-quslab*costhetas-qvslab*sinthetas,'bl')

    return g


def getARAxis(g, quslab, qvslab, mask):
    '''Find AR axis from AR region mask

    Args:
        g (networkx.DiGraph): directed planar graph constructed from AR mask
                              and flows. See maskToGraph().
        quslab (cdms.TransientVariable): 2D map of u-flux.
        qvslab (cdms.TransientVariable): 2D map of v-flux.
        mask (ndarray): 2D binary map showing the location of an AR with 1s.

    Returns:
        path (ndarray): Nx2 array storing the AR axis coordinate indices in
                        (y, x) format.
        axismask (ndarray): 2D binary map with same shape as <mask>, with
                            grid cells corresponding to coordinates in <path>
                            set to 1s.
    '''

    nodes=list(g.nodes())

    #---------------Find boundary nodes---------------
    edge=mask-morphology.binary_erosion(mask)

    gy,gx=np.gradient(np.array(mask))
    inedge=(gx*quslab+gy*qvslab)*edge
    inedgecoor=np.where(inedge>0)
    inedgecoor=zip(inedgecoor[0],inedgecoor[1])
    inedgecoor=list(set(inedgecoor).intersection(nodes))

    outedgecoor=np.where(inedge<0)
    outedgecoor=zip(outedgecoor[0],outedgecoor[1])
    outedgecoor=list(set(outedgecoor).intersection(nodes))

    n1=len(inedgecoor)
    n2=len(outedgecoor)

    # when mask is at edge of the map. Rarely happens.
    if n1==0:
        inedgecoor=nodes
        n1=len(inedgecoor)
    if n2==0:
        outedgecoor=nodes
        n2=len(outedgecoor)

    dists=np.zeros((n1,n2))

    def sumDists(path,attr,g):
        '''Sum edge distances along a path'''
        s=0
        for ii in range(len(path)-1):
            if NX_VERSION=='2':
                sii=g[path[ii]][path[ii+1]][attr]
            else:
                sii=g.edge[path[ii]][path[ii+1]][attr]

            # penalize sharp turns. Doesn't make big difference but notably
            # slower
            '''
            if ii+2<len(path):
                pii1=(lats[path[ii][0]], lons[path[ii][1]])
                pii2=(lats[path[ii+1][0]], lons[path[ii+1][1]])
                pii3=(lats[path[ii+2][0]], lons[path[ii+2][1]])
                theta1=computeTheta(pii1,pii2)
                theta2=computeTheta(pii2,pii3)
                dtheta=theta1.dot(theta2)
                dtheta=abs(dtheta)**1

                #if ii==0:
                    #dtheta_old=1.

                #dtheta=np.mean([dtheta,dtheta_old])
                #sii=sii*dtheta
                #dtheta_old=dtheta
            '''
            s+=sii

        return s

    #---------------Find "longest" path---------------
    for ii in range(n1):
        eii=inedgecoor[ii]
        pathsii=nx.single_source_dijkstra_path(g,eii,weight='weight')
        pathsii=dict([(kk,vv) for kk,vv in pathsii.items() if kk in outedgecoor])
        if len(pathsii)>0:
            distdict=dict([(kk, sumDists(vv,'ivt',g)) for kk,vv in pathsii.items()])
            nodeii=sorted(distdict,key=distdict.get)[-1]
            distii=distdict[nodeii]
            dists[ii,outedgecoor.index(nodeii)]=distii

    if np.max(dists)==0:
        # this may happen when a mask is touching the map edges, and inedgecoor
        # outedgecoor can't be linked by a path. Very rarely happen, but damn
        # annoying. A fallback solution is to use an undirected graph linking
        # the most inward and most outward pixels.
        mostin=np.unravel_index(np.argmax(inedge), mask.shape)
        mostout=np.unravel_index(np.argmin(inedge), mask.shape)
        g_und=g.to_undirected()
        try:
            path=nx.dijkstra_path(g_und,mostin,mostout,weight='weight')
        except:
            # if it still can't find a path, make a full connected network
            g_full=maskToGraph(mask, quslab, qvslab, np.ones(mask.shape),
                    np.ones(mask.shape), -np.inf)
            path=nx.dijkstra_path(g_full,mostin,mostout,weight='weight')
    else:
        maxidx=np.argmax(dists)
        yidx,xidx=np.unravel_index(maxidx,(n1,n2))
        path=nx.dijkstra_path(g,inedgecoor[yidx],outedgecoor[xidx],weight='weight')

    # get a mask for axis
    axismask=np.zeros(mask.shape)
    for (y,x) in path:
        axismask[y,x]=1

    path=np.array(path)

    return path, axismask


def cropMask(mask, edge=4):
    '''Cut out a bounding box around mask==1 areas

    Args:
        mask (ndarray): 2D binary map showing the location of an AR with 1s.
        edge (int): number of pixels as edge at 4 sides.

    Returns:
        mask[y1:y2,x1:x2] (ndarray): a sub region cut from <mask> surrouding
                                      regions with value=1.
        (yy,xx): y-, x- indices of the box of the cut region. Can later by
                 used in applyCropIdx(new_slab, (yy,xx)) to crop out the same
                 region from a new array <new_slab>.
    '''

    yidx,xidx=np.where(mask==1)
    if len(yidx)==0:
        raise Exception("mask empty")

    y1=np.min(yidx)
    y2=np.max(yidx)
    x1=np.min(xidx)
    x2=np.max(xidx)

    y1=max(0,y1-edge)
    y2=min(mask.shape[0],y2+edge)
    x1=max(0,x1-edge)
    x2=min(mask.shape[1],x2+edge)
    xx=np.arange(x1,x2)
    yy=np.arange(y1,y2)

    return mask[y1:y2,x1:x2], (yy,xx)


def applyCropIdx(slab, cropidx):
    '''Cut out a bounding box from given 2d slab given corner indices

    Args:
        slab (ndarray): 2D array to cut a box from.
        cropidx (tuple): (y, x) coordinate indices, output from cropMask().

    Returns:
        cropslab (ndarray): 2D sub array cut from <slab> using <cropidx> as
                            boundary indices.
    '''

    cropslab=np.array(slab)[np.ix_(*cropidx)]
    try:
        croplat=slab.getLatitude()[:][cropidx[0]]
        croplon=slab.getLongitude()[:][cropidx[1]]

        croplat=cdms.createAxis(croplat)
        croplat.designateLatitude()
        croplat.id='y'
        croplat.units='degree'
        croplat.name='latitude'

        croplon=cdms.createAxis(croplon)
        croplon.designateLongitude()
        croplon.id='x'
        croplon.units='degree'
        croplon.name='longitude'

        cropslab=MV.array(cropslab)
        cropslab.setAxis(0,croplat)
        cropslab.setAxis(1,croplon)
    except:
        pass

    return cropslab


def insertCropSlab(shape, cropslab, cropidx, axislist=None):
    '''Insert the cropped sub-array back to a larger empty slab

    Args:
        shape (tuple): (n, m) size of the larger slab.
        cropslab (ndarray): 2D array to insert.
        cropidx (tuple): (y, x) coordinate indices, output from cropMask(),
                         defines where <cropslab> will be inserted into.
    Keyword Args:
        axislist (list or None): if list, a list of cdms.TransientAxis objs.

    Returns:
        result (ndarray): 2D slab with shape (n, m), an empty array with a
                          box at <cropidx> replaced with data from <cropslab>.
                          Optionally, axes information is added if <axistlist>
                          is not None, making it an TransientVariable.
    '''

    result=np.zeros(shape)
    result[np.ix_(*cropidx)]=cropslab

    if axislist is not None:
        result=MV.array(result)
        result.setAxisList(axislist)

    return result


def partPeaks(cropmask, cropidx, orislab, max_ph_ratio):
    '''Separate local maxima by topographical prominence

    Args:
        cropmask (ndarray): 2D binary array, defines regions of local maxima.
        cropidx (tuple): (y, x) coordinate indices, output from cropMask().
        orislab (ndarray): 2D array, giving magnitude/height/intensity values
                           defining the topography.
        max_ph_ratio (float): maximum peak/height ratio. Local peaks with
                              a peak/height ratio larger than this value is
                              treated as an independent peak.

    Returns:
        result (ndarray): 2D binary array, similar as the input <cropmask>
                          but with connected peaks (if any) separated so that
                          each connected region (with 1s) denotes an
                          independent local maximum.
    '''

    cropslab=applyCropIdx(orislab,cropidx)
    if 0 in cropidx[0] or 0 in cropidx[1] or orislab.shape[0]-1 in\
            cropidx[0] or orislab.shape[1]-1 in cropidx[1]:
        include_edge=True
    else:
        include_edge=False

    # compute prominences
    peaks,peakid,peakpro,peakparents=pp2d.getProminence(cropslab*cropmask,
            10.,include_edge=include_edge,centroid_num_to_center=1,verbose=False)

    peakheights=(peakpro>0)*cropslab*cropmask
    ratios=cropmask*peakpro/peakheights

    # take maxima whose prominence/height ratio> max_ph_ratio
    localpeaks=np.where(ratios>max_ph_ratio)
    localpeaks=zip(localpeaks[0],localpeaks[1])

    mask1=np.zeros(cropmask.shape)    # modified mask
    # residual mask, the union of complimentary masks. A complimentary mask
    # is the sea level mask that separates a peak from its parent. Note that
    # peaks' sea levels are not necessarily at the same height.
    resmask=np.zeros(cropmask.shape)

    def breakPeaks(yidx,xidx,localpeaks,col):
        '''Separate the contour of a peak from its parent by iteratively
        rising the sea level
        '''
        labels=morphology.label(cropslab*cropmask>col)
        dropthis=False
        while True:
            #plabels=[labels[yjj,xjj] for yjj,xjj in localpeaks]
            plabels=[labels[yjj,xjj] for yjj,xjj in localpeaks if labels[yjj, xjj]==labels[yidx, xidx]]
            #if len(set(plabels))==len(localpeaks):
            if len(plabels)==1:
                break
            col+=5.
            if col>cropslab[yidx,xidx]:
                dropthis=True
                col-=5.
                break
            labels=morphology.label(cropslab*cropmask>col)

        tmpmask=np.zeros(cropmask.shape)
        if not dropthis:
            tmpmask[yidx,xidx]=1
            tmpmask=morphology.reconstruction(tmpmask,cropslab>col,'dilation')

        return tmpmask,col


    if len(localpeaks)==1:
        mask1=cropmask
    else:
        # sort by prominence/height ratios
        ratios=[ratios[int(yjj), int(xjj)] for yjj,xjj in localpeaks]
        heights=[peakheights[int(yjj), int(xjj)] for yjj, xjj in localpeaks]
        sortidx=np.argsort(ratios)
        ratios.sort()
        localpeaks=[localpeaks[idjj] for idjj in sortidx]
        '''
        localpeakids=[peakid[yjj,xjj] for yjj,xjj in localpeaks]
        localpeakcols=[peaks[idjj]['col_level'] for idjj in localpeakids]
        c_col=np.min(localpeakcols)-10

        labs=morphology.label(cropslab*cropmask>c_col)
        while True:
            if c_col>np.max(cropslab*cropmask):
                break
            plabels=[labs[yjj,xjj] for yjj,xjj in localpeaks]
            if len(set(plabels))==len(localpeaks):
                break
            c_col+=5.
            labs=morphology.label(cropslab*cropmask>c_col)
        mask1=morphology.reconstruction(localpeaks_map, cropslab>c_col, 'dilation')
        '''

        for yjj,xjj in localpeaks:
            yjj=int(yjj)
            xjj=int(xjj)
            idjj=peakid[yjj,xjj]
            coljj=peaks[idjj]['col_level']
            if peaks[idjj]['parent']==0 and peakheights[yjj,xjj]==np.max(heights):
                # the heighest peak
                tmpmask=np.zeros(cropslab.shape)
                tmpmask[yjj,xjj]=1
                tmpmask=morphology.reconstruction(tmpmask,
                        (cropmask-resmask)>0,'dilation')
                mask1=mask1+tmpmask
            else:
                # separate local peaks
                tmpmask,coljj2=breakPeaks(yjj,xjj,localpeaks,coljj)
                # if childrens overlap, may not need this anymore
                if (tmpmask+mask1).max()>1:
                    tmpmask2=np.zeros(cropmask.shape)
                    tmpmask2[yjj,xjj]=1
                    tmpmask=morphology.reconstruction(tmpmask2,tmpmask-tmpmask*resmask,'dilation')
                mask1=mask1+tmpmask
                resmask=np.where((resmask==1) | (cropslab<coljj2),1,0)

    result=insertCropSlab(orislab.shape,mask1,cropidx)

    return result


def getARData(slab, quslab, qvslab, anoslab, quano, qvano, areas,
        mask_list, axis_list, timestr, param_dict):
    '''Fetch AR related data

    Args:
        slab (cdms.TransientVariable): (n * m) 2D array of IVT, in kg/m/s.
        quslab (cdms.TransientVariable): (n * m) 2D array of u-flux, in kg/m/s.
        qvslab (cdms.TransientVariable): (n * m) 2D array of v-flux, in kg/m/s.
        anoslab (cdms.TransientVariable): (n * m) 2D array of IVT anomalies,
                                        in kg/m/s.
        quano (cdms.TransientVariable): (n * m) 2D array of u-flux anomalies,
                                        in kg/m/s.
        qvano (cdms.TransientVariable): (n * m) 2D array of v-flux anomalies,
                                        in kg/m/s.
        areas (cdms.TransientVariable): (n * m) 2D grid cell area slab, in km^2.
        mask_list (list): list of 2D binary masks, each with the same shape as
                      <anoslab> etc., and with 1s denoting the location of a
                      found AR.
        axis_list (list): list of AR axis coordinates. Each coordinate is defined
                     as a Nx2 ndarray storing (y, x) indices of the axis
                     (indices defined in the matrix of corresponding mask
                     in <masks>.)
        timestr (str): string of time snap.
        param_dict (dict): parameter dict defined in Global preamble.

    Returns:
        labels (cdms.TransientVariable): (n * m) 2D int map showing all ARs
                                       at current time. Each AR is labeled by
                                       an int label, starting from 1. Background
                                       is filled with 0s.
        angles (cdms.TransientVariable): (n * m) 2D map showing orientation
                                       differences between AR axes and fluxes,
                                       for all ARs. In degrees.
        crossfluxes (cdms.TransientVariable): (n * m) 2D map showing cross-
                                           sectional fluxes in all ARs.
                                           In kg/m/s.
        anocrossflux (cdms.TransientVariable): similar as <crossfluxes> but for
                                             anomalous fluxes (corresponding
                                             to <anoslab>).
        df (pandas.DataFrame): AR record table. Each row is an AR, see code
                               below for columns.
    '''

    max_isoq=param_dict['max_isoq']
    min_length=param_dict['min_length']
    min_length_hard=param_dict['min_length_hard']
    rdp_thres=param_dict['rdp_thres']
    min_area=param_dict['min_area']

    lonax=slab.getLongitude()   # NOTE: max > 360
    latax=slab.getLatitude()

    # prepare outputs
    labels=MV.zeros(slab.shape)
    angles=MV.zeros(slab.shape)
    crossfluxes=MV.zeros(slab.shape)
    results={}

    #-----------------Loop through ARs-----------------
    for ii in range(len(mask_list)):

        maskii=mask_list[ii]

        # region properties, in pixel units
        rpii=measure.regionprops(maskii, intensity_image=np.array(slab))[0]

        # get centroid
        centroidy,centroidx=rpii.weighted_centroid
        centroidy=latax[int(centroidy)]
        centroidx=lonax[int(centroidx)]

        # get axis coordinate array
        skelii=axis_list[ii]
        latsii=latax[skelii[:,0]]
        lonsii=lonax[skelii[:,1]]
        axisii=np.c_[latsii,lonsii]

        # segment axis using rdp
        axis_rdpii=np.array(rdp.rdpGC(axisii.tolist(),rdp_thres))  # lat,lon

        # area
        areaii=(maskii*areas).sum()  # km^2

        # compute length
        lens=funcs.greatCircle(axis_rdpii[:-1,0], axis_rdpii[:-1,1],
                axis_rdpii[1:,0], axis_rdpii[1:,1])/1e3
        lenii=lens.sum() #km

        # skip if too small and too short
        if areaii<min_area or lenii<min_length_hard:
            continue

        # mean width
        widthii=areaii/lenii # km

        # mask contour
        contii=funcs.getBinContour(maskii,lonax,latax)

        # isoperimetric quotient
        isoquoii=4*np.pi*rpii.area/rpii.perimeter**2

        # length/width ratio
        ratioii=lenii/widthii

        # mean strength
        slabii=MV.masked_where(maskii==0,slab)
        strengthii=cdutil.averager(slabii,axis='xy',
                weights=['generate','generate'])

        # strength std
        strengthstdii=float(stats.std(slabii,axis='xy'))

        # anomaly strength
        anoslabii=MV.masked_where(maskii==0,anoslab)
        anostrengthii=cdutil.averager(anoslabii,axis='xy',
                weights=['generate','generate'])

        # max strength
        max_strengthii=float(MV.max(slabii))

        # compute angles and cross-section flux of total flux
        cropmask,cropidx=cropMask(maskii)
        #cropskelii=skelii-np.array([cropidx[0].min(), cropidx[1].min()])

        cropu=applyCropIdx(quslab,cropidx)
        cropv=applyCropIdx(qvslab,cropidx)

        anglesii,anglesmeanii,crossfluxii,seg_thetasii=crossSectionFlux(
                cropmask, cropu, cropv, axis_rdpii)

        # create plots
        #if isplot:
            #pass
            #plotARCrosssectionFlux(cropmask, cropu, cropv, cropskelii, axis_rdpii,
                    #'%s AR-%d' %(timestr, ii+1), shift_lon, anglesii, anglesmeanii,
                    #crossfluxii, seg_thetasii, outputdir)

        # insert crop back to the big map
        anglesii=insertCropSlab(maskii.shape,anglesii,cropidx,
                slab.getAxisList())
        anglesii=MV.where(maskii==1,anglesii,0)

        crossfluxii=insertCropSlab(maskii.shape,crossfluxii,cropidx,
                slab.getAxisList())
        crossfluxii=MV.where(maskii==1,crossfluxii,0)

        # mean meridional flux
        cropv=applyCropIdx(qvslab,cropidx)
        cropv=MV.masked_where(cropmask==0,cropv)
        qvmeanii=cdutil.averager(cropv,axis='xy',weights=['generate',\
                'generate'])

        # is candidate a strict AR
        is_relaxedii=False
        if isoquoii>max_isoq or ratioii<2:
            is_relaxedii=True
        if lenii<min_length:
            is_relaxedii=True
        if qvmeanii<=0:
            is_relaxedii=True

        labels=labels+maskii*(ii+1)
        angles=angles+anglesii
        crossfluxes=crossfluxes+crossfluxii

        results[ii+1]={
                'id': ii+1,
                'time':timestr,
                'contour_y': contii.vertices[:,1],
                'contour_x': contii.vertices[:,0],
                'centroid_y': centroidy,
                'centroid_x': centroidx,
                'axis_y':axisii[:,0],
                'axis_x':axisii[:,1],
                'axis_rdp_y':axis_rdpii[:,0],
                'axis_rdp_x':axis_rdpii[:,1],
                'area': areaii,
                'length': lenii,
                'width': widthii,
                'iso_quotient':isoquoii,
                'LW_ratio':ratioii,
                'strength':strengthii,
                'strength_ano':anostrengthii,
                'strength_std':strengthstdii,
                'max_strength':max_strengthii,
                'mean_angle': float(anglesmeanii),
                'is_relaxed':is_relaxedii,
                'qv_mean':qvmeanii
                }

    labels.setAxisList(slab.getAxisList())
    angles.setAxisList(slab.getAxisList())
    crossfluxes.setAxisList(slab.getAxisList())

    labels.id='labels'
    labels.long_name='AR labels'
    labels.standard_name=labels.long_name
    labels.title=labels.long_name
    labels.units=''

    angles.id='angles'
    angles.long_name='AR moisture flux orientation difference'
    angles.standard_name=angles.long_name
    angles.title=angles.long_name
    angles.units='degree'

    crossfluxes.id='ivt_cross'
    crossfluxes.long_name='AR total cross sectional moisture flux'
    crossfluxes.standard_name=crossfluxes.long_name
    crossfluxes.title=crossfluxes.long_name
    crossfluxes.units=getattr(slab, 'units', '')

    keys=['id', 'time', 'contour_y', 'contour_x', 'centroid_y', 'centroid_x',
            'axis_y', 'axis_x', 'axis_rdp_y', 'axis_rdp_x',
            'area', 'length', 'width', 'iso_quotient', 'LW_ratio',
            'strength', 'strength_ano', 'strength_std', 'max_strength',
            'mean_angle', 'is_relaxed', 'qv_mean']

    df=pd.DataFrame(results).T
    if len(df)>0:
        df=df[keys]


    return labels,angles,crossfluxes,df


def uvDecomp(u0, v0, i1, i2):
    '''Decompose background-transient components of u-, v- fluxes

    Args:
        u0 (cdms.TransientVariable): nd array of total u-flux.
        v0 (cdms.TransientVariable): nd array of total v-flux.
        i1 (cdms.TransientVariable): nd array of the reconstruction component
                                   of IVT.
        i2 (cdms.TransientVariable): nd array of the anomalous component
                                   of IVT (i2 = IVT - i1).
    Returns:
        u1 (cdms.TransientVariable): nd array of the u-flux component
                                   corresponding to <i1>, i.e. the background
                                   component.
        v1 (cdms.TransientVariable): nd array of the v-flux component
                                   corresponding to <i1>, i.e. the background
                                   component.
        u2 (cdms.TransientVariable): nd array of the u-flux component
                                   corresponding to <i2>, i.e. the transient
                                   component.
        v2 (cdms.TransientVariable): nd array of the v-flux component
                                   corresponding to <i2>, i.e. the transient
                                   component.
    '''

    i0=i1+i2
    v1=v0*i1/i0
    v2=v0*i2/i0
    u1=u0*i1/i0
    u2=u0*i2/i0

    return u1,u2,v1,v2


def save2DF(result_dict):
    '''Save AR records to a pandas DataFrame

    Args:
        result_dict (dict): key: time str in 'yyyy-mm-dd hh:00'
                            value: pandas dataframe. See getARData().

    Returns:
        result_df (pandas.DataFrame): AR record table containing records
                                      from multiple time steps sorted by
                                      time.
    '''

    for ii,kk in enumerate(result_dict.keys()):
        vv=result_dict[kk]
        if ii==0:
            result_df=vv
        else:
            result_df=pd.concat([result_df,vv],axis=0,ignore_index=True)

    result_df['time']=pd.to_datetime(result_df.time)
    result_df=result_df.sort_values(by='time')

    return result_df


def plotAR(ardf, ax, bmap):
    '''Helper function to plot the regions and axes of ARs

    Args:
        ardf (pandas.DataFrame): table containing AR records.
        ax (matplotlib axis): axis to plot onto.
        bmap (Basemap obj): defining the geo map.
    '''

    for ii in range(len(ardf)):

        vv=ardf.iloc[ii]
        isrelaxkk=vv['is_relaxed']

        # plot contour
        px=vv['contour_x']
        py=vv['contour_y']

        px,py=bmap(px,py)
        linewidth=1.5 if isrelaxkk else 1.5
        linestyle=':' if isrelaxkk else '-'
        ax.plot(px,py,color='k',linestyle=linestyle,linewidth=linewidth)

        # plot axis
        px=vv['axis_x']
        py=vv['axis_y']

        px,py=bmap(px,py)
        ax.plot(px,py,'g:',linewidth=2.0)

        # plot cross flux text
        '''
        lenkk=vv['length']
        areakk=vv['area']
        widthkk=vv['width']
        cx=float(vv['centroid_x'])%360
        cy=float(vv['centroid_y'])
        cx,cy=bmap(cx,cy)

        strkk=r'ID=%d, $R=%.0f$' %(ii+1,np.sqrt(areakk/3.14)) + '\n'+\
                r'$L = %d km$' %lenkk +'\n'+\
                r'$W = %d km$' %widthkk

        ax.annotate(strkk,xy=(cx,cy),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=8,
                bbox=dict(facecolor='white',alpha=0.5))
        '''

    return


def getNormalVectors(point_list, idx):
    '''Get the normal vector and the tagent vector to the plane dividing
    2 sections along the AR axis.

    Args:
        point_list (list): list of (lat, lon) coordinates.
        idx (int): index of the point in <point_list>, denoting the point
                   in question.

    Returns:
        normi (tuple): the (x, y, z) Cartesian coordinate of the unit
                       normal vector, at the point denoted by <idx>,
                       on the Earth surface. This is the normal vector to
                       the plane spanned by the vector Theta and P.
                       Where P is the vector pointing to the point in
                       question (point_list[idx]), and Theta is the
                       tangent vector evenly dividing the angle formed by
                       <P,P1>, and <P,P2>. Where P1, P2 are 2 points on
                       both side of P.
        thetai (tuple): the (x, y, z) Cartesian coordinate of the tangent
                        vector Theta above.
    '''
    pi=point_list[idx]
    pic=spherical2Cart(*pi)
    theta1=computeTheta(pi,point_list[idx-1])
    theta2=computeTheta(pi,point_list[idx+1])
    thetai=(theta1+theta2)/2.

    normi=np.cross(pic,thetai)
    normi=normi/np.linalg.norm(normi)
    thetai=thetai/np.linalg.norm(thetai)

    return normi,thetai


def crossSectionFlux(mask, quslab, qvslab, axis_rdp):
    '''Compute setion-wise orientation differences and cross-section fluxes
    in an AR

    Args:
        mask (ndarray): CROPPED (see cropMask and applyCropIdx) 2D binary map
                        showing the location of an AR with 1s.
        quslab (cdms.TransientVariable): CROPPED (n * m) 2D array of u-flux,
                                       in kg/m/s.
        qvslab (cdms.TransientVariable): CROPPED (n * m) 2D array of v-flux,
                                       in kg/m/s.
        axis_rdp (ndarray): Nx2 array storing the (lat, lon) coordinates of
                           rdp-simplified AR axis.

    Returns:

        angles (TransientVariable): 2D map with the same shape as <mask>,
                                    showing section-wise orientation
                                    differences between horizontal flux (as
                                    in <quslab>, <qvslab>) and the AR axis of
                                    that section. In degrees. Regions outside
                                    of AR (0s in <mask>) are masked.

        anglesmean (float): area-weighted averaged of <angles> inside <mask>.

        crossflux (TransientVariable): 2D map with the same shape as <mask>,
                                       the section-wise cross-section fluxes
                                       in the AR, defined as the projection
                                       of fluxes onto the AR axis, i.e. flux
                                       multiplied by the cos of <angles>.

        seg_thetas (list): list of (x, y, z) Cartesian coordinates of the
                           tangent vectors along section boundaries.
    '''
    # get coordinates
    axislist=quslab.getAxisList()
    lats=quslab.getLatitude()[:]
    lons=quslab.getLongitude()[:]
    lonss,latss=np.meshgrid(lons,lats)

    # convert to cartesian coordinates
    carts=spherical2Cart(latss,lonss)
    vs=wind2Cart(quslab,qvslab,latss,lonss)
    vsnorm=np.linalg.norm(vs,axis=0)
    vsnorm=vs/vsnorm[None,:,:]

    # loop through segments to get orientation differences
    nsegs=len(axis_rdp)-1
    seg_thetas=[]
    angles=np.zeros(mask.shape)

    for ii in range(nsegs):

        pic=spherical2Cart(*axis_rdp[ii])
        pi1c=spherical2Cart(*axis_rdp[ii+1])

        if ii==0:
            setL=1.
            thetai=0 # dummy place holder
        else:
            # get evenly dividing angle theta and normal vector to theta
            normi,thetai=getNormalVectors(axis_rdp, ii)
            # dot products between normal vector and grid coordinates
            dotsi=(normi[:,None,None]*carts).sum(axis=0)
            setL=np.where(dotsi*(normi.dot(pi1c))>=0,1,0)

        if ii==nsegs-1:
            setR=1.
            thetai=0 # dummy place holder
        else:
            normi1,thetai=getNormalVectors(axis_rdp, ii+1)
            dotsi1=(normi1[:,None,None]*carts).sum(axis=0)
            setR=np.where(dotsi1*(normi1.dot(pic))>0,1,0)

        segii=setL*setR*mask
        seg_thetas.append(thetai)

        # sel the correct region if shape too curvy
        segregii=measure.label(segii)
        if segregii.max()>1:
            piidx=[funcs.findIndex(axis_rdp[ii][0],lats),\
                    funcs.findIndex(axis_rdp[ii][1],lons)]
            for jj in range(segregii.max()):
                segjj=np.where(segregii==jj+1,1,0)
                if segjj[piidx[0],piidx[1]]==1:
                    segii=segjj
                    break

        # mean orientation of AR axis segment
        meanori=np.cross(pic,pi1c)
        meanori=meanori/np.linalg.norm(meanori)

        # orientation of flux vectors
        '''
        fluxori=np.cross(carts,vs,axisa=0,axisb=0) # output: ny,nx,3
        norms=np.linalg.norm(fluxori,axis=2)
        fluxori=fluxori/(1e-6+norms[:,:,None])

        # get angles as arccos of the dot product of meanori and fluxori
        anglesii=(meanori[None,None,:]*fluxori).sum(axis=-1)
        anglesii=anglesii*segii
        '''

        # get sin(angle of flux vector and axis segment plane)
        # sign of sin() is: >0: flux vector aligns with meanori, and it
        # is pointing towards the cold side (according to thermal wind).
        # sin() <0: flux vector aligns against meanori, and it is pointing
        # towards the warm side.
        anglesii=(meanori[:,None,None]*vsnorm).sum(axis=0)
        #anglesii=np.sqrt(1-cos_alphaii**2)*np.where(cos_alphaii<0,-1,1)
        anglesii=anglesii*segii

        angles=angles+anglesii


    # compute cross section flux
    angles=np.array(angles)
    cos_angles=np.sqrt(1-angles**2)
    crossflux_c=cos_angles[None,:,:]*vs
    # convert to local tangent winds
    crossflux_u,crossflux_v=cart2Wind(crossflux_c,latss,lonss)
    crossflux=np.sqrt(crossflux_u**2+crossflux_v**2)
    crossflux=MV.masked_where(mask==0,crossflux)
    crossflux.setAxisList(axislist)

    # convert cos to angle in degrees
    angles=np.arcsin(angles)/np.pi*180
    angles=MV.masked_where(mask==0,angles)
    angles.setAxisList(axislist)

    anglesmean=cdutil.averager(angles,axis='xy',weights=['generate','generate'])

    return angles,anglesmean,crossflux,seg_thetas


def findARAxis(quslab, qvslab, armask_list, costhetas, sinthetas, param_dict,
        verbose=True):
    '''Find AR axis

    Args:
        quslab (cdms.TransientVariable): (n * m) 2D u-flux slab, in kg/m/s.
        qvslab (cdms.TransientVariable): (n * m) 2D v-flux slab, in kg/m/s.
        armask_list (list): list of 2D binary masks, each with the same shape
            as <quslab> etc., and with 1s denoting the location of a found AR.
        costhetas (cdms.TransientVariable): (n * m) 2D slab of grid cell shape:
                                          cos=dx/sqrt(dx^2+dy^2).
        sinthetas (cdms.TransientVariable): (n * m) 2D slab of grid cell shape:
                                          sin=dy/sqrt(dx^2+dy^2).
        param_dict (dict): parameter dict defined in Global preamble.

    Returns:
        axes (list): list of AR axis coordinates. Each coordinate is defined
                     as a Nx2 ndarray storing (y, x) indices of the axis
                     (indices defined in the matrix of corresponding mask
                     in <armask_list>.)
        axismask (ndarray): 2D binary mask showing all axes in <axes> merged
                            into one map.

    New in v2.0.
    '''

    edge_eps=param_dict['edge_eps']

    #-----------------Prepare outputs-----------------
    axismask=np.zeros(quslab.shape)
    axes=[]

    #--------------------Find axes--------------------
    for maskii in armask_list:

        #----------Convert mask to directed graph----------
        gii=maskToGraph(maskii,quslab,qvslab,costhetas,sinthetas,edge_eps)

        #--------------Get AR axis from graph--------------
        axisarrii,axismaskii=getARAxis(gii,quslab,qvslab,maskii)
        axes.append(axisarrii)
        axismask=axismask+axismaskii

    return axes, axismask


def prepareMeta(lats, lons, times, ntime, nlat, nlon,
        ref_time='days since 1900-01-01', verbose=True):
    '''Prepare metadata for AR detection function calls

    Args:
        lats (ndarray): 1D, latitude coordinates, the length needs to equal
            <nlat>.
        lons (ndarray): 1D, longitude coordinates, the length needs to equal
            <nlon>.
        times (list or array): time stamps of the input data as a list of strings,
            e.g. ['2007-01-01 06:00:00', '2007-01-01 12:00', ...].
            Needs to have the a length of <ntime>.

    Keyword Args:
        ref_time (str): reference time point to create time axis.

    Returns:
        timeax (cdms2.axis.TransientAxis): a time axis obj created from strings
            in <times>.
        areamap (ndarray): grid cell areas in km^2, with shape (<nlat> x <nlon>).
        costhetas (ndarray): ratios of dx/sqrt(dx^2 + dy^2) for all grid cells.
            with shape (<nlat> x <nlon>).
        sinthetas (ndarray): ratios of dy/sqrt(dx^2 + dy^2) for all grid cells.
            with shape (<nlat> x <nlon>).

    New in v2.0.
    '''

    #-------------------Check inputs-------------------
    lats=np.asarray(lats).squeeze()
    lons=np.asarray(lons).squeeze()
    if np.ndim(lats) != 1:
        raise Exception("<lats> needs to be an 1-D array.")
    if np.ndim(lons) != 1:
        raise Exception("<lons> needs to be an 1-D array.")
    if len(lats) != nlat:
        raise Exception("Length of <lats> doesn't equal <nlat>.")
    if len(lons) != nlon:
        raise Exception("Length of <lons> doesn't equal <nlon>.")

    #------Make sure lats and lons are increasing------
    # NOTE: important to make sure lat is increasing
    lats=np.sort(lats)
    lons=np.sort(lons)

    #------------------Get time axis------------------
    timeax=funcs.getTimeAxis(times, ntime, ref_time).asComponentTime()

    #-------------Compute grid geometries-------------
    dxs=funcs.dLongitude2(lats, lons, R=6371)
    dys=funcs.dLatitude2(lats, lons, R=6371)
    areamap=dxs*dys # km^2
    costhetas=dxs/np.sqrt(dxs**2+dys**2)
    sinthetas=dys/np.sqrt(dxs**2+dys**2)

    if verbose:
        print('\n# <prepareMeta>: Metadata created.')

    return timeax, areamap, costhetas, sinthetas, lats, lons


def _findARs(anoslab, areas, param_dict):
    '''Find ARs from THR results at a time snap

    Args:
        anoslab (cdms.TransientVariable): (n * m) 2D anomalous IVT slab, in kg/m/s.
        areas (cdms.TransientVariable): (n * m) 2D grid cell area slab, in km^2.
        param_dict (dict): parameter dict controlling the detection process.

    Returns:
        masks (list): list of 2D binary masks, each with the same shape as
                      <anoslab>, and with 1s denoting the location of a
                      found AR.
        armask (ndarray): 2D binary mask showing all ARs in <masks> merged into
                         one map.

    If no ARs are found, <masks> will be []. <armask> will be zeros.
    '''

    # fetch parameters
    thres_low=param_dict['thres_low']
    min_area=param_dict['min_area']
    max_area=param_dict['max_area']
    min_lat=param_dict['min_lat']
    max_lat=param_dict['max_lat']
    single_dome=param_dict['single_dome']
    max_ph_ratio=param_dict['max_ph_ratio']
    max_isoq_hard=param_dict['max_isoq_hard']
    fill_radius=param_dict['fill_radius']
    latax=anoslab.getLatitude()

    def paddedClosing(mask, ele, pad):
        # pad
        padmask=np.pad(mask, (pad,pad), mode='constant', constant_values=0)
        # closing
        padmask=morphology.closing(padmask, selem=ele)
        # trim
        padmask=padmask[slice(pad,-pad), slice(pad,-pad)]
        return padmask

    mask0=np.where(anoslab>thres_low,1,0)
    if single_dome:
        # NOTE: if using single_dome, should not filter max_area here.
        mask0=areaFilt(mask0,areas,min_area,np.inf)
    else:
        mask0=areaFilt(mask0,areas,min_area,max_area)

    # prepare outputs
    masks=[]
    armask=np.zeros(mask0.shape)

    if mask0.max()==0:
        return masks, armask

    #---------------Separate grouped peaks---------------
    if single_dome:
        labels=measure.label(mask0,connectivity=2)
        mask1=np.zeros(mask0.shape)

        for ii in range(labels.max()):
            maskii=np.where(labels==ii+1,1,0)

            #-------------Skip if latitude too low or too high---------
            latsii=np.where(maskii==1)[0]
            latmaxii=latax[np.max(latsii)]
            if latmaxii<min_lat:
                continue
            latminii=latax[np.min(latsii)]
            if latminii>max_lat:
                continue

            cropmask,cropidx=cropMask(maskii)
            maskii2=partPeaks(cropmask,cropidx,anoslab,max_ph_ratio)
            # should I revert to maskii if the peak separation results in
            # a large area loss?
            mask1=mask1+maskii2
    else:
        mask1=mask0

    mask1=areaFilt(mask1,areas,min_area,max_area)

    if mask1.max()==0:
        return masks, armask

    #-------Latitude and iso_quotient filtering-------
    labels=measure.label(mask1,connectivity=2)
    mask1=np.zeros(mask0.shape)

    for ii in range(labels.max()):
        maskii=np.where(labels==ii+1,1,0)

        #-------------Skip if latitude too low or too high---------
        rpii=measure.regionprops(maskii, intensity_image=np.array(anoslab))[0]
        centroidy,centroidx=rpii.weighted_centroid
        centroidy=latax[int(centroidy)]
        min_lat_idx=np.argmin(np.abs(latax[:]-min_lat))

        if (centroidy<=min_lat and\
                maskii[:min_lat_idx].sum()/float(maskii.sum())>=0.5)\
                or centroidy>=max_lat:
            continue

        # filter by isoperimetric quotient
        isoquoii=4*np.pi*rpii.area/rpii.perimeter**2

        if isoquoii>=max_isoq_hard:
            continue

        mask1=mask1+maskii

    if mask1.max()==0:
        return masks, armask

    #--------Fill some small holes and make the contour smoother--------
    labels=measure.label(mask1,connectivity=2)
    filldisk=morphology.disk(fill_radius)

    for ii in range(labels.max()):
        maskii=np.where(labels==ii+1,1,0)
        maskii=paddedClosing(maskii, filldisk, fill_radius)
        masks.append(maskii)
        armask=armask+maskii

    return masks, armask


def findARs(ivt, ivtrec, ivtano, qu, qv, lats, lons, param_dict,
        times=None, ref_time='days since 1900-01-01',
        verbose=True):
    '''Find ARs from THR results, get all results in one go.

    Args:
        ivt (TransientVariable): 3D or 4D input IVT data, with dimensions
            (time, lat, lon) or (time, level, lat, lon).
        ivtrec (TransientVariable): 3D or 4D array, the reconstruction
            component from the THR process.
        ivtano (TransientVariable): 3D or 4D array, the difference between
            input <ivt> and <ivtrec>.
        qu (TransientVariable): 3D or 4D array, zonal component of
            integrated moisture flux.
        qv (TransientVariable): 3D or 4D array, meridional component of
            integrated moisture flux.
        lats (ndarray): 1D, latitude coordinates, the length needs to be the
            same as the lat dimension of <ivt>.
        lons (ndarray): 1D, longitude coordinates, the length needs to be the
            same as the lon dimension of <ivt>.
        param_dict (dict): parameter dict controlling the detection process.
    Keyword Args:
        times (list or array): time stamps of the input data as a list of strings,
            e.g. ['2007-01-01 06:00:00', '2007-01-01 12:00', ...].
            Needs to have the same length as the time dimension of <ivt>.
            If None, default to create a dummy 6-hourly time axis, using
            <ref_time> as start, with a length as the time dimension of <ivt>.
        ref_time (str): reference time point to create dummy time axis, if
            no time stamps are given in <times>.
    Returns:
        time_idx (list): indices of the time dimension when any AR is found.
        labels_all (TransientVariable): 3D array, with dimension
            (time, lat, lon). At each time slice, a unique int label is assign
            to each detected AR at that time, and the AR region is filled
            out with the label value in the (lat, lon) map.
        angles_all (TransientVariable): 3D array showing orientation
            differences between AR axes and fluxes, for all ARs. In degrees.
        crossfluxes_all (TransientVariable): 3D array showing cross-
            sectional fluxes in all ARs. In kg/m/s.
        result_df (DataFrame): AR record table. Each row is an AR, see code
            in getARData() for columns.
    See also:
        findARsGen(): generator version, yields results at time points
            separately.
    '''

    time_idx=[]
    result_dict={}
    labels_all=[]
    angles_all=[]
    crossfluxes_all=[]

    finder_gen = findARsGen(ivt, ivtrec, ivtano, qu, qv, lats, lons,
            param_dict, times=times, ref_time=ref_time, verbose=verbose)
    next(finder_gen)  # prime the generator to prepare metadata

    for (tidx, timett, label, angle, cross, ardf) in finder_gen:

        time_idx.append(tidx)
        labels_all.append(label)
        angles_all.append(angle)
        crossfluxes_all.append(cross)
        result_dict[timett]=ardf

    labels_all=MV.concatenate(labels_all, axis=0)
    angles_all=MV.concatenate(angles_all, axis=0)
    crossfluxes_all=MV.concatenate(crossfluxes_all, axis=0)

    labels_all.id='labels'
    labels_all.long_name='AR labels'
    labels_all.standard_name=labels_all.long_name
    labels_all.title=labels_all.long_name
    labels_all.units=''

    angles_all.id='angles'
    angles_all.long_name='AR moisture flux orientation difference'
    angles_all.standard_name=angles_all.long_name
    angles_all.title=angles_all.long_name
    angles_all.units='degree'

    crossfluxes_all.id='ivt_cross'
    crossfluxes_all.long_name='AR total cross sectional moisture flux'
    crossfluxes_all.standard_name=crossfluxes_all.long_name
    crossfluxes_all.title=crossfluxes_all.long_name
    crossfluxes_all.units=getattr(ivt, 'units', '')

    result_df=save2DF(result_dict)

    return time_idx, labels_all, angles_all, crossfluxes_all, result_df


def findARsGen(ivt, ivtrec, ivtano, qu, qv, lats, lons, param_dict,
        times=None, ref_time='days since 1900-01-01',
        verbose=True):
    '''Find ARs from THR results, generator version

    Args:
        ivt (TransientVariable): 3D or 4D input IVT data, with dimensions
            (time, lat, lon) or (time, level, lat, lon).
        ivtrec (TransientVariable): 3D or 4D array, the reconstruction
            component from the THR process.
        ivtano (TransientVariable): 3D or 4D array, the difference between
            input <ivt> and <ivtrec>.
        qu (TransientVariable): 3D or 4D array, zonal component of
            integrated moisture flux.
        qv (TransientVariable): 3D or 4D array, meridional component of
            integrated moisture flux.
        lats (ndarray): 1D, latitude coordinates, the length needs to be the
            same as the lat dimension of <ivt>.
        lons (ndarray): 1D, longitude coordinates, the length needs to be the
            same as the lon dimension of <ivt>.
        param_dict (dict): parameter dict controlling the detection process.

    Keyword Args:
        times (list or array): time stamps of the input data as a list of strings,
            e.g. ['2007-01-01 06:00:00', '2007-01-01 12:00', ...].
            Needs to have the same length as the time dimension of <ivt>.
            If None, default to create a dummy 6-hourly time axis, using
            <ref_time> as start, with a length as the time dimension of <ivt>.
        ref_time (str): reference time point to create dummy time axis, if
            no time stamps are given in <times>.

    Returns:
        ii (int): index of the time dimension when any AR is found.
        timett_str (str): time when any AR is found, in string format.
        labels (TransientVariable): 2D array, with dimension
            (lat, lon). A unique int label is assign
            to each detected AR at the time, and the AR region is filled
            out with the label value in the (lat, lon) map.
        angles (TransientVariable): 2D array showing orientation
            differences between AR axes and fluxes, for all ARs. In degrees.
        crossfluxes (TransientVariable): 2D array showing cross-
            sectional fluxes in all ARs. In kg/m/s.
        ardf (DataFrame): AR record table. Each row is an AR, see code
            in getARData() for columns.
    See also:
        findARs(): collect and return all results in one go.
    New in v2.0.
    '''

    def squeezeTo3D(vv):
        if np.ndim(vv) not in [3, 4]:
            raise Exception("Input <ivt>, <ivtrec>, <ivtano>, <qu> and <qv> should be 3D or 4D.")
        if np.ndim(vv)==4 and vv.shape[0]!=1:
            vv=vv(squeeze=1)
        elif np.ndim(vv)==4 and vv.shape[0]==1:
            vv=vv[:,0,:,:]
        elif np.ndim(vv)==3 and vv.shape[0]==1:
            pass
        # NOTE: important to make sure lat is increasing
        vv=funcs.increasingLatitude(vv)
        return vv

    #-----------------Get coordinate metadata-----------------
    # squeeze to 3D
    ivt=squeezeTo3D(ivt)
    ivtrec=squeezeTo3D(ivtrec)
    ivtano=squeezeTo3D(ivtano)
    qu=squeezeTo3D(qu)
    qv=squeezeTo3D(qv)

    timeax, areamap, costhetas, sinthetas, lats, lons = prepareMeta(
            lats, lons, times, ivt.shape[0], ivt.shape[1], ivt.shape[2],
            ref_time=ref_time, verbose=verbose)
    yield # this bare yield prepares the generator by advancing it to the 1st yield

    #######################################################################
    #                          Start processing                           #
    #######################################################################

    #----------------Loop through time----------------
    for ii, timett in enumerate(timeax):

        timett_str='%d-%02d-%02d %02d:00' %(timett.year,timett.month,\
                timett.day,timett.hour)

        if verbose:
            print('\n# <findARsGen>: Processing time: %s' %timett_str)

        slab=ivt[ii]
        slabano=ivtano[ii]
        slabrec=ivtrec[ii]
        quslab=qu[ii]
        qvslab=qv[ii]

        # find ARs
        mask_list,armask=_findARs(slabano, areamap, param_dict)

        # skip if none
        if armask.sum()==0:
            continue

        # find AR axis
        axis_list, axismask=findARAxis(quslab, qvslab, mask_list, costhetas,
                sinthetas, param_dict)

        # decompose background-transient
        qurec,quano,qvrec,qvano=uvDecomp(quslab,qvslab,slabrec,slabano)

        # fetch AR related data
        labels,angles,crossfluxes,ardf=getARData(
                slab,quslab,qvslab,
                slabano,quano,qvano,
                areamap,
                mask_list,axis_list,timett_str,param_dict)

        if verbose:
            print('# <findARsGen>: NO. of ARs found =  %d' %len(ardf))

        # prepare nc output
        timeaxii=cdms.createAxis([timett.torel('days since 1900-1-1').value])
        timeaxii.designateTime()
        timeaxii.id='time'
        timeaxii.units='days since 1900-1-1'

        labels=funcs.addExtraAxis(labels,timeaxii)
        angles=funcs.addExtraAxis(angles,timeaxii)
        crossfluxes=funcs.addExtraAxis(crossfluxes,timeaxii)

        yield ii, timett_str, labels, angles, crossfluxes, ardf


