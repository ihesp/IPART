Detect AR appearances from THR output
=====================================

.. _detect_ars:

Definition of AR occurrence
###########################


An AR occurrence at a given time point is defined using these following rules:

1. A connected region in the IVT anomaly field (:math:`I - \delta(I)`, 
   computed in the section ":ref:`compute_thr`") where its values is greater than 0.
2. The centroid (weighted by IVT values of the grid cells) of the region is north of :math:`20 ^{\circ} N`,
   and south of :math:`80 ^{\circ}`, i.e. we are only interested in mid-latitude systems.
3. The region's area has to be within :math:`50 -- 1800 \times 10^4 km^2`.
4. The region's `isoperimeteric quotient <https://en.wikipedia.org/wiki/Isoperimetric_inequality>`_ :math:`q = \frac{4 \pi A}{L^2} \ge 0.6`. This is to filter out circular features like tropical cyclones."
5. After the computation of this AR candidates axis (see :ref:`compute_axis`) and the effective width (defined as area/length ratio), the length has to be :math:`\ge\, 1500 km`, and length/width ratio has to be :math:`\ge \,2`.


.. _detect_params:

Input data
##########

These are the input data required for AR occurrence detection:

* u- and v- components of the integrated vapor fluxes (:math:`F_u` and :math:`F_v`).
* IVT (as :math:`\sqrt{F_u^2 + F_v^2}`), see :ref:`get_ivt`.
* The output from the THR process: the reconstruction component (:math:`\delta(I)`) and the anomaly
  component (:math:`I - \delta(I)`). See :ref:`compute_thr`.


Additional inputs:

* latitude, longitude and time axis. See :ref:`metadata`.
* detection parameters, see below.

::

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

            # max prominence/height ratio of a local peak. Only used when SINGLE_DOME=True
            'max_ph_ratio': 0.4,

            # minimal proportion of flux component in a direction to total flux to
            # allow edge building in that direction
            'edge_eps': 0.4
            }


Usage in Python scripts
#######################

The following snippet shows the detection function calls:
::

        import MV2 as MV
        from utils import funcs
        from river_tracker1 import findARs
        from river_tracker1_funcs import uvDecomp

        timeax=ivt.getTime().asComponentTime()
        latax=qu.getLatitude()
        lonax=qu.getLongitude()

        dxs=funcs.dLongitude(qu,R=6371)
        dys=funcs.dLatitude(qu,R=6371)
        areamap=dxs*dys # km^2
        costhetas=dxs/MV.sqrt(dxs**2+dys**2)
        sinthetas=dys/MV.sqrt(dxs**2+dys**2)

        #----------------Loop through time----------------
        for ii, timett in enumerate(timeax):

            timett_str='%d-%02d-%02d %02d:00' %(timett.year,timett.month,\
                timett.day,timett.hour)

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

where

* ``ivt`` is the IVT data, in ``(time, level, latitude, longitude)`` or ``(time, latitude, longitude)``.
* ``ivtrec`` is :math:`\delta(I)`, and ``ivtano`` is :math:`I-\delta(I)`, see :ref:`compute_thr`.
* ``qu``: is :math:`F_u`, and ``qv`` is :math:`F_v`.
* ``PARAM_DICT`` is the parameter dictionary as defined above.

After this process, one can optionally call the
``river_tracker1_funcs.getARData()`` function to obtain more AR-related
attributes, including length, width, area, mean IVT values etc..
::

    from river_tracker1_funcs import getARData

    labels, angles, crossfluxes, ardf = getARData(
        slab,quslab,qvslab,
        slabano,quano,qvano,
        areamap,
        mask_list, axis_list, timett_str, PARAM_DICT, 80,
        False, OUTPUTDIR)



Example output
##############

The resultant detected ARs can be visualized using the following snippet:
::

    import matplotlib.pyplot as plt
    from utils import plot
    from river_tracker1_funcs import plotAR

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
    plotAR(ardf,ax,bmap)
    figure.show()


One example output figure is shown below:

.. figure:: ar_1984-01-04_06:00.png
    :width: 700px
    :align: center
    :figclass: align-center

    (a) The IVT field in kg/m/s at 1984-01-04 06:00 UTC over the North
    Hemisphere. (b) the IVT reconstruction field at the same time point. (c)
    the IVT anomaly field from the THR process at the same time point. In all
    three subplots, the detected ARs are outlined in black contour. The AR axes
    are drawn in green dashed lines.


Notebook example
################

An example of this process is given in this `notebook <https://github.com/ihesp/AR_tracker/notebooks/3 detect ARs.ipynb>`_.



