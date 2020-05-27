---
title: 'AR tracker: A Python package for detecting and tracking atmospheric rivers using image-processing'
tags:
  - Python
  - meteorology
  - hydrological cycle
  - image processing
authors:
  - name: Guangzhi XU
    orcid: 0000-0001-5877-5608
    affiliation: 1
  - name: Xiaohui MA
    affiliation: "1, 3"
  - name: Ping CHANG
    affiliation: "2, 3"
affiliations:
 - name: Key Laboratory of Physical Oceanography, Institute for Advanced Ocean Studies, Ocean University of China and Qingdao National Laboratory for Marine Science and Technology, Qingdao, China
   index: 1
 - name: Department of Oceanography and Department of Atmospheric Sciences, Texas A&M University, College Station, Texas, USA
   index: 2
 - name: The International Laboratory for High-Resolution Earth System Prediction, Texas A&M University, College Station, Texas, USA
   index: 3

date: 02 April 2020
bibliography: paper.bib

---

# Summary

An atmospheric river (AR) in the field of meteorology/climatology
refers to enhanced water vapor content in the lower troposphere.  The
term was coined as an analogy to the terrestrial rivers in a sense
that when viewed from satellite imagery or large scale atmospheric
observation, they appear as narrow and elongated vapor filaments,
representing transient intensified horizontal moisture fluxes
(e.g. @Gimeno2014, @Dettinger2011). A typical atmospheric river
can carry 7-15 times the water in the Mississippi River
[@Ralph2011], and at any time in winter, there are four to five such
systems in the Northern Hemisphere alone [@Zhu1998], accounting for
$80-90 \,\%$ of the total north-south integrated vapor transport
[@Guan2015; @Zhu1998].  Its dual hydrological role, both as a fresh
water source for some water-stressed areas [@Dettinger2011;
@Dettinger2013; @Rutz2012] and as a potential trigger for floods
[@Lavers2012; @Lavers2013; @Neiman2008; @Moore2012], has granted it
increasing attentions among the research community.  And their
long-term change in a warming climate also stands as pressing research
question. However, an important prerequisite to answer such questions
is a robust and consistent detection method. As meteorologists and
climatologist often deal with observational or simulation data in
large sizes, an algorithmic method can ensure better efficiency,
consistency and objectivity compared with human
identification.

In many existing applications, a magnitude thresholding
approach is used. For instance, @Ralph2004, @Neiman2008,
@Hagos2015 and @Dettinger2011 identified ARs by first locating
regions where the Integrated Water Vapor (IWV) is greater than $20\,
mm$.  A $250 \, kg/m/s$ threshold on the Integrated Vapor Transport
(IVT) was used by @Rutz2014 and @Rutz2015.
However, an implicit assumption
with this magnitude thresholding approach is that the atmospheric
moisture level stays unchanged throughout the analysis period.  Such
an assumption may not be fully justifiable under a warming climate as
the atmospheric moisture level is expected to increase.

In this package we propose a suite of new detection/tracking
to help solve the above difficulties.
Through a
systematic analysis using 7 years of Reanalysis data [@Dee2011], we
have found that the proposed detection algorithm has reduced sensitivity to
parameters and data resolution.  Long-lived ARs spanning multiple
days, having cross-continent or cross-basin tracks can be more reliably
traced through their tropical/sub-tropical origins to high-latitude
landfalls.

The *AR_tracker* package includes a collection of Python functions/classes
designed for an analysis workflow covering the detection of ARs,
the simplification of the AR's geographical location,
to the subsequent tracking through time.
The algorithms are implemented using the Python
programming language as a wrapper to some well-established numeric packages
including *numpy* ,*scikit-image* and *networkx* etc..
The input and output data use the NetCDF format, an industry standard in
the geoscience field. Optional graphical outputs can also be saved,
making it suitable for production usage and educational purposes as well.
A series of [Jupyter notebooks](https://github.com/ihesp/AR_tracker/tree/master/notebooks)
are also included to help guide the users through the entire workflow and quickly
build their own production scripts.


# Example use case

The AR detection algorithm is inspired and modified from the image
processing algorithm *Top-hap by Reconstruction (THR)*, which consists
of subtracting from the original image a *greyscale reconstruction by
dilation* image.  Some more details of the THR algorithm and its
applications can be found in this work of [@Vincent1993].

In the context of AR detection, the greyscale image in question is the
non-negative IVT distribution, denoted as $I$.  The greyscale
reconstruction by dilation component corresponds to the background IVT
component, denoted as $\delta(I)$.  The difference $I - \delta(I)$
gives the transient IVT component, from which AR candidates are
searched. \autoref{fig:thr} shows this decomposition process.  It
could be seen that after this separation of back/transient components,
it comes trivial to locate AR-like features.

![(a) The IVT field in kg/m/s at 1984-01-26 00:00 UTC over the North Hemisphere. (b) the IVT reconstruction field ($\delta(I)$) at the same time point. (c) the IVT anomaly field ($I-\delta(I)$) from the THR process at the same time point.\label{fig:thr}](fig3.png)

After locating ARs at various time steps, an single curve is sought for
each AR as a summary of its location. A directed planar graph model
is used in this process, and weighted Dijkstra path searching
algorithm is used to find this "AR axis". Further details can be found
in the [documentation page](https://ar-tracker.readthedocs.io/en/latest/Find-AR-axis.html).


Lastly, a [modified Hausdorff distance definition](https://ar-tracker.readthedocs.io/en/latest/Track-ARs.html) is used as inter-AR
distance estimate, and an exclusive nearest neighbor approach is used to link
ARs at consecutive time points. \autoref{fig:track} shows an
example of this tracking process. A black-to-yellow color scheme is used to
indicate this AR's evolution.

![Locations of a track labelled "198424" found in year 1984. Black to yellow color scheme indicates the evolution.\label{fig:track}](ar_track_198424.png)


# External libraries used

Manipulation of the NetCDF data is achieved using the *CDAT* [@CDAT],
*numpy* [@numpy] and *scipy* [@scipy] packages.
Detection process utilizes the image-processing package *scikit-image* [@scikit-image].
AR axis finding process utilizes the *networkx* package [@networkx].
Generated outputs are further manipulated with *pandas* [@pandas] and
displayed using *matplotlib* [@matplotlib].


# Acknowledgements

This work is supported by the National Key Research and Development Program of
China with grant NO. 2017YFC1404000.
This work is also supported by the National Science Foundation of China (NSFC) with
grant No. 41490644 (41490640). This is a collaborative project between the
Ocean University of China (OUC), Texas A\&M University (TAMU) and the National
Center for Atmospheric Research (NCAR) and completed through the International
Laboratory for High Resolution Earth System Prediction (iHESP)- a collaboration
by the Qingdao National Laboratory for Marine Science and Technology
Development Center, Texas A&M University, and the National Center for
Atmospheric Research.  We also received support from the
NSFC grant \#41776013, and National Key Research and Development Program
of China grant 2017YFC1404100 and 2017YFC1404101.

# References
