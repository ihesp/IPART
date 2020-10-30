---
title: 'IPART: A Python Package for Image-Processing based Atmospheric River Tracking'
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

# Background

An atmospheric river (AR) in the field of meteorology/climatology
refers to enhanced water vapor content in the lower troposphere.  The
term was coined as an analogy to terrestrial rivers in a sense
that when viewed from satellite imagery or large scale atmospheric
observation, they appear as narrow and elongated vapor filaments,
representing transient intensified horizontal moisture fluxes
[e.g., @Gimeno2014; @Dettinger2011]. A typical atmospheric river
can carry 7-15 times the water in the Mississippi River
[@Ralph2011], and at any time in winter, there are four to five such
systems in the Northern Hemisphere alone [@Zhu1998], accounting for
80-90% of the total north-south integrated vapor transport
[@Guan2015; @Zhu1998].  Its dual hydrological role, both as a fresh
water source for some water-stressed areas [@Dettinger2011;
@Dettinger2013; @Rutz2012] and as a potential trigger for floods
[@Lavers2012; @Lavers2013; @Neiman2008; @Moore2012], has granted it
increasing attention among the research community.  Their
long-term change in a warming climate also stands as a pressing research
question. However, an important prerequisite to answer such questions
is a robust and consistent detection method. As meteorologists and
climatologists often deal with observational or simulation data in
large sizes, an algorithmic method can ensure better efficiency,
consistency and objectivity compared with human
identification.

In many existing applications, a magnitude thresholding approach is used. For
instance, @Ralph2004, @Neiman2008, @Hagos2015 and @Dettinger2011 identified ARs
by first locating regions where the Integrated Water Vapor (IWV) is greater
than 20 mm. A 250 kg/(m $\cdot$ s) threshold on the Integrated Vapor Transport
(IVT) was used by @Rutz2014 and @Rutz2015.  However, an implicit assumption
with this magnitude thresholding approach is that the atmospheric moisture
level stays unchanged throughout the analysis period.  Such an assumption may
not be fully justifiable under a warming climate as the atmospheric moisture
level is expected to increase.

# Summary

In this package we propose a suite of new detection/tracking algorithms to help
solve the above difficulties.  Through a systematic analysis using seven years of
reanalysis data [@Dee2011], we have found that the proposed detection algorithm
has reduced sensitivity to parameters [@Xu2020b]. Long-lived ARs
spanning multiple days, having cross-continent or cross-basin tracks, can be
more reliably traced through their tropical/sub-tropical origins to
high-latitude landfalls. As the research on ARs matures, new AR
detection/tracking methods are being developed, and the inter-comparisons of
various AR detection/tracking methods are carried out by, for instance,
the Atmospheric River Tracking Method Intercomparison Project [ARTMIP, @Rutz2019b;@Shields2018].
Using the terminology of ARTMIP [@Shields2018], the proposed method is a
"tracking" (Lagrangian approach) type, with length and shape geometrical
requirements. It imposes no threshold on IVT/IWV, but instead imposes absolute
thresholds on the spatio-temporal scale of AR-like systems. The detected ARs
can be optionally time-stitched to identify coherent AR objects.  We have
performed some systematic comparisons with two magnitude thresholding based AR
detection methods included in ARTMIP, and the proposed method displays better
correspondence between North Hemisphere AR tracks and storm tracks, better
identification of the strong mid-latitude AR-related moisture transports, and
longer AR track durations. The detailed comparison analysis is given in @Xu2020, and
a more thorough description of the detection/tracking methods is given in @Xu2020b.

`IPART` is therefore intended for researchers and students who are interested
in the field of atmospheric river studies in the present day climate or future
projections.

The `IPART` package includes a collection of Python functions and classes designed
for an analysis workflow covering the detection of ARs, the simplification of
the AR's geographical location, to the subsequent tracking through time.  The
algorithms are implemented using the Python programming language as a wrapper
to some well-established numeric packages including `numpy`, `scikit-image` and
`networkx` etc.  The input and output data use the NetCDF format, an industry
standard in the geoscience field. Optional graphical outputs can also be saved,
making it suitable for production usage and educational purposes as well.  A
series of [Jupyter notebooks](https://github.com/ihesp/IPART/tree/master/notebooks) are also
included to help guide the users through the entire workflow, and some example
scripts are provided as templates to help the user quickly build their own
production scripts.



# Example use case

The AR detection algorithm is inspired and modified from the image
processing algorithm *Top-hap by Reconstruction (THR)*, which consists
of subtracting from the original image a *greyscale reconstruction by
dilation* image [@Vincent1993].

In the context of AR detection, the greyscale image in question is the
non-negative IVT distribution, denoted as $I$.  The greyscale
reconstruction by dilation component corresponds to the background IVT
component, denoted as $\delta(I)$.  The difference $I - \delta(I)$
gives the transient IVT component, from which AR candidates are
searched. \autoref{fig:thr} shows this decomposition process.  It
could be seen that after this separation of background/transient components,
it becomes trivial to locate AR-like features.

![(a) The IVT field in kg/(m $\cdot$ s) at 1984-01-26 00:00 UTC over the North Hemisphere. (b) the IVT reconstruction field ($\delta(I)$) at the same time point. (c) the IVT anomaly field ($I-\delta(I)$) from the THR process at the same time point.\label{fig:thr}](fig3.png)

After locating ARs at various time steps, a single curve is sought for
each AR as a summary of its location. A directed planar graph model
is used in this process, and weighted Dijkstra path searching
algorithm is used to find this "AR axis". Further details can be found
in the [documentation page](https://ar-tracker.readthedocs.io/en/latest/Find-AR-axis.html).


Lastly, a [modified Hausdorff distance definition](https://ar-tracker.readthedocs.io/en/latest/Track-ARs.html) is used as an inter-AR
distance estimate, and an exclusive nearest neighbor approach is used to link
ARs at consecutive time points. \autoref{fig:track} shows an
example of this tracking process.

![Locations of a track labelled "198424" found in year 1984. A color scheme
of black to yellow through purple indicates the evolution,
where black curves represent the AR at earlier times and yellow curves at
later times.\label{fig:track}](ar_track_198424.png)


# External libraries used

Manipulation of the NetCDF data is achieved using the Python interface of the
`NetCDF` software [@netcdf], `numpy` [@numpy] and `scipy`
[@scipy] packages.  The detection process utilizes the image-processing package
`scikit-image` [@scikit-image].  The AR axis finding process utilizes the
`networkx` package [@networkx].  Generated outputs are further manipulated with
`pandas` [@pandas] and displayed using `matplotlib` [@matplotlib].


# Acknowledgements

This work is supported by the National Key Research and Development Program of
China with grant NO. 2017YFC1404000.
This work is also supported by the National Science Foundation of China (NSFC) with
grant No. 41490644 (41490640). This is a collaborative project between the
Ocean University of China (OUC), Texas A\&M University (TAMU) and the National
Center for Atmospheric Research (NCAR) and completed through the International
Laboratory for High Resolution Earth System Prediction (iHESP) - a collaboration
by the Qingdao National Laboratory for Marine Science and Technology
Development Center, Texas A&M University, and NCAR.  We also received support from the
NSFC grant \#41776013, and National Key Research and Development Program
of China grant 2017YFC1404100 and 2017YFC1404101.

# References
