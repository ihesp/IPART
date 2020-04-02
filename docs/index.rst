.. AR tracker documentation master file, created by
   sphinx-quickstart on Sun Mar 29 17:13:48 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AR tracker's documentation
==========================

.. contents:: Table of Contents
  :local:


Introduction
############

**AR tracker** is a Python package for automated Atmospheric River (AR)
detection, axis finding and AR tracking from gridded Integrated Vapor Transport (IVT) data, for
instance `Reanalysis
datasets <https://www.esrl.noaa.gov/psd/data/gridded/reanalysis/>`_.

An overview of what ARs are can be found in this review paper: `Atmospheric
rivers: a mini-review <https://doi.org/10.3389/feart.2014.00002>`_.

Different from commonly used AR detection methods that rely on thresholding on
the IVT magnitude, this package includes a method inspired by an image
processing technique -- `Top-hap by reconstruction (THR) <https://ieeexplore.ieee.org/document/217222>`_.

Below is an example output figure:

.. figure:: ar_1984-01-04_06:00.png
    :width: 700px
    :align: center
    :figclass: align-center

    (a) The IVT field in kg/m/s at 1984-01-04 06:00 UTC over the North
    Hemisphere. (b) the IVT reconstruction field at the same time point. (c)
    the IVT anomaly field from the THR process at the same time point. In all
    three subplots, the detected ARs are outlined in black contour. The AR axes
    are drawn in green dashed lines.


Installation
############

This package at the moment does not support installation via ``pip`` or ``conda``.
Users are advised to build a working Python virtual environment with required `Dependencies`_, and obtain
the code of this pacakge from the `github page <https://github.com/ihesp/AR_tracker>`_. 




Recommend building the Python environment using `Anaconda <https://www.anaconda.com/distribution/>`_.

After Anaconda installation, create a working environment:
::

    conda create -n ar python=3.7
    conda activate ar


This creates a Python environment named "ar" with Python version 3.7.

Then install the above listed dependencies, e.g.
::

    conda install numpy


This installs the current latest version of ``numpy``. Most likely the latest versions will work, in case of compatibility issues, consider forcing a given version of a package, e.g.
::

    conda install matplotlib=2.2.3


For installation of ``CDAT``, checkout the `installation guides <https://github.com/CDAT/cdat/wiki/Install>`_. This is likely the most difficult package to install, and consider leaving it to the end. To verify the CDAT installation, in a python session:
::

    import cdms2
    import MV2
    import cdutil


If nothing prints out, the installation is successful. In case of errors, also consider their `partial installation instructions <https://github.com/CDAT/cdat/wiki/Additional-Installation-Configurations>`_. Only the ``cdms2`` and ``cdutil`` modules are needed, the ``vcs`` module is not required.

Lastly, get a copy of this repository, either by doing a ``git clone``:
::

    git clone git@github.com:ihesp/AR_tracker.git


or download a ``zip`` archive of the repository by clicking on the **Clone or download** button, then extract the contents of the ``zip`` file into a folder named ``AR_tracker``.


Dependencies
############

* OS: Linux or MacOS. Windows is likely not going to work, due to the lack of a Windows version of the CDAT package.
* Python2.7 or Python3.7.
* CDAT `Climatic Data Analysis Tool <https://github.com/CDAT/cdat>`_.
* numpy (developed in 1.16.5 in py2, tested 1.18.1 in py3)
* scipy (developed in 1.2.1 in py2, tested 1.4.1 in py3)
* matplotlib (2.2.3 for both py2 and py3, having `issues <https://github.com/matplotlib/matplotlib/issues/12820>`_ with 3.1.3.
* basemap (developed in 1.2.0 in py2, tested 1.2.0 in py3)
* pandas (developed in 0.23.4 in py2, tested 1.0.3 in py3)
* networkx (developed in 1.11 and 2.2 in py2, tested 2.4 in py3)
* scikit-image (developed in 0.14.2 in py2, tested 0.16.2 in py3)



Main functionalities
####################

There are four main functionalities provided by the package that collectively constitute a
specific workflow for the automated AR detection/tracking task:

1. Perform THR computation on input data.
2. Detect ARs from the outputs from the previous step, and at the same time,
3. Identify the AR axis.
4. Track ARs detected at individual time steps to form tracks.


More details regarding these steps are provided in separate pages below.

Applications on example data can be found in a series of example notebooks at
`github repository <https://github.com/ihesp/AR_tracker/blob/master/notebooks/Index.ipynb>`_.


.. toctree::
   :maxdepth: 2
   :caption: The automated AR detect/tracking Workflow:

   Data-preparation
   Compute-THR
   Detect-ARs
   Find-AR-axis
   Track-ARs



Main scripts/modules
####################
   

.. toctree::
   :maxdepth: 1

   compute_thr_singlefile.py <compute_thr_singlefile>
   compute_thr_multifile.py <compute_thr_multifile>
   river_tracker1.py <river_tracker1>
   river_tracker2.py <river_tracker2>
   river_tracker1_funcs.py <river_tracker1_funcs>
   utils.funcs (selected parts) <funcs>
   utils.plot (selected parts) <plot>




Github and Contact
##################

The code of this package is hosted at https://github.com/ihesp/AR_tracker.

For any queries, please contact xugzhi1987@gmail.com.


License
#######

:ref:`license`



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
