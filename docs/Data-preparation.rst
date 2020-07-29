Data preparation
================

netcdf data
###########

Source data are the u- and v- components of the vertically integrated vapor fluxes, in a rectangular
grid.

These are usually computed as:

.. math::

    \left\{\begin{matrix}
    F_u & = \int_{TOA}^{P_s} \frac{u q }{g} dP \\
    F_v & = \int_{TOA}^{P_s} \frac{v q }{g} dP
    \end{matrix}\right.

where

* :math:`F_u` (:math:`F_v`) is the zonal (meridional) component of the integrated flux, both
  in :math:`kg/m/s`.
* :math:`u` (:math:`v`): is the zonal (meridional) wind speed (in :math:`m/s`) at a given level, and
  :math:`q` is the specific humidity (in `kg/kg`) at the same level.
* :math:`dP` is the pressure increment, in :math:`Pa`, and :math:`g` is acceleration by gravity.
* :math:`TOA` can be substituted with a sensibly high level, e.g. :math:`300 hPa`.

No strict requirement for the spatial or temporal resolution of input data is imposed, however, for
better results, one should use something greater than :math:`\sim 1.5 ^{\circ}` latitude/longitude.
6-hourly temporal resolution is a standard for `Reanalysis
datasets <https://www.esrl.noaa.gov/psd/data/gridded/reanalysis/>`_, daily would probably work, but
one should do some parameter adjustments in such a case.

Data are supposed to saved in `netcdf files <https://www.unidata.ucar.edu/software/netcdf/docs/index.html>`_.

.. _metadata:

Metadata
########


.. note:: the user is responsible for making sure that the data are saved in the following rank order:

::

    (time, level, latitude, longitude)

or::

    (time, latitude, longitude)

The ``level`` dimension is optional. As data are vertical integrals, the length
of the level dimension, if present, should be 1.

The user also needs to provide the time, latitude and longitude axes values.
This is because these temporal and geographical information is used in the computation.

To test the sanity of your input data, run this script against the netcdf data file:
::

    cd /path/to/IPART/folder/you/cloned
    python test_data.py /path/to/your/netcdf/file 'id_of_variable'

For instance

::

    cd ~/Downloads/IPART
    python test_data.py ~/datasets/erai/uflux_file.nc 'uflux'

The expected outputs would look like this:

::

        ### <test_data>: Read in file:
         /home/guangzhi/scripts/IPART/notebooks/uflux_s_6_1984_Jan.nc
        ##################################################
        Variable time axis:
        ##################################################
        ### Description of slab ###
          id: time
          shape: (124,)
          filename: None
          missing_value: None
          comments: None
          grid_name: None
          grid_type: None
          long_name: time
          units: hours since 1900-01-01 00:00:00.0
          standard_name: None
          Order: []

        ### End of description ###

        None
        Time axis values:
        [datetime.datetime(1984, 1, 1, 0, 0) datetime.datetime(1984, 1, 1, 6, 0)
         datetime.datetime(1984, 1, 1, 12, 0) datetime.datetime(1984, 1, 1, 18, 0)
         datetime.datetime(1984, 1, 2, 0, 0) datetime.datetime(1984, 1, 2, 6, 0)
         datetime.datetime(1984, 1, 2, 12, 0) datetime.datetime(1984, 1, 2, 18, 0)
         datetime.datetime(1984, 1, 29, 0, 0) datetime.datetime(1984, 1, 29, 6, 0)
         ...
         datetime.datetime(1984, 1, 31, 12, 0)
         datetime.datetime(1984, 1, 31, 18, 0)]

        ##################################################
        Variable latitude axis:
        ##################################################
        ### Description of slab ###
          id: latitude
          shape: (94,)
          filename: None
          missing_value: None
          comments: None
          grid_name: None
          grid_type: None
          long_name: latitude
          units: degrees_north
          standard_name: None
          Order: []

        ### End of description ###

        None
        Latitude axis values:
        [10.   10.75 11.5  12.25 13.   13.75 14.5  15.25 16.   16.75 17.5  18.25
         19.   19.75 20.5  21.25 22.   22.75 23.5  24.25 25.   25.75 26.5  27.25
         28.   28.75 29.5  30.25 31.   31.75 32.5  33.25 34.   34.75 35.5  36.25
         37.   37.75 38.5  39.25 40.   40.75 41.5  42.25 43.   43.75 44.5  45.25
         46.   46.75 47.5  48.25 49.   49.75 50.5  51.25 52.   52.75 53.5  54.25
         55.   55.75 56.5  57.25 58.   58.75 59.5  60.25 61.   61.75 62.5  63.25
         64.   64.75 65.5  66.25 67.   67.75 68.5  69.25 70.   70.75 71.5  72.25
         73.   73.75 74.5  75.25 76.   76.75 77.5  78.25 79.   79.75]

        ##################################################
        Variable longitude axis:
        ##################################################
        ### Description of slab ###
          id: longitude
          shape: (480,)
          filename: None
          missing_value: None
          comments: None
          grid_name: None
          grid_type: None
          long_name: longitude
          units: degrees_east
          standard_name: None
          Order: []

        ### End of description ###

        None
        Longitude axis values:
        [-180.   -179.25 -178.5  -177.75 -177.   -176.25 -175.5  -174.75 -174.
         -173.25 -172.5  -171.75 -171.   -170.25 -169.5  -168.75 -168.   -167.25
         -166.5  -165.75 -165.   -164.25 -163.5  -162.75 -162.   -161.25 -160.5
         ...
          164.25  165.    165.75  166.5   167.25  168.    168.75  169.5   170.25
          171.    171.75  172.5   173.25  174.    174.75  175.5   176.25  177.
          177.75  178.5   179.25]

        Data have unit of "kg m**-1 s**-1"


Pay some attention to the values listed in the **latitude** and **longitude**
axes blocks, to make sure the values make physical sense. For high resolution
data, the input variable may have a fairly large size, e.g. a longitude axis of
length 720 (if your data have a resolution of :math:`0.5 \times 0.5 ^{\circ}`).
If the longitude axis reports a largest value of 720, it is probably reporting
the size of the longitude dimension, rather than the actual longitude label (as
the maximum possible longitude label should be 360). In such cases, the user
should take some extra steps to make sure that the data have proper metadata
associated with them.



.. _get_ivt:

Compute IVT
###########


With :math:`F_u` and :math:`F_v`, compute the IVT as

.. math::
    IVT = \sqrt{F_u^2 + F_v^2}


This is trivial to achieve, you can use the ``compute_ivt.py`` script provided
in the package for this computation.


