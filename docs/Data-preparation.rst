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

The user should also make sure that the dimensions (time, latitude and longitude)
have proper axes associated with them. This is because these temporal and geographical information
is used in the computation.

To test the sanity of your input data, run this script against the netcdf data file:
::

    cd /path/to/AR_tracker/folder/you/cloned
    python test_data.py /path/to/your/netcdf/file 'id_of_variable'

For instance

::

    cd ~/Downloads/AR_tracker
    python test_data.py ~/datasets/erai/uflux_file.nc 'uflux'

The expected outputs would look like this:

::

        ### <test_data>: Read in file:
         /home/guangzhi/datasets/erai/erai_qflux/uflux_m1-60_6_1984_cln.nc
        ##################################################
        Variable time axis:
        ##################################################
        [1984-1-1 0:0:0.0, 1984-1-1 6:0:0.0, 1984-1-1 12:0:0.0, 1984-1-1
        18:0:0.0, 1984-1-2 0:0:0.0, 1984-1-2 6:0:0.0, 1984-1-2 12:0:0.0,
        ...
        1984-12-30 12:0:0.0, 1984-12-30 18:0:0.0, 1984-12-31 0:0:0.0,
        1984-12-31 6:0:0.0, 1984-12-31 12:0:0.0, 1984-12-31 18:0:0.0]

        ##################################################
        Variable latitude axis:
        ##################################################
           id: latitude
           Designated a latitude axis.
           units:  degrees_north
           Length: 241
           First:  -90.0
           Last:   90.0
           Other axis attributes:
              long_name: latitude
              axis: Y
           Python id:  0x7f7b68896250


        ##################################################
        Variable longitude axis:
        ##################################################
           id: longitude
           Designated a longitude axis.
           units:  degrees_east
           Length: 480
           First:  0.0
           Last:   359.25
           Other axis attributes:
              modulo: [360.]
              long_name: longitude
              axis: X
              topology: circular
           Python id:  0x7f7b688962d0


        Data have unit of "kg m**-1 s**-1"


Pay some attention to the values listed in the **latitude** and **longitude** axes
blocks, to make sure the values make physical sense. For high resolution data,
the input variable may have a fairly large size, e.g. a longitude axis of length 720 (if
your data have a resolution of :math:`0.5 \times 0.5 ^{\circ}`. If the longitude
axis reports a ``Last: 720``, it is probably reporting the size of the longitude
dimension, rather than the actual longitude label, as the maximum possible
longitude label should be 360. In such cases, the user should take some extra steps
to make sure that the data have proper metadata associated with them.



.. _get_ivt:

Get IVT
#######


With :math:`F_u` and :math:`F_v`, compute the IVT as

.. math::
    IVT = \sqrt{F_u^2 + F_v^2}


This is trivial to achieve, you can use the ``compute_ivt.py`` script provided in the package for this computation.


