'''Test input netcdf data.
'''


#--------Import modules-------------------------
from __future__ import print_function
import os
import unittest
import numpy as np
from ipart.utils import funcs


class TestDataMetaData(unittest.TestCase):
    """Test input netcdf data"""


    def setUp(self):
        thisdir=os.path.dirname(__file__)
        fixture_dir=os.path.join(thisdir, 'fixtures')
        abpath_in=os.path.join(fixture_dir,'uflux_vflux_ivt_test.nc')

        uflux=funcs.readVar(abpath_in, 'uflux')
        vflux=funcs.readVar(abpath_in, 'vflux')
        ivt=funcs.readVar(abpath_in, 'ivt')

        self.uflux=uflux
        self.vflux=vflux
        self.ivt=ivt

    def test_dataShape(self):
        self.assertTupleEqual(self.uflux.shape, (25,1,93,213))
        self.assertTupleEqual(self.vflux.shape, (25,1,93,213))
        self.assertTupleEqual(self.ivt.shape, (25,1,93,213))

    def test_timeAxis(self):
        timeax=self.uflux.getTime()
        timeax_cp=timeax.asComponentTime()
        t1=timeax[0]
        t2=timeax[-1]
        tc1=timeax_cp[0]
        tc2=timeax_cp[-1]
        self.assertEqual(len(timeax), 25, "Lenght of time axis wrong.")
        self.assertEqual(t1, 30711.0, "Start time point not correct.")
        self.assertEqual(t2, 30717.0, "End time point not right.")
        self.assertEqual(str(tc1), '1984-2-1 0:0:0.0', "Start time point not correct.")
        self.assertEqual(str(tc2), '1984-2-7 0:0:0.0', "End time point not correct.")

    def test_latAxis(self):
        latax=self.uflux.getLatitude()
        self.assertIsNotNone(latax, "No latitude axis found.")
        self.assertEqual(len(latax), 93, "Latitude length wrong.")
        self.assertAlmostEqual(latax[0], 10.5, 5, "First value in latitude wrong.")
        self.assertAlmostEqual(latax[-1], 79.5, 5, "Last value in latitude wrong.")

    def test_lonAxis(self):
        lonax=self.uflux.getLongitude()
        self.assertIsNotNone(lonax, "No longitude axis found.")
        self.assertEqual(len(lonax), 213, "Longitude length wrong.")
        self.assertAlmostEqual(lonax[0], 100.5, 5, "First value in longitude wrong.")
        self.assertAlmostEqual(lonax[-1], 259.5, 5, "Last value in longitude wrong.")

    def test_dlat(self):
        dlat=funcs.dLatitude(self.uflux)
        self.assertTrue(np.allclose(dlat, 83396.1949, rtol=1e-3), "dLatitude() wrong.")

    def test_dlon(self):
        dlon=funcs.dLongitude(self.uflux)
        self.assertAlmostEqual(np.min(dlon), 15197.74930, 4, "dLongitude() wrong.")
        self.assertAlmostEqual(np.max(dlon), 81999.71815, 4, "dLongitude() wrong.")


if __name__=='__main__':

    unittest.main()


