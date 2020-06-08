'''Test THR process.
'''


#--------Import modules-------------------------
from __future__ import print_function
import os
import unittest
import numpy as np
from ipart.utils import funcs
from ipart.thr import THR


class TestTHR(unittest.TestCase):
    """Test input netcdf data"""


    def setUp(self):
        thisdir=os.path.dirname(__file__)
        fixture_dir=os.path.join(thisdir, 'fixtures')
        abpath_in=os.path.join(fixture_dir,'uflux_vflux_ivt_test.nc')

        ivt=funcs.readVar(abpath_in, 'ivt')
        self.ivt=ivt

    def test_THR(self):
        ivt,ivtrec,ivtano=THR(self.ivt, [3,3,3], verbose=False)
        self.assertTupleEqual(self.ivt.shape, ivt.shape, "Output <ivt> shape worng.")
        self.assertTupleEqual(self.ivt.shape, ivtrec.shape, "Output <ivtrec> shape worng.")
        self.assertTupleEqual(self.ivt.shape, ivtano.shape, "Output <ivtano> shape worng.")

        #------Test the number of pixels >0 in ivtano------
        nums=[np.sum(ivtano[ii]>0) for ii in range(len(ivtano))]
        ans=[2201, 1608, 1544, 1475, 1532, 1562, 1622, 1402, 1316, 1317, 1402,
                1394, 1482, 1358, 1318, 1078, 1090, 1152, 1064, 944, 1080,
                1285, 1386, 1399, 1848]

        self.assertTrue(np.all(nums==ans), "Output <ivtano> values wrong.")



if __name__=='__main__':

    unittest.main()



