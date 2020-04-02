'''Test THR process.
'''


#--------Import modules-------------------------
from __future__ import print_function
import os
import unittest
import numpy as np
from utils import funcs
from compute_thr_singlefile import filterData


class TestTHR(unittest.TestCase):
    """Test input netcdf data"""


    def setUp(self):
        thisdir=os.path.dirname(__file__)
        fixture_dir=os.path.join(thisdir, 'fixtures')
        abpath_in=os.path.join(fixture_dir,'uflux_vflux_ivt_test.nc')

        ivt=funcs.readVar(abpath_in, 'ivt')
        self.ivt=ivt

    def test_THR(self):
        ivt,ivtrec,ivtano=filterData(self.ivt, [3,3,3], verbose=False)
        self.assertTupleEqual(self.ivt.shape, ivt.shape, "Output <ivt> shape worng.")
        self.assertTupleEqual(self.ivt.shape, ivtrec.shape, "Output <ivtrec> shape worng.")
        self.assertTupleEqual(self.ivt.shape, ivtano.shape, "Output <ivtano> shape worng.")

        #------Test the number of pixels >0 in ivtano------
        nums=[np.sum(ivtano[ii]>0) for ii in range(len(ivtano))]
        ans=[1966, 1424, 1353, 1260, 1357, 1323, 1412, 1220, 1140, 1150, 1195, 1191,
            1275, 1171, 1153, 915, 935, 942, 848, 743, 918, 1130, 1212, 1274, 1660]

        self.assertTrue(np.all(nums==ans), "Output <ivtano> values wrong.")



if __name__=='__main__':

    unittest.main()



