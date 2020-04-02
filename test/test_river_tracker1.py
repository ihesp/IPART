'''Test detection in river_tracker1.
'''


#--------Import modules-------------------------
from __future__ import print_function
import os
import unittest
import MV2 as MV
from utils import funcs
from compute_thr_singlefile import filterData
from river_tracker1 import findARs
from river_tracker1_funcs import uvDecomp


class TestDetect(unittest.TestCase):
    """Test input netcdf data"""


    def setUp(self):
        thisdir=os.path.dirname(__file__)
        fixture_dir=os.path.join(thisdir, 'fixtures')
        abpath_in=os.path.join(fixture_dir,'uflux_vflux_ivt_test.nc')

        self.uflux=funcs.readVar(abpath_in, 'uflux')
        self.vflux=funcs.readVar(abpath_in, 'vflux')
        ivt=funcs.readVar(abpath_in, 'ivt')

        _,ivtrec,ivtano=filterData(ivt, [3,3,3], verbose=False)
        self.ivt=ivt
        self.ivtrec=ivtrec
        self.ivtano=ivtano

        dxs=funcs.dLongitude(self.uflux,R=6371)
        dys=funcs.dLatitude(self.uflux,R=6371)
        self.areamap=dxs*dys # km^2
        self.costhetas=dxs/MV.sqrt(dxs**2+dys**2)
        self.sinthetas=dys/MV.sqrt(dxs**2+dys**2)

        self.param_dict={
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
            'fill_radius': max(1,int(4*0.75/0.75)),
            # max prominence/height ratio of a local peak. Only used when SINGLE_DOME=True
            'max_ph_ratio': 0.4,
            # minimal proportion of flux component in a direction to total flux to
            # allow edge building in that direction
            'edge_eps': 0.4
            }

    def test_THR(self):

        idx=10
        quslab=self.uflux[idx](squeeze=1)
        qvslab=self.vflux[idx](squeeze=1)
        slabrec=self.ivtrec[idx](squeeze=1)
        slabano=self.ivtano[idx](squeeze=1)

        qurec,quano,qvrec,qvano=uvDecomp(quslab,qvslab,slabrec,slabano)

        # find ARs
        mask_list,axis_list,armask,axismask=findARs(slabano,quano,qvano,
                self.areamap,self.costhetas,self.sinthetas,self.param_dict)

        self.assertEqual(len(mask_list), 2, "Wrong number of ARs found.")
        self.assertEqual(len(axis_list), 2, "Wrong number of ARs found.")
        self.assertAlmostEqual(armask.sum(), 719, delta=50, msg="AR masks not right.")
        self.assertAlmostEqual(axismask.sum(), 73, delta=7, msg="AR axes not right.")

    def test_THR2(self):

        idx=15
        quslab=self.uflux[idx](squeeze=1)
        qvslab=self.vflux[idx](squeeze=1)
        slabrec=self.ivtrec[idx](squeeze=1)
        slabano=self.ivtano[idx](squeeze=1)

        qurec,quano,qvrec,qvano=uvDecomp(quslab,qvslab,slabrec,slabano)

        # find ARs
        mask_list,axis_list,armask,axismask=findARs(slabano,quano,qvano,
                self.areamap,self.costhetas,self.sinthetas,self.param_dict)

        self.assertEqual(len(mask_list), 1, "Wrong number of ARs found.")
        self.assertEqual(len(axis_list), 1, "Wrong number of ARs found.")
        self.assertAlmostEqual(armask.sum(), 276, delta=30, msg="AR masks not right.")
        self.assertAlmostEqual(axismask.sum(), 50, delta=5, msg="AR axes not right.")

if __name__=='__main__':

    unittest.main()




