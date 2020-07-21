'''Test detection in river_tracker1.
'''


#--------Import modules-------------------------
from __future__ import print_function
import os
import numpy as np
import unittest
from ipart.utils import funcs2 as funcs
from ipart.thr import THR
from ipart.AR_detector2 import findARs, _findARs, prepareMeta


class TestDetect(unittest.TestCase):
    """Test input netcdf data"""


    def setUp(self):
        thisdir=os.path.dirname(__file__)
        fixture_dir=os.path.join(thisdir, 'fixtures')
        abpath_in=os.path.join(fixture_dir,'uflux_vflux_ivt_test.nc')

        uflux=funcs.readNC(abpath_in, 'uflux')
        vflux=funcs.readNC(abpath_in, 'vflux')
        ivt=funcs.readNC(abpath_in, 'ivt')

        _,ivtrec,ivtano=THR(ivt, [3,3,3], verbose=False)
        self.ivt=ivt.squeeze()
        self.ivtrec=ivtrec.squeeze()
        self.ivtano=ivtano.squeeze()
        self.uflux=uflux.squeeze()
        self.vflux=vflux.squeeze()

        self.latax=ivt.getLatitude()
        self.lonax=ivt.getLongitude()
        timeax=ivt.getTime()
        self.timeax=timeax

        self.param_dict={
            # kg/m/s, define AR candidates as regions >= than this anomalous ivt.
            'thres_low' : 1,
            # km^2, drop AR candidates smaller than this area.
            'min_area': 50*1e4,
            # km^2, drop AR candidates larger than this area.
            'max_area': 1800*1e4,
            # float, isoperimetric quotient. ARs larger than this (more circular in shape) is treated as relaxed.
            #'max_isoq': 0.6,
            # float, isoperimetric quotient. ARs larger than this is discarded.
            #'max_isoq_hard': 0.7,
            'min_LW': 2,
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
            'single_dome': False,
            'max_ph_ratio': 0.4,
            # minimal proportion of flux component in a direction to total flux to
            # allow edge building in that direction
            'edge_eps': 0.4,
            'zonal_cyclic': False,
            }

    def test_Metadata(self):

        timeax, areamap, costhetas, sinthetas, lats, lons, reso = prepareMeta(
                self.latax, self.lonax, self.timeax,
                self.ivt.shape[0], self.ivt.shape[1], self.ivt.shape[2])

        t1=timeax[0]
        t1_str='%04d-%02d-%02d %02d:%02d:%02d' %(t1.year, t1.month, t1.day,
                t1.hour, t1.minute, t1.second)
        t2=timeax[10]
        t2_str='%04d-%02d-%02d %02d:%02d:%02d' %(t2.year, t2.month, t2.day,
                t2.hour, t2.minute, t2.second)

        self.assertEqual(t1_str, '1984-02-01 00:00:00', "Time axis wrong.")
        self.assertEqual(t2_str, '1984-02-03 12:00:00', "Time axis wrong.")

        self.assertAlmostEqual(np.min(areamap), 1267.4345, 3,
                "Min value in areamap wrong.")
        self.assertAlmostEqual(np.max(areamap), 6838.4643, 3,
                "Max value in areamap wrong.")
        self.assertAlmostEqual(np.min(costhetas), 0.17928, 3,
                "Min value in costhetas wrong.")
        self.assertAlmostEqual(np.max(costhetas), 0.70111, 3,
                "Max value in costhetas wrong.")

    def test_findARs(self):

        # find ARs
        time_idx, labels, angles, crossfluxes, result_df = findARs(
                self.ivt.data, self.ivtrec.data, self.ivtano.data,
                self.uflux.data, self.vflux.data, self.latax, self.lonax,
                self.timeax, **self.param_dict)

        self.assertEqual(len(result_df), 51, "Wrong number of ARs found.")
        self.assertTrue(np.all(time_idx==np.arange(25)), msg="time_idx wrong.")
        self.assertEqual(labels.data.sum(), 24870, "Labels sum doesn't match.")

    def test_findARsinner(self):

        idx=15
        slabano=self.ivtano.data[idx].squeeze()

        timeax, areamap, costhetas, sinthetas, lats, lons, reso = prepareMeta(
                self.latax, self.lonax, self.timeax,
                self.ivt.shape[0], self.ivt.shape[1], self.ivt.shape[2])

        # find ARs
        mask_list,armask=_findARs(slabano, lats, areamap, self.param_dict)

        self.assertEqual(len(mask_list), 1, "Wrong number of ARs found.")
        self.assertAlmostEqual(armask.sum(), 276, delta=30,
                msg="AR masks not right.")

if __name__=='__main__':

    unittest.main()




