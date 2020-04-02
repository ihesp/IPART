'''Test tracking in river_tracker2.
'''


#--------Import modules-------------------------
from __future__ import print_function
import os
import unittest
import numpy as np
import pandas as pd
from compute_thr_singlefile import filterData
from river_tracker2 import convarray, trackARs2, filterTracks


class TestDetect(unittest.TestCase):
    """Test input netcdf data"""


    def setUp(self):

        thisdir=os.path.dirname(__file__)
        fixture_dir=os.path.join(thisdir, 'fixtures')
        abpath_in=os.path.join(fixture_dir,'sample_record.csv')

        convkeys=['contour_y', 'contour_x',
                'axis_y', 'axis_x', 'axis_rdp_y', 'axis_rdp_x']

        converters=dict([(keyii, convarray) for keyii in convkeys])

        dtypes={'id': 'int', 'time': 'str',
                'area': np.float64, 'length': np.float64, 'width': np.float64,
                'iso_quotient': np.float64, 'LW_ratio': np.float64,
                'strength': np.float64, 'strength_ano': np.float64,
                'strength_std': np.float64,
                'mean_angle': np.float64, 'is_relaxed': 'bool'}

        self.ardf=pd.read_csv(abpath_in,dtype=dtypes,converters=converters)

        self.track_params={}
        self.track_params['time_gap_allow']=6
        self.track_params['num_anchors']=7
        self.track_params['track_scheme']='simple'  # 'simple' | 'full'
        self.track_params['max_dist_allow']=1200  # km
        self.track_params['min_duration']=24
        self.track_params['min_nonrelax']=1


    def test_track(self):

        self.assertTupleEqual(self.ardf.shape, (99,22), "Input <ardf> shape wrong.")

        #-------------------Track ars-------------------
        track_list=trackARs2(self.ardf, self.track_params, isplot=False, verbose=False)
        self.assertEqual(len(track_list), 19, "Number of tracks wrong.")

        #------------------Filter tracks------------------
        track_list=filterTracks(track_list, self.track_params)
        self.assertEqual(len(track_list), 4, "Number of filtered tracks wrong.")

        nums=[len(track_list[ii].data) for ii in range(len(track_list))]
        ans=(12, 21, 11, 11)
        self.assertTupleEqual(tuple(nums), ans, "Track lengths wrong.")


if __name__=='__main__':

    unittest.main()





