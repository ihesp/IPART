'''Test installation of ipart.
'''


#--------Import modules-------------------------
from __future__ import print_function
import unittest

class TestIPRT(unittest.TestCase):

    def test_ipart(self):

        from ipart.thr import THR, rotatingTHR
        from ipart.utils import funcs, plot
        from ipart.AR_detector import _findARs, findARs, findARsGen,\
                prepareMeta, findARAxis
        from ipart.AR_tracer import AR, trackARs, filterTracks

if __name__=='__main__':

    unittest.main()


