import unittest
from pyeegtools import *
import numpy as np
import pylab as pl

class Test_imageproc(unittest.TestCase):

    def setUp(self):
        pass

    def test_disttransform(self):
        z=np.arange(1,100,dtype=np.int)
        self.assertRaises(ValueError, disttransform_deadreckoning, z)
        z2=np.zeros((100,100),dtype=np.int32)
        z2[20,20]=1;
        dt=disttransform_deadreckoning(z2)
        pl.imshow(dt)
        pl.show()


if __name__ == '__main__':
    improcsuite = unittest.TestLoader().loadTestsFromTestCase(Test_imageproc)

    suite=unittest.TestSuite([improcsuite])
    unittest.TextTestRunner(verbosity=2).run(suite)
