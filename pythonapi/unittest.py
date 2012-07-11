import unittest
import numpy as np

class Test_imageproc(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        pass

    def test_disttransform():




if __name__ == '__main__':
    improcsuite = unittest.TestLoader().loadTestsFromTestCase(Test_imageproc)

    suite=unittest.TestSuite([improcsuite])
    unittest.TextTestRunner(verbosity=2).run(suite)
