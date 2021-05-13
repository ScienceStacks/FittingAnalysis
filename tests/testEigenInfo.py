import src.sympyUtil as su
from src.eigenInfo import EigenInfo
import src.constants as cn

import numpy as np
import pandas as pd
import sympy
import unittest


IGNORE_TEST = False
IS_PLOT = False
VARIABLES = "k0 k1 k2"
su.addSymbols(VARIABLES, dct=globals())
SUBS = {k0: 1, k1: 2, k2: 3}


#############################
# Tests
#############################
class TestEigenInfo(unittest.TestCase):

    def setUp(self):
        su.addSymbols(VARIABLES, dct=globals())
        self.aMat = sympy.Matrix([ 
                              [   0,           0,   0,  0, 0],
                              [   0,           0,   0,  0, 1],
                              [   1,  -(k0 + k2),   0,  0, 0], 
                              [   0,          k0, -k1,  0, 0], 
                              [   0,          k2,  k1, -1, 0],
                            ])
        aMat = su.evaluate(self.aMat, subs=SUBS, isNumpy=False)
        self.eigenInfo = EigenInfo(aMat)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        import pdb; pdb.set_trace()
        


if __name__ == '__main__':
  unittest.main()
