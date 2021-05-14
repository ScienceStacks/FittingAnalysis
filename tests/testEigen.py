import src.sympyUtil as su
from src.eigenCollection import EigenCollection
import src.constants as cn

import numpy as np
import pandas as pd
import sympy
import unittest


IGNORE_TEST = True
IS_PLOT = True
VARIABLES = "k0 k1 k2"
su.addSymbols(VARIABLES, dct=globals())
# No. eigenvalues: 2, Algebraic multiplicity: 1, Geometric multiplicity: 1
FULL_MAT = sympy.Matrix([
      [1, 2],
      [2, 1],
      ])
# No. eigenvalues: 1, Algebraic multiplicity: 2, Geometric multiplicity: 1
DEFICIENT_MAT = sympy.Matrix([
      [1, 0],
      [2, 1],
      ])
SUBS = {k0: 1, k1: 2, k2: 3}


#############################
# Tests
#############################
class TestEigenCollection(unittest.TestCase):

    def setUp(self):
        su.addSymbols(VARIABLES, dct=globals())
        self.aMat = FULL_MAT.copy()
        self.eigenCollection = EigenCollection(self.aMat)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(len(self.eigenCollection.eigenInfos), 2)
        self.assertEqual(self.eigenCollection.eigenInfos[0].val, -1.0)
        self.assertEqual(self.eigenCollection.eigenInfos[1].val, 3.0)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        eigen = EigenCollection(DEFICIENT_MAT)
        self.assertEqual(len(eigen.eigenInfos), 1)
        eigenInfo = eigen.eigenInfos[0]
        self.assertEqual(eigenInfo.val, 1.0)
        self.assertEqual(eigenInfo.mul, 2)


class TestEigenInfo(unittest.TestCase):

    def setUp(self):
        su.addSymbols(VARIABLES, dct=globals())
        aMat = FULL_MAT.copy()
        eigenCollection = EigenCollection(self.aMat)
        self.eigenInfo = eigenCollection[0]

    def testConstructor(self)
        # TESTING
        import pdb; pdb.set_trace()


if __name__ == '__main__':
  unittest.main()
