import src.constants as cn
from src.LTIModel import LTIModel
import src.sympyUtil as su

import numpy as np
import pandas as pd
import sympy
import unittest


IGNORE_TEST = False
IS_PLOT = False


#############################
# Tests
#############################
class TestLTIModel(unittest.TestCase):

    def setUp(self):
        su.addSymbols("aMat mat k1 k2 k3 t X_0")
        aMat = sympy.Matrix( [[0, 0, 0 ], [k1, -k2, 0], [0, k2, -k3]])
        initialVec = sympy.Matrix([X_0, 0, 0])
        self.model = LTIModel(aMat, initialVec)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(self.model.aMat.rows, 3)



if __name__ == '__main__':
  unittest.main()
