import src.constants as cn
from src.LTIModel import LTIModel
import src.sympyUtil as ut

import numpy as np
import pandas as pd
import sympy
import unittest


IGNORE_TEST = True
IS_PLOT = True


#############################
# Tests
#############################
class TestLTIModel(unittest.TestCase):

    def setUp(self):
Amat = sympy.Matrix( [[0, 0, 0 ], [k1, -k2, 0], [0, k2, -k3]])

    def testConstructor(self):
        if IGNORE_TEST:
            return



if __name__ == '__main__':
  unittest.main()
