import src.sympyUtil as ut
import src.constants as cn

import numpy as np
import pandas as pd
import sympy
import unittest


IGNORE_TEST = False
IS_PLOT = False
VARIABLES = "X Y Z"




#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

    def setUp(self):
        names = VARIABLES.split(" ")
        dct = globals()
        for name in names:
            dct[name] = sympy.Symbol(name)

    def testSubstitute(self):
        if IGNORE_TEST:
            return
        Y = ut.substitute(2*X + 1, subs={X: Z})
        self.assertTrue("Z" in str(Y))

    def testEvaluate(self):
        if IGNORE_TEST:
            return
        val = ut.evaluate(2*X + 1, subs={X: 2})
        self.assertEqual(val, 5)

    def testEvaluate2(self):
        if IGNORE_TEST:
            return
        expr = sympy.Matrix( [2*Z, Z**2])
        val = ut.evaluate(expr, subs={Z: 2})
        self.assertEqual(expr.rows, np.shape(val)[0])
        self.assertEqual(expr.cols, np.shape(val)[1])


if __name__ == '__main__':
  unittest.main()
