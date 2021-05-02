import src.sympyUtil as ut
import src.constants as cn

import numpy as np
import pandas as pd
import sympy
import unittest


IGNORE_TEST = True
IS_PLOT = True
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
        pass

    def testSubstitute(self):
        if IGNORE_TEST:
            return
        Y = ut.substitute(2*X + 1, subs={X: Z})
        self.assertTrue("Z" in str(Y))

    def testEvaluate(self):
        # TESTING
        val = ut.evaluate(2*X + 1, subs={X: 2})
        self.assertEqual(val, 5)

    def testEvaluate2(self):
        if IGNORE_TEST:
            return
        val = ut.evaluate(2*X + 1, subs={X: 2})
        self.assertTrue(np.isclose(val, 5))

          


if __name__ == '__main__':
  unittest.main()
