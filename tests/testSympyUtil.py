import src.sympyUtil as su
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
        su.addSymbols(VARIABLES)

    def testAddSymbols(self):
        if IGNORE_TEST:
            return
        names = ["xx", "yy"]
        su.addSymbols(" ".join(names))
        for name in names:
            self.assertTrue(name in globals().keys())
            expr = "isinstance(%s, sympy.Symbol)" % name
            self.assertTrue(eval(expr))

    def testAddSymbols2(self):
        if IGNORE_TEST:
            return
        names = ["xx", "yy"]
        su.addSymbols(" ".join(names), dct=globals())  # Ensure explicit opt works
        for name in names:
            self.assertTrue(name in globals().keys())
            expr = "isinstance(%s, sympy.Symbol)" % name
            self.assertTrue(eval(expr))

    def testRemoveSymbols(self):
        if IGNORE_TEST:
            return
        names = ["xx", "yy"]
        su.addSymbols(" ".join(names))
        su.removeSymbols(" ".join(names))
        for name in names:
            self.assertFalse(name in globals().keys())
            
    def testSubstitute(self):
        if IGNORE_TEST:
            return
        Y = su.substitute(2*X + 1, subs={X: Z})
        self.assertTrue("Z" in str(Y))

    def testEvaluate(self):
        if IGNORE_TEST:
            return
        val = su.evaluate(2*X + 1, subs={X: 2})
        self.assertEqual(val, 5)

    def testEvaluate2(self):
        if IGNORE_TEST:
            return
        expr = sympy.Matrix( [2*Z, Z**2])
        val = su.evaluate(expr, subs={Z: 2})
        self.assertEqual(expr.rows, np.shape(val)[0])
        self.assertEqual(expr.cols, np.shape(val)[1])


if __name__ == '__main__':
  unittest.main()
