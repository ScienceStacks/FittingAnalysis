import src.constants as cn
from src.LTIModel import LTIModel
import src.sympyUtil as su

import numpy as np
import pandas as pd
import sympy
import unittest


IGNORE_TEST = True
IS_PLOT = True
SYMBOLS = "X_0 X_1 x y z k0 k1 k2 t"
su.addSymbols(SYMBOLS, dct=globals())
SUBS = {k0: 1, k1: 2, k2: 3, t: 100}


#############################
# Tests
#############################
class TestLTIModel(unittest.TestCase):

    def setUp(self):
        su.addSymbols(SYMBOLS, dct=globals())
        aMat = sympy.Matrix([ [   0,    0,          0,  0],
                              [   1,  -(k0 + k2),   0,  0], 
                              [   0,          k0, -k1,  0], 
                              [   0,          k2,  k1, -1],
                            ])
        initialVec = sympy.Matrix([1, 0, 0, 0])
        rVec = sympy.Matrix([k0, 0, -k0])
        rVec = None
        self.model = LTIModel(aMat, initialVec, rVec=rVec)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(self.model.aMat.rows, 4)

    def testSolve1(self):
        # TESTING
        # Homogeneous equation with initial values
        subs = dict(SUBS)
        resultVec = self.model.solve(subs=subs)
        resultVec = resultVec.subs(subs)
        resultVec = sympy.simplify(resultVec)
        vals = [1, 0.25, 0.125, 1.0]
        for pos in range(len(resultVec)):
            num1 = float(resultVec[pos])
            num2 = float(vals[pos])
            self.assertTrue(np.isclose(num1, num2))

    def testEvaluate1(self):
        if IGNORE_TEST:
            return
        solutionVec = self.model.solve()
        resultVec = su.evaluate(solutionVec, subs=SUBS)
        vals = [1, 0.25, 0.125, 1.0]
        for pos in range(len(resultVec)):
            num1 = float(resultVec[pos][0])
            num2 = float(vals[pos])
            self.assertTrue(np.isclose(num1, num2))

    def testEvaluate2(self):
        # TESTING
        # Geometric multiplicity < algebraic multiplicity
        subs = dict(SUBS)
        subs[k1] = 1
        solutionVec = self.model.solve(subs=subs)
        import pdb; pdb.set_trace()
        resultVec = su.evaluate(solutionVec, subs=subs)
        vals = [1, 0.25, 0.125, 1.0]
        for pos in range(len(resultVec)):
            num1 = float(resultVec[pos][0])
            num2 = float(vals[pos])
            self.assertTrue(np.isclose(num1, num2))

    # TODO: Need better test matrix
    def testMkEigenVectors(self):
        return
        if IGNORE_TEST:
            return
        mat = 3*sympy.eye(2)
        infos = su.getEigenInfo(mat)
        eigenInfo = infos[0]
        eigenInfo.vecs.reverse()
        eigenInfo.vecs.pop()
        result = self.model._mkEigenVectors(mat, eigenInfo)

    def testPlot(self):
        if IGNORE_TEST:
            return
        subs = dict(SUBS)
        del subs[t]
        self.model.plot(0, 10, 100, subs, isPlot=IS_PLOT, ylabel="values")



if __name__ == '__main__':
  unittest.main()
