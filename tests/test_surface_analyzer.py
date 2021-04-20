"""
Tests for Kinetic Law
"""
from src.surface_analyzer import SurfaceAnalyzer
import src.constants as cn

import matplotlib
matplotlib.use("Tkagg")
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
MODEL = """

J1: $X0 -> x; k1*X0
J2: x -> $X1; k2*x

X0 = 1
x = 0
k1 = 1
k2 = 1
"""
MODEL2 = """

J1: $X0 -> x; k1*X0
J1a: $X0 -> x2; k1*X0
J2: x -> $X1; k2*x
J2a: x2 -> $X1; k2*x2

X0 = 1
x = 0
k1 = 1
k2 = 1
"""
PARAMETER_DCT = {"k1": 1, "k2": 1}


#############################
# Tests
#############################
class TestSurfaceAnalyzer(unittest.TestCase):

  def setUp(self):
    self.analyzer = SurfaceAnalyzer(MODEL, PARAMETER_DCT)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    def test(analyzer):
      self.assertTrue(isinstance(analyzer.baseArr, np.ndarray))
      self.assertEqual(np.ndim(analyzer.baseArr), 1)
      self.assertTrue(isinstance(analyzer.nrmseNrm, np.ndarray))
    #
    test(self.analyzer)
    test(SurfaceAnalyzer(MODEL2, PARAMETER_DCT))

  def testSimulate(self):
    if IGNORE_TEST:
      return
    arr = self.analyzer.simulate(PARAMETER_DCT)[:, 1]
    diffSq = sum((self.analyzer.baseArr - arr)**2)
    self.assertTrue(np.isclose(diffSq, 0))

  def testGetFlatValues(self):
    if IGNORE_TEST:
      return
    analyzer = SurfaceAnalyzer(MODEL2, PARAMETER_DCT)
    values = analyzer._getFlatValues(PARAMETER_DCT)
    self.assertEqual(np.ndim(values), 1)

  def testCalcNrmse(self):
    if IGNORE_TEST:
      return
    def test(model):
      analyzer = SurfaceAnalyzer(model, PARAMETER_DCT)
      values = analyzer._getFlatValues(PARAMETER_DCT)
      rssq = np.sum((values - analyzer.baseArr)**2)
      self.assertTrue(np.isclose(rssq, 0))
    #
    test(MODEL)
    test(MODEL2)

  def testMkParameterRange(self):
    if IGNORE_TEST:
      return
    def test(size):
      initialValue = 2
      changeFrc = 0.5
      arr = self.analyzer._mkParameterRange(initialValue, changeFrc, size)
      if size % 2 == 0:
        self.assertEqual(len(arr), size+1)
      else:
        self.assertEqual(len(arr), size)
      expected = initialValue*(1-changeFrc)
      self.assertTrue(np.isclose(arr[0], expected))
    #
    test(10)
    test(11)
    
  def testMkFactorialDesign(self):
    if IGNORE_TEST:
      return
    analyzer = SurfaceAnalyzer(MODEL2, PARAMETER_DCT)
    numLevel = 5
    changeFrc = 1.0
    parameterDct = {"k1": 1, "k2": 2, "k3": 3}
    designDF = analyzer._mkFactorialDesign(parameterDct, 1.0, numLevel)
    # Correct length
    self.assertEqual(len(designDF), numLevel**len(parameterDct))
    # First values are 0
    trues = [v == 0 for v in designDF.loc[0,:]]
    self.assertTrue(all(trues))
    # Middle value is the original value of the parameter
    midIdx = (len(designDF) - 1) // 2
    trues = [designDF.loc[midIdx, k] == parameterDct[k]
        for k in parameterDct.keys()]
    self.assertTrue(all(trues))

  def testRunExperiments(self):
    if IGNORE_TEST:
      return
    numLevel = 5
    maxFrc = 0.5
    self.analyzer.runExperiments(maxFrc, numLevel)
    self.assertTrue(isinstance(self.analyzer.simDF, pd.DataFrame))
    columns = list(PARAMETER_DCT.keys())
    columns.append(cn.NRMSE)
    diff = set(columns).symmetric_difference(self.analyzer.simDF.columns)
    self.assertEqual(len(diff), 0)

  def testPlotSurface(self):
    if IGNORE_TEST:
      return
    with self.assertRaises(ValueError):
      self.analyzer.plotSurface()
    #
    numLevel = 10
    for maxFrc in [1.0, 0.1, 0.01]:
      self.analyzer.runExperiments(maxFrc, numLevel)
      self.analyzer.plotSurface(isPlot=IS_PLOT)


if __name__ == '__main__':
  unittest.main()

