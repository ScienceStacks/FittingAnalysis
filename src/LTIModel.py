""" Model of linear time-invariant systems.

SUffix conventions
  *Arr - numpy array
  *Mat - sympy.Matrix N X N
  *Vec - sympy.Matrix N X 1
  *s - list

"""
import src.constants as cn
import src.sympyUtil as su

import matplotlib.pyplot as plt
import numpy as np
import sympy

X = "x"
t = sympy.Symbol("t")

class LTIModel():

    def __init__(self, aMat, initialVec, rVec=None):
        """
        Parameters
        ----------
        aMat: sympy.Matrix (N X N)
            A marix
        initialVec: sympy.Matrix (N X 1)
            Initial values
        rVec: sympy.Matrix (N X 1)
            r matrix in the differential equation
        """
        self.aMat = aMat
        self.numRow = self.aMat.rows
        self.initialVec = initialVec
        self.rVec = rVec
        # Results
        self.solutionVec = None
        self.evaluatedSolutionVec = None

    def solve(self):
        """
        Solves the LTI system symbolically.
        Updates self.solutionVec.
        
        Returns
        -------
        sympy.Matrix (N X 1)
        """
        vecs = []
        # Construct the fundamental matrix
        results = aMat.eigenvects()
        for result in results:
            eigen = result[0]
            term = sympy.exp(eigen*t)
            vecs.append(result[2][0] * term)
        self.solutionVec = sympy.Matrix(vecs)
        self.solutionVec = self.solutionVec.reshape(numRow, numRow)
        self.solutionVec = self.solutionVec.transpose()
        # Find the coefficients for the eigenvecs
        smat = self.solutionVec.subs(t, 0) # evaluate at time 0
        # Solve the linear system
        system = smat, initialVec
        mkVector(X, numRow)
        result = sympy.linsolve(system, *x) # returns a finite set, not a vector
        coefVec = sympy.Matrix(result.args[0])
        fundamentalMat = self.solutionVec*sympy.diag(*coefVec)
        # Solution
        if rVec is not None:
            dCoefMat = fundamentalMat.inv()*rVec
            coefMat = sympy.integrate(dCoefMat, t)
            self.solutionVec = fundamentalMat * dCoefMat
        else:
            self.solutionVec = fundamentalMat * sympy.ones(numRow, cols=1)
        return self.solutionVec

    def evaluate(self, subs={t: 0}):
        """
        Returns a numerical solution.
        
        Returns
        -------
        numpy.ndarray (N X 1)
        """
        if self.solutionVec is None:
            _ = self.solve()
        return su.evaluate(self.solutionVec, subs=subs)
