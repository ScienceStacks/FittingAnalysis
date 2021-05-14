"""Collects information for eigenvalue and vectors.
Key properties:
    value - eigenvalue
    vectors - eigenvectors
    algebraicMultipliciaty - algebraic multiplicity
""" 

import src.constants as cn
import src.sympyUtil as su

import collections
import numpy as np
import sympy

SMALL_VALUE = 1e-8
t = sympy.Symbol(cn.SYM_T)


class EigenInfo():
    # Information about one eigenvalue and its eigenvectors

    def __init__(self, val, vecs, mul):
        """
        Parameters
        ----------
        val: float
            Eigenvalue
        vecs: list-sympy.Matrix
            eigenvectors
        mul: int
            algebraic multiplicity
        """
        self.val = val
        self.vecs = vecs
        self.mul = mul

    # TODO: Not sure that I'm using the correct vectors for those added
    def addVectors(self, numvec):
        """
        Adds the specified number of eigenvectors using the first eigenvector.
        Updates self.vecs

        Parameters
        ----------
        num: int
            number of eigenvectors to add
        """
        newVecs = []
        for num in range(1, numvec+1):
            term = self.vecs[0]
            for termIdx in range(1, num+1):
                term = term + t**termIdx / sympy.factorial(termIdx) * self.vecs[0]
            newVecs.append(term)
        self.vecs.extend(newVecs)
        
    


class EigenCollection():
    # Container for all EigenInfo for a matrix

    def __init__(self, mat):
        """
        Parameters
        ----------
        mat: sympy.Matrix
        """
        def simplify(v):
            newV = self._roundToZero(v)
            return su.expressionToNumber(newV.evalf())
        #
        eigenInfos = []  # Container for eigenInfos
        self.eigenvalDct = {simplify(k): v for k, v in mat.eigenvals().items()}
        # Create the raw EigenInfo
        for eigenInfo in mat.eigenvects():
            eigenvalue = simplify(eigenInfo[0])
            algebraicMultiplicity = self.eigenvalDct[eigenvalue]
            vecs = [self._vectorRoundToZero(v) for v in eigenInfo[2]]
            vecs = [v.evalf() for v in vecs]
            eigenInfos.append(EigenInfo(
                  val=eigenvalue,
                  vecs=vecs,
                  mul=algebraicMultiplicity))
        # TODO:
        # Prune: Sort by Eigenvalue. Then do pairwise _merge.
        self.eigenInfos = sorted(eigenInfos, key=lambda e: np.abs(e.val))
        for idx, eigenInfo in enumerate(self.eigenInfos[:-1]):
            if np.isclose(eigenInfo.val, self.eigenInfos[idx+1].val):
                raise RuntimeError("Duplicate EigenInfo")

    # TODO: Is this needed?
    @staticmethod
    def _merge(eigenInfo1, eigenInfo2):
        """
        Merges the two EigenInfo if they have the same eigenvalue.

        Parameters
        ----------
        eigneEntry1: EigenInfo
        eigneEntry2: EigenInfo
        
        Returns
        -------
        list-EigenInfo
        """
        if isZero(eigenInfo1.val - eigenInfo2.val):
            algebraicMultiplicity = eigenInfo1.mul + eigenInfo2.mul
            vecs = list(eigenInfo1.vecs)
            vecs.extend(eigenInfo2.vecs)
            eigenInfo = [EigenInfo(val=eigenInfo1.val,
                mul=algebraicMultiplicity,
                vecs = eliminateDuplicateVectors(vecs))
                ]
        else:
            eigenInfo = [eigenInfo1, eigenInfo2]
        return eigenInfo

    # TODO: Is this needed?
    @staticmethod
    def eliminateDuplicateVectors(vecs):
        """
        Eliminates vectors that have the same values.

        Parameters
        ----------
        vecs: list-sympy.Matrix
        
        Returns
        -------
        vecs: list-sympy.Matrix
        """
        results = []
        for idx in range(len(vecs) - 1):
            curVec = vecs[idx]
            results.append(curVec)  # Assume it differs from other vectors
            for vec in vecs[idx+1:]:
                if curVec.rows != vec.rows:
                    continue
                if su.isVecZero(curVec - vec):
                    _ = results.pop()  # Remove curVec
                    break
   
    @staticmethod 
    def _vectorRoundToZero(vec):
        if vec.cols > 1:
            RuntimeError("Can only handle vectors.")
        newValues = [EigenCollection._roundToZero(v) for v in vec]
        return sympy.Matrix(newValues)
    
    @staticmethod 
    def _roundToZero(v):
        if "is_symbol" in dir(v):
            if not v.is_Number:
                return v
        if np.abs(v) < SMALL_VALUE:
            return 0
        return v
