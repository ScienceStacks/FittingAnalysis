"""Collects information for eigenvalue and vectors.
Key properties:
    value - eigenvalue
    vectors - eigenvectors
    algebraicMultipliciaty - algebraic multiplicity
""" 

import src.sympyUtil as su

import collections
import numpy as np
import sympy

SMALL_VALUE = 1e-8

# Characterization for a single eigenvalue
# val - eigenvalue
# vecs - eigenvectors
# mul - algebraic multiplicity
EigenEntry = collections.namedtuple("EigenEntry", "val, vecs, mul")


class EigenInfo():
    # Container for all EigenEntry for a matrix

    def __init__(self, mat):
        """
        Parameters
        ----------
        mat: sympy.Matrix
        """
        eigenEntries = []  # Container for eigenEntries
        self.eigenvalDct = {self._roundToZero(k): v for k, v
               in mat.eigenvals().items()}
        # Create the raw EigenEntry
        for entry in mat.eigenvects():
            eigenvalue = self._roundToZero(entry[0])
            algebraicMultiplicity = self.eigenvalDct[eigenvalue]
            vecs = [self._vectorRoundToZero(v) for v in entry[2]]  # Eigenvectors
            algebraicMultiplicity = len(vecs)
            eigenEntries.append(EigenEntry(
                  val=eigenvalue,
                  vecs=vecs,
                  mul=algebraicMultiplicity))
        # TODO:
        # Prune: Sort by Eigenvalue. Then do pairwise _merge.
        self.eigenEntries = sorted(eigenEntries, key=lambda e: np.abs(e.val))
        import pdb; pdb.set_trace()

    @staticmethod
    def _merge(eigenEntry1, eigenEntry2):
        """
        Merges the two EigenEntry if they have the same eigenvalue.

        Parameters
        ----------
        eigneEntry1: EigenEntry
        eigneEntry2: EigenEntry
        
        Returns
        -------
        list-EigenEntry
        """
        if isZero(eigenEntry1.val - eigenEntry2.val):
            algebraicMultiplicity = eigenEntry1.mul + eigenEntry2.mul
            vecs = list(eigenEntry1.vecs)
            vecs.extend(eigenEntry2.vecs)
            eigenEntry = [EigenEntry(val=eigenEntry1.val,
                mul=algebraicMultiplicity,
                vecs = eliminateDuplicateVectors(vecs))
                ]
        else:
            eigenEntry = [eigenEntry1, eigenEntry2]
        return eigenEntry

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
        newValues = [EigenInfo._roundToZero(v) for v in vec]
        return sympy.Matrix(newValues)
    
    @staticmethod 
    def _roundToZero(v):
        if "is_symbol" in dir(v):
            if not v.is_Number:
                return v
        if np.abs(v) < SMALL_VALUE:
            return 0
        return v
