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
import time

X = "x"
t = sympy.Symbol("t")
IS_TIMER = False


class Timer():

    def __init__(self, name):
        self.name = name
        if IS_TIMER:
            self.startTime = time.time()

    def start(self):
        if IS_TIMER:
            self.startTime = time.time()

    def print(self, name=None):
        if not IS_TIMER:
            return
        if name is None:
            name = self.name
        elapsed = time.time() - self.startTime
        print("***%s: %2.4f" % (name, elapsed))
        self.start()


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

    # TODO: TEST
    @staticmethod
    def _mkEigenVectors(aMat, eigenInfo):
        """
        Creates an orthogonal set of eigenvectors to handle
        situations where the algebraic multiplicity is greater than
        geometric multiplicity.
  
        From "3.7: Multiple Eigenvalues", Libre Tests.

        Parameters
        ----------
        aMat: sympy.Matrix (N X N)
        eigenInfo: EigenInfo
        
        Returns
        -------
        list-(sympy.Matrix N X 1)
        """
        if eigenInfo.mul == len(eigenInfo.vecs):
            return eigenInfo.vecs
        if len(eigenInfo.vecs) != 1:
            raise RuntimeError("Cannot handle this case.")
        # Construct a set of orthogonal eigenvectors
        vecs = list(eigenInfo.vecs)
        curVec = vecs[0]
        numRow = curVec.rows
        dummyVec = su.mkVector("x", numRow)
        dummySymbols = [v for v in dummyVec]
        newVecs = []
        for _ in range(eigenInfo.mul - len(vecs)):
            mat = aMat - eigenInfo.val*sympy.eye(numRow)
            system = mat, curVec
            result = sympy.linsolve(system, *dummyVec)
            lst = su.flatten(result)
            # Handle case of multiple solutions
            subs = {s: 1 for s in lst if s in dummySymbols}
            newVec = sympy.Matrix(lst)
            newVec = su.substitute(newVec, subs=subs)
            newVecs.append(newVec)
            curVec = newVec
        # Construct linear independent vector
        eigenVecs = list(newVecs[0])
        numVec = len(newVecs)
        for idx in range(numVec):
            n = numVec - idx - 1
            newEigenVec = 1.0/np.math.factorial(n) * (t**n) * newVecs[idx]
            eigenVecs.append(newEigenVec)
        return eigenVecs

    def solve(self, subs={}):
        """
        Solves the LTI system symbolically.
        Updates self.solutionVec.
        
        Returns
        -------
        sympy.Matrix (N X 1)
        """
        # TODO: Handle algebraich multiplicy > 2
        # TODO: Handle imaginary eigenvalues
        # FIXME: fundamental matrix does not include the constants and so
        #        I'm incorrectly calculating Yp
        def vectorRoundToZero(vec):
            if vec.cols > 1:
                RuntimeError("Can only handle vectors.")
            newValues = [roundToZero(v) for v in vec]
            return sympy.Matrix(newValues)
        #
        def roundToZero(v):
            if np.abs(v) < SMALL_VALUE:
                return 0
            return v
        #
        timer = Timer("solve")
        vecs = []
        # Do evaluations
        aMat = su.evaluate(self.aMat, isNumpy=False, subs=subs)
        initialVec = su.evaluate(self.initialVec, isNumpy=False, subs=subs)
        # Find the complete set of eigenvectors, handling cases
        # in which the algebraic multiplicity > geometric multiplicity
        # TODO: combine together EigenInfo with same eigenValue
        eigenInfos = su.getEigenInfo(aMat)
        for eigenInfo in eigenInfos:
            import pdb; pdb.set_trace()
            eigenVecs = self._mkEigenVectors(aMat, eigenInfo)
            term = sympy.exp(eigenInfo.val * t)
            vecs.extend([v * term for v in eigenVecs])
        timer.print(name="solve_1")
        # Construct the fundamental matrix
        fundamentalMat= sympy.Matrix(vecs)
        fundamentalMat = fundamentalMat.reshape(self.numRow, self.numRow)
        fundamentalMat = sympy.simplify(fundamentalMat.transpose())
        timer.print(name="solve_2")
        # Find the coefficients for the homogeneous system
        # TODO: Instead of finding the inverse, solve the linear system
        t0Mat = sympy.simplify(fundamentalMat.subs(t, 0)) # evaluate at time 0
        import pdb; pdb.set_trace()
        t0InvMat = t0Mat.inv()
        coefHomogeneousVec = sympy.simplify(t0InvMat*initialVec)
        timer.print(name="solve_3")
        # Particular solution 
        if self.rVec is not None:
            rVec = su.evaluate(self.rVec, isNumpy=False, subs=subs)
            invfundMat = fundamentalMat.inv()
            timer.print(name="solve_4a")
            dCoefMat = invfundMat*rVec
            timer.print(name="solve_4b")
            coefMat = sympy.integrate(dCoefMat, t)
            timer.print(name="solve_4c")
            particularVec = fundamentalMat * dCoefMat
            timer.print(name="solve_4")
        else:
            particularVec = sympy.zeros(self.numRow, 1)
            timer.print(name="solve_5")
        # Full solution
        self.solutionVec =  particularVec + fundamentalMat*coefHomogeneousVec
        timer.print(name="solve_6")
        self.solutionVec = sympy.simplify(self.solutionVec)
        timer.print(name="solve_7")
        return self.solutionVec

    def evaluate(self, subs={}, **kwargs):
        """
        Returns a numerical solution.

        Parameters
        ----------
        kwargs: dict
            dictionary of symbol value assignments
        
        Returns
        -------
        numpy.ndarray (N X 1)
        """
        if self.solutionVec is None:
            _ = self.solve(subs=subs)
        return su.evaluate(self.solutionVec, subs=subs, **kwargs)

    def plot(self, startTime, endTime, numPoint, subs, isPlot=True,
          ylabel="", title=""):
        """
        Does time series plot. Should not substitute 't'.

        Parameters
        ----------
        startTime: float
        endTime: float
        numPoint: int
        subs: dict
            arguments for evaluate
        isPlot: bool
            display the plot
        ylabel: str
        title: str
        """
        # Create the plot data
        vec = self.evaluate(subs=subs, isNumpy=False)
        delta = (endTime - startTime)/numPoint
        tVals = [startTime + n*delta for n in range(numPoint)]
        arrs = [vec.subs(t, tval) for tval in tVals]
        # Plot the results
        labels = ["x_%d" % n for n in range(len(vec))]
        fig, ax = plt.subplots(1)
        for idx in range(len(vec)):
            yvals = [arr[idx] for arr in arrs]
            ax.plot(tVals, yvals)
        ax.legend(labels)
        ax.set_xlabel("time")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if isPlot:
            plt.show()
        

