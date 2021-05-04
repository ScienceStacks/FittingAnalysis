""" Utilities for sympy. ***

SUffix conventions
  *Arr - numpy array
  *Mat - sympy.Matrix N X N
  *Vec - sympy.Matrix N X 1
  *s - list

"""

import inspect
import matplotlib.pyplot as plt
import numpy as np
import sympy


def addSymbols(symbolStr, dct=None):
    """
    Adds symbols to the dictionary.

    Parameters
    ----------
    symbolStr: str
    dct: dict
        default: globals() of caller
    """
    if dct is None:
        frame = inspect.currentframe()
        dct = frame.f_back.f_globals
    symbols = symbolStr.split(" ")
    for idx, symbol in enumerate(symbols):
        dct[symbol] = sympy.Symbol(symbol)

def removeSymbols(symbolStr, dct=None):
    """
    Removes symbols from the dictionary.

    Parameters
    ----------
    symbolStr: str
    dct: dict
    """
    if dct is None:
        frame = inspect.currentframe()
        dct = frame.f_back.f_globals
    symbols = symbolStr.split(" ")
    for symbol in symbols:
        del dct[symbol]

def substitute(expression, subs={}):
    """
    Substitutes into the expression.

    Parameters
    ----------
    subs: dict
       key: sympy.symbol
       value: number
    
    Returns
    -------
    sympy.expression
    """
    expr = expression.copy()
    for key, value in subs.items():
        expr = expr.subs(key, value)
    return sympy.simplify(expr)

def evaluate(expression, isNumpy=True, **kwargs):
    """
    Evaluates the solution for the substitutions provided.

    Parameters
    ----------
    expression: sympy.Add
    isNumpy: bool
        return float or ndarray of float
    kwargs: dict
        keyword arguments for substitute
    
    Returns
    -------
    float/np.ndarray
    """
    expr = substitute(expression, **kwargs)
    val = expr.evalf()
    if isNumpy:
        if "rows" in dir(expression):
            result = np.array(val)
        else:
            result = float(val)
    else:
        result = val
    return result

def mkVector(name, numRow):
    """
    Constructs a vector of symbols.

    Parameters
    ----------
    name: str
        root name for elements of the vector
    numRow: int
    
    Returns
    -------
    sympy.Matrix numRow X 1
    """
    # Create the solution vector. The resulting vector is in the global name space.
    symbols = ["%s_%d" % (name, n) for n in range(numRow)]
    symbolStr = " ".join(symbols)
    addSymbols(symbolStr)
    return sympy.Matrix([ [s] for s in symbols])

