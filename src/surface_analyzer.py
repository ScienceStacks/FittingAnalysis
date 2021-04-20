#!/usr/bin/env python

"""
Analyzes fitting surfaces for reaction networks.
"""
import src.constants as cn

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import tellurium as te
import pandas as pd
import seaborn as sn


class SurfaceAnalyzer():

  def __init__(self, model, trueParameterDct,
      startTime=0, endTime=None, numPoint=None):
    """
    Parameters
    ----------
    model: str/ExtendedAntimony
      parameter values are the true values
    trueParameterDct: dict
      True values of parameters for the model
      key: parameter name
      value: value of parameter
    
    Returns
    -------
    """
    self.startTime = startTime
    self.endTime = endTime
    self.numPoint = numPoint
    self.trueParameterDct = trueParameterDct
    if isinstance(model, str):
      self.antimonyModel = model
      self.roadrunnerModel = te.loada(self.antimonyModel)
    else:
      self.roadrunnerModel = model
      self.antimonyModel = self.roadrunnerModel.getAntimony()
    self.baseArr = self._getFlatValues(self.trueParameterDct)
    self.nrmseNrm = np.sqrt(np.sum(self.baseArr**2))
    self.simDF = None # Results of last simulation experiments
    self.simTable = None  # 2 parameter pivot table created in plot

  def simulate(self, parameterDct):
    """
    Runs the simulation with the parameter values in the Dct
    
    Parameters
    ----------
    parameterDct
        key: parameterName
        value: value of parameter
        
    Returns
    -------
    NamedArray
    """
    self.roadrunnerModel.reset()
    for parameterName, value in parameterDct.items():
      self.roadrunnerModel[parameterName] = value
    data = self.roadrunnerModel.simulate(self.startTime, self.endTime,
        self.numPoint)
    return data

  def _getFlatValues(self, parameterDct):
    """
    Provides one dimensional array of simulation results.
    
    Parameters
    ----------
    parameterDct: dict
    model: str
        Antimony string
    """
    data = self.simulate(parameterDct)
    flattenedArr = (data[:, 1:]).flatten()
    return flattenedArr

  def calcNrmse(self, parameterDct):
    """
    Calculates the R^2 value for the difference between the aseline and parameters.
    
    Parameters
    ----------
    parameterDct: dict
        changed parameters
         
    Returns
    -------
    float
        root of the mean square error normalized by the variance of the base
    """
    newArr = self._getFlatValues(parameterDct)
    residualsArr = self.baseArr - newArr
    nrmse = np.sqrt(np.sum(residualsArr**2))/self.nrmseNrm
    return nrmse

  @staticmethod
  def _mkParameterRange(initialValue, maxFrc, numPoint):
    """
    Creates an array of values centered on the initial value
    
    Parameters
    ----------
    initialValue: float
    maxFrc: float
        maximum fractional change
    numPoint: int
    
    Result
    ------
    np.array
    """
    if numPoint % 2 == 0:
      numPoint += 1
    maxChange = initialValue*maxFrc
    return np.linspace(start=initialValue-maxChange,
                    stop=initialValue+maxChange,
                    num=numPoint)

  @staticmethod  
  def _mkFactorialDesign(parameterDct, maxFrc, numLevel):
    """
    Creates a factorial design for the parameters. Levels are determined
    by fractional changes in the initial value based on the maximum
    fractional change, either positive or negative.
    
    Parameter
    ---------
    parameterDct: dict
    maxFrc: float
        maximum change fraction in positive and negative direction
    numLevel: int
        total number of points in an axis
    
    Returns
    -------
    pd.DataFrame
       Columns:
           parameter name: contains parameter value
           rsq: rsq w.r.t. a change of 0
    """
    # Calculate number of points in one direction
    if numLevel // 2 == 0:
      numHalf = numLevel // 2
    else:
      numHalf = (numLevel - 1) // 2
    #
    incrFrc = maxFrc / numHalf
    posFrcs = np.array([n*incrFrc for n in range(1, numHalf+1)])
    negFrcs = - posFrcs
    negFrcs = np.sort(negFrcs)
    frcs = np.concatenate([negFrcs, np.array([0]), posFrcs])
    # Construct parameter values
    parameterNames = list(parameterDct.keys())
    levelDct = {}
    for parameterName in parameterNames:
      levelDct[parameterName] = (1 + frcs)*parameterDct[parameterName]
    # Construct the parameter columns
    designDF = pd.DataFrame({parameterNames[0]: levelDct[parameterNames[0]]})
    for parameterName in parameterNames[1:]:
      values = levelDct[parameterName]
      dfs = []
      for value in values:
          newDF = designDF.copy()
          newDF[parameterName] = value
          dfs.append(newDF)
      designDF = pd.concat(dfs)
    designDF.index = range(len(designDF))
    # Run the simulations
    return designDF
        
  def runExperiments(self, maxFrc, numLevel):
      """
      Creates an experimental design such that each parameter has the same
      maximum fractional change (both positive and negative) and the same
      number of levels. The normalized RMSE is calculated for each
      combination of parameter values.
      
      Parameters
      ----------
      maxFrac: float
        maximum fractional change in a parameter value
      numLevel: int
          
      Returns
      -------
      pd.DataFrame
          Columns: parameter names, cn.NRMSE
      """
      designDF = self._mkFactorialDesign(self.trueParameterDct, maxFrc, numLevel)
      #
      nrmses = []
      for idx, rowSer in designDF.iterrows():
          parameterDct = rowSer.to_dict()
          nrmse = self.calcNrmse(parameterDct)
          nrmses.append(nrmse)
      #
      self.simDF = designDF.copy()
      self.simDF[cn.NRMSE] = nrmses
        
  def plotSurface(self, ax=None, fig=None, scale=1.0,
      isColorBar=True, isPlot=True, title=""):
    """
    Plots the results of simulation experiments. Updats self.simTable
    
    Parameters
    ----------
    scale: float in [0, 1]
      Scale the range of values displayed
    isColorBar: bool
      include the color bar
    """
    if self.simDF is None:
      raise ValueError("Must first run simulation.")
    # Handle axes and figure
    if (ax is None) or (fig is None):
      fig, ax = plt.subplots(1)
    # Create a table for the contour
    columns = list(self.simDF.columns)
    rowName = columns[0]
    colName = columns[1]
    df = self.simDF.copy()
    df = df.astype(float)
    self.simTable = df.pivot_table(values=[cn.NRMSE], columns=[colName],
        index=[rowName])
    # Construct the contour
    my_cmap = sn.cubehelix_palette(start=3, rot=0.6, reverse=False, as_cmap=True)
    xv = np.array(list(set(self.simDF[rowName])))
    xv.sort()
    yv = np.array(list(set(self.simDF[colName])))
    yv.sort()
    zv = self.simTable.to_numpy()
    maxLevel = scale*self.simDF.max().max()
    levels = [0.01*n*maxLevel for n in range(100)]
    contourPlot = ax.contourf(xv, yv, zv, levels=levels, cmap=my_cmap)
    if isColorBar:
      cbar = plt.colorbar(contourPlot)
      cbar.ax.set_ylabel('Normalized RMSE')
    ax.set_xlabel(rowName)
    ax.set_ylabel(colName)
    # Mark the center
    ax.yaxis.grid(True, linewidth=0.2, linestyle='-', color='0.05')
    ax.xaxis.grid(True, linewidth=0.2, linestyle='-', color='0.05')
    ax.set_title(title)
    # Mark the center
    xIdx = (len(xv) - 1) // 2
    yIdx = (len(yv) - 1) // 2
    ax.scatter(xv[xIdx], yv[yIdx], color="red")
    if isPlot:
      plt.show()
