# resources: https://www.datatechnotes.com/2021/11/scattered-data-spline-fitting-example.html

import numpy as np
import pandas as pd

# libraries for smooth splining
from scipy import interpolate


def max_delta_by_spline(timepoints, indata: pd.DataFrame) -> tuple:
    """<description>

    Parameters
    ----------
    timepoints : str

    indata : pd.DataFrame
        matrix of Z scores for sequence

    Returns
    -------
    tuple
        Contains maximum Z score, the difference (delta) in actual from
        predicted Z scores, the spline values for timepoint1 (x) and
        timepoint2 (y)
    """
    maxZ = np.apply_over_axes(np.max, indata.loc[:, [timepoints[0], timepoints[1]]], 1)

    # perform smoothing spline prediction
    y = indata.loc[:, timepoints[0]].to_numpy()
    x = indata.loc[:, timepoints[1]].to_numpy()
    # tentative magic number 5 (knots) came from tutorial linked above
    smooth_spline = spline(5, y)
    deltaZ = y - smooth_spline(x)

    # convert maxZ and deltaZ to pandas Series - allows for association of
    # values with peptides
    maxZ = pd.Series(data=[num for num in maxZ], index=indata.index)
    deltaZ = pd.Series(data=deltaZ, index=indata.index)
    
    return (maxZ, deltaZ, smooth_spline(x), smooth_spline(y))


def spline(knots, y):
    """<description>

    Parameters
    ----------
    knots : int

    y : float array

    Returns
    -------
    BSpline
        Used to make predictions
    """
    x = range(0, len(y))
    x_new = np.linspace(0, 1, knots+2)[1:-1]
    q_knots = np.quantile(x, x_new)
    # smoothing condition `s` from smooth.spline() in original R code
    t, c, k = interpolate.splrep(x, y, t=q_knots, s=0.788458)
    return interpolate.BSpline(t, c, k)
