# resources: https://www.datatechnotes.com/2021/11/scattered-data-spline-fitting-example.html

import numpy as np
import pandas as pd

# libraries for smooth splining
from scipy import interpolate


def max_delta_by_spline(timepoint1, timepoint2, indata: pd.DataFrame) -> list:
    # specific code for smooth spline
    # TODO: probably gets wrapped in a dedicated function called "spline" - most likely
    knot_numbers = 5 # default value from tutorial
    y = indata.loc[:, timepoint2].to_numpy()
    x = range(0, len(y))
    x_new = np.linspace(0, 1, knot_numbers+2)[1:-1]

    # giving column at timepoint1 because R code passes indata[[timepoint1]] to predict()
    # tutorial x is indata[[timepoint1]] at the moment
    q_knots = np.quantile(x, x_new)

    # smoothing condition from smooth.spline() in original R code
    t, c, k = interpolate.splrep(x, y, t=q_knots, s=0.788458)
    yfit = interpolate.BSpline(t, c, k)(indata[[timepoint1]])
    print("yfit: %s" % yfit) # TODO: remove when finished troubleshooting for expected test values

    # code which is the general idea of max delta by splining
    maxZ = np.apply_over_axes(np.max, indata[[timepoint1, timepoint2]], 1)
    return list
