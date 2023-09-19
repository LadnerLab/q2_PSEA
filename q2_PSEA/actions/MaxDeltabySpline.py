# resources: https://rpy2.github.io/doc/latest/html/generated_rst/pandas.html
#            https://stackoverflow.com/questions/51167555/how-to-predict-using-rpy2
#            https://scikit-learn.org/stable/auto_examples/linear_model/plot_polynomial_interpolation.html

import numpy as np
import pandas as pd

# libraries for smooth splining
from scipy import interpolate


def max_delta_by_spline(timepoint1, timepoint2, indata: pd.DataFrame) -> list:
    # specific code for smooth spline
    # TODO: probably gets wrapped in a dedicated function called "spline" - most likely
    y = indata.loc[:, timepoint2].to_numpy()
    x = range(0, len(y))
    knot_numbers = 5
    x_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
    print("x_new: %s" % x_new)

    # giving column at timepoint1 because R code passes indata[[timepoint1]] to predict()
    # tutorial x is indata[[timepoint1]] at the moment
    q_knots = np.quantile(x, x_new)

    # smoothing condition from smooth.spline() in original R code
    t, c, k = interpolate.splrep(x, y, t=q_knots, s=0.788458)
    yfit = interpolate.BSpline(t, c, k)(indata[[timepoint1]])
    print("yfit: %s" % yfit)
    # print("indata col=timepoint1: %s" % len(indata.loc[:, timepoint1].to_numpy()))
    # print("indata col=timepoint2: %s" % len(indata.loc[:, timepoint2].to_numpy()))

    # code which is the general idea of max delta by splining
    maxZ = np.apply_over_axes(np.max, indata[[timepoint1, timepoint2]], 1)
    return list
