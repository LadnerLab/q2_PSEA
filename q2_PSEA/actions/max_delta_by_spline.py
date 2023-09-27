# resources: https://www.datatechnotes.com/2021/11/scattered-data-spline-fitting-example.html

import numpy as np
import pandas as pd

# libraries for smooth splining
from scipy import interpolate


def max_delta_by_spline(timepoint1, timepoint2, indata: pd.DataFrame) -> list:
    # 5 is the default value from above tutorial
    print("yfits: %s" % spline(5, indata.loc[:, timepoint2].to_numpy()))

    # code which is the general idea of max delta by splining
    maxZ = np.apply_over_axes(np.max, indata[[timepoint1, timepoint2]], 1)
    return list


def spline(knots, y):
    x = range(0, len(y))
    x_new = np.linspace(0, 1, knots+2)[1:-1]
    q_knots = np.quantile(x, x_new)
    # smoothing condition (s) from smooth.spline() in original R code
    t, c, k = interpolate.splrep(x, y, t=q_knots, s=0.788458)
    yfit = interpolate.BSpline(t, c, k)(x)
    return yfit
