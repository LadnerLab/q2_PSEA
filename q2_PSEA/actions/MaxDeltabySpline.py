# resources: https://rpy2.github.io/doc/latest/html/generated_rst/pandas.html
#            https://stackoverflow.com/questions/51167555/how-to-predict-using-rpy2

import numpy as np
import pandas as pd

def max_delta_by_spline(timepoint1, timepoint2, indata: pd.DataFrame) -> list:
    # TODO: make sure maxZ is what is expected/needed
    maxZ = np.apply_over_axes(np.max, indata[[timepoint1, timepoint2]], 1)
    # deltaZ = indata[[timepoint2]].apply(lambda val: )
    return list

def smooth_spline(x, y):
    pass

def predict():
    pass
