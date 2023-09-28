import numpy as np
import pandas as pd

from math import pow, log
from scipy import interpolate


def make_psea_table(ctx, data, timepoints):
    """

    Parameters
    ----------

    Returns
    -------
    """
    # TODO: optimize dataframe creation by using pd.read_csv()
    scores = []
    columns = []
    indexes = []
    with open(data, "r") as data_fh:
        # read matrix into memory
        lines = data_fh.readlines()

        # grab column names from header
        columns = [name for name in lines[0].replace("\n", "").split("\t")]
        columns.pop(0)  # remove "Sequence name"
        lines.pop(0)  # remove header line

        for line in lines:
            split_line = line.replace("\n", "").split("\t")
            indexes.append(split_line[0])  # add peptide to index list
            split_line.pop(0)  # remove peptide name
            scores.append(split_line)  # add data line

    base = 2
    offset = 3
    power = pow(base, offset)

    datain = pd.DataFrame(
        data=scores, index=indexes, columns=columns,
        dtype=float
    )  # TODO: check how much precision is needed (float, double (32, 64))?
    data1 = datain.apply(lambda row: power + row, axis=0)
    data1 = data1.apply(
        lambda row: row.apply(lambda val: 1 if val < 1 else val),
        axis=0
    )
    # might need to catch an exception for no columns
    data2 = data1.drop(columns=[])

    pass_data = data2.apply(
        lambda row: row.apply(lambda val: log(val, base) - offset)
    )

    # TODO: user should be able to specify
    timepoint1 = "070060_D360.Pro_PV2T"
    timepoint2 = "070060_D540.Pro_PV2T"

    maxDelta = max_delta_by_spline([timepoint1, timepoint2], pass_data)
    maxZ = maxDelta[0]
    deltaZ = maxDelta[1]
    print(f"MaxZ: {maxZ}")
    print(f"DeltaZ: {deltaZ}")
    # TODO: add an option for the user to dictate the threshold
    # table = psea(maxZ, deltaZ, threshold=0.75)
    # table = psea(maxZ, deltaZ, 1.00, "../../example/input.tsv")


def max_delta_by_spline(timepoints, data: pd.DataFrame) -> tuple:
    """<description>

    Parameters
    ----------
    timepoints : str

    data : pd.DataFrame
        matrix of Z scores for sequence

    Returns
    -------
    tuple
        Contains maximum Z score, the difference (delta) in actual from
        predicted Z scores, the spline values for timepoint1 (x) and
        timepoint2 (y)
    """
    maxZ = np.apply_over_axes(
        np.max,
        data.loc[:, [timepoints[0], timepoints[1]]],
        1
    )

    # perform smoothing spline prediction
    y = data.loc[:, timepoints[0]].to_numpy()
    x = data.loc[:, timepoints[1]].to_numpy()
    # tentative magic number 5 (knots) came from tutorial linked above
    smooth_spline = spline(5, y)
    deltaZ = y - smooth_spline(x)

    # convert maxZ and deltaZ to pandas Series - allows for association of
    # values with peptides
    maxZ = pd.Series(data=[num for num in maxZ], index=data.index)
    deltaZ = pd.Series(data=deltaZ, index=data.index)

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
