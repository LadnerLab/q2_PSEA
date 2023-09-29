import gseapy as gp
import numpy as np
import pandas as pd

from math import pow, log
from scipy import interpolate


def make_psea_table():
    base = 2
    offset = 3
    power = pow(base, offset)

    # table will also be used for GSEA processing
    data = read_table("../../example/IM0031_PV2T_25nt_raw_2mm_i1mm_Z-HDI75.tsv")
    # data1 = data.apply(lambda row: power + row, axis=0)
    # data1 = data1.apply(
    #     lambda row: row.apply(lambda val: 1 if val < 1 else val),
    #     axis=0
    # )
    # # might need to catch an exception for no columns
    # data2 = data1.drop(columns=[])

    # data = data2.apply(
    #     lambda row: row.apply(lambda val: log(val, base) - offset)
    # )

    # maxDelta = max_delta_by_spline(timepoint1, timepoint2, data)
    # maxZ = maxDelta[0]
    # deltaZ = maxDelta[1]
    # # TODO: add an option for the user to dictate the threshold
    # #table = psea(maxZ, deltaZ, 0.75)
    # table = psea(
    #     maxZ, deltaZ, 1.00,
    #     "../../example/input.tsv", "../../example/PV2species.csv"
    # ) # for testing purposes

    gs_res = gp.gsea(
        data=data,
        gene_sets="../../example/input.gmt",
        cls="../../example/Pro_PV2T.cls",
        permutation_type="gene_set",
        threads=4
    )

    print(f"GSEA head: {gs_res.res2d.head()}")


def max_delta_by_spline(timepoint1, timepoint2, indata: pd.DataFrame) -> tuple:
    """

    Parameters
    ----------
    timepoint1 : str
        x values for which y values will be predicted 

    timepoint2 : str
        values onto which a smooth cubic spline will be fit

    indata : pd.DataFrame
        matrix of Z scores for sequence

    Returns
    -------
    tuple
        Contains maximum Z score, the difference (delta) in actual from
        predicted Z scores, the spline values for timepoint1 (x) and
        timepoint2 (y)
    """
    maxZ = np.apply_over_axes(np.max, indata.loc[:, [timepoint1, timepoint2]], 1)

    # perform smoothing spline prediction
    y = indata.loc[:, timepoint2].to_numpy()
    x = indata.loc[:, timepoint1].to_numpy()
    # tentative magic number 5 (knots) came from tutorial linked above
    smooth_spline = spline(5, y)
    deltaZ = y - smooth_spline(x)

    # convert maxZ and deltaZ to pandas Series - allows for association of
    # values with peptides
    maxZ = pd.Series(data=[num for num in maxZ], index=indata.index)
    deltaZ = pd.Series(data=deltaZ, index=indata.index)
    
    return (maxZ, deltaZ, smooth_spline(x), smooth_spline(y))


def psea(
        data: pd.DataFrame,
        maxZ: pd.Series,
        deltaZ: pd.Series,
        threshold: float,
        input: str, # rename this file?
        species: str
):
    """

    Parameters
    ----------
    maxZ : pd.Series

    deltaZ : pd.Series

    threshold : float

    input : pd.DataFrame

    species : pd.DataFrame

    Returns
    -------
    """
    # # grab indexes where condition is true
    # maxZ_above_thresh = np.where(maxZ > threshold)
    # deltaZ_not_zero = np.where(deltaZ != 0)
    # # create gene list
    # gene_list = deltaZ.iloc[np.intersect1d(maxZ_above_thresh, deltaZ_not_zero)].sort_values(ascending=False)

    # # read and rename columns
    # # TODO: make sure we don't want a header
    # term_gene = pd.read_csv(input, sep="\t", header=None)


def read_table(file_path):
    data = []
    columns = []
    indexes = []
    with open("../../example/IM0031_PV2T_25nt_raw_2mm_i1mm_Z-HDI75.tsv", "r") \
            as data_fh:
        # read matrix into memory
        lines = data_fh.readlines()

        # grab column names from header
        columns = [name for name in lines[0].replace("\n", "").split("\t")]
        columns.pop(0) # remove "Sequence name"
        lines.pop(0) # remove header line

        for line in lines:
            split_line = line.replace("\n", "").split("\t")
            indexes.append(split_line[0]) # add peptide to index list
            split_line.pop(0) # remove peptide name
            data.append(split_line) # add data line

    return pd.DataFrame(
        data=data,
        index=indexes,
        columns=columns,
        dtype=float
    )


def spline(knots, y):
    """<description>

    Parameters
    ----------
    knots : int

    y : float array

    Returns
    -------
    BSpline
    """
    x = range(0, len(y))
    x_new = np.linspace(0, 1, knots+2)[1:-1]
    q_knots = np.quantile(x, x_new)
    # smoothing condition `s` from smooth.spline() in original R code
    t, c, k = interpolate.splrep(x, y, t=q_knots, s=0.788458)
    return interpolate.BSpline(t, c, k)


make_psea_table()  # remove after plugin code has been full setup
