import gseapy as gp
import numpy as np
import pandas as pd

from math import pow, log
from scipy import interpolate


def make_psea_table(
    ctx,
    scores_file,
    timepoints_file,
    pairs_file,
    gene_sets_file,
    cls,
    min_size,
    max_size,
    threads
):
    # TODO: tentative code for reading and formatting time points and pairs
    #     1) ensure pairs file MUST be specified
    #     2) ensure timepoints file MUST be specified
    # collect sample pairs
    with open(pairs_file, "r") as fh:
        lines = [line.replace("\n", "") for line in fh.readlines()]
        pairs = [tuple(lines[i:i+2]) for i in range(0, len(lines), 2)]
    # collect timepoints
    # TODO: make sure timepoints file MUST be specified
    with open(timepoints_file, "r") as fh:
        timepoints = fh.readlines()[0].strip().split()

    print(f"Pairs: {pairs}")
    print(f"Timepoints: {timepoints}")

    processed_scores = load_scores(scores_file, pairs)
    processed_scores.to_csv("processed_scores.tsv", sep="\t")

    for pair in pairs:
        maxDelta = max_delta_by_spline(processed_scores, pair)
        maxZ = maxDelta[0]
        deltaZ = maxDelta[1]
        # probably need to extract other returned information
        print(f"MaxZ: {maxZ}\n")
        print(f"DeltaZ: {deltaZ}")

    # TODO: make sure this happens for each pair
    # table = psea(
    #     scores=processed_scores,
    #     maxZ=maxZ,
    #     deltaZ=deltaZ,
    #     threshold=1.00,
    #     gene_sets_file=gene_sets_file,
    #     cls=cls,
    #     min_size=min_size,
    #     max_size=max_size,
    #     threads=threads
    # )
    # table.to_csv("gsea_table.tsv", sep="\t")

    # # visualization
    # repScatters_tsv = ctx.get_action("ps-plot", "repScatters_tsv")

    # # for testing purposes
    # return repScatters_tsv(
    #     source=,
    #     user_space_pairs=,  # might be timepoints?
    #     pn_filepath=None,
    #     plot_log=False,
    #     zscore_filepath=scores_file,
    #     col_sum_filepath=None,
    #     facet_charts=False,
    #     xy_threshold=None, # definitely need for highlighting outliers
    # )


def max_delta_by_spline(data, pair) -> tuple:
    """Finds the maximum value between two samples, and calculates the
    difference in Z score for each peptide

    Parameters
    ----------
    data : pd.DataFrame
        Matrix of Z scores for sequence

    pair : Tuple
        Tuple pair of samples for which to run max spline

    Returns
    -------
    tuple
        Contains maximum Z score, the difference (delta) in actual from
        predicted Z scores, the spline values for x and y
    """
    maxZ = np.apply_over_axes(np.max, data.loc[:, pair], 1)

    # perform smoothing spline prediction
    y = data.loc[:, pair[0]].to_numpy()
    x = data.loc[:, pair[1]].to_numpy()
    # tentative magic number 5 (knots) came from tutorial linked above
    smooth_spline = spline(5, y)
    deltaZ = y - smooth_spline(x)

    # convert maxZ and deltaZ to pandas Series - allows for association of
    # values with peptides
    maxZ = pd.Series(data=[num for num in maxZ], index=data.index)
    deltaZ = pd.Series(data=deltaZ, index=data.index)
    
    return (maxZ, deltaZ, smooth_spline(x), smooth_spline(y))


def psea(
    # TODO: figure out if this is actually needed for process
    scores: pd.DataFrame,
    maxZ: pd.Series,
    deltaZ: pd.Series,
    threshold: float,
    gene_sets_file: str,
    cls: str,
    min_size: int = 15,
    max_size: int = 500,
    threads: int = 1
):
    """

    Parameters
    ----------
    maxZ : pd.Series

    deltaZ : pd.Series

    threshold : float

    gene_sets : str

    Returns
    -------
    """
    # TODO: figure out how to integrate into gsea function
    # grab indexes where condition is true
    maxZ_above_thresh = np.where(maxZ > threshold)
    deltaZ_not_zero = np.where(deltaZ != 0)
    # create gene list
    gene_list = deltaZ.iloc[
        np.intersect1d(maxZ_above_thresh, deltaZ_not_zero)
    ].sort_values(ascending=False)
    gene_list.to_frame()
    gene_list.to_csv("gene_list.tsv", sep="\t")

    # print(f"MaxZ > {threshold}: {maxZ_above_thresh}")
    # print(f"DeltaZ not 0: {deltaZ_not_zero}")
    # print(f"Gene list:\n{gene_list}\n")

    # TODO: make sure to add other parameters
    return gp.gsea(
        data=gene_list,
        gene_sets=gene_sets_file,
        cls=cls,
        permutation_type="gene_set",  # TODO: force gene set?
        min_size=min_size,
        max_size=max_size,
        threads=threads
    )


def load_scores(file_path, pairs):
    """Reads Z score matrix

    Returns a Pandas DataFrame of processed Z scores
    """
    base = 2
    offset = 3
    power = pow(base, offset)

    data = []
    columns = []
    indexes = []
    with open("../../example/IM0031_PV2T_25nt_raw_2mm_i1mm_Z-HDI75.tsv", "r") \
            as data_fh:
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
            data.append(split_line)  # add data line

    scores = pd.DataFrame(
        data=data,
        index=indexes,
        columns=columns,
        dtype=float
    )

    data1 = scores.apply(lambda row: power + row, axis=0)
    data1 = data1.apply(
        lambda row: row.apply(lambda val: 1 if val < 1 else val),
        axis=0
    )
    # might need to catch an exception for no columns

    # extrapolate list of column names from pairs list
    pairs_list = []
    for pair in pairs:
        for sample in pair:
            pairs_list.append(sample)
    pairs_list = list(np.unique(pairs_list))
    # grab columns from data1
    data2 = data1.loc[:, pairs_list]

    return data2.apply(
        lambda row: row.apply(lambda val: log(val, base) - offset)
    )


def spline(knots, y):
    """Creates spline object from which a prediction can be made based on given
    y-values

    Parameters
    ----------
    knots : int

    y : float array

    Returns
    -------
    BSpline
        Spline object from which predictions can be made given a set of
        x-values
    """
    x = range(0, len(y))
    x_new = np.linspace(0, 1, knots+2)[1:-1]
    q_knots = np.quantile(x, x_new)
    # smoothing condition `s` from smooth.spline() in original R code
    t, c, k = interpolate.splrep(x, y, t=q_knots, s=0.788458)
    return interpolate.BSpline(t, c, k)
