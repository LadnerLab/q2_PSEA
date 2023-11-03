import gseapy as gp
import numpy as np
import pandas as pd
import qiime2

from math import pow, log
from scipy import interpolate
from q2_pepsirf.format_types import (
    PepsirfContingencyTSVFormat, PepsirfInfoSNPNFormat
)


def generate_metadata(replicates):
    base_reps = []

    replicates.sort()

    for replicate in replicates:
        base_seq_name = replicate.split("_")[2]
        base_reps.append(base_seq_name)

    meta_series = pd.Series(data=base_reps, index=replicates)
    meta_series.index.name = "sample-id"
    meta_series.name = "source"
    print(f"\n\nMeta-series:\n{meta_series}\n\n")

    return qiime2.metadata.CategoricalMetadataColumn(meta_series)


def make_psea_table(
        ctx,
        scores_file,
        timepoints_file,
        pairs_file,
        gene_sets_file,
        threshold,
        min_size,
        max_size,
        threads=4,
        pepsirf_binary="pepsirf"
):
    norm = ctx.get_action("pepsirf", "norm")
    zenrich = ctx.get_action("ps-plot", "zenrich")

    # collect zscores -> is there a better place to do this?
    scores = pd.read_csv(scores_file, sep="\t", index_col=0)
    # collect pairs
    with open(pairs_file, "r") as fh:
        pairs = [
            tuple(line.replace("\n", "").split("\t"))
            for line in fh.readlines()
        ]
    # collect timepoints
    with open(timepoints_file, "r") as fh:
        timepoints = fh.readlines()[0].strip().split()

    # process scores
    processed_scores = process_scores(scores, pairs)
    # save to disk for zenrich plot creation
    processed_scores.to_csv("py_processed_scores.tsv", sep="\t", index=True)

    # # run psea for pairs
    # table_num = 1
    # for pair in pairs:
    #     # run PSEA operation for current pair
    #     spline_tup = max_delta_by_spline(processed_scores, pair)
    #     maxZ = spline_tup[0]
    #     deltaZ = spline_tup[1]
    #     # spline_x = spline_tup[2]
    #     # spline_y = spline_tup[3]

    #     table = psea(
    #         maxZ=maxZ,
    #         deltaZ=deltaZ,
    #         gene_sets_file=gene_sets_file,
    #         threshold=threshold,
    #         min_size=min_size,
    #         max_size=max_size,
    #         threads=threads,
    #         permutation_num=10000,
    #         outdir="table_outdir"
    #     )
    #     print(f"PSEA result table:\n{table.res2d}")
    #     table_num += 1

    # import score data as artifacts
    # scores_artifact = ctx.make_artifact(
    #     type="FeatureTable[Zscore]",
    #     view=scores_file,
    #     view_type=PepsirfContingencyTSVFormat
    # )
    # processed_artifact = ctx.make_artifact(
    #     type="FeatureTable[RawCounts]",
    #     view="processed_scores.tsv",
    #     view_type=PepsirfContingencyTSVFormat
    # )
    # TODO: ensure formatting does not cause some unforeseen problem
    # highlight_probes = ctx.make_artifact(
    #     type="InfoSNPN",
    #     view="lead_genes.tsv",
    #     view_type=PepsirfInfoSNPNFormat
    # )

    # generate zenrich plot
    # col_sum, = norm(
    #     peptide_scores=processed_artifact,
    #     normalize_approach="col_sum",  # TODO: check if user should decide
    #     pepsirf_binary=pepsirf_binary
    # )
    # col_sum.view(PepsirfContingencyTSVFormat).save("norm_processed_scores.tsv")
    # zenrich_plot = zenrich(
    #     data=col_sum,
    #     zscores=scores_artifact,
    #     source=generate_metadata([rep for pair in pairs for rep in pair]),
    #     highlight_probes=highlight_probes,
    #     negative_controls=None
    # )

    # TODO: remove temp "processed_scores.tsv" file
    # return zenrich_plot
    return qiime2.sdk.Visualization


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
    maxZ: pd.Series,
    deltaZ: pd.Series,
    gene_sets_file: str,
    threshold: float,
    permutation_num: int = 0,
    min_size: int = 15,
    max_size: int = 500,
    threads: int = 1,
    outdir: str = None
):
    """Finds the intersection between peptides in `maxZ` with a value greater
    than `threshold` and those in `deltaZ` which have a value not equal to 0

    Parameters
    ----------
    maxZ : pd.Series

    deltaZ : pd.Series

    gene_sets_file : str

    threshold : float

    min_size : int

    max_size : int

    threads : int

    Returns
    -------
    SingleSampleGSEA
        Containes the resulting table with information associating a sample
        name (Name) to a Term, and the normalized and raw enrichment scores
        for each Term
    """
    # grab indexes where condition is true
    maxZ_above_thresh = np.where(maxZ > threshold)
    deltaZ_not_zero = np.where(deltaZ != 0)

    # create gene list will
    gene_list = deltaZ.iloc[
        np.intersect1d(maxZ_above_thresh, deltaZ_not_zero)
    ].sort_values(ascending=False)

    # TODO:
    # 1) ask if `seed` needs to be set, and to what?
    # 2) ask if `scale` should be True
    # 3) ask if we want to pass output director to ssgsea
    # 4) see about if plotting feature is useful and generates the plots we
    # want (`no_plot`)
    return gp.ssgsea(
        data=gene_list,
        gene_sets=gene_sets_file,
        outdir=outdir,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,  # TODO: keep in mind this is here
        weight=threshold,
        threads=threads
    )


# TODO: maybe this is the best place to also load zscores to reduce memory
# usage
def process_scores(scores, pairs) -> pd.DataFrame:
    """Grabs replicates specified `pairs` from scores matrix and processes
    those remaining scores

    Returns a Pandas DataFrame of processed Z scores
    """
    base = 2
    offset = 3
    power = pow(base, offset)
    # collect unique replicates from pairs
    reps_list = []
    for pair in pairs:
        for rep in pair:
            reps_list.append(rep)
    reps_list = list(np.unique(reps_list))
    # exclude unused replicates
    processed_scores = scores.loc[:, reps_list]
    # process scores
    processed_scores = scores.apply(lambda row: power + row, axis=0)
    processed_scores = processed_scores.apply(
        lambda row: row.apply(lambda val: 1 if val < 1 else val),
        axis=0
    )

    return processed_scores.apply(
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
