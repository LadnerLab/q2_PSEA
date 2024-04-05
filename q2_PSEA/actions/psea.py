import numpy as np
import os
import pandas as pd
import rpy2.robjects as ro
import qiime2
import q2_PSEA.utils as utils
import tempfile

from math import isnan, log, pow
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from scipy import interpolate
from q2_pepsirf.format_types import PepsirfContingencyTSVFormat
from q2_PSEA.actions.r_functions import INTERNAL


cluster_profiler = importr("clusterProfiler")
pandas2ri.activate()


def make_psea_table(
        ctx,
        scores_file,
        pairs_file,
        peptide_sets_file,
        threshold,
        p_val_thresh=0.05,
        es_thresh=0.4,
        species_taxa_file="",
        min_size=15,
        max_size=2000,
        permutation_num=10000,  # as per original PSEA code
        table_dir="./psea_table_outdir",
        # TODO: will expand and notion of an R- or Python-spline will be
        # infered by the plug-in
        spline=["r", "py"],
        threads=4,
        pepsirf_binary="pepsirf"
):    
    volcano = ctx.get_action("ps-plot", "volcano")
    zscatter = ctx.get_action("ps-plot", "zscatter")

    if not os.path.exists(table_dir):
        os.mkdir(table_dir)
    else:
        print(
            f"Warning: the directory '{table_dir}' already exists; files may"
            "be overwritten!"
        )

    with open(pairs_file, "r") as fh:
        pairs = [
            tuple(line.replace("\n", "").split("\t"))
            for line in fh.readlines()
        ]
    scores = pd.read_csv(scores_file, sep="\t", index_col=0)
    processed_scores = process_scores(scores, pairs)

    with tempfile.TemporaryDirectory() as tempdir:
        rread_gmt = ro.r["read.gmt"]
        processed_scores.to_csv(f"{tempdir}/proc_scores.tsv", sep="\t")
        processed_scores = utils.remove_peptides_in_gmt_format(
            processed_scores, peptide_sets_file
        )
        peptide_sets = rread_gmt(peptide_sets_file)

        titles = []
        taxa_access = "species_name"
        p_val_thresholds = []
        used_pairs = []
        pair_spline_dict = {}

        if not species_taxa_file:
            taxa_access = "ID"

        for pair in pairs:
            print(f"Working on pair ({pair[0]}, {pair[1]})...")
        
            spline_tup = r_max_delta_by_spline(
                processed_scores,
                pair
            )
            maxZ = spline_tup[0]
            deltaZ = spline_tup[1]
            spline_x = np.array(spline_tup[2]["x"])
            spline_y = np.array(spline_tup[2]["y"])
            pair_spline_dict[pair[0]] = pd.Series(spline_x)
            pair_spline_dict[pair[1]] = pd.Series(spline_y)
            used_pairs.append(pair)

            table = INTERNAL.psea(
                maxZ,
                deltaZ,
                peptide_sets,
                species_taxa_file,
                threshold,
                permutation_num,
                min_size,
                max_size
            )
            with (ro.default_converter + pandas2ri.converter).context():
                table = ro.conversion.get_conversion().rpy2py(table)
            prefix = f"{pair[0]}~{pair[1]}"
            table.to_csv(
                f"{table_dir}/{prefix}_psea_table.tsv",
                sep="\t", index=False
            )

            taxa = table.loc[:, taxa_access].to_list()
            p_val_thresholds.append(p_val_thresh / len(taxa))

            titles.append(prefix)

            pd.DataFrame(used_pairs).to_csv(
                f"{tempdir}/used_pairs.tsv", sep="\t",
                header=False, index=False
            )
            pd.DataFrame(pair_spline_dict).to_csv(
                f"{tempdir}/timepoint_spline_values.tsv", sep="\t", index=False
            )

        processed_scores_art = ctx.make_artifact(
            type="FeatureTable[Zscore]",
            view=f"{tempdir}/proc_scores.tsv",
            view_type=PepsirfContingencyTSVFormat
        )
    
        scatter_plot, = zscatter(
            zscores=processed_scores_art,
            pairs_file=f"{tempdir}/used_pairs.tsv",
            spline_file=f"{tempdir}/timepoint_spline_values.tsv",
            highlight_data=table_dir,
            highlight_thresholds=p_val_thresholds,
            species_taxa_file=species_taxa_file
        )

        volcano_plot, = volcano(
            xy_dir=table_dir,
            xy_access=["NES", "p.adjust"],
            taxa_access=taxa_access,
            x_threshold=es_thresh,
            y_thresholds=p_val_thresholds,
            xy_labels=["Enrichment score", "Adjusted p-values"],
            titles=titles
        )

    return scatter_plot, volcano_plot


def py_max_delta_by_spline(data, timepoints) -> tuple:
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
    maxZ = np.apply_over_axes(np.max, data.loc[:, timepoints], 1)

    y = data.loc[:, timepoints[0]].to_numpy()
    x = data.loc[:, timepoints[1]].to_numpy()
    # tentative magic number 5 knots came from tutorial linked above
    smooth_spline = spline(5, y)
    deltaZ = y - smooth_spline(x)

    maxZ = pd.Series(
        data=[num for elem in maxZ for num in elem],
        index=data.index
    )
    deltaZ = pd.Series(data=deltaZ, index=data.index)
    return (maxZ, deltaZ)


def r_max_delta_by_spline(data, timepoints) -> tuple:
    """Uses R functions to find the maximum Z score between two points in time,
    and to calculate spline for time points and the delta

    Parameters
    ----------
    data : pd.DataFrame
        Matrix of Z scores for sequence

    pair : Tuple
        Tuple pair of samples for which to run max spline

    Returns
    -------
    tuple
        Contains maximum Z score and the difference (delta) in actual from
        predicted Z scores
    """
    rapply = ro.r["apply"]
    rsmooth_spline = ro.r["smooth.spline"]

    with (ro.default_converter + pandas2ri.converter).context():
        maxZ = rapply(data.loc[:, timepoints], 1, "max")
        # TODO: `spline` object ends up coming back as an OrdDict; when passed
        # to another function, an error is thrown explaining that py2rpy is not
        # defined for rpy2.rlike.containers.OrdDict; this must be revisted to
        # fix having to do the spline operation twice
        spline = rsmooth_spline(
            data.loc[:, timepoints[0]], data.loc[:, timepoints[1]]
        )
        spline = dict(spline)
        deltaZ = INTERNAL.delta_by_spline(
            data.loc[:, timepoints[0]], data.loc[:, timepoints[1]]
        )

    maxZ = pd.Series(data=maxZ, index=data.index.to_list())
    deltaZ = pd.Series(data=deltaZ, index=data.index.to_list())

    return (maxZ, deltaZ, spline)


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

    processed_scores = processed_scores.apply(lambda row: power + row, axis=0)
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
