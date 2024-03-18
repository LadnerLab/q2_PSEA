import gseapy as gp
import numpy as np
import os
import pandas as pd
import rpy2.robjects as ro
import qiime2
import q2_PSEA.utils as utils
import tempfile

from math import isnan, log, pow
from rpy2.robjects import pandas2ri
from scipy import interpolate
from q2_pepsirf.format_types import PepsirfContingencyTSVFormat
from q2_PSEA.actions.r_functions import INTERNAL


pandas2ri.activate()


def make_psea_table(
        ctx,
        scores_file,
        timepoints_file,
        pairs_file,
        peptide_sets_file,
        threshold,
        p_val_thresh=float("nan"),
        es_thresh=0.4,
        species_taxa_file=None,
        min_size=15,
        max_size=2000,
        permutation_num=10000,  # as per original PSEA code
        table_dir="./psea_table_dir",
        # True by default since Python implementation is still being developed
        r_ctrl=True,
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
    with open(timepoints_file, "r") as fh:
        timepoints = fh.readlines()[0].strip().split()
    scores = pd.read_csv(scores_file, sep="\t", index_col=0)
    processed_scores = process_scores(scores, pairs)

    with tempfile.TemporaryDirectory() as tempdir:
        processed_scores.to_csv(f"{tempdir}/proc_scores.tsv", sep="\t")

        if r_ctrl:
            processed_scores = utils.remove_peptides_in_csv_format(
                processed_scores, peptide_sets_file
            )

            titles = []
            taxa_access = "species_name"
            p_val_thresholds = []
            used_pairs = []
            pair_spline_dict = {}
            for pair in pairs:
                if (timepoints[0] not in pair[0]
                    or timepoints[1] not in pair[1]):
                    continue

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
                    peptide_sets_file,
                    species_taxa_file,
                    threshold,
                    permutation_num,
                    min_size,
                    max_size
                )
                with (ro.default_converter + pandas2ri.converter).context():
                    table = ro.conversion.get_conversion().rpy2py(table)
                # TODO: make consistent the column names from both R and Py
                prefix = f"{pair[0]}~{pair[1]}"
                table.to_csv(
                    f"{table_dir}/{prefix}_psea_table.tsv",
                    sep="\t", index=False
                )

                if species_taxa_file:
                    taxa = table.loc[:, "species_name"].to_list()
                else:
                    taxa = table.loc[:, "ID"].to_list()
                    taxa_access = "ID"

                if isnan(p_val_thresh):
                    p_val_thresholds.append(0.05 / len(taxa))

                titles.append(prefix)

            pd.DataFrame(used_pairs).to_csv(
                f"{tempdir}/used_pairs.tsv", sep="\t",
                header=False, index=False
            )
            pd.DataFrame(pair_spline_dict).to_csv(
                f"{tempdir}/timepoint_spline_values.tsv", sep="\t", index=False
            )
        else:
            print(
                "The '--p-r-ctrl' parameter has been unset, please set the"
                " parameter to 'True' as PSEA using Python is still being"
                " worked on."
            )
            # if out_table_name == "":
            #     out_table_name = "py_table.tsv"

            # processed_scores = remove_peptides_in_gmt_format(
            #     processed_scores, peptide_sets_file
            # )

            # # TODO: reimplement loop to cover all defined pairs
            # # run PSEA operation for current pair
            # spline_tup = py_max_delta_by_spline(
            #     processed_scores,
            #     ["070060_D360.Pro_PV2T", "070060_D540.Pro_PV2T"]
            # )
            # maxZ = spline_tup[0]
            # deltaZ = spline_tup[1]

            # table = psea(
            #     maxZ=maxZ,
            #     deltaZ=deltaZ,
            #     peptide_sets_file=peptide_sets_file,
            #     species_taxa_file=species_taxa_file,
            #     threshold=threshold,
            #     min_size=min_size,
            #     max_size=max_size,
            #     permutation_num=permutation_num,
            #     threads=threads,
            #     outdir="table_outdir"
            # )
            # table.to_csv(out_table_name, sep="\t", index=False)

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


def psea(
    maxZ: pd.Series,
    deltaZ: pd.Series,
    peptide_sets_file: str,
    species_taxa_file: str,
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

    peptide_sets_file : str

    species_taxa_file : str

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

    # create peptide list will
    peptide_list = deltaZ.iloc[
        np.intersect1d(maxZ_above_thresh, deltaZ_not_zero)
    ].sort_values(ascending=False)

    # TODO:
    # 1) ask if `seed` needs to be set, and to what?
    # 2) ask if `scale` should be True
    # 3) ask if we want to pass output director to ssgsea
    # 4) see about if plotting feature is useful and generates the plots we
    # want (`no_plot`)
    res = gp.ssgsea(
        data=peptide_list,
        gene_sets=peptide_sets_file,
        sample_norm_method="custom",
        correl_norm_type="rank",
        outdir=outdir,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,  # TODO: keep in mind this is here
        # weight_score_type=1, # for earlier versions
        weight=1,
        threads=threads
    )

    # check that species names are important
    if species_taxa_file is not None:
        # grab result table and include "Name" column
        pre_res = res.res2d.loc[:, [
            "Name", "Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"
        ]]
        # grab species name to taxanomic ID mapping
        species_tax_map = pd.read_csv(
            species_taxa_file,
            sep="\t",
            header=None,
            index_col=0
        )
        # match names in PSEA result table with taxanomic IDs
        names = species_tax_map.index.to_list()
        for i in range(species_tax_map.iloc[:, 0].size):
            id = species_tax_map.iloc[i, 0]
            pre_res["Name"][pre_res.Term == str(id)] = names[i]
    # otherwise, assume species names are not important
    else:
        pre_res = res.res2d.loc[:, [
            "Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"
        ]]

    # collect names of peptides with max z scores above threshold
    maxZ_peptides = maxZ.index.to_list()
    peps_above_thresh = []
    for i in list(maxZ_above_thresh[0]):
        peps_above_thresh.append(maxZ_peptides[i])
    # collect peptides which have max z scores greater than `threshold`
    tested_peps = res.res2d.loc[:, "Lead_genes"].apply(
        # joins peptides in a semicolon (;) delimited string
        lambda peps: ";".join(
            [pep for pep in np.intersect1d(peps.split(";"), peps_above_thresh)]
        )
    )
    # append tested peptides `pre_res`
    pre_res.insert(len(pre_res.columns), "Tested Peptides", tested_peps.values)

    return pre_res


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
