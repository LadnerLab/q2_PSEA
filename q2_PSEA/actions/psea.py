import gseapy as gp
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import qiime2

from math import pow, log
from rpy2.robjects import pandas2ri
from scipy import interpolate


def generate_metadata(replicates):
    base_reps = []

    replicates.sort()

    for replicate in replicates:
        base_seq_name = replicate.split("_")[2]
        base_reps.append(base_seq_name)

    meta_series = pd.Series(data=base_reps, index=replicates)
    meta_series.index.name = "sample-id"
    meta_series.name = "source"

    return qiime2.metadata.CategoricalMetadataColumn(meta_series)


def make_psea_table(
        ctx,
        scores_file,
        timepoints_file,
        pairs_file,
        peptide_sets_file,
        threshold,  # TODO: ask for suggested default value
        species_tax_file=None,
        min_size=15,
        max_size=2000,
        permutation_num=0,  # TODO: check if our user would be required to pass
        threads=4,
        r_ctrl=False,
        pepsirf_binary="pepsirf"
):
    # collect zscores -> is there a better place to do this?
    scores = pd.read_csv(scores_file, sep="\t", index_col=0)
    # collect pairs
    with open(pairs_file, "r") as fh:
        pairs = [
            tuple(line.replace("\n", "").split("\t"))
            for line in fh.readlines()
        ]
    # collect timepoints
    # with open(timepoints_file, "r") as fh:
    #     timepoints = fh.readlines()[0].strip().split()
    # process scores
    processed_scores = process_scores(scores, pairs, peptide_sets_file)

    # TODO: reimplement loop to cover all defined pairs
    # run PSEA operation for current pair
    spline_tup = max_delta_by_spline(
        processed_scores,
        ["070060_D360.Pro_PV2T", "070060_D540.Pro_PV2T"],
        r_ctrl
    )
    maxZ = spline_tup[0]
    deltaZ = spline_tup[1]
    # spline_x = spline_tup[2]
    # spline_y = spline_tup[3]

    table = psea(
        maxZ=maxZ,
        deltaZ=deltaZ,
        peptide_sets_file=peptide_sets_file,
        species_tax_file=species_tax_file,
        threshold=threshold,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        threads=threads,
        outdir="table_outdir"
    )
    table.to_csv("table.tsv", sep="\t", index=False)

    return qiime2.sdk.Result


def max_delta_by_spline(data, timepoints, r_ctrl) -> tuple:
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
    if r_ctrl:
        rapply = ro.r["apply"]
        rspline = ro.r["smooth.spline"]
        rpredict = ro.r["predict"]

        data_dict = {
            timepoints[0]:
                ro.FloatVector(data.loc[:, timepoints[0]].to_list()),
            timepoints[1]:
                ro.FloatVector(data.loc[:, timepoints[1]].to_list()),
        }
        rdata = ro.DataFrame(data_dict)
        rdata.rownames = data.index.to_list()
        print(f"Data as R DataFrame:\n{rdata}")
        # maxZ = rapply(r_from_pd_data, 1, max)
    else:
        maxZ = np.apply_over_axes(np.max, data.loc[:, timepoints], 1)

        # perform smoothing spline prediction
        y = data.loc[:, timepoints[0]].to_numpy()
        x = data.loc[:, timepoints[1]].to_numpy()
        # tentative magic number 5 knots came from tutorial linked above
        smooth_spline = spline(5, y)
        deltaZ = y - smooth_spline(x)

        # convert maxZ and deltaZ to pandas Series - allows for association of
        # values with peptides
        maxZ = pd.Series(
            data=[num for elem in maxZ for num in elem],
            index=data.index
        )
        deltaZ = pd.Series(data=deltaZ, index=data.index)

    print("Return from max_delta_by_spline...")
    return (maxZ, deltaZ, smooth_spline(x), smooth_spline(y))


def psea(
    maxZ: pd.Series,
    deltaZ: pd.Series,
    peptide_sets_file: str,
    species_tax_file: str,
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

    species_tax_file : str

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
        outdir=outdir,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,  # TODO: keep in mind this is here
        # weight=threshold,  # TODO: verify equivocation
        threads=threads
    )

    # check that species names are important
    if species_tax_file is not None:
        # grab result table and include "Name" column
        pre_res = res.res2d.loc[:, [
            "Name", "Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"
        ]]
        # grab species name to taxanomic ID mapping
        species_tax_map = pd.read_csv(
            species_tax_file,
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
def process_scores(scores, pairs, peptide_sets_file) -> pd.DataFrame:
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
    processed_scores = processed_scores.apply(lambda row: power + row, axis=0)
    processed_scores = processed_scores.apply(
        lambda row: row.apply(lambda val: 1 if val < 1 else val),
        axis=0
    )
    # remove peptides not present in peptide set file
    pep_list = []
    # TODO: maybe I can pull this info out and pass to ssgsea instead of the
    # file name
    with open(peptide_sets_file, "r") as fh:
        lines = [line.replace("\n", "").split("\t") for line in fh.readlines()]
        # remove tax IDs
        for line in lines:
            line.pop(0)
            for pep in line:
                pep_list.append(pep)
    # remove unspecified peptides from processed score matrix
    pep_list = processed_scores.index.difference(pep_list)
    processed_scores = processed_scores.drop(index=pep_list)
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
