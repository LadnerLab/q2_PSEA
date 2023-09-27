import gseapy as gp
import numpy as np
import pandas as pd

def psea(
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
    # grab indexes where condition is true
    maxZ_above_thresh = np.where(maxZ > threshold)
    deltaZ_not_zero = np.where(deltaZ != 0)
    # create gene list
    gene_list = deltaZ.iloc[np.intersect1d(maxZ_above_thresh, deltaZ_not_zero)].sort_values(ascending=False)

    # read and rename columns
    # TODO: make sure we don't want a header
    term_gene = pd.read_csv(input, sep="\t", header=None)

    enr = gp.enrichr(gene_list=gene_list, gene_sets="KEGG_2016", outdir=None)
    print(f"Enrichr: {enr}")
