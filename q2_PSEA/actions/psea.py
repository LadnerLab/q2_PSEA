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
    term_gene = pd.read_csv(input, sep="\t")
