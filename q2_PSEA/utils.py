import pandas as pd
import qiime2


def generate_metadata(replicates):
    """

    Parameters
    ----------

    Returns
    -------
    """
    base_reps = []

    replicates.sort()

    for replicate in replicates:
        base_seq_name = replicate.split("_")[2]
        base_reps.append(base_seq_name)

    meta_series = pd.Series(data=base_reps, index=replicates)
    meta_series.index.name = "sample-id"
    meta_series.name = "source"

    return qiime2.metadata.CategoricalMetadataColumn(meta_series)


def make_metadata(df, length):
    """Given a Pandas DataFrame, and the length of columns, returns Qiime2
    Metadata
    """
    indexes = []
    for i in range(length):
        indexes.append(f"{i}")
    df.index = indexes
    df.index.name = "sample-id"
    return qiime2.Metadata(df)


def save_taxa_leading_peps_file(
        taxa_peps_filepath,
        taxa,
        leading_peps
    ) -> None:
    """

    Parameters
    ----------

    Returns
    -------
    """
    with open(taxa_peps_filepath, "w") as fh:
        for i in range(len(taxa)):
            fh.write(
                taxa[i] + "\t" + leading_peps[i].replace("/", "\t") + "\n"
            )


def remove_peptides(scores, peptide_sets_file, r_ctrl) -> pd.DataFrame:
    """Removes peptides not present in peptide sets file from a matrix of
    Zscores - this is a helper function to control the process based on the
    user's choice to use Python or R for their analysis

    Returns
    -------
    pd.DataFrame
        Contains remaining peptides which were found in the peptide sets file
    """
    if r_ctrl:
        return remove_peptides_in_csv_format(scores, peptide_sets_file)
    else:
        return remove_peptides_in_gmt_format(scores, peptide_sets_file)


def remove_peptides_in_gmt_format(scores, peptide_sets_file) -> pd.DataFrame:
    """Removes peptides not present in GMT formatted file from a matrix of
    Zscores

    Returns
    -------
    pd.DataFrame
        Contains remaining peptides which were found in the peptide sets file
    """
    pep_list = []
    # TODO: maybe I can pull this info out and pass to ssgsea instead of the
    # file name
    with open(peptide_sets_file, "r") as fh:
        lines = [line.replace("\n", "").split("\t") for line in fh.readlines()]
        for line in lines:
            line.pop(0)
            for pep in line:
                pep_list.append(pep)
    pep_list = scores.index.difference(pep_list)
    return scores.drop(index=pep_list)


def remove_peptides_in_csv_format(scores, peptide_sets_file) -> pd.DataFrame:
    """Removes peptides not present in CSV formatted file from a matrix of
    Zscores

    Returns
    -------
    pd.DataFrame
        Contains remaining peptides which were found in the peptide sets file
    """
    pep_list = []
    with open(peptide_sets_file, "r") as fh:
        lines = [line.replace("\n", "").split(",") for line in fh.readlines()]
        lines.pop(0)
        for line in lines:
            pep_list.append(line[0])
    pep_list = scores.index.difference(pep_list)
    return scores.drop(index=pep_list)
