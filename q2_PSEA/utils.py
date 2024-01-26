import pandas as pd


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
