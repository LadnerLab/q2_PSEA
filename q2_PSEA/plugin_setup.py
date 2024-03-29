#! /usr/bin/env python

import q2_PSEA


from q2_PSEA.actions.psea import make_psea_table
from qiime2.plugin import (
    Bool, Float, Int, Plugin, Str, Visualization
)


# q2-PSEA plugin object
plugin = Plugin(
    "psea", version=q2_PSEA.__version__,
    website="https://github.com/LadnerLab/q2-PSEA.git",
    description="Qiime2 Plugin for PSEA."  # TODO: get a description
)


# register make_psea_table function
plugin.pipelines.register_function(
    function=make_psea_table,
    inputs={},
    input_descriptions=None,
    parameters={
        "scores_file": Str,
        "pairs_file": Str,
        "peptide_sets_file": Str,
        "species_taxa_file": Str,
        "threshold": Float,
        "p_val_thresh": Float,
        "es_thresh": Float,
        "min_size": Int,
        "max_size": Int,
        "permutation_num": Int,
        "table_dir": Str,
        "r_ctrl": Bool,
        "threads": Int,
        "pepsirf_binary": Str
    },
    parameter_descriptions={
        "scores_file": "Name of Z score matrix file.",
        "pairs_file": "Name of tab-delimited file containing pairs of"
            " sample names.",
        "peptide_sets_file": "Name of file containing information about"
            " species and their peptides. If `--p-r-ctrl` is False, the file"
            " should be in GMT format; otherwise, the format should be CSV."
            " Please refer to 'input.gmt' and 'input.csv' files in the"
            " 'examples' directory for example files.",
        "species_taxa_file": "Name of tab-delimited file containing species"
            " name and taxanomy ID associations.",
        "threshold": "Minimum Z score a peptide must maintain to be"
            " considered in Gene Set Enrichment Analysis.",
        "p_val_thresh": "Specifies the value adjusted p-values must meet to be"
            " considered for highlighting in volcano and scatter plots. If a"
            " threshold is not provided, then the theshold is set to 0.05 / N"
            " - where N is the number of taxa being tested.",
        "es_thresh": "Specifies the value ",
        "min_size": "Minimum allowed number of peptides from peptide set also"
            " the data set.",
        "max_size": "Maximum allowed number of peptides from peptide set also"
            " the data set.",
        "permutation_num": "Number of permutations. Minimal possible nominal"
            " p-value is about 1/perm.",
        "table_dir": "Directory where resulting PSEA tables will be stored.",
        "r_ctrl": "Specifies to run PSEA using Python or R functions. If set"
            " to True, then R functions will be used.",
        "threads": "Number of threads with which to run ssGSEA operation.",
        "pepsirf_binary": "Path to pepsirf binary."
    },
    outputs=[("scatter_plot", Visualization), ("volcano_plot", Visualization)],
    output_descriptions={
        "scatter_plot": "Name of plot file visualization comparison between"
            " two samples. This plot includes the smooth spline fit to the"
            " given data and highlights the leading edge peptides for all"
            " significant taxa.",
        "volcano_plot": "Name of plot file visualization comparison between"
            " enrichment scores (ES) and p-values."
    },
    name="Make PSEA Table",
    description="Qiime2 plug-in which provides a Python wrapper around R"
        " functions to perform Peptide Set Enrichment Analysis. A Z score"
        " scatter plot for each sample replicate pair complete with a spline"
        " and highlighting of significant taxa. As well as a volcano plot"
        " showing relationship between adjusted p-values and normalized"
        " enrichment scores (NES) with highlighting of points which pass"
        " provided thresholds."
)
