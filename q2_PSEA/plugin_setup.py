#! /usr/bin/env python

import q2_PSEA

from q2_PSEA.actions.psea import make_psea_table
from qiime2.plugin import (
    Float, Int, Plugin, Str, Visualization
)

# q2-PSEA plugin object
plugin = Plugin(
    "psea", version=q2_PSEA.__version__,
    website="https://github.com/LadnerLab/q2-PSEA.git",
    description="Qiime2 Plugin for PSEA." # TODO: get a description
)


# register make_psea_table function
plugin.pipelines.register_function(
    function=make_psea_table,
    inputs={},
    input_descriptions=None,
    parameters={
        "scores_file": Str,
        "timepoints_file": Str,
        "pairs_file": Str,
        "peptide_sets_file": Str,
        "species_tax_file": Str,
        "threshold": Float,
        "min_size": Int,
        "max_size": Int,
        "permutation_num": Int,
        "threads": Int,
        "pepsirf_binary": Str
    },
    parameter_descriptions={
        "scores_file": "Name of Z score matrix file.",
        "timepoints_file": "Name of tab-delimited file containing information"
            " referencing the time a sample was taken.",
        "pairs_file": "Name of tab-delimited file containing pairs of"
            " sample names.",
        "peptide_sets_file": "Name of GMT formatted file containing"
            " information about species and their peptides.",
        "species_tax_file": "Name of tab-delimited file containing species"
            " name and taxanomy ID associations.",
        "threshold" : "",
        "min_size": "Minimum allowed number of peptides from peptide set also"
            " the data set.",
        "max_size": "Maximum allowed number of peptides from peptide set also"
            " the data set.",
        "permutation_num": "Number of permutations. Minimal possible nominal"
            " p-value is about 1/perm.",
        "threads": "Number of threads with which to run ssGSEA operation.",
        "pepsirf_binary": "Path to pepsirf binary."
    },
    outputs=[("zenrich_plot", Visualization)],
    output_descriptions={
        "zenrich_plot": "Name of plot file visualization comparison between"
            " two samples. This plot includes the smooth spline fit to the given"
            " data and highlights the leading edge peptides for all"
            " significant taxa."
    },
    name="Make PSEA Table",
    description="" # TODO: get a description
)
