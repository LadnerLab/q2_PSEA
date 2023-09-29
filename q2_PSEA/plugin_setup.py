#! /usr/bin/env python

import q2_PSEA

from q2_PSEA.actions.psea import make_psea_table
from qiime2.plugin import (
    Plugin, Str, Visualization
)

# q2-PSEA plugin object
plugin = Plugin(
    "PSEA", version=q2_PSEA.__version__,
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
        "gene_set": Str,
        "timepoints": Str
    },
    parameter_descriptions={
        "scores_file": "Name of Z score matrix file.",
        "gene_set": "Name of file containing information about viruses and"
            " their peptides.",
        "timepoints": "Name of tab-delimited file containing sequence names to"
            " compare via smooth splining."
    },
    outputs=[("zscore_scatter", Visualization)],
    output_descriptions={
        "zscore_scatter": "Name of plot file visualization comparison between"
            " two samples. This plot includes the smooth spline fit to the given"
            " data and highlights the leading edge peptides for all"
            " significant taxa."
    },
    name="Make PSEA Table",
    description="" # TODO: get a description
)
