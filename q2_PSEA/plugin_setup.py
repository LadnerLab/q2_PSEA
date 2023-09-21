#! /usr/bin/env python

import q2_PSEA

from q2_PSEA.actions.PSEAcode import psea
# TODO: tenative import until it is known what modules are required
from qiime2.plugin import (
    Plugin, Str, Visualization
)

# q2-PSEA plugin object
plugin = Plugin(
    "q2-PSEA", version=q2_PSEA.__version__,
    website="https://github.com/LadnerLab/q2-PSEA.git",
    description="Qiime2 Plugin for PSEA." # TODO: get a description
)

# register PSEA function
plugin.pipelines.register_function(
    function=psea,
    inputs={},
    input_descriptions=None,
    parameters={
        "data": Str,
        "timepoints": Str
    },
    parameter_descriptions={
        "data": "Name of Z score matrix file.",
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
    name="PSEA", # TODO: verify name
    description="**ADD DESCRIPTION**" # TODO: get a description
)
