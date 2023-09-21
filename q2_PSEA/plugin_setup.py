#! /usr/bin/env python

import q2_PSEA

from q2_PSEA.actions.PSEAcode import PSEA
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
    function=PSEA,
    inputs={},
    input_descriptions=None,
    parameters={
        "data_fn": Str,
        "exclude": Str
    },
    parameter_descriptions={
        # TODO: verify the file is of Z scores
        "data_fn": "Name of Z score matrix file",
        # TODO: verify we want a TSV and that it identifies replicates
        "exclude": "Name of tab-delimited file containing name of replicates"
            " to exclude"
    },
    outputs=[
        ("zscore_scatter", Visualization)
    ],
    output_descriptions=None,
    name="PSEA", # TODO: verify name
    description="" # TODO: get a description
)
