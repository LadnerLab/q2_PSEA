#! /usr/bin/env python

import q2_PSEA

from q2_PSEA.actions.make_psea_table import make_psea_table
from q2_PSEA.actions.max_delta_by_spline import max_delta_by_spline
from qiime2.plugin import (
    Plugin, Str, Visualization
)

# q2-PSEA plugin object
plugin = Plugin(
    "PSEA", version=q2_PSEA.__version__,
    website="https://github.com/LadnerLab/q2-PSEA.git",
    description="Qiime2 Plugin for PSEA." # TODO: get a description
)

shared_parameters={
    "data": Str,
    "timepoints": Str
}
shared_parameter_descriptions={
    "data": "Name of Z score matrix file.",
    "timepoints": "Name of tab-delimited file containing sequence names to"
        " compare via smooth splining."
}

shared_outputs=[("zscore_scatter", Visualization)]
shared_output_descriptions={
    "zscore_scatter": "Name of plot file visualization comparison between"
        " two samples. This plot includes the smooth spline fit to the given"
        " data and highlights the leading edge peptides for all"
        " significant taxa."
}


# register make_psea_table function
plugin.pipelines.register_function(
    function=make_psea_table,
    inputs={},
    input_descriptions=None,
    parameters=shared_parameters,
    parameter_descriptions=shared_parameter_descriptions,
    outputs=shared_outputs,
    output_descriptions=shared_output_descriptions,
    name="Make PSEA Table",
    description="" # TODO: get a description
)

plugin.methods.register_function(
    function=max_delta_by_spline,
    inputs={},
    input_descriptions=None,
    parameters=shared_parameters,
    parameter_descriptions=shared_parameter_descriptions,
    outputs=shared_outputs,
    output_descriptions=shared_output_descriptions,
    name="Max Delta by Spline",
    description="" # TODO: get a description
)
