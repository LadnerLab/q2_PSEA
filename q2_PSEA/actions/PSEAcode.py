import pandas as pd

from math import pow, log


# TODO: make sure this is not supposed to be lower-case
def PSEA(ctx, data_fn, exclude):
    max_delta_by_spline = ctx.get_action("")
    
    # TODO: optimize dataframe creation by using pd.read_csv()
    scores = []
    columns = []
    indexes = []
    with open(data_fn, "r") as data_fh:
        # read matrix into memory
        lines = data_fh.readlines()

        # grab column names from header
        columns = [name for name in lines[0].replace("\n", "").split("\t")]
        columns.pop(0) # remove "Sequence name"
        lines.pop(0) # remove header line

        for line in lines:
            split_line = line.replace("\n", "").split("\t")
            indexes.append(split_line[0]) # add peptide to index list
            split_line.pop(0) # remove peptide name
            scores.append(split_line) # add data line

    base = 2
    offset = 3
    power = pow(base, offset)

    datain = pd.DataFrame(
        data=scores, index=indexes, columns=columns,
        dtype=float
    ) # TODO: check how much precision is needed (float, double (32, 64))?
    data1 = datain.apply(lambda row: power + row, axis=0)
    data1 = data1.apply(
        lambda row: row.apply(lambda val: 1 if val < 1 else val),
        axis=0
    )
    # might need to catch an exception for no columns
    data2 = data1.drop(columns=exclude)

    data = data2.apply(
        lambda row: row.apply(lambda val: log(val, base) - offset)
    )

    # TODO: user should be able to specify
    input = pd.read_csv("../../example/input.csv")
    s = pd.read_csv("../../example/PV2species.csv", header=None)

    # TODO: user should be able to specify
    timepoint1 = "070060_D360.Pro_PV2T"
    timepoint2 = "070060_D540.Pro_PV2T"

    maxDelta = max_delta_by_spline(timepoint1, timepoint2, data)
