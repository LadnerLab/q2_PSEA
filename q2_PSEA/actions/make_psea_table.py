from max_delta_by_spline import max_delta_by_spline
from psea import psea
import pandas as pd

from math import pow, log

# TODO:
# 1. allow user be able to specify a matrix file
# 2. allow user to specify sequences to exclude

# wrap in function called psae
# driver function for PSEA program

# TODO: optimize by using pd.read_csv()
# read data from file
data = []
columns = []
indexes = []
with open("../../example/IM0031_PV2T_25nt_raw_2mm_i1mm_Z-HDI75.tsv", "r") \
        as data_fh:
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
        data.append(split_line) # add data line

# sequences to exclude from data2
toexclude = ["070208_D540.Pro_PV2T","041058_D180.Pro_PV2T","041073_D180.Pro_PV2T","041119_D180.Pro_PV2T"]

base = 2
offset = 3
power = pow(base, offset)

datain = pd.DataFrame(
    data=data, index=indexes, columns=columns,
    dtype=float
) # TODO: check how much precision is needed (float, double (32, 64))?
data1 = datain.apply(lambda row: power + row, axis=0)
data1 = data1.apply(
    lambda row: row.apply(lambda val: 1 if val < 1 else val),
    axis=0
)
# might need to catch an exception for no columns
data2 = data1.drop(columns=toexclude)

data = data2.apply(
    lambda row: row.apply(lambda val: log(val, base) - offset)
)

timepoint1 = "070060_D360.Pro_PV2T"
timepoint2 = "070060_D540.Pro_PV2T"

maxDelta = max_delta_by_spline(timepoint1, timepoint2, data)
maxZ = maxDelta[0]
deltaZ = maxDelta[1]
# TODO: add an option for the user to dictate the threshold
#table = psea(maxZ, deltaZ, 0.75)
table = psea(
    maxZ, deltaZ, 1.00,
    "../../example/input.tsv", "../../example/PV2species.csv"
) # for testing purposes
