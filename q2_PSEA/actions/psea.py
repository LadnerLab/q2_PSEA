import numpy as np
import os
import pandas as pd
import rpy2.robjects as ro
import q2_PSEA.actions.splines as splines
import qiime2
import q2_PSEA.utils as utils
import tempfile

import time
import concurrent.futures
import multiprocessing

from math import isnan, log, pow
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from q2_pepsirf.format_types import PepsirfContingencyTSVFormat
from q2_PSEA.actions.r_functions import INTERNAL


pandas2ri.activate()


def make_psea_table(
        ctx,
        scores_file,
        pairs_file,
        peptide_sets_file,
        threshold,
        p_val_thresh=0.05,
        nes_thresh=1,
        species_taxa_file="",
        min_size=15,
        max_size=2000,
        permutation_num=10000,  # as per original PSEA code
        spline_type="r-smooth",
        degree=3,
        dof=None,
        table_dir="./psea_table_outdir",
        pepsirf_binary="pepsirf",
        iterative_analysis=True,
        iter_tables_dir="",
        max_workers=None,
        event_summary=True,
        taxa_matrix_out="taxa_matrix.tsv",
        positive_nes_vaes_out="Positive_NES_VAEs.txt",
        negative_nes_vaes_out="Negative_NES_VAEs.txt"
):    
    start_time = time.perf_counter()

    volcano = ctx.get_action("ps-plot", "volcano")
    zscatter = ctx.get_action("ps-plot", "zscatter")

    assert spline_type in splines.SPLINE_TYPES, \
        f"'{spline_type}' is not a valid spline method!"
    assert not os.path.exists(table_dir), \
        f"'{table_dir}' already exists! Please move or remove this directory."
    if iterative_analysis:
        assert ".gmt" in peptide_sets_file.lower(), \
            "You are running iterative analysis without a GMT peptide sets file."
    assert max_workers <= multiprocessing.cpu_count(), \
        f"Max workers excedes {multiprocessing.cpu_count()}, the number of CPUs on your machine."

    os.mkdir(table_dir)

    with open(pairs_file, "r") as fh:
        pairs = [
            tuple(line.replace("\n", "").split("\t"))
            for line in fh.readlines()
        ]
    scores = pd.read_csv(scores_file, sep="\t", index_col=0)
    processed_scores = process_scores(scores, pairs)

    scores_file_split = scores_file.rsplit("/", 1)
    if len(scores_file_split) > 1:
        processed_scores_file = f"transformed_{scores_file_split[1]}"
    else:
        processed_scores_file = f"transformed_{scores_file_split[0]}"

    # temporary directory to hold iterative analysis tables
    with tempfile.TemporaryDirectory() as temp_peptide_sets_dir:
        if iterative_analysis:
            if iter_tables_dir:
                if not os.path.exists(iter_tables_dir):
                    os.mkdir(iter_tables_dir)
                else:
                    print(
                        f"\nWarning: the directory '{iter_tables_dir}' already exists; files may"
                        " be overwritten!"
                    )

            # run iterative peptide analysis, writes to temp_peptide_sets_dir files
            pair_pep_sets_file_dict = run_iterative_peptide_analysis(
                pairs=pairs,
                processed_scores=processed_scores,
                og_peptide_sets_file=peptide_sets_file,
                species_taxa_file=species_taxa_file,
                threshold=threshold,
                permutation_num=permutation_num,
                min_size=min_size,
                max_size=max_size,
                spline_type=spline_type,
                degree=degree,
                dof=dof,
                p_val_thresh=p_val_thresh,
                nes_thresh=nes_thresh,
                peptide_sets_out_dir = temp_peptide_sets_dir,
                iter_tables_dir=iter_tables_dir,
                max_workers=max_workers
            )
        else:
            # each pair will have the same gmt file
            pair_pep_sets_file_dict = dict.fromkeys(pairs, peptide_sets_file)

        with tempfile.TemporaryDirectory() as tempdir:
            event_matrix = dict()
            empty_pair_row = [0] * len(pairs)
            pos_nes_count_dict = dict()
            neg_nes_count_dict = dict()
            zero_nes_count_dict = dict()

            processed_scores.to_csv(processed_scores_file, sep="\t")

            titles = []
            taxa_access = "species_name"
            used_pairs = []
            pair_spline_dict = { "x": list(), "y": list(), "pair": list() }


            if not dof:
                dof = ro.NULL
            if not species_taxa_file:
                taxa_access = "ID"
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
                pair_futures = [executor.submit(create_fgsea_table_for_pair,
                                pair,
                                processed_scores,
                                pair_pep_sets_file_dict[ pair ],
                                species_taxa_file,
                                threshold,
                                permutation_num,
                                min_size,
                                max_size,
                                spline_type,
                                degree,
                                dof,
                                p_val_thresh,
                                nes_thresh,
                                False,
                                table_dir
                                ) for pair in pairs]

                for future in concurrent.futures.as_completed(pair_futures):
                    result = future.result()
                    x = result[0]
                    yfit = result[1]
                    table_prefix = result[2]
                    pair = result[3]

                    pair_spline_dict["x"].extend(x.tolist())
                    pair_spline_dict["y"].extend(yfit.tolist())
                    pair_spline_dict["pair"].extend([table_prefix] * len(x))
                    used_pairs.append(pair)

                    # taxa = table.loc[:, taxa_access].to_list()

                    titles.append(table_prefix)

                    if event_summary:
                        # populate event matrix with species that are significant this pair
                        tableDf = pd.read_csv(f"{table_dir}/{table_prefix}_psea_table.tsv", sep="\t")
                        for i, row in tableDf.iterrows():
                            taxa = row[taxa_access]
                            if taxa not in event_matrix.keys():
                                event_matrix[taxa] = empty_pair_row.copy()

                            if row["p.adjust"] < p_val_thresh and np.absolute(row["NES"]) > nes_thresh:
                                event_matrix[taxa][pairs.index(pair)] = 1

                                if row["NES"] > 0:
                                    if taxa not in pos_nes_count_dict.keys():
                                        pos_nes_count_dict[taxa] = 0
                                    pos_nes_count_dict[taxa] += 1
                                elif row["NES"] < 0:
                                    if taxa not in neg_nes_count_dict.keys():
                                        neg_nes_count_dict[taxa] = 0
                                    neg_nes_count_dict[taxa] += 1
                                # output message if nes is 0
                                else:
                                    if taxa not in zero_nes_count_dict.keys():
                                        zero_nes_count_dict[taxa] = 0
                                    zero_nes_count_dict[taxa] += 1

            if event_summary:
                event_matrix_df = pd.DataFrame.from_dict(event_matrix)
                event_matrix_df.index = [f"{pair[0]}~{pair[1]}" for pair in pairs]
                # sort rows and columns
                event_matrix_df.sort_index(inplace=True)
                sorted_cols = sorted(event_matrix_df.columns.tolist(), key=lambda col: event_matrix_df[col].sum(), reverse=True)
                event_matrix_df = event_matrix_df[sorted_cols]
                event_matrix_df.to_csv(taxa_matrix_out, sep="\t")

                pos_nes_count_dict = {k:v for k, v in sorted(pos_nes_count_dict.items(), key=lambda item: item[1], reverse=True)}
                neg_nes_count_dict = {k:v for k, v in sorted(neg_nes_count_dict.items(), key=lambda item: item[1], reverse=True)}
                # create column sums for positive and negative NES
                with open(positive_nes_vaes_out, "w") as pos_file:
                    for taxa in pos_nes_count_dict.keys():
                        if pos_nes_count_dict[taxa] > 1:
                            pos_file.write(f"{pos_nes_count_dict[taxa]} events for {taxa}\n")
                        else:
                            pos_file.write(f"{pos_nes_count_dict[taxa]} event for {taxa}\n")
                with open(negative_nes_vaes_out, "w") as neg_file:
                    for taxa in neg_nes_count_dict.keys():
                        if neg_nes_count_dict[taxa] > 1:
                            neg_file.write(f"{neg_nes_count_dict[taxa]} events for {taxa}\n")
                        else:
                            neg_file.write(f"{neg_nes_count_dict[taxa]} event for {taxa}\n")
                if len(zero_nes_count_dict) > 0:
                    print("\n")
                    for taxa in zero_nes_count_dict.keys():
                        if zero_nes_count_dict[taxa] > 1:
                            print(f"{zero_nes_count_dict[taxa]} events for {taxa}, which has an NES of 0")
                        else:
                            print(f"{zero_nes_count_dict[taxa]} event for {taxa}, which has an NES of 0")

            pd.DataFrame(used_pairs).to_csv(
                f"{tempdir}/used_pairs.tsv", sep="\t",
                header=False, index=False
            )
            pd.DataFrame(pair_spline_dict).to_csv(
                f"{tempdir}/spline_data.tsv", sep="\t", index=False
            )

            processed_scores_art = ctx.make_artifact(
                type="FeatureTable[Zscore]",
                view=processed_scores_file,
                view_type=PepsirfContingencyTSVFormat
            )
        
            scatter_plot, = zscatter(
                zscores=processed_scores_art,
                pairs_file=f"{tempdir}/used_pairs.tsv",
                spline_file=f"{tempdir}/spline_data.tsv",
                p_val_access="p.adjust",
                le_peps_access="core_enrichment",
                taxa_access=taxa_access,
                highlight_data=table_dir,
                highlight_threshold=p_val_thresh
            )

            volcano_plot, = volcano(
                xy_dir=table_dir,
                xy_access=["NES", "p.adjust"],
                taxa_access=taxa_access,
                x_threshold=nes_thresh,
                y_threshold=p_val_thresh,
                xy_labels=["Enrichment score", "Adjusted p-values"],
                titles=titles
            )

    end_time = time.perf_counter()

    print(f"\nFinished in {round(end_time-start_time, 2)} seconds")
    
    return scatter_plot, volcano_plot


def create_fgsea_table_for_pair(
    pair,
    processed_scores,
    pep_sets_file,
    species_taxa_file,
    threshold,
    permutation_num,
    min_size,
    max_size,
    spline_type,
    degree,
    dof,
    p_val_thresh,
    nes_thresh,
    iteration,
    table_dir=""
    ):
    print(f"Working on pair ({pair[0]}, {pair[1]})...")

    table_prefix = f"{pair[0]}~{pair[1]}"

    # each pair has a different peptide set file
    processed_scores, peptide_sets = utils.remove_peptides(
        processed_scores, pep_sets_file
    )

    data_sorted = processed_scores.loc[:, pair].sort_values(by=pair[0])
    x = data_sorted.loc[:, pair[0]].to_numpy()
    y = data_sorted.loc[:, pair[1]].to_numpy()

    # TODO: optimize with a dictionary, if possible
    if spline_type == "py-smooth":
        yfit = splines.smooth_spline(x, y)
    elif spline_type == "cubic":
        yfit = splines.R_SPLINES.cubic_spline(x, y, degree, dof)
    else:
        yfit = splines.R_SPLINES.smooth_spline(x, y)

    maxZ = np.apply_over_axes(
        np.max,
        data_sorted.loc[:, pair],
        1
    )
    maxZ = pd.Series(
        data=[num for elem in maxZ for num in elem],
        index=data_sorted.index
    )
    deltaZ = pd.Series(
        data=y - yfit, index=data_sorted.index
    )

    table = INTERNAL.psea(
        maxZ,
        deltaZ,
        peptide_sets,
        species_taxa_file,
        threshold,
        permutation_num,
        min_size,
        max_size
    )
    with (ro.default_converter + pandas2ri.converter).context():
        table = ro.conversion.get_conversion().rpy2py(table)

    if iteration:
        return table

    # Note: I had to write table here instead of main because of error with rpy2
    table.to_csv(
        f"{table_dir}/{table_prefix}_psea_table.tsv",
        sep="\t", index=False
                    )
    return x, yfit, table_prefix, pair


def process_scores(scores, pairs) -> pd.DataFrame:
    """Grabs replicates specified `pairs` from scores matrix and processes
    those remaining scores

    Returns a Pandas DataFrame of processed Z scores
    """
    base = 2
    offset = 3
    power = pow(base, offset)
    # collect unique replicates from pairs
    reps_list = []
    for pair in pairs:
        for rep in pair:
            reps_list.append(rep)
    reps_list = list(np.unique(reps_list))
    # exclude unused replicates
    processed_scores = scores.loc[:, reps_list]

    processed_scores = processed_scores.apply(lambda row: power + row, axis=0)
    processed_scores = processed_scores.apply(
        lambda row: row.apply(lambda val: 1 if val < 1 else val),
        axis=0
    )
    return processed_scores.apply(
        lambda row: row.apply(lambda val: log(val, base) - offset)
    )


def run_iterative_peptide_analysis(
    pairs,
    processed_scores,
    og_peptide_sets_file,
    species_taxa_file,
    threshold,
    permutation_num,
    min_size,
    max_size,
    spline_type,
    degree,
    dof,
    p_val_thresh,
    nes_thresh,
    peptide_sets_out_dir,
    iter_tables_dir,
    max_workers
    ) -> dict:

    iteration_num = 1
    
    # initialize gmt dict
    gmt_dict = dict()
    with open(og_peptide_sets_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split("\t")
            # species : set of peptides
            gmt_dict[ line[0] ] = set( line[2:] )

    # each pair should have its own gmt
    pair_gmt_dict = dict()
    # keep track of which pairs are fully processed
    sig_species_found_dict = dict()
    # keep a dict for each pair and it's tested species
    tested_species_dict = dict()
    for pair in pairs:
        pair_gmt_dict[ pair ] = gmt_dict.copy()
        sig_species_found_dict[ pair ] = True
        tested_species_dict[ pair ] = set()

    # keep a dict for the output gmt file of each pair
    pair_sets_filename_dict = dict()

    # loop until no other significant peptides were found
    while any(sig_species_found_dict.values()):

        print("\nIteration:", iteration_num)

        if iter_tables_dir:
            iter_out_dir = f"{iter_tables_dir}/Iteration_{iteration_num}"
            if not os.path.exists(iter_out_dir):
                os.mkdir(iter_out_dir)
        else:
            iter_out_dir = ""

        # -------------------------------
        # note: rpy2 is not compatible with multithreading, only multiprocessing
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            pair_futures = [executor.submit(run_iterative_process_single_pair,
                            pair, 
                            tested_species_dict[pair], 
                            pair_gmt_dict[pair], 
                            processed_scores,
                            species_taxa_file,
                            threshold,
                            permutation_num,
                            min_size,
                            max_size,
                            spline_type,
                            degree,
                            dof,
                            p_val_thresh,
                            nes_thresh,
                            peptide_sets_out_dir,
                            iter_out_dir
                            ) for pair in pairs if sig_species_found_dict[pair]]

            # concurrent.futures.wait(pair_futures, timeout=None, return_when=concurrent.futures.ALL_COMPLETED)

            for future in concurrent.futures.as_completed(pair_futures):
                res_tup = future.result()
                pair = res_tup[4]
                sig_species_found_dict[pair] = res_tup[0]
                pair_sets_filename_dict[pair] = res_tup[1]
                tested_species_dict[pair] = res_tup[2]
                pair_gmt_dict[pair]= res_tup[3]        
        # -------------------------------

        iteration_num += 1

    print("\nEnd of Iterative Peptide Analysis\n")
    return pair_sets_filename_dict


def run_iterative_process_single_pair(
    pair, 
    tested_species, 
    gmt_dict, 
    processed_scores,
    species_taxa_file,
    threshold,
    permutation_num,
    min_size,
    max_size,
    spline_type,
    degree,
    dof,
    p_val_thresh,
    nes_thresh,
    peptide_sets_out_dir,
    iter_out_dir
    ):
    if not dof:
        dof = ro.NULL

    sig_species_found = False

    # create gmt file for this pair (filtered gmt from prev iteration)
    pair_sets_filename = f"{peptide_sets_out_dir}/{pair[0]}_{pair[1]}".replace(".","-") + ".gmt"

    write_gmt_from_dict(pair_sets_filename, gmt_dict)

    table = create_fgsea_table_for_pair(
                                        pair=pair,
                                        processed_scores=processed_scores,
                                        pep_sets_file=pair_sets_filename,
                                        species_taxa_file=species_taxa_file,
                                        threshold=threshold,
                                        permutation_num=permutation_num,
                                        min_size=min_size,
                                        max_size=max_size,
                                        spline_type=spline_type,
                                        degree=degree,
                                        dof=dof,
                                        p_val_thresh=p_val_thresh,
                                        nes_thresh=nes_thresh,
                                        iteration = True
                                        )

    # sort the table by ascending p-value (lowest on top)
    table.sort_values(by=["p.adjust"], ascending=True)

    if iter_out_dir:
        table.to_csv(f"{iter_out_dir}/{pair}.tsv", sep="\t")

    # iterate through each row
    for index, row in table.iterrows():

        # test for significant species that has not already been used for this pair
        if row["p.adjust"] < p_val_thresh and np.absolute(row["NES"]) > nes_thresh \
                                    and row["ID"] not in tested_species:

            print(f"Found {row['species_name']} in {pair} to be significant")

            # set sig species found to true for this pair
            sig_species_found = True

            tested_species.add(row["ID"])

            # take out all tested peptides from peptide set in gmt for all other species in the gmt
            all_tested_peps = set(row["all_tested_peptides"].split("/"))
            for gmt_species in gmt_dict.keys():
                if gmt_species != row['ID']:
                    gmt_dict[ gmt_species ] = gmt_dict[ gmt_species ] - all_tested_peps

            # only get top significant species
            break

    return (sig_species_found, pair_sets_filename, tested_species, gmt_dict, pair)


def write_gmt_from_dict(outfile_name, gmt_dict)->None:
    with open( outfile_name, "w" ) as gmt_file:
        for species in gmt_dict.keys():
            gmt_file.write(f"{species}\t\t")

            for peptide in gmt_dict[ species ]:
                gmt_file.write(f"{peptide}\t")

            gmt_file.write("\n")
