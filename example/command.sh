qiime psea make-psea-table --p-scores-file example/IM0031_PV2T_25nt_raw_2mm_i1mm_Z-HDI75.tsv \
--p-pairs-file example/pairs.tsv \
--p-peptide-sets-file example/input.csv \
--p-species-tax-file example/species_tax.tsv \
--p-threshold 0.750000 \
--p-min-size 3 \
--p-max-size 5000 \
--p-permutation-num 10000 \
--p-threads 4 \
--output-dir testing
