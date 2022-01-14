#!/bin/sh

search_dir=/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/testing/

input_csvs=$(find $search_dir -type f -name "*.csv")

for entry in $input_csvs
do
  echo "$entry"
  gzip -v ${entry}
  echo "__________"
done

