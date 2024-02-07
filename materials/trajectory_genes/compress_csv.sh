#!/bin/sh
search_dir=/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/cis_trans_genes_lineages/

input_csvs=$(find $search_dir -type f -name "*.csv")

for entry in $input_csvs
do
  echo "$entry"
  gzip -v ${entry}
  echo "__________"
done

