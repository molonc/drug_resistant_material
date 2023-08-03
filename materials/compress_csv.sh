#!/bin/sh

search_dir=/home/htran/Projects/farhia_project/drug_resistant_material/materials/Supplementary_Tables/

input_csvs=$(find $search_dir -type f -name "*.csv")

for entry in $input_csvs
do
  echo "$entry"
  gzip -v ${entry}
  echo "__________"
done

