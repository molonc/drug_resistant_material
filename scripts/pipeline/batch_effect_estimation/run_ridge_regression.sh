#!/bin/sh

configPath=/home/htran/Projects/farhia_project/rnaseq/pipeline/utils
# input_dir=/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation
input_dir=/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/BE_mtx_v2


datatag=SA609
input_fn="${input_dir}/${datatag}_norm.csv.gz"
output_fn="${input_dir}/${datatag}_corrected.csv.gz"
meta_fn="${input_dir}/${datatag}_meta_info.csv"

# datatag=SA535
# input_fn="${input_dir}/${datatag}_cisplatin/${datatag}_norm.csv.gz"
# output_fn="${input_dir}/${datatag}_cisplatin/${datatag}_corrected.csv.gz"
# meta_fn="${input_dir}/${datatag}_cisplatin/${datatag}_meta_info.csv"


# datatag=SA1035
# input_fn="${input_dir}/${datatag}/${datatag}_norm.csv.gz"
# output_fn="${input_dir}/${datatag}/${datatag}_corrected.csv.gz"
# meta_fn="${input_dir}/${datatag}/${datatag}_meta_info.csv"

echo "Batch effect correction using ridge regression"
now=$(date +"%T")
echo "Starting time : $now"

python $configPath/ridge_regression_batch_effect_removal.py --input_fn $input_fn --meta_fn $meta_fn --output_fn $output_fn --datatag $datatag --batch batch --confounder clone

now=$(date +"%T")
echo "End time : $now"
echo "Completed!"


# datatag=SA609
# input_fn="${input_dir}/${datatag}/${datatag}_norm.csv.gz"
# output_fn="${input_dir}/${datatag}/${datatag}_corrected.csv.gz"
# meta_fn="${input_dir}/${datatag}/${datatag}_meta_info.csv"
# 
# echo "Batch effect correction using ridge regression"
# now=$(date +"%T")
# echo "Starting time : $now"
# 
# python $configPath/ridge_regression_batch_effect_removal.py --input_fn $input_fn --meta_fn $meta_fn --output_fn $input_fn --datatag $datatag --batch batch --confounder clone
# 
# now=$(date +"%T")
# echo "End time : $now"
# echo "Completed!"
