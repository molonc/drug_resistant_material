source(paste0("/home/htran/Projects/farhia_project/rscript/pipeline/utils/normalize_utils.R"))
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation_v1/')
output_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/')

datatag <- 'SA609'
sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce.rds')

# sce_fn <- paste0(input_dir,'rnaseq_v6/',datatag,'-v6/total_sce_treated.rds')
sce <- readRDS(sce_fn)
dim(sce)

print("SCTransform normalization")
normalize_SCTransform(sce, input_dir, output_dir, return_data=F)