Rscript run_clonealign.R \
  --sce /cellassign/fitness-scrna/mirela/data/scrna_human_normed/SA609X5XB03231.rdata \
  --cnv /cellassign/fitness-scrna/mirela/data/dlp/pdx_sample_names/SA609X5XB03231.csv \
  --clone_prevs /cellassign/fitness-scrna/mirela/results/new-cut-v1/outputs/parse_cnv/clone_prevalences/SA001X5.csv \
  --samples SA609X5XB03231 \
  --n_extra_genes 200 \
  --max_copy_number 6 \
  --conda_env r-tensorflow \
  --conda_path /home/rstudio/miniconda/bin/conda \
  --clones_to_ignore "None" \
  --mean_counts 0.05 \
  --outfname clonealign_result_SA609X5XB03231.rds


