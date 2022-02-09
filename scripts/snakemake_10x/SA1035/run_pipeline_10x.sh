#!/bin/sh
# SA919
# Tyler corrupt tree encoding version
# echo "Download data"

configPath=/home/htran/Projects/farhia_project/rscript/snakemake_10x

echo "Download data"

snakemake -s $configPath/download_Snakefile_human --configfile $configPath/snakemake_config.yaml -j 8 2> download.log


# snakemake -s $configPath/download_Snakefile_mouse --configfile $configPath/snakemake_config.yaml -j 8 2> download.log



# echo "Complete"
