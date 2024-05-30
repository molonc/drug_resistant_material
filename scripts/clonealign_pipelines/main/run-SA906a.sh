snakemake -s Snakefile-SA906a -j 10 -k -p --use-singularity --singularity-args "--home $HOME -B /cellassign/:/cellassign" --singularity-prefix /cellassign/fitness-scrna/docker/image/ $1

#MA: added --home $HOME for singularity to use this as a mount point

# If the run was interrupted, type snakemake --unlock and then run again 


# -p --use-singularity --singularity-args "-B /juno/work/shah/alzhang:/juno/work/shah/alzhang" --singularity-prefix /juno/work/shah/alzhang/projects/cellassign-paper/.singularity --keep-going
