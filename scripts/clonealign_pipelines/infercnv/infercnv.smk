workflow_name = 'infercnv'

output_dir = os.path.join(config['outdir'], workflow_name)
scratch_dir = os.path.join(config['scratchdir'], workflow_name)
log_dir = os.path.join(config['logdir'], workflow_name)
workspace_dir = config['workspaces'][workflow_name]

infercnv_figs = expand('{outdir}/{{sample}}/{{sample}}_all_genes.png'.format(
            outdir=output_dir
            ), sample=samples)



rule infercnv:
    input:
        sce="{outdir}/preprocess/sce_annotated/{{sample}}.rds".format(
            outdir=config['outdir']
            ),
        gtex="../../data/external/recount2/rse_gene_breast.Rdata",
    output:
        rds="{outdir}/{{sample}}.rds".format(outdir=scratch_dir),
        infercnv_dir= directory(os.path.join(scratch_dir, "{sample}")),
    log:
        '{outdir}/infercnv/{{sample}}.log'.format(
            outdir=log_dir
        ),
    script:
        "inferCNV.R"


rule plot_cnv:
    input:
        rds="{outdir}/{{sample}}.rds".format(outdir=scratch_dir),        
        clonealign_fit='{outdir}/align_clones/clonealign_fit/{{sample}}.rds'.format(
            outdir=config['outdir']
        ),
    output:
        png_all_genes = '{outdir}/{{sample}}/{{sample}}_all_genes.png'.format(
            outdir=output_dir
            ),
        png_clonealign_only = '{outdir}/{{sample}}/{{sample}}_clonealign_genes.png'.format(
            outdir=output_dir
            ),
    script:
        "plotCNV.R"
