
workflow_name = 'de'

output_dir = os.path.join(config['outdir'], workflow_name)
scratch_dir = os.path.join(config['scratchdir'], workflow_name)
log_dir = os.path.join(config['logdir'], workflow_name)
workspace_dir = config['workspaces'][workflow_name]

print(workspace_dir)

#de_figs = expand("{output_dir}/{{sample}}_heatmap.png".format(output_dir=output_dir),
#                  sample=samples)
de_rds = expand("{output_dir}/{{sample}}_de.rds".format(output_dir=output_dir),
                 sample=samples)
de_html = expand("{output_dir}/{{sample}}_de_figures.html".format(output_dir=output_dir),
                 sample=samples)


rule differential_expression:
    input:
        sce='{outdir}/sce_annotated/{{sample}}.rds'.format(
                outdir=os.path.join(config['outdir'], 'preprocess')
         ),
        clonealign_fit='{outdir}/clonealign_fit/{{sample}}.rds'.format(
            outdir=os.path.join(config['outdir'], 'align_clones')
        ),
        gs_hallmark=config['genesets']['hallmark'],
        gs_go=config['genesets']['go'],
        gs_oncogenic=config['genesets']['oncogenic'],
        clone_cn=lambda wildcards: '{outdir}/gene_clone_cn/{group}.tsv'.format(
            outdir=os.path.join(config['outdir'], 'parse_cnv'),
            group=get_sample_group(wildcards.sample)
        ),        
    output:
        rds = "{output_dir}/{{sample}}_de.rds".format(output_dir=output_dir),
    params:
        interesting_clones=lambda wildcards: config['differential_expression']['clones'][wildcards.sample],
        sample='{sample}',
        workspace=workspace_dir,
    # log:
    #     '{outdir}/de/{{sample}}.log'.format(
    #         outdir=log_dir
    #     ),
    script:
        '{params.workspace}/R/run_de.R'

rule de_figures:
    input:
        rds = "{output_dir}/{{sample}}_de.rds".format(output_dir=output_dir),
    params:
        sample = '{sample}',
        workspace=workspace_dir,
    output:
        "{output_dir}/{{sample}}_de_figures.html".format(output_dir=output_dir),
    script:
        '{params.workspace}/R/de-figs.Rmd'        





