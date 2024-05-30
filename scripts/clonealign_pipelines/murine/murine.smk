
workflow_name = 'murine'

print("Sub-pipeline murine")

output_dir = os.path.join(config['outdir'], workflow_name)
scratch_dir = os.path.join(config['scratchdir'], workflow_name)
log_dir = os.path.join(config['logdir'], workflow_name)
workspace_dir = config['workspaces'][workflow_name]
# MA: adding these
human_dir = config['human_dir']
mouse_dir = config['mouse_dir']

rule identify_murine:
    input:
#        sce_human=lambda wildcards: "../../data/scrna2/" + \
#                        config['scrnaseq_data']['sce_library_id'][wildcards.sample] + \
#                        ".rdata",
#        sce_mouse=lambda wildcards: "../../data/scrna-mouse/" + \
#                        config['scrnaseq_data']['sce_library_id'][wildcards.sample] + \
#                        ".rdata",
        sce_human=lambda wildcards: human_dir + "/" + \
                        config['scrnaseq_data']['sce_library_id'][wildcards.sample] + \
                        ".rdata",
        sce_mouse=lambda wildcards: mouse_dir + "/" + \
                        config['scrnaseq_data']['sce_library_id'][wildcards.sample] + \
                        ".rdata",
    output:
        sce='{outdir}/sces/{{sample}}.rds'.format(
            outdir=output_dir
        ),
        figure='{outdir}/figs/{{sample}}.png'.format(
            outdir=output_dir
        ),
        csv='{outdir}/counts/{{sample}}.csv'.format(
            outdir=output_dir
        ),
    params:
        workspace=workspace_dir,
    log:
        '{outdir}/murine/{{sample}}.log'.format(
            outdir=log_dir
        ),
    shell:
        'Rscript {params.workspace}/R/identify-murine.R '
        '--sce_human {input.sce_human} '
        '--sce_mouse {input.sce_mouse} '
        '--sample {wildcards.sample} '
        '--output_figure {output.figure} '
        '--output_csv {output.csv} '
        '--output_sce {output.sce} '
        '>& {log}'
