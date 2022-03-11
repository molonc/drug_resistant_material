## Question: need a statistical test to get significant pathway related to a small set of genes
## fgsea is not fit for a small set of genes


library(gprofiler2)
# pathway_fn: a hallmark genes symbol gmt file
custom_id <- gprofiler2::upload_GMT_file(pathway_fn)
#gp__6TUj_rpgo_WS8 is custom id return from web after uploading


## fdr option
gostres <- gprofiler2::gost(list(genes_use), organism = 'gp__6TUj_rpgo_WS8', 
                            correction_method='fdr')
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
# HALLMARK_APOPTOSIS
# HALLMARK_INTERFERON_ALPHA_RESPONSE
# HALLMARK_INTERFERON_GAMMA_RESPONSE
# HALLMARK_UV_RESPONSE_DN



## "gSCS" (synonyms: "analytical", "g_SCS") option
gostres <- gprofiler2::gost(list(genes_use), organism = 'gp__6TUj_rpgo_WS8', 
                            correction_method='gSCS')
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
# HALLMARK_APOPTOSIS
# HALLMARK_INTERFERON_GAMMA_RESPONSE
# HALLMARK_UV_RESPONSE_DN


## bonferroni option
gostres <- gprofiler2::gost(list(genes_use), organism = 'gp__6TUj_rpgo_WS8', 
                            correction_method='bonferroni')
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION


## @Kieran: in my bootstrap testing function, the result return HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
## similar to output of bonferroni above

## Results if you use gprofiler web interface with default option gSCS (similar to gSCS output above)
# https://biit.cs.ut.ee/gprofiler/gost
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
# HALLMARK_APOPTOSIS
# HALLMARK_INTERFERON_GAMMA_RESPONSE
# HALLMARK_UV_RESPONSE_DN