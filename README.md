
# Drug Resistant Manuscript Analysis

This repository contains all materials and scripts used in manuscript: 
```
Farhia Kabeer, Hoa Tran, Mirela Andronescu, Gurdeep Singh, Hakwoo Lee, Sohrab Salehi, Justina Biele, Jazmine Brimhall, David Gee, Viviana Cerda, Ciara O'Flanagan, Teresa Algara, Takako Kono, Sean Beatty, Elena Zaikova, Daniel Lai, Eric Lee, Richard Moore, Andrew J. Mungall, IMAXT Consortium,  Marc William, Andrew Roth, Kieran R. Campbell, Sohrab P. Shah, Samuel Aparicio.

Single-cell decoding of drug induced transcriptomic reprogramming in triple negative breast cancers
Genome Biology June 2024. 

```

- [Overview](#overview)
- [Materials](#materials)
  - [Uploaded data link](#uploaded-data-link)
  - [Meta data](#meta-data)
  - [Fitness results](#fitness-coefficient-results)
  - [Phylogenetic tree results](#phylogenetic-tree-results)
  - [Reference gene sets](#reference-gene-sets)
  - [Clone alignment results](#clone-alignment-results)
  - [Tissue screening](#tissue-screening)
  - [Cis trans genes](#cis-trans-genes)
  - [Main figures](#main-figures)
  - [Supplementary figures](#supplementary-figures)
  - [Supplementary tables](#supplementary-tables)
  - [Pseudotime analysis results](#pseudotime-analysis-results)
 
- [Scripts](#scripts)
  - [Preprocessing functions](#preprocessing-functions)
  - [Clonealign pipelines](#clonealign-pipelines)
  - [Differential expression analysis](#differential-expression-analysis)
  - [Gene type identification](#gene-type-identification)
  - [Treatment cycles](#treatment-cycles)
  - [Pseudotime analysis script](#pseudotime-analysis-script)
- [Citation](#citation)
- [Licence](#licence)


## Overview

The encoding of resistance states in cancer reflects the contributions of genomic and non-genomic variation, however identifying the potential contributions of each has remained problematic. Here we show that **clonal and non-clonal transcriptional dynamics** of TNBC tumours serially exposed to platinum can separate different clonal responses. Pathway analysis shows that cis and trans transcripts converge on EMT and cytokine signaling states associated with resistance. We observe that copy number clones with strong genotype associated fitness under platinum may become fixed in their states, resulting in minimal transcriptional reversion on **drug withdrawal**. In contrast clones with weaker fitness undergo non-genomic transcriptional plasticity. Together the data show that copy number mediated and copy number independent processes contribute to chemotherapeutic drug resistance.



## Materials

### Uploaded data link
- scRNA-seq cellranger alignment libraries are at: [Uploaded Data URL](https://ega-archive.org/studies/EGAS00001007242)


### Meta data
Noted: patient ID: Pt1, Pt2, Pt3, Pt4, Pt5, Pt6 correspond to the series id SA501, SA530, SA604, SA609, SA535, SA1035 in fitness previously published paper. 
- Library infos, drug treatment status and time series - passages are noted at [materials/metadata_drug_resistance](https://github.com/molonc/drug_resistant_material/tree/main/materials/metadata_drug_resistance/).
We'll add SA501, SA530, SA604 here too.


### Fitness coefficient results
*Clonal fitness inferred from time-series modelling of single-cell cancer genomes* [published paper](http://dx.doi.org/10.1038/s41586-021-03648-3)

- Fitness coefficient file [materials/fitness_paper_DLP/SUPP_Table2_fitness_coefficients.csv.gz](https://github.com/molonc/drug_resistant_material/tree/main/materials/fitness_paper_DLP/SUPP_Table2_fitness_coefficients.csv.gz)
- Fitness coefficient file with full details of all results from fitness paper [materials/fitness_paper_DLP/master_file_fitness_materials_373358_2_data_set_3595534_qnqbt5_results.xlsx](https://github.com/molonc/drug_resistant_material/tree/main/materials/fitness_paper_DLP/master_file_fitness_materials_373358_2_data_set_3595534_qnqbt5_results.xlsx)


### Phylogenetic tree results

- The inferred sitka phylogeny trees for patient Pt4, Pt5, Pt6 (corresponding to SA609, SA535, SA1035) are downloaded from previously [published paper](http://dx.doi.org/10.1038/s41586-021-03648-3), see directory [materials/DLP_trees/](https://github.com/molonc/drug_resistant_material/tree/main/materials/DLP_trees/)
- The inferred sitka phylogeny trees for patient Pt1, Pt2, Pt3 (corresponding to SA501, SA530, SA604) are generated in this manuscript using sitka Bayesian inference algorithm from [sitka published paper](https://peercommunityjournal.org/articles/10.24072/pcjournal.292/), [sitka git repo](https://github.com/UBC-Stat-ML/sitkatree/) and results newick trees are at directory [materials/DLP_trees/](https://github.com/molonc/drug_resistant_material/tree/main/materials/DLP_trees/)


### Reference gene sets
Reference genes sets are at [materials/reference_gene_sets/](https://github.com/molonc/drug_resistant_material/tree/main/materials/reference_gene_sets/)
- COSMIC cancer-related reference gene set
- Curated cisplatin reference gene set 


### Clone alignment results
Cell clone assignment results
- Clonealign results for each sample are at [materials/clonealign_results/clonealign/](https://github.com/molonc/drug_resistant_material/tree/main/materials/clonealign_results/clonealign/)
- Input data for Main Figure 2 are at [materials/clonealign_results](https://github.com/molonc/drug_resistant_material/tree/main/materials/clonealign_results)

### Tissue screening
- TMA score results for TMA microarray tissue screening are at file  [materials/metadata_drug_resistance/TMA20-001%20TMA_FK3%20with%20scores.xls](https://github.com/molonc/drug_resistant_material/blob/main/materials/metadata_drug_resistance/TMA20-001%20TMA_FK3%20with%20scores.xls) and from previous [published paper](http://dx.doi.org/10.1038/s41586-021-03648-3)

### Cis trans genes 
- DE analysis are at [materials/cis_trans/](https://github.com/molonc/drug_resistant_material/tree/main/materials/cis_trans/) 
- Pathway analysis results are at [materials/pathway](https://github.com/molonc/drug_resistant_material/tree/main/materials/pathway)

### Main Figures 
- Figure files at [materials/main_figures/](https://github.com/molonc/drug_resistant_material/tree/main/materials/main_figures/) 
- UMAP files at [materials/umap_figs/](https://github.com/molonc/drug_resistant_material/tree/main/materials/umap_figs/) 

### Supplementary figures 
- Supplementary figure files at [materials/supplementary_figures/](https://github.com/molonc/drug_resistant_material/tree/main/materials/supplementary_figures/) 

### Supplementary tables
- Supplementary tables at [materials/supplementary_tables/](https://github.com/molonc/drug_resistant_material/tree/main/materials/supplementary_tables/) 

### Pseudotime analysis results
- Significant genes across pseudotime analysis are at [materials/trajectory_genes/](https://github.com/molonc/drug_resistant_material/tree/main/materials/trajectory_genes/)


## Scripts

### Preprocessing functions
- Preprocessing scripts are at [scripts/pipeline/utils/](https://github.com/molonc/drug_resistant_material/tree/main/scripts/pipeline/utils/)

### Clonealign pipelines
- Snakemake pipelines and R scripts at [scripts/clonealign_pipelines](https://github.com/molonc/drug_resistant_material/tree/main/scripts/clonealign_pipelines)
### Differential expression analysis
- Scripts are at [scripts/de_edgeR/](https://github.com/molonc/drug_resistant_material/tree/main/scripts/de_edgeR/)
### Gene type identification
- Scripts are at [scripts/cis_trans/](https://github.com/molonc/drug_resistant_material/tree/main/scripts/cis_trans/)
### Treatment cycles
- Scripts for Figures 4 and 5 at [scripts/treatment_cycles/](https://github.com/molonc/drug_resistant_material/tree/main/scripts/treatment_cycles/)
### Pseudotime analysis script
- Scripts for Figure 6, and Supplementary Figures S13, S14 are at [scripts/trajectory_analysis/](https://github.com/molonc/drug_resistant_material/tree/main/scripts/trajectory_analysis/)



## Citation
```
Farhia Kabeer, Hoa Tran, Mirela Andronescu, Gurdeep Singh, Hakwoo Lee, Sohrab Salehi, Justina Biele, Jazmine Brimhall, David Gee, Viviana Cerda, Ciara O'Flanagan, Teresa Algara, Takako Kono, Sean Beatty, Elena Zaikova, Daniel Lai, Eric Lee, Richard Moore, Andrew J. Mungall, IMAXT Consortium,  Marc William, Andrew Roth, Kieran R. Campbell, Sohrab P. Shah, Samuel Aparicio.

Single-cell decoding of drug induced transcriptomic reprogramming in triple negative breast cancers
Genome Biology June 2024. 

```

## Licence
[Apache licence v2.0](https://github.com/molonc/drug_resistant_material/blob/main/LICENCE)

