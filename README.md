
# Drug Resistant Manuscript Analysis
In preparing for publication...

This repo contains all materials and scripts used in manuscript: 
```
Longitudinal tracking of drug induced transcriptomic reprogramming in triple negative breast cancers
Aparicio Lab 2023 
```

- [Overview](#overview)
- [Materials](#materials)
  - [scRNAseq meta data](#scRNAseq-meta-data)
  - [Uploaded data link](#uploaded-data-link)
  - [Fitness results](#fitness-coefficient-results)
  - [Phylogenetic tree results](#Phylogenetic-tree-results)
  - [Reference gene sets](#reference-gene-sets)
 
- [Scripts](#scripts)
  - [Dataset preview](#dataset-preview)
  - [Full run](#full-run)
- [Citation](#citation)


## Overview

The encoding of resistance states in cancer reflects the contributions of genomic and non-genomic variation, however identifying the potential contributions of each has remained problematic. Here we show that **clonal and non-clonal transcriptional dynamics** of TNBC tumours serially exposed to platinum can separate different clonal responses. Pathway analysis shows that cis and trans transcripts converge on EMT and cytokine signaling states associated with resistance. We observe that copy number clones with strong genotype associated fitness under platinum may become fixed in their states, resulting in minimal transcriptional reversion on **drug withdrawal**. In contrast clones with weaker fitness undergo non-genomic transcriptional plasticity. Together the data show that copy number mediated and copy number independent processes contribute to chemotherapeutic drug resistance.
The details of the method are described in the [pre-print](https://www.biorxiv.org/content/uploadme) 



## Materials
Here are all materials that are used in this manuscript. 

### Uploaded data link
- scRNA-seq cellranger alignment libraries are at: [Uploaded Data URL](https://ega-archive.org/studies/EGAS00001007242)


### scRNAseq meta data
Noted: patient ID: Pt1, Pt2, Pt3, Pt4, Pt5, Pt6 correponding to the series id SA501, SA530, SA604, SA609, SA535, SA1035 in fitness previously published paper. 
- Library infos, drug treatment status and time series - passages are noted at [directory](https://github.com/molonc/drug_resistant_material/tree/main/materials/metadata_drug_resistance/). 
- Add SA501, SA530, SA604 here TODO


### Fitness coefficient results
*Clonal fitness inferred from time-series modelling of single-cell cancer genomes* [published paper](http://dx.doi.org/10.1038/s41586-021-03648-3)

- Fitness coefficient [file](https://github.com/molonc/drug_resistant_material/tree/main/materials/fitness_paper_DLP/SUPP_Table2_fitness_coefficients.csv.gz)
- Fitness coefficient with full details of all results from fitness paper[file](https://github.com/molonc/drug_resistant_material/tree/main/materials/fitness_paper_DLP/master_file_fitness_materials_373358_2_data_set_3595534_qnqbt5_results.xlsx)
- 


### Phylogenetic tree results
- Patient 4 - Pt4 (SA609)[file](https://github.com/molonc/drug_resistant_material/tree/main/materials/fitness_paper_DLP/) data is from previous [published paper](http://dx.doi.org/10.1038/s41586-021-03648-3)
- Patient 5 - Pt5 (SA535)[file](https://github.com/molonc/drug_resistant_material/tree/main/materials/fitness_paper_DLP/master_file_fitness_materials_373358_2_data_set_3595534_qnqbt5_results.xlsx) data is from previous [published paper](http://dx.doi.org/10.1038/s41586-021-03648-3)
- Cell clone assignment

### Reference gene sets
Reference genes sets are at [directory](https://github.com/molonc/drug_resistant_material/tree/main/materials/biodatabase/)
- COSMIC cancer-related reference gene set
- Curated cisplatin reference gene set 


### CloneAlign results
- Clone align results for each sample are at [directory](https://github.com/molonc/drug_resistant_material/tree/main/materials/clonealign_results/clonealign/)
- Input data for Main Figure 2 are at [directory](https://github.com/molonc/drug_resistant_material/tree/main/materials/clonealign_plot)

### TMA microarray tissue screening data
- TMA score results at [file](https://github.com/molonc/drug_resistant_material/blob/main/materials/metadata_drug_resistance/TMA20-001%20TMA_FK3%20with%20scores.xls) and from previous [published paper](http://dx.doi.org/10.1038/s41586-021-03648-3)


### Cis trans genes 
- DE analysis are at [directory](https://github.com/molonc/drug_resistant_material/tree/main/materials/cis_trans/) 
- Pathway analysis results are at [directory]https://github.com/molonc/drug_resistant_material/tree/main/materials/pathway

### Main Figure 1, 2, 3 
- UMAP files at at [directory](https://github.com/molonc/drug_resistant_material/tree/main/materials/umap_figs/) 


### Pseudotime analysis results
- Significant genes across pseudotime analysis are at [directory] (https://github.com/molonc/drug_resistant_material/tree/main/materials/trajectory_genes/)


## Scripts


### Clonealign execution script 
### DE analysis
### Cis, trans gene identification
### Pseudotime analysis script



## Citation
```
Farhia Kabeer, Hoa Tran, Mirela Andronescu, Gurdeep Singh, Hakwoo Lee, Sohrab Salehi, Beixi Wang, Justina Biele, Jazmine Brimhall, David Gee, Viviana Cerda, Ciara O'Flanagan, Teresa Algara, Takako Kono, Sean Beatty, Elena Zaikova, Daniel Lai, Eric Lee, Richard Moore, Andrew J. Mungall, IMAXT Consortium,  Marc William, Andrew Roth, Kieran R. Campbell, Sohrab P. Shah, Samuel Aparicio.

Longitudinal tracking of drug induced transcriptomic reprogramming in triple negative breast cancers
Aparicio Lab 2023. 

```


