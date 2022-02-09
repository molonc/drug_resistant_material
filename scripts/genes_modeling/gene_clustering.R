# Preparing fission data for lesson

# BiocManager::install("fission")
library(fission)
library(DESeq2)
library(magrittr)
library(tibble)

data("fission")

# modify colnames of the data to be a bit less confusing
colnames(fission) <- colData(fission)$id

# Create DESeq dataset - simple design, no interaction
dds <- DESeqDataSet(fission, design = ~ strain + minute)

# Remove very low count genes - at least half the samples with more than 5 reads
genes_to_keep <- rowSums((counts(dds) > 5)) > 18
dds <- dds[genes_to_keep, ]


#
# Run tests ----
#
# Run DESeq model
dds <- DESeq(dds)

# Make a contrast between first and other time points for WT
test_result <- lapply(c("15", "30", "60", "120", "180"), function(t){
  res <- results(dds, contrast = c("minute", t, "0"), tidy = TRUE, lfcThreshold = 1) %>%
    as_tibble()
  res$comparison <- as.numeric(t)
  return(res)
})
test_result <- do.call("rbind", test_result)

# Rename first column
names(test_result)[1] <- "gene"


#
# Gene counts ----
#
# Extract raw counts
raw_cts <- counts(dds, normalized = TRUE)

# Applying normalization to data (ignoring design) - for clustering, etc.
norm_cts <- vst(dds, blind = TRUE) %>% assay()

# Get sample information
sample_info <- colData(dds) %>% as.data.frame() %>% as_tibble(rownames = "sample")
sample_info <- sample_info[, c("sample", "strain", "minute", "replicate")]  # retain only a few columns


#
# Save objects ----
#
# save(raw_cts, norm_cts, sample_info, test_result,
#      file = "data/fission_data.RData", compress = "bzip2")
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/testing_only/'
# dir.create(save_dir)

raw_cts %>%
  as_tibble(rownames = "gene") %>%
  readr::write_csv(paste0(save_dir,"counts_raw.csv"))

norm_cts %>%
  as_tibble(rownames = "gene") %>%
  readr::write_csv(paste0(save_dir,"counts_transformed.csv"))


readr::write_csv(sample_info, paste0(save_dir,"sample_info.csv"))
readr::write_csv(test_result, paste0(save_dir,"test_result.csv"))



##### setup ####

# load packages
library(tidyverse)

# read the data
trans_cts <- read_csv(paste0(save_dir,"counts_transformed.csv"))
sample_info <- read_csv(paste0(save_dir,"sample_info.csv"))
test_result <- read_csv(paste0(save_dir,"test_result.csv"))

dim(trans_cts)
View(trans_cts[1:3,1:3])
##### get counts for candidate genes ####

# set of candidate genes for clustering
candidate_genes <- test_result %>% 
  filter(padj < 0.01) %>%    # filter table
  pull(gene) %>%             # extract the gene column as a vector
  unique()                   # retain only unique values

length(candidate_genes)
trans_cts_mean_test <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = wt_0_r1:mut_180_r3, names_to = "sample", values_to = "cts")

dim(trans_cts_mean_test)
View(trans_cts_mean_test[1:4,1:3])
length(unique(trans_cts_mean_test$gene))
# Summarise counts 
trans_cts_mean <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = wt_0_r1:mut_180_r3, names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>% 
  # filter to retain only genes of interest
  filter(gene %in% candidate_genes) %>% 
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene, strain, minute) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()

length(unique(trans_cts_mean$gene))


# Create a matrix
hclust_matrix <- trans_cts %>% 
  select(-gene) %>% 
  as.matrix()

# assign rownames
rownames(hclust_matrix) <- trans_cts$gene

hclust_matrix <- hclust_matrix[candidate_genes, ]
hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()


gene_dist <- dist(hclust_matrix)

gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 10, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

cutree(gene_hclust, k = 5)

gene_cluster <- cutree(gene_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)

head(gene_cluster)

trans_cts_cluster <- trans_cts_mean %>% 
  inner_join(gene_cluster, by = "gene")

head(trans_cts_cluster)

trans_cts_cluster %>% 
  ggplot(aes(minute, mean_cts_scaled)) +
  geom_line(aes(group = gene)) +
  facet_grid(rows = vars(strain), cols = vars(cluster))

trans_cts_cluster %>% 
  ggplot(aes(minute, mean_cts_scaled)) +
  geom_line(aes(group = gene), alpha = 0.3) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(strain), cols = vars(cluster))




