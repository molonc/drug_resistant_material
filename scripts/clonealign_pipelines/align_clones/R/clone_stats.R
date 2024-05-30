
library(tidyverse)
library(cowplot)
library(forcats)
library(yaml)

theme_set(theme_cowplot(font_size = 11))

clone_prevs <- read_csv(snakemake@input[['clone_prevs']])
# clone_prevs <- read_csv("/cellassign/fitness-scrna/results/v1/outputs/parse_cnv/clone_prevalences/SA609.csv")

# metadata <- read_yaml("/cellassign/fitness-scrna/config/metadata/sample_metadata.yaml")
metadata <- read_yaml(snakemake@input[['metadata']])

threshold <- 0.5
# sample <- "SA609X3XB01584"
sample <- snakemake@params[['sample']]
passage <- metadata[[ sample ]]$timepoint

# clonealign_fit <- readRDS("/cellassign/fitness-scrna/results/v1/outputs/align_clones/clonealign_fit/SA609X3XB01584.rds")
clonealign_fit <- readRDS(snakemake@input[['clonealign_fit']])

## DLP clonal prevalences

tp <- gsub("X", "", clone_prevs$time)
p_rna <- gsub("X", "", passage)

dlp_passage <- clone_prevs$time[which.min(abs(as.numeric(tp) - as.numeric(p_rna)))]

clone_prevs <- filter(clone_prevs, time == dlp_passage)

clone_prevs <- mutate(clone_prevs, freq = n / sum(n))


## RNA clone prevalences

clone_probs <- clonealign_fit$ml_params$clone_probs

#' Assign clones
tag_clone <- function(x) {
  clones <- NULL
  x <- sort(x, decreasing = TRUE)
  
  cs <- cumsum(x)
  n_over <- sum(cs > threshold)
  
  if(n_over == length(x)) {
    clones <- names(x[1])
  } else if(n_over == (length(x) - 1)) {
    clones <- names(x[1:2])
  } else {
    clones <- "U"
  }
  
  clones <- sort(clones)
  
  clones <- paste0(clones, collapse="-")
  
  clones
}

clones <- apply(clone_probs, 1, tag_clone)
clone_tbl <- table(clones)

rna_df <- as_tibble(clone_tbl) %>% 
  mutate(freq = n / sum(n))

rna_df <- filter(rna_df, freq > 0.005)


## Make plots

ggplot(clone_prevs, aes(x =fct_reorder(cluster, 1/(freq)), y = freq)) +
  geom_bar(stat = 'identity') +
  labs(x = "Clone", y = "Prevalence", subtitle = paste("From DLP passage", dlp_passage),
       title = sample)

dlp_plot <- last_plot()

ggplot(rna_df, aes(x = fct_reorder(clones, 1/(freq)), y = freq)) +
  geom_bar(stat = 'identity') +
  labs(x = "Clone", y = "Prevalence", subtitle = "From clonealign", title = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

rna_plot <- last_plot()

n_clones_dna <- nrow(clone_prevs)
n_clones_rna <- nrow(rna_df)

plot_grid(dlp_plot, rna_plot,
          rel_widths = c(n_clones_dna + 5, n_clones_rna + 5))

ggsave(snakemake@output[['figure']], width = 0.5 * (n_clones_rna + n_clones_dna), height = 4)


data <- list(
  sample = snakemake@params[['sample']],
  dlp_clone_prevs = clone_prevs,
  rna_clone_prevs = rna_df
)

saveRDS(data, snakemake@output[['data']])

