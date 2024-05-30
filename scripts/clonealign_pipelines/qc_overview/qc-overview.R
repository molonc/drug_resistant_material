

library(SingleCellExperiment)
library(dplyr)
library(tidyr)

p <- "../../results/v1/outputs/preprocess/sce_annotated/"

fs <- dir(p, full.names = TRUE)

sces <- lapply(fs, readRDS)

col_names <- c("id", "total_counts", "total_features_by_counts", 
               "pct_counts_mitochondrial", "pct_counts_ribosomal",
               "scrublet_scores", "series")

col_data <- lapply(sces, function(sce) {
  dplyr::select(as.data.frame(colData(sce)), one_of(col_names))
}) %>% 
  bind_rows() %>% 
  as_tibble()

col_data <- gather(col_data, metric, value, -id, -series)

ggplot(col_data, aes(x = id, y = value, fill = series)) +
  geom_boxplot(outlier.size = 0.1, size = .2) +
  facet_wrap(~ metric, scales = "free_x") +
  cowplot::theme_cowplot(font_size = 11) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip()
