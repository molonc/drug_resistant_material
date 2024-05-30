
suppressPackageStartupMessages({
  library(tidyverse)
})

prev_path <- "/cellassign/fitness-scrna/results/new-cut-v1//outputs/align_clones/clonal_prevalences//"

prev_data <- dir(prev_path, full.names = TRUE, patter = "rds")

prevs <- lapply(prev_data, readRDS)

samples <- sapply(prevs, `[[`, 'sample')
names(prevs) <- samples
samples <- sort(samples)

prevs <- prevs[samples]


sort_item <- function(item) {
  rna <- item$dlp_clone_prevs %>% 
    rename(clone = cluster) %>% 
    select(-time) %>% 
    mutate(modality = "scWGS (DLP)", sample = item$sample)
  
  dna <- item$rna_clone_prevs %>% 
    rename(clone = clones) %>% 
    mutate(modality = "scRNA-seq (10x)", sample = item$sample)
  
  bind_rows(rna, dna)
}

df <- lapply(prevs, sort_item) %>% 
  bind_rows() %>% 
  rename(n_cells = n) %>% 
  ungroup()

df <- group_by(df, modality, sample) %>% 
  mutate(freq = n_cells / sum(n_cells)) %>% 
  ungroup()

all_dfs <- list()

## Outrageous R code
append_row <- function(x) {
  if(grepl("-|_", x$clone)) {
    orig_clone <- x$clone
    clones <- strsplit(x$clone, "-|_")[[1]]
    for(c in clones) {
      x$clone <- c
      x$uncertainty <- "Assignment uncertain"
      x$link <- orig_clone
      all_dfs[[length(all_dfs) + 1]] <<- x
    }
  } else {
    x$link <- NA
    x$uncertainty <- "Assignment certain"
    all_dfs[[length(all_dfs) + 1]] <<- x
  }
}

for(i in seq_len(nrow(df))) {
  append_row(df[i,])
}

df_all <- bind_rows(all_dfs)


## Make segment file

needs_link <- filter(df_all, !is.na(link), modality == "scRNA-seq (10x)") %>% 
  count(modality, sample, link)

x <- needs_link[5,]


## Sum weird double

df_all <- group_by(df_all, modality, clone, sample, uncertainty) %>% 
  summarize(freq = sum(freq)) %>% 
  ungroup()



# df_test <- filter(df_all, sample %in% c("SA906_p11a", "SA609X10XB02454"))


ggplot(arrange(df_all, desc(freq)), aes(x = clone, y = sample)) +
  geom_point(aes(size = freq, colour = uncertainty, fill = uncertainty), shape = 21) +
  facet_wrap(~ modality, scales = "free_x") +
  # scale_fill_brewer(palette = "Set1", guide = FALSE) +
  scale_fill_manual(values = c("Assignment certain"="grey50", "Assignment uncertain"="grey90"),
                    name = "clonealign assignment") +
  cowplot::theme_cowplot(font_size = 11) +
  scale_colour_manual(values = c("Assignment certain"="white", "Assignment uncertain"="black"),
                      name = "clonealign assignment") +
  scale_alpha_manual(values = c("Assignment certain"=1, "Assignment uncertain"=0.8),
                     name = "clonealign assignment") +
  labs(x = "Clone", y = "Sample") +
  scale_size(name = "Relative clonal prevalence") +
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = "bold"))
  
 ggsave("../../../tmpfigs/overall-clonal-prevalence.png", width = 7.5, height = 3)


