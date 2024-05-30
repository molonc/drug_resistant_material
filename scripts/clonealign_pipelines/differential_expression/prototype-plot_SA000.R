library(tidyverse)
library(glue)
library(cowplot)


source("de-utils.R")

## Get prevalences ---

prevalence_dir <- "/cellassign/fitness-scrna/results/new-cut-v1/outputs/align_clones/clonal_prevalences/"

fix_prev <- function(l) {
  inner_join(
    mutate(l$dlp_clone_prevs, sample = l$sample) %>% 
      select(clone = cluster, freq, sample),
    mutate(l$rna_clone_prevs, sample = l$sample) %>% 
      select(clone = clones, freq, sample),
    by = c("clone", "sample"),
    suffix = c("_DNA", "_RNA")
  )
}

df_prev <- bind_rows(
  fix_prev(readRDS(file.path(prevalence_dir, "SA609X4XB003083.rds"))),
  fix_prev(readRDS(file.path(prevalence_dir, "SA609X5XB03230.rds"))),
  fix_prev(readRDS(file.path(prevalence_dir, "SA609X6XB03404.rds")))
)

filter(df_prev, freq_DNA > 0.01, freq_RNA > 0.01) %>% 
  ggplot(aes(x = freq_DNA, y = freq_RNA, colour = clone)) +
#  geom_point() +
  facet_wrap(~ sample, ncol = 1) +
  geom_text(aes(label = clone)) +
  scale_colour_manual(values = get_cluster_colours(), guide=F) +
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold')) +
  labs(x = "Clonal prevalence DLP", y = "Clonal prevalence CloneAlign")

prevalence_plot <- last_plot()


## Finish prevalences

input_dir <- "/cellassign/fitness-scrna/results/new-cut-v1/outputs/differential_expression/output_figures/"

files <- dir(input_dir, pattern = "rds$")
files <- files[!grepl("image", files)]

plots <- lapply(files, function(f) readRDS(glue("{input_dir}/{f}")))

names(plots) <- sapply(strsplit(files, ".", fixed=T), `[`, 1)


upper_grid <- plot_grid(plots$SA609X4XB003083_tsne_by_clone +
                          labs(subtitle = "SA000 X4") +
                          theme(legend.position = "bottom"), 
                        plots$SA609X5XB03230_tsne_by_clone +
                          labs(subtitle = "SA000 X5") +
                          theme(legend.position = "bottom"), 
                        plots$SA609X6XB03404_tsne_by_clone +
                          labs(subtitle = "SA000 X6") +
                          theme(legend.position = "bottom"),
                        prevalence_plot,
                        rel_widths = c(1,1,1,0.6),
                        labels = "AUTO",
                        nrow = 1)

middle_grid <- plot_grid(plots$SA609X4XB003083_tsne_gene_MDM4 +
                           labs(subtitle = "SA000 X4 MDM4") +
                           theme(legend.position = "bottom"),
                         plots$SA609X5XB03230_tsne_gene_MDM4 +
                           labs(subtitle = "SA000 X5 MDM4") +
                           theme(legend.position = "bottom"),
                         plots$SA609X6XB03404_tsne_gene_MDM4 +
                           labs(subtitle = "SA000 X6 MDM4") +
                           theme(legend.position = "bottom"),
                         rel_widths = c(1,1,1),
                         labels = "AUTO",
                         nrow=1
                         #NULL,
                         #rel_widths = c(1, 2.6),
                         #labels = c("E", "xseq stuff")
                         )

bottom_grid <- plot_grid(#plots$SA609X10XB02454_gex_cnv_plot_e_vs_g  +
                         #  labs(subtitle = "SA609 X10"),
  
                         plots$SA609X4XB003083_gex_cnv_plot_b_vs_a  +
                           labs(subtitle = "SA000 X4"),
                         plots$SA609X5XB03230_gex_cnv_plot_b_vs_a  +
                           labs(subtitle = "SA000 X5"),                         
                         plots$SA609X6XB03404_gex_cnv_plot_b_vs_a  +
                           labs(subtitle = "SA000 X6"),
                         labels = c("B", "A"),
                         ncol = 1,
                         nrow = 3)

main_fig <- plot_grid(
  upper_grid,
  middle_grid,
  bottom_grid,
  ncol = 1,
  rel_heights = c(1, 1, 3.2)
)

ggsave("prototype_SA000.png", width = 15, height = 20)
