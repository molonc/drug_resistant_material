library(ggplot2)
library(ggrepel)

d <- "/cellassign/fitness-scrna/results/new-cut-v1/outputs/differential_expression/output_figures/"

fig_input <- dir(d, pattern = 'cnv_plot')
fig_output <- gsub(".rds", ".png", fig_input)

# adding these manually in
descriptions <- c(
  "SA609X10XB02454",
  "SA609X10XB02454",
  "SA906_p50b",
  "SA906_p50b",
  "SA906_p57a"
)

figs <- lapply(fig_input, function(f) readRDS(file.path(d, f)))

for(i in seq_along(figs)) {
  f <- figs[[i]] +
    theme(plot.title = element_text()) +
    ggtitle(descriptions[i])
  #print(f)
  ggsave(file.path("gex_plot_output_tmp", fig_output[i]), plot=f, width=14, height = 7)
}
