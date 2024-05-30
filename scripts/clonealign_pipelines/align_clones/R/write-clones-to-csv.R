
library(clonealign)
library(readr)

ca <- readRDS(snakemake@input[['clonealign_fit']])

ca <- clonealign:::recompute_clone_assignment(ca, 0.5)

df <- ca$clone_fit
df$clone <- ca$clone

write_csv(df, snakemake@output[['clonealign_fit']])
