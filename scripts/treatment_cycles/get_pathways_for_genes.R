gmt <- data.frame(read.csv("h.all.v7.0.symbols.gmt",sep='\t',header = FALSE))
gmt[,2] <- NULL

library(reshape2)

# Specify id.vars: the variables to keep but not split apart on
pgenes <- melt(gmt, id.vars=c("V1"))
pgenes$variable <- NULL
colnames(pgenes) <- c("pathway","gene_symbol")

mygenes <- read.csv("SuppTable8_pseudotime_genes_modules_Pt4_Pt5_Pt6.csv")
mygenes <- mygenes[mygenes$patient=="Pt4" & mygenes$gene_type_module=="Module6",]

table(pgenes[pgenes$gene_symbol %in% mygenes$gene_symbol,]$pathway)