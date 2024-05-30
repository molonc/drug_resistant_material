#Summary plots
library(tidyverse)
library(glue)
library(cowplot)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(grid)
library("ggplotify")


source("de-utils.R")
input_dir <- "."

#####################################
# Put all the volcano plots together
#####################################

# These files were obtained by differential-expression-by-clone.Rmd
files <- dir(input_dir, pattern = "_differential_expression.rds$")
i <- 1
p <- vector('list',length(files))
for (file in files) {
  tt <- readRDS(file)
  df_annot <- top_n(tt, 20, -log10(FDR))[1:20,]
  file <- gsub("star", "*", file)
  sp <- strsplit(file,"_")
  p[[i]] <- ggplot(tt, aes(x = logFC, y = -log10(FDR))) +
    geom_point() +
    geom_text_repel(data = df_annot, aes(label = gene_symbol)) +
    labs(x = paste0("log FC, ", sp[[1]][1], " ", sp[[1]][2], " vs. ", sp[[1]][3], " ", sp[[1]][4]))
  i <- i + 1
}
#grid.arrange(p[[1]], p[[2]], ncol=2)
png("SA609_volcano.png", width = 10, height = 15, units="in", res=300, pointsize=4)
#plot_grid(p[[1]],p[[2]],p[[3]],p[[4]], labels="AUTO")
do.call("grid.arrange", c(p, ncol=2))
#ggsave("SA609_volcano.png", width = 15, height = 25)
dev.off()

#####################################
# Put all the track plots together 
#####################################

files <- dir(input_dir, pattern = "_gex_cnv_plot.rds$")
i <- 1
p <- vector('list',length(files))
#p <- vector('list',2)
for (file in files) {
  p[[i]] <- readRDS(file)
  i <- i+1
}
png("SA609_track_plots.png", width = 15, height = 45, units="in", res=100, pointsize=2)
#plot_grid(p[[1]],p[[2]],p[[3]],p[[4]], labels="AUTO")
do.call("grid.arrange", c(p, ncol=1))
#ggsave("SA609_volcano.png", width = 15, height = 25)
dev.off()

#####################################
# Make a correlation plot between the proportions of DLP clones and TENX clones
#####################################
## May 2: this moved to a separate script in fitness_material/scripts/de

# 8 Apr 2020 - updating plot
# dlpdir <- "/cellassign/fitness-scrna/mirela/data/dlp/all_dlp_inputs/"
# #dlpdir <- "/cellassign/fitness-scrna/mirela/data/dlp/pdx_sample_names/"
# 
# 
# #cadir <- "/cellassign/fitness-scrna/mirela/results/new-cut-v1/outputs/align_clones/clonealign_fit_csv/"
# cadir <- "/cellassign/fitness-scrna/mirela/results/new-cut-v2/outputs/align_clones/clonealign_fit_csv/"
# 
# #samples <- c("SA609X10XB02454", "SA609X3XB01584", "SA609X4XB003083", "SA609X4XB03080", "SA609X5XB03230", "SA609X5XB03231", "SA609X7XB03510")
# #samples <- c("SA609X10XB02454", "SA609X3XB01584", "SA609X4XB003083", "SA609X4XB03080", "SA609X5XB03230", "SA609X7XB03510")
# #samples <- c("SA609X4XB003083", "SA609X4XB03080", "SA609X5XB03230", "SA609X5XB03231", "SA609X7XB03510")
# samples <- c("SA906_p57a","SA906_p50b","SA609X10XB02454", "SA609X3XB01584")
# 
# points <- data.frame(prop.dlp=c(), prop.ca=c(), sample=c(), clone=c(), stringsAsFactors = FALSE)
# 
# adjust_clone <- function(ca) {
#   ca <- ca[ca$clone != "unassigned",]  
#   for (i in 1:nrow(ca)) {
#     #print(ca[i,])
#     myclones <- unlist(strsplit(paste0("",ca$clone[i]), "_"))
#     #print(myclones)
#     levels(ca$clone) <- c(levels(ca$clone), myclones)
#     ca$clone[i] <- sample(myclones,1)
#   }
#   return(ca)
# }
# 
# 
# for (sample in samples) {
#   dlp <- read.csv(file=paste0(dlpdir,sample,".csv"), header=TRUE)
#   ca <- read.csv(file=paste0(cadir,sample,".csv"), header=TRUE)
#   ca <- adjust_clone(ca)
#   clusters <- unique(dlp$cluster)
#   for (cluster in clusters) {
#     prop.dlp <- nrow(dlp[dlp$cluster==cluster,])/nrow(dlp)
#     prop.ca <- nrow(ca[ca$clone==cluster,])/nrow(ca)
#     print(cluster)
#     print(prop.dlp)
#     print(prop.ca)
#     points <- rbind(points,data.frame(prop.dlp,prop.ca,sample,cluster))
#   }
# }
# 
# # adding the 1 samples
# #levels(points$sample) <- c(levels(points$sample), "SA609X6XB03404", "SA609X7XB03505", "SA609X6XB03401")
# #points <- rbind(points,c(1.,1.,"SA609X6XB03404","A"))
# #points <- rbind(points,c(1.,1.,"SA609X7XB03505","A"))
# #points <- rbind(points,c(1.,1.,"SA609X6XB03401","A"))
# 
# # rename the sample and clone names
# 
# levels(points$sample)[match("SA609X10XB02454",levels(points$sample))] <- "TNBC PDX X10"
# levels(points$sample)[match("SA609X3XB01584",levels(points$sample))] <- "TNBC PDX X3"
# levels(points$sample)[match("SA906_p57a",levels(points$sample))] <- "p53-/-a X57"
# levels(points$sample)[match("SA906_p50b",levels(points$sample))] <- "p53-/-b X50"
# # put the levels in the desired order
# points$sample <- factor(points$sample, levels <- c("p53-/-a X57","p53-/-b X50","TNBC PDX X3","TNBC PDX X10"))
# 
# # levels(points$sample)[match("SA609X10XB02454",levels(points$sample))] <- "SA609 X10"
# # levels(points$sample)[match("SA609X3XB01584",levels(points$sample))] <- "SA609 X3"
# # levels(points$sample)[match("SA609X4XB003083",levels(points$sample))] <- "SA000 X4"
# # levels(points$sample)[match("SA609X5XB03230",levels(points$sample))] <- "SA000 X5"
# # levels(points$sample)[match("SA609X6XB03404",levels(points$sample))] <- "SA000 X6"
# # levels(points$sample)[match("SA609X7XB03505",levels(points$sample))] <- "SA000 X7"
# # levels(points$sample)[match("SA609X4XB03080",levels(points$sample))] <- "SA001 X4"
# # levels(points$sample)[match("SA609X5XB03231",levels(points$sample))] <- "SA001 X5"
# # levels(points$sample)[match("SA609X6XB03401",levels(points$sample))] <- "SA001 X6"
# # levels(points$sample)[match("SA609X7XB03510",levels(points$sample))] <- "SA001 X7"
# 
# #levels(points$cluster)[match("A",levels(points$cluster))] <- paste("A"^"*")
# 
# points$prop.dlp <- as.numeric(points$prop.dlp)
# points$prop.ca <- as.numeric(points$prop.ca)
# 
# ll <- LETTERS[1:10]; names(ll) <- ll
# 
# #install.packages("wesanderson")
# library(wesanderson)
# 
# #ggplot(points, aes(x=points$prop.dlp,y=points$prop.ca), xlabel="Proportion DLP data", color = factor(points$cluster), shape=factor(points$sample), cex=2.5)
# png(paste0("allsamples_correlation.png"), width = 4, height = 4, units="in", res=300, pointsize=10)
# pal <- wes_palette("Moonrise2", 4, type = "discrete")
# ggplot(points, aes(x = prop.dlp, y = prop.ca, color = sample)) +
# # qplot(points$prop.dlp, points$prop.ca, color = factor(points$sample), cex=6) + 
#   #scale_shape_manual(values=ll) +
#   scale_color_manual(values = pal) + 
#   labs(shape = "Clone", color="Sample") + 
#   labs(x="Clonal proportion DLP", y="Clonal proportion Clonealign") + 
#   guides(size=FALSE, shape=FALSE) + 
#   ggtitle(paste0("Pearson correlation ", round(cor(as.numeric(points$prop.dlp),as.numeric(points$prop.ca)),digits=2))) + 
#   theme(axis.text.x=element_blank(), axis.text.y=element_text()) +
#   cowplot::theme_cowplot() + 
#   theme(legend.position = c(0.1, 0.8)) +
#   geom_point(size = 5) + 
#   #guides(colour = guide_legend(override.aes = list(size=5))) +
#   geom_abline(intercept =0 , slope = 1)
# dev.off()



ggplot(points, aes(x = prop.dlp, y = prop.ca, colour = factor(points$sample))) +
  #  geom_point() +
  #facet_wrap(~ sample, ncol = 1) +
  #geom_text(aes(label = points$clone)) +
  #scale_colour_manual(values = get_cluster_colours(), guide=F) +
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold')) +
  labs(x = "Clonal prevalence DLP", y = "Clonal prevalence CloneAlign")


#####################################
# Plot trajectories for some genes
#####################################

#interesting_genes <- c("TOP2A", "UBE2C", "ARHGAP15", "TCF4", "FOXC1", "NUF2")
# get all the sces I need, qc them, give them proper clone names and then bind them together

interesting_genes <- c()
gene_labels <- c()

# Add housekeeping genes
# interesting_genes <- c(interesting_genes, "GAPDH", "SNRPD3", "PFN1", "RPL36", "RPL8", "HNRNPM")

# Add previously selected interesting genes
interesting_genes <- c(interesting_genes, "TOP2A", "UBE2C", "ARHGAP15", "TCF4", "FOXC1", "NUF2", "MYC", "MDM4", "USP51", "SHISA2", "GPR12", "WASF3", "RUNX1T1", "TCF12")
# genes that didn't work "CAPRIN2", "IGKV1D-8", "CYP39A1", "DAB2IP", "AQP7", "CDK8", "RNF6", "TMEM205", "RFWD3", "PRF1"
gene_labels <-       c(gene_labels,       "TOP2A", "UBE2C", "ARHGAP15", "TCF4", "FOXC1", "NUF2", "MYC CN", "MDM4 CN", "USP51", "SHISA2", "GPR12", "WASF3", "RUNX1T1 CN", "TCF12 CN")


# Add the newest genes from Fatemeh SA000
interesting_genes <- c(interesting_genes, "TMEM132B", "TP53",  "ZNF804A", "ABLIM3", "MROH2B", "PABPC1", "ARC", "USP51")
# Not found "RP11-248J23.6",  "CAPRIN2","GTF3C1", "LINC00661", "PTPRF", "GTF3C2", "IGKV1D-8", "CYP39A1", "STEAP2-AS1", "SCRIB", "DAB2IP", "AQP7", 
gene_labels <- c(gene_labels, "TMEM132B SNV_SA000_A*_B", "TP53 SNV_SA000_A*_B",  "ZNF804A SNV_SA000_A*_B", "ABLIM3 SNV_SA000_A*_B", "MROH2B SNV_SA000_A*_B", "PABPC1 SNV_SA000_A*_B", "ARC SNV_SA000_A*", "USP51 SNV_SA000_A*_B")


interesting_genes <- c("MYC", "TP53", "RUNX1T1", "MDM4", "ARC", "FOXC1", "TCF4", "TOP2A")
gene_labels <- c("MYC (CN)", "TP53 (SNV clones A* and B)", "RUNX1T1 (CN)", "MDM4 (CN)", "ARC (SNV clone A*)", "FOXC1", "TCF4", "TOP2A")

# Add the newest genes from Fatemeh SA001
# interesting_genes <- c(interesting_genes, "HPSE2")
# Not found HPSE2
# gene_labels <- c(gene_labels, "HPSE2 SNV_SA001_A*")

# housekeeping genes
#interesting_genes <- c("GAPDH", "SNRPD3", "PFN1", "RPL36", "RPL8", "HNRNPM", "CKS1B")
# "SRRM1", "HNRNPA2B1", "YWHAB", "GDI2", "CSDE1", "C14ORF2", "NCL", "ARF1", "TARDBP", "STOML2", "RPS5", "THRAP3", "POLR2E", "SRSF3", 

#input_dir <- "../../mirela/results/new-cut-v1/outputs/preprocess/sce_annotated/"
input_dir <- "../../mirela/data/scrna_human_normed/"

prepare_sce <- function(id, label) {
  input_sce <- paste0(input_dir, id, ".rdata")
  sce <- readRDS(input_sce)
  #sce <- fixsce(sce)
  # could select for clone here if I want
  #sce_qc <- sce[,colData(sce)$Barcode %in% mybarcodes2]
  sce_qc <- sce
  sce_qc$clone <- label
  #sce_qc <- sce_qc[!grepl("^RP[L|S]|^MT-|^FOS|^JUN|^HSP", rowData(sce_qc)$Symbol)]
  #sce_qc <- sce_qc[, sce_qc$total_features_by_counts > 3500]
  return(sce_qc)
}

sce000X4 <- prepare_sce("SA609X4XB003083", "SA000X4")
sce000X5 <- prepare_sce("SA609X5XB03230", "SA000X5")
sce000X6 <- prepare_sce("SA609X6XB03404", "SA000X6")
sce000X7 <- prepare_sce("SA609X7XB03505", "SA000X7")

sce001X5 <- prepare_sce("SA609X5XB03231", "SA001X5")
sce001X6 <- prepare_sce("SA609X6XB03401", "SA001X6")
sce001X7 <- prepare_sce("SA609X7XB03510", "SA001X7")

sce609X3 <- prepare_sce("SA609X3XB01584", "SA609X3")
sce609X4 <- prepare_sce("SA609X4XB03080", "SA609X4")
sce609X5 <- prepare_sce("SA609X5XB03223", "SA609X5")
sce609X6 <- prepare_sce("SA609X6XB03447", "SA609X6")
sce609X7 <- prepare_sce("SA609X7XB03554", "SA609X7")
#sce10 <- prepare_sce("SA609X10XB02454", "SA609X10")

rowData(sce000X4) <- rowData(sce000X5) <- rowData(sce000X6) <- rowData(sce000X7) <- NULL
rowData(sce001X5) <- rowData(sce001X6) <- rowData(sce001X7) <- NULL
rowData(sce609X3) <- rowData(sce609X4) <- rowData(sce609X5) <- rowData(sce001X6) <- rowData(sce001X7) <- NULL

sce000 <- cbind(sce000X4, sce000X5, sce000X6, sce000X7)
sce001 <- cbind(sce001X5, sce001X6, sce001X7)
sce609 <- cbind(sce609X3, sce609X4, sce609X5, sce609X6, sce609X7)

#png(paste0("SA609_gene_plots_violins.png"), width = 6, height = 10, units="in", res=150, pointsize=6)
#sce <- cbind(sce1, sce2, sce3, sce4, sce5, sce6, sce7, sce8)
#plotExpression(sce, x = "clone", features = interesting_genes, exprs_values = "counts", ncol=1, show_median=TRUE)
#dev.off()


# TO add error bars, see example here
# http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
#png(paste0("SA609_gene_plots_lines.png"), width = 7, height = 15, units="in", res=150, pointsize=6)
#n_col <- 3
png(paste0("SA609_gene_plots_lines.png"), width = 6, height = 8, units="in", res=150, pointsize=6)
n_col <- 2
par(oma = c(4, 1, 1, 1))
par( mfrow = c( ceiling(length(interesting_genes)/n_col), n_col ), oma=c(4,4,0,0) )
i <- 1
for (gene in interesting_genes) {
  y1 <- c(mean(logcounts(sce000X4)[gene,]), mean(logcounts(sce000X5)[gene,]), mean(logcounts(sce000X6)[gene,]), mean(logcounts(sce000X7)[gene,]))
  y2 <- c(mean(logcounts(sce001X5)[gene,]), mean(logcounts(sce001X6)[gene,]), mean(logcounts(sce001X7)[gene,]))
  y3 <- c(mean(logcounts(sce609X3)[gene,]), mean(logcounts(sce609X4)[gene,]), mean(logcounts(sce609X5)[gene,]), mean(logcounts(sce609X6)[gene,]), mean(logcounts(sce609X7)[gene,]))
  xlabels <- c("X3", "X4", "X5", "X6", "X7")
  #plot(c(2:5),y1,col="coral3",type="o",main=gene_labels[i], xaxt="n",lwd=4, ylab="Mean logcounts", xlab="", xlim=c(1,5), ylim=c(min(y1,y2,y3),max(y1,y2,y3)), cex.lab=3.5, cex.axis=2.5, cex.main=3.5, cex.sub=3.5)
  plot(c(2:5),y1,col="coral3",type="o",main=gene_labels[i], xaxt="n",lwd=4, ylab="", xlab="", xlim=c(1,5), ylim=c(min(y1,y2,y3),max(y1,y2,y3)), cex.lab=3.5, cex.axis=3.0, cex.main=3.5, cex.sub=3.5)
  lines(c(3:5),y2,col="cyan4",type="o",lwd=4)
  lines(c(1:5),y3,col="darkolivegreen4",type="o",lwd=4, lty=2)
  axis(1, at=1:5, labels=c("X3", "X4", "X5", "X6", "X7"),  cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5)
  grid()
  print(gene)
  print(y1)
  print(y2)
  print(y3)
  i <- i + 1
}
# overlay another empty plot over so I can do a common legend
# http://dr-k-lo.blogspot.com/2014/03/the-simplest-way-to-plot-legend-outside.html
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt="n")
legend("bottom", c("Rx", "Rx-vacation", "Untreated"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", lwd=5, col = c("coral3","cyan4","darkolivegreen4"), cex = 4)
#mtext("Timepoint",side=1,line=0,outer=TRUE) #,cex=1.3)
#mtext("Expression (mean counts)",side=2,line=0,outer=TRUE) #,cex=1.3,las=0)
dev.off()











# doesn't work 
png(paste0("SA609_gene_plots_lines.png"), width = 7, height = 4, units="in", res=150, pointsize=6)
p <- vector('list',length(interesting_genes))
i <- 1
for (gene in interesting_genes) {
  y1 <- c(mean(counts(sce1)[gene,]), mean(counts(sce2)[gene,]), mean(counts(sce3)[gene,]), mean(counts(sce4)[gene,]))
  y2 <- c(mean(counts(sce5)[gene,]), mean(counts(sce6)[gene,]), mean(counts(sce7)[gene,]), mean(counts(sce8)[gene,]))
  xlabels <- c("X4", "X5", "X6", "X7")
  p[[i]] <- qplot(x,y1,col="coral3",type="o",main=gene, xaxt="n",lwd=4, cex=1.5, xlab="", ylab="", ylim=c(min(y1,y2),max(y1,y2)))
  lines(x,y2,col="darkcyan",type="o",lwd=4)
  axis(1, at=1:4, labels=c("X4", "X5", "X6", "X7"))
  grid()
  i <- i + 1
}
grid.arrange(p, ncol=3)
dev.off()

#set.seed(123L)
#sce_qc <- runPCA(sce_qc, ncomponents = 20)
#sce_qc <- runTSNE(sce_qc)


p000 <- plotExpression(sce000, x = "clone", features = interesting_genes[1], exprs_values = "counts")
p001 <- plotExpression(sce001, x = "clone", features = interesting_genes[1], exprs_values = "counts")
#png(paste0("SA609_gene_plots_",interesting_genes[1],".png"), width = 15, height = 15, units="in", res=100, pointsize=2)
png(paste0("SA609_gene_plots_",interesting_genes[1],".png"))
par( mfrow = c( 2, 1 ) )
#plot(p000)
#plot(p001)
plotExpression(sce000, x = "clone", features = interesting_genes[1], exprs_values = "counts")
plotExpression(sce001, x = "clone", features = interesting_genes[1], exprs_values = "counts")
dev.off()

p <- vector('list',2)
p[[1]] <- p000
p[[2]] <- p001
png(paste0("SA609_gene_plots_",interesting_genes[1],".png"))
par( mfrow = c( 2, 1 ) )
do.call("grid.arrange", c(p, ncol=1))
dev.off()
#####################################
# Plot tSNE plots by gene
#####################################
genes_to_plot <- c("MYC", "TMEM132B", "CAPRIN2", "TP53", "IGKV1D-8", "CYP39A1", "ARC", "DAB2IP", "AQP7", "USP51", "CDK8", "SHISA2", "RNF6", "GPR12", "WASF3", "TMEM205")
files <- dir(input_dir, pattern = "rds$")
files <- files[!grepl("image", files)]

plots <- lapply(files, function(f) readRDS(glue("{input_dir}/{f}")))

names(plots) <- sapply(strsplit(files, ".", fixed=T), `[`, 1)

ggsave("summary_tsne_genes.png", width = 15, height = 20)