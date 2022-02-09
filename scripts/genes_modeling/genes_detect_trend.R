## MA started in June 29, 2020
## addapted from plot_expression_distribution_heatmap.R
# /cellassign/fitness-scrna/mirela-drug-res/extranalysis/scripts on krcds, 
# so you can find the inputs relative to this path
library(ggplot2)
library(SingleCellExperiment)
library(wesanderson)
library(RColorBrewer)
library(ggrepel)
### NOW doing for dataset SA1035
#dataset <- "SA1035"
#dataset <- "SA535"
dataset <- "SA609"
# If true, it doesn't matter how untreated changed
untreated_any <- TRUE
version <- "v4"
firstdir <- paste0("../results/",dataset, "-", version)
outdir <- paste0("../results/",dataset, "-", version,"/monotonically_changed_genes")
outdir_tenx <- paste0("../results/",dataset, "-", version,"/tenx")
covthr <- 10
rm_mito_ribo <- "Yes"
all_sce_file <- paste0(outdir_tenx, "/all", dataset, "-", version, "_covthr", covthr, ".rds")
dir.create(firstdir, showWarnings = TRUE)
dir.create(outdir, showWarnings = TRUE)
dir.create(outdir_tenx, showWarnings = TRUE)
# Load the data and keep only the genes that have data in more than "threshold" cells
load_data <- function(libpath, condition) {
  print(paste0("Loading ",libpath))
  sce <- readRDS(libpath)
  rowData(sce) <- rowData(sce)[,c("ID","Symbol")]
  # MA Jul 30, 2020 Adding this below because otherwise Symbol can be duplicated
  print("Rewriting row names")
  rowData(sce)$Symbol <- rownames(sce)
  sce$condition <- condition
  m <- as.matrix(logcounts(sce))
  m[m==0] <- NA
  sce <- sce[rowSums(!is.na(m)) >= covthr,]
  if (rm_mito_ribo == "Yes") {
    sce <- sce[!grepl("^RP[L|S]|^MT-|^FOS|^JUN|^HSP", rowData(sce)$Symbol)]  
  }
  print(dim(sce))
  return(sce)
}
## First create the sceall file
if (file.exists(all_sce_file)) {
  print("Loading all sce file")
  sce <- readRDS(all_sce_file)
} else {
  print("all sce file not found, creating ...")
  dir <- paste0("../../results/",dataset, "-", version,"/outputs/preprocess/sce_twice_scran/")
  ext <- ".rdata"  
  ## NOTE: remove the timepoints, X4 U is now just "U" because for 535 it is a later timepoint.
  if (dataset == "SA1035") {
    sce_u <- load_data(paste0(dir, "SA1035X4XB02879", ext), "U")
    sce_uu <- load_data(paste0(dir, "SA1035X5XB03021", ext), "UU")
    sce_uuu <- load_data(paste0(dir, "SA1035X6XB03216", ext), "UUU")
    sce_uuuu <- load_data(paste0(dir, "SA1035X7XB03502", ext), "UUUU")
    sce_uuuuu <- load_data(paste0(dir, "SA1035X8XB03631", ext), "UUUUU")  
    sce_ut <- load_data(paste0(dir, "SA1035X5XB03015", ext),"UT")
    sce_utt <- load_data(paste0(dir, "SA1035X6XB03211", ext),"UTT")
    sce_uttt <- load_data(paste0(dir, "SA1035X7XB03338", ext),"UTTT")
    sce_utttt <- load_data(paste0(dir, "SA1035X8XB03425", ext),"UTTTT")
    sce_utu <- load_data(paste0(dir, "SA1035X6XB03209", ext), "UTU")
    sce_uttu <- load_data(paste0(dir, "SA1035X7XB03340", ext), "UTTU")
    sce_utttu <- load_data(paste0(dir, "SA1035X8XB03420", ext), "UTTTU")
  }  else if (dataset == "SA535") {
    sce_u <- load_data(paste0(dir, "SA535X5XB02895", ext), "U")
    sce_uu <- load_data(paste0(dir, "SA535X6XB03099", ext), "UU")
    sce_uuu <- load_data(paste0(dir, "SA535X7XB03448", ext), "UUU")
    sce_uuuu <- load_data(paste0(dir, "SA535X8XB03664", ext), "UUUU")
    ## NNote: we don't have X9, making it X8 for now
    sce_uuuuu <- load_data(paste0(dir, "SA535X8XB03664", ext), "UUUUU")  
    sce_ut <- load_data(paste0(dir, "SA535X6XB03101", ext),"UT")
    sce_utt <- load_data(paste0(dir, "SA535X7XB03304", ext),"UTT")
    sce_uttt <- load_data(paste0(dir, "SA535X8XB03431", ext),"UTTT")
    sce_utttt <- load_data(paste0(dir, "SA535X9XB03617", ext),"UTTTT")
    sce_utu <- load_data(paste0(dir, "SA535X7XB03305", ext), "UTU")
    sce_uttu <- load_data(paste0(dir, "SA535X8XB03434", ext), "UTTU")
    sce_utttu <- load_data(paste0(dir, "SA535X9XB03616", ext), "UTTTU")
  } else if (dataset == "SA609") {
    sce_u <- load_data(paste0(dir, "SA609X3XB01584", ext), "U")
    sce_uu <- load_data(paste0(dir, "SA609X4XB03080", ext), "UU")
    sce_uuu <- load_data(paste0(dir, "SA609X5XB03223", ext), "UUU")
    sce_uuuu <- load_data(paste0(dir, "SA609X6XB03447", ext), "UUUU")
    sce_uuuuu <- load_data(paste0(dir, "SA609X7XB03554", ext), "UUUUU")  
    sce_ut <- load_data(paste0(dir, "SA609X4XB003083", ext),"UT")
    sce_utt <- load_data(paste0(dir, "SA609X5XB03230", ext),"UTT")
    sce_uttt <- load_data(paste0(dir, "SA609X6XB03404", ext),"UTTT")
    sce_utttt <- load_data(paste0(dir, "SA609X7XB03505", ext),"UTTTT")
    sce_utu <- load_data(paste0(dir, "SA609X5XB03231", ext), "UTU")
    sce_uttu <- load_data(paste0(dir, "SA609X6XB03401", ext), "UTTU")
    sce_utttu <- load_data(paste0(dir, "SA609X7XB03510", ext), "UTTTU")    
  }
  # I need to get the common genes
  common_genes <- Reduce(intersect, list(
    rowData(sce_u)$Symbol, 
    rowData(sce_uu)$Symbol, 
    rowData(sce_uuu)$Symbol,
    rowData(sce_uuuu)$Symbol,
    rowData(sce_uuuuu)$Symbol,
    rowData(sce_ut)$Symbol, 
    rowData(sce_utt)$Symbol,
    rowData(sce_uttt)$Symbol,
    rowData(sce_utttt)$Symbol,
    rowData(sce_utu)$Symbol,
    rowData(sce_uttu)$Symbol,
    rowData(sce_utttu)$Symbol
  ))
  print(paste0(length(common_genes), " common genes with at least ", covthr, " cells in each library"))
  sce_u <- sce_u[common_genes,]
  sce_uu <- sce_uu[common_genes,]
  sce_uuu <- sce_uuu[common_genes,]      
  sce_uuuu <- sce_uuuu[common_genes,]  
  sce_uuuuu <- sce_uuuuu[common_genes,]  
  sce_ut <- sce_ut[common_genes,]
  sce_utt <- sce_utt[common_genes,]      
  sce_uttt <- sce_uttt[common_genes,]  
  sce_utttt <- sce_utttt[common_genes,]        
  sce_utu <- sce_utu[common_genes,]      
  sce_uttu <- sce_uttu[common_genes,]  
  sce_utttu <- sce_utttu[common_genes,]        
  sce <- do.call("cbind", list(sce_u, sce_uu, sce_uuu, sce_uuuu, sce_uuuuu, 
                               sce_ut, sce_utt, sce_uttt, sce_utttt, 
                               sce_utu, sce_uttu, sce_utttu))           
  saveRDS(sce, all_sce_file)  
}
#### GOT HERE, I have to produce these de files first
### Select the genes that are significantly differentially expressed from X6 to X7 and from X7 to X8 in treated
if (dataset == "SA1035") {
  de1  <- read.csv("../results/SA1035-v4/de/comp6-SA1035X7XB03338-SA1035X6XB03211_logfc_results.csv", header=TRUE)
  de2  <- read.csv("../results/SA1035-v4/de/comp7-SA1035X8XB03425-SA1035X7XB03338_logfc_results.csv", header=TRUE)
} else if (dataset == "SA535") {
  de1  <- read.csv("../results/SA535-v4/de/comp6-SA1035X7XB03338-SA1035X6XB03211_logfc_results.csv", header=TRUE)
  de2  <- read.csv("../results/SA535-v4/de/comp7-SA1035X8XB03425-SA1035X7XB03338_logfc_results.csv", header=TRUE)
} else if (dataset == "SA609") {
  de1  <- read.csv("../../../mirela/pipelines/differential_expression/de/comp31-SA609X6XB03404-R-SA609X5XB03230-R_logfc_results.csv", header=TRUE)
  de2  <- read.csv("../../../mirela/pipelines/differential_expression/de/comp10-SA609X7XB03505-R-SA609X6XB03404-R_logfc_results.csv", header=TRUE)  
}
genes1 <- de1[de1$FDR < 0.01,"gene_symbol"]
genes2 <- de2[de2$FDR < 0.01,"gene_symbol"]
common <- intersect(genes1,genes2)
sceall <- sce[rowData(sce)$Symbol %in% common,]
saveRDS(sceall, file=paste0(outdir, "drug_res_samples_with_treatedDEgenes.rds"))
#genedf <- data.frame(name=rownames(sceall), row.names = rownames(sceall), slopeU=NA, slopeT=NA, 
#                     monotonicT=NA, pathway=NA, line=NA)
genedf <- NULL
for (cond in unique(sceall$condition)) {
  print(cond)
  data <- as.matrix(logcounts(sceall[,sceall$condition==cond]))
  data[data==0] <- NA
  if (grepl("T$", cond) == TRUE) {
    line <- "Rx"
  } else if (grepl("T", cond) == TRUE) {
    line <- "Rx-H"
  } else {
    line <- "Un"
  }  
  timepoint <- sceall[,sceall$condition==cond]$timepoint[1]
  time <- as.numeric(gsub("X","",timepoint))
  genedf <- rbind(genedf, data.frame(name=rownames(data), condition=cond, timepoint=timepoint, 
                                     slope=NA, monotonic=NA, treatment=line, time=time,
                                     mean=rowMeans(data, na.rm=TRUE), sd=rowSds(data, na.rm=TRUE)))
}