#file <- "../../materials/comparisons/scrande_SA535_4_SA535_UUT_A_G_UUU_G_logfc_results.csv"
#file <- "../../materials/comparisons/scrande_SA609_10_SA609_UUUU_H_UUUU_C_logfc_results.csv"
file <- "../../materials/comparisons/scrande_SA535_61_SA535_UUT_A_UUU_G_logfc_results.csv"
data <- read.csv(file)
data <- data[data$FDR < 0.01 & (data$logFC > 0.5 | data$logFC < -0.5),]

get_cnv2 <- function(file) {
  cnv <- read.csv(file=file, header=TRUE)
  row.names(cnv) <- cnv$ensembl_gene_id
  # add the symbol gene name
  symbs <- read.csv("Symbol_ensembl.csv")
  row.names(symbs) <- symbs$Ensembl
  cnv$Symbol <- symbs[rownames(cnv),]$Symbol
  return(cnv)
}

#pt <- get_cnv2(file="../../materials/dlp_cnv/SA609_cnv_mat.csv.gz")
pt <- get_cnv2(file="../../materials/dlp_cnv/SA535_cnv_mat.csv")

#pt <- get_cnv2(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/mapped_wholedata_SA535.csv.gz")

pt$ct <- ifelse(pt[["A"]]==pt[["G"]], "trans", "cis")

## keep only cis
pt <- pt[pt$ct=="cis",]

m <- merge(data,pt, by.x="gene_symbol", by.y="Symbol")
m <- m[!is.na(m$ct),]
m$CNG <- m$A-m$G


ntrans <- nrow(pt[pt$A == pt$G,])
ncis <- nrow(pt[pt$A != pt$G,])