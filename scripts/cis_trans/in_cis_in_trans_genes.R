library(DESeq2)
# In-cis, in-trans genes
dds <- makeExampleDESeqDataSet()
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

View(head(dds$condition))
class(dds)

class(dds$condition)
head(dds@assays$counts)


# see vignette for suggestions on generating
# count tables from RNA-Seq data
cnts <- matrix(rnbinom(n=4000, mu=100, size=1/0.5), ncol=40)
dim(cnts)
ncol(cnts) 
cnts[1:3,1:3]
cond <- factor(rep(c("A","B","C","D"), c(10,15,10,5)))
cond <- model.matrix(~0+cond)
# https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
nrow(cond)
colData_mtx <- data.frame(cond)
nrow(colData_mtx)
colnames(colData_mtx)
length(cond)
# object construction
 
dds <- DESeqDataSetFromMatrix(cnts, colData_mtx, cond)

dds <- DESeqDataSetFromMatrix(cnts, colData_mtx, ~ cond)

# countData <- cnts
# colData <- colData_mtx
# 

# standard analysis
dds <- DESeq(dds)
res <- results(dds,contrast=list(c("cond1"),c("cond2","cond3")))
res <- results(dds,contrast=list(c("cond2"),c("cond1","cond3")))
res

# moderated log2 fold changes
# resultsNames(dds)
# resLFC <- lfcShrink(dds, coef=2, type="apeglm")

# an alternate analysis: likelihood ratio test
ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)

dds <- makeExampleDESeqDataSet(n=100,m=6)
dds$condition <- factor( c( "A","A","B","B","C","C") )
dds <- DESeq(dds)
res = results(dds, contrast=c("condition","C","A"))
res <- res[order(res$padj),]
kable(res[1:5,-(3:4)])



cnts <- matrix(rnbinom(n=10000, mu=100, size=1/0.5), ncol=200)
data <- as.data.frame(cnts)
dim(data)
colnames(data) <- paste0('C',rep(c(1:ncol(cnts)),1))
rownames(data) <- paste0('G',rep(c(1:nrow(cnts)),1))
threshold <- 0.05
outputPrefix <- '/home/htran/storage/datasets/drug_resistance/rna_results/Mirela_output/gene_regression/deseq2_'
dir.create(outputPrefix)
# commandArgs()[1] = HTSeq count table
# commandArgs()[2] = prefix of output files
# commandArgs()[3] = num. of replicates
# commandArgs()[4] = threshold


#----- input main arguments -----#
args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1]
outputPrefix <- args[2]
NUM <- as.numeric(args[3])
threshold <- as.numeric(args[4])
data <- read.table(inputFile, header = T, row.names = 1)
outputNormalizedCount <- paste(outputPrefix, ".normalized_count", sep = "")
outputNormalizedMean <- paste(outputPrefix, ".normalized_mean", sep = "")
outputTestResults <- paste(outputPrefix, ".results", sep = "")
outputSignResults <- paste(outputPrefix, ".significant", sep = "")
#--------------------------------#


#----- input HTSeq count data -----#
nbreplicates <- c(50,100,30,20)
obs_clones <- c("A","B","C","D")
names(nbreplicates) <- obs_clones
cond <- factor(rep(obs_clones, nbreplicates))
colData = data.frame(condition = cond, type = "single-end", row.names = colnames(data))

pair_clones <- combinations(length(obs_clones),2,obs_clones)

# colData = data.frame(condition = rep(c("A", "B", "C", "D"), 
# each = NUM), type = "single-end", row.names = colnames(data))

# check whether this function do normalization
dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)
#----------------------------------#


#----- perform statistical test -----#
library(parallel)
cores_use <- 8
nbcores <- detectCores()
if(cores_use > nbcores){
  cores_use <- nbcores
}
p <- MulticoreParam(cores_use)
dds <- DESeq(dds, fitType = "local", parallel = TRUE, BPPARAM = p)  # TO DO: parallel computation
# AB <- results(dds, contrast = c("condition", "B", "A"))
# AC <- results(dds, contrast = c("condition", "C", "A"))
# AD <- results(dds, contrast = c("condition", "D", "A"))
# BC <- results(dds, contrast = c("condition", "C", "B"))
# BD <- results(dds, contrast = c("condition", "D", "B"))
# CD <- results(dds, contrast = c("condition", "D", "C"))

DE_results <- list()
for(i in rep(1:nrow(pair_clones),1)){
  lb <- paste0(pair_clones[i,1],pair_clones[i,2])
  res <- results(dds, contrast = c("condition", pair_clones[i,2], pair_clones[i,1]), 
                 parallel = TRUE, BPPARAM = p)
  rownames(res) <- rownames(data)
  
  res_trans <- data.frame(transform(res, sample1 = pair_clones[i,1], sample2 = pair_clones[i,2]), 
                   row.names = paste0(rownames(res), "_", lb))
  DE_results[[i]] <- res_trans
}
length(DE_results)
dim(res_trans)
all <- do.call(rbind, DE_results)
print(dim(all))
# FDR test, BH adjusted test
q <- p.adjust(all$pval, method = "BH")
result <- transform(all, qval = q)
print(summary(result$qval))
#------------------------------------#
# rownames(AB) <- rownames(data)
# rownames(AC) <- rownames(data)
# rownames(AD) <- rownames(data)
# rownames(BC) <- rownames(data)
# rownames(BD) <- rownames(data)
# rownames(CD) <- rownames(data)

#----- calc q-values -----#
# ab <- data.frame(transform(AB, sample1 = "A", sample2 = "B"))
# dim(ac)


# ab <- data.frame(transform(AB, sample1 = "A", sample2 = "B"), row.names = paste(rownames(AB), "_AB", sep = ""))
# ac <- data.frame(transform(AC, sample1 = "A", sample2 = "C"), row.names = paste(rownames(AC), "_AC", sep = ""))
# ad <- data.frame(transform(AD, sample1 = "A", sample2 = "D"), row.names = paste(rownames(AD), "_AD", sep = ""))
# bc <- data.frame(transform(BC, sample1 = "B", sample2 = "C"), row.names = paste(rownames(BC), "_BC", sep = ""))
# bd <- data.frame(transform(BD, sample1 = "B", sample2 = "D"), row.names = paste(rownames(BD), "_BD", sep = ""))
# cd <- data.frame(transform(CD, sample1 = "C", sample2 = "D"), row.names = paste(rownames(CD), "_CD", sep = ""))
# all <- rbind(ab, ac, ad, bc, bd, cd)
# q <- p.adjust(all$pval, method = "BH")
# result <- transform(all, qval = q)
#-------------------------#




#----- calc mean expression levels for each condition -----#
normalizedCount <- as.data.frame(counts(dds, normalized = TRUE))
dim(normalizedCount)
# get obs genes only: normalizedCount <- normalizedCount[obs_genes,]

# gene_cn_summarized <- mclapply2(rownames(normalizedCount), function(g) {
#   dplyr::filter(gene_cn, ensembl_gene_id == g)  %>%
#     dplyr::group_by(ensembl_gene_id, cell_id, clone) %>%
#     dplyr::summarise(percent_overlap=sum(percent_overlap),
#                      n_parts=n(), .groups = 'drop') %>% 
#     ungroup()
# }, mc.cores = cores_use) %>% bind_rows()


ls_val <- list()
tmp <- rep(0, times = nrow(normalizedCount))
for(c in obs_clones){
  ls_val[[c]] <- tmp
}
normalizedMean <- data.frame(matrix(unlist(ls_val), ncol=length(ls_val), byrow=F),
                             row.names = rownames(normalizedCount))
colnames(normalizedMean) <- names(ls_val)
dim(normalizedMean)

cells_clones <- list()
for(cn in obs_clones){
  cells_use <- rownames(colData[colData$condition==cn,])
  cells_clones[[cn]] <- cells_use
}

library(dplyr)

gene_mean <- mclapply(rownames(normalizedMean), function(g) {
  g_tmp <- normalizedMean[g,]
  for(cn in obs_clones){
    SUM = 0
    for(c in cells_clones[[cn]]){
      SUM = SUM + normalizedCount[g, c]
    }
    g_tmp[g,cn] <- SUM / length(cells_clones[[cn]])
  }
  g_tmp
}, mc.cores = cores_use) %>% bind_rows()

# print(dim(gene_mean))
# gene_mean[1:3,]
# normalizedMean[1:3,]
# normalizedMean <- data.frame(A = rep(0, times = nrow(normalizedCount)),
#                              B = rep(0, times = nrow(normalizedCount)),
#                              C = rep(0, times = nrow(normalizedCount)),
#                              D = rep(0, times = nrow(normalizedCount)),
#                              row.names = rownames(normalizedCount))
# 
# for(i in 1 : nrow(normalizedCount)){
#   for(j in 0 : (length(nbreplicates)-1)){
#     SUM = 0
#     for(k in 1 : nbreplicates[j+1]){
#       SUM = SUM + normalizedCount[i, (j * nbreplicates[j+1] + k)]
#     }
#     normalizedMean[i, (j + 1)] <- SUM / nbreplicates[j+1]
#   }
# }
# #----------------------------------------------------------#
# normalizedMean[1:2,]

#----- output results -----#
out <- file(outputNormalizedCount, "w")
write.table(normalizedCount, out, col.names = F, row.names = T, append = F, quote = F, sep = "\t")
close(out)

out <- file(outputNormalizedMean, "w")
write.table(normalizedMean, out, col.names = F, row.names = T, append = F, quote = F, sep = "\t")
close(out)

out <- file(outputTestResults, "w")
write.table(result, out, col.names = T, row.names = T, quote = F, sep = "\t")
close(out)

out <- file(outputSignResults, "w")
write.table(subset(result, qval < threshold), out, col.names = T, row.names = T, quote = F, sep = "\t")
close(out)
#--------------------------#

t <- subset(result, qval < threshold)
t

summary(result$qval)


# TO DO

# + Get data for each series, with clone label, treatment status, timepoint 
# Get variable genes, filtering
# + Test algo 
# + Get significant genes list
# + save to file
# + Get segment data and remove bad genes only, convert to ensemble genes id
# + get median genotype and variant genes 
# + get in-cis, in-trans genes







