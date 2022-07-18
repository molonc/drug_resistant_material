# Summarize how many genes changed over treatment and how many reversed

stepFC <- 0.25
FDRthr <- 0.01
maxFC <- 1.5

source("ggplot_utils.R")

clone <- "aware"
#clone <- "unaware"


### get the pathway genes
gmt <- data.frame(read.csv("h.all.v7.0.symbols.gmt",sep='\t',header = FALSE))
gmt[,2] <- NULL

library(reshape2)

# Specify id.vars: the variables to keep but not split apart on
pgenes <- melt(gmt, id.vars=c("V1"))
pgenes$variable <- NULL
colnames(pgenes) <- c("pathway","gene_symbol")
# 4385 unique genes

add_summary_data <- function (df, label, series, geneset=NULL, file, ignore_cis_trans=FALSE) {
  t1 <- read.csv(file)
  
  t1 <- t1[t1$FDR <= FDRthr,]
  
  if(!is.null(geneset)) {
    t1 <- t1[t1$gene_symbol %in% geneset,]
  }

  
  #maxFC <- max(abs(t1$logFC))
  
  for (logFC in seq(0,maxFC,stepFC)) {
    n <- length(t1[abs(t1$logFC) >= logFC,]$logFC)
    df <- rbind(df,data.frame(label=label, type="all", series=series, number=n, logFCthr=logFC))
  }
  
  if (!ignore_cis_trans) {
    ## cis
    for (logFC in seq(0,maxFC,stepFC)) {
      n <- length(t1[abs(t1$logFC) >= logFC & t1$cnchange!=0,]$logFC)
      df <- rbind(df,data.frame(label=label, type="in cis", series=series, number=n, logFCthr=logFC))
    }
    
    ## trans
    for (logFC in seq(0,maxFC,stepFC)) {
      n <- length(t1[abs(t1$logFC) >= logFC & t1$cnchange==0,]$logFC)
      df <- rbind(df,data.frame(label=label, type="in trans", series=series, number=n, logFCthr=logFC))
    }  
  }

  return(df)
}


## NOTE: these are clone aware, perhaps better to use no clones

df <- NULL
### to get the files from Hoa
### For SA609 and SA535 I am using clone aware here
### I used to have geneset=pgenes$gene_symbol, but now we want all the genes, not just the ones from the pathways
df <- add_summary_data(df, label="X4 Rx:A UnRx:H", series="Pt4", file="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_1_SA609_UT_R_UU_H_logfc_results2.csv")
df <- add_summary_data(df, label="X5 Rx:A UnRx:H", series="Pt4", file="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_2_SA609_UTT_R_UUU_H_logfc_results2.csv")
#df <- add_summary_data(df, label="2RxH", series="SA609", geneset=pgenes$gene_symbol, file="../results/SA609-v6/comps2/scrande_SA609_22_SA609_UTU_B_UUU_H_logfc_results2.csv")
df <- add_summary_data(df, label="X6 Rx:A UnRx:H", series="Pt4", file="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_3_SA609_UTTT_R_UUUU_H_logfc_results2.csv")
#df <- add_summary_data(df, label="3RxH", series="SA609", geneset=pgenes$gene_symbol, file="../results/SA609-v6/comps2/scrande_SA609_23_SA609_UTTU_R_UUUU_H_logfc_results2.csv")
df <- add_summary_data(df, label="X7 Rx:A UnRx:H", series="Pt4", file="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_4_SA609_UTTTT_R_UUUUU_H_logfc_results2.csv")
#df <- add_summary_data(df, label="4RxH", series="SA609", geneset=pgenes$gene_symbol, file="../results/SA609-v6/comps2/scrande_SA609_24_SA609_UTTTU_R_UUUUU_H_logfc_results2.csv")

#df <- add_summary_data(df, label="1 Rx 1 UnRx", series="SA609", file="../results/SA609-v6/comps2/scrande_SA609_22_SA609_UTU_B_UUU_H_logfc_results2.csv")

# removing next two lines because they they are not completely clone aware
#df <- add_summary_data(df, label="X6 Rx:AG UnRx:G", series="Pt5", file="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_4_SA535_UUT_A_G_UUU_G_logfc_results2.csv")
#df <- add_summary_data(df, label="X7 Rx:AEJ UnRx:G", series="Pt5", file="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_5_SA535_UUTT_A_E_J_UUUU_G_logfc_results2.csv")
df <- add_summary_data(df, label="X8 Rx:A UnRx:G", series="Pt5", file="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_1_SA535_UUTTT_A_UUUUU_G_logfc_results2.csv")
#df <- add_summary_data(df, label="4 Rx", series="SA535", geneset=pgenes$gene_symbol, file="../results/SA535-v7/comps2/scrande_SA535_2_SA535_UUTTTT_A_UUUUUU_G_logfc_results2.csv")
df <- add_summary_data(df, label="X9+1 Rx:A UnRx:G", series="Pt5", file="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_3_SA535_UUTTTTT_A_UUUUUU_G_logfc_results2.csv")

#  SA1035
df <- add_summary_data(df, label="X6 Rx:G UnRx:G", series="Pt6", file="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_3_SA1035_UTT_G_UUU_G_logfc_results2.csv")
df <- add_summary_data(df, label="X7 Rx:H UnRx:E", series="Pt6", file="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_1_SA1035_UTTT_H_UUUU_E_logfc_results2.csv")
df <- add_summary_data(df, label="X8 Rx:H UnRx:E", series="Pt6", file="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_2_SA1035_UTTTT_H_UUUUU_E_logfc_results2.csv")


#df <- add_summary_data(df, label="1Rx", series="Pt6", geneset=pgenes$gene_symbol, file="../results/SA1035-v6/comps2/scran-nocl-SA1035-1_SA1035_UT_UU-SA1035-SA1035X5XB03015-SA1035X5XB03021_logfc_results2.csv")
#df <- add_summary_data(df, label="2Rx", series="Pt6", geneset=pgenes$gene_symbol, file="../results/SA1035-v6/comps2/scran-nocl-SA1035-2_SA1035_UTT_UUU-SA1035-SA1035X6XB03211-SA1035X6XB03216_logfc_results2.csv")
#df <- add_summary_data(df, label="3Rx", series="Pt6", geneset=pgenes$gene_symbol, file="../results/SA1035-v6/comps2/scran-nocl-SA1035-3_SA1035_UTTT_UUUU-SA1035-SA1035X7XB03338-SA1035X7XB03502_logfc_results2.csv")
#df <- add_summary_data(df, label="4Rx", series="Pt6", geneset=pgenes$gene_symbol, file="../results/SA1035-v6/comps2/scran-nocl-SA1035-4_SA1035_UTTTT_UUUUU-SA1035-SA1035X8XB03425-SA1035X8XB03631_logfc_results2.csv")




df$series <- factor(df$series, level=c("Pt4","Pt5","Pt6"))

library(ggplot2)
# Basic line plot with points
ggplot(data=df, aes(x=logFCthr, y=number, group=label)) +
  geom_line(aes(color=label), size=2)+
  geom_point(aes(color=label), size=3) +
  labs(y = "Number of DEGs", x = "|log2FC| threshold") + 
  facet_grid(series ~ type,scales = "free") 



#df$series <- factor(df$series, level=c("SA609","SA535","SA1035"))

library(ggplot2)
# Basic line plot with points
ggplot(data=df, aes(x=logFCthr, y=number, group=label)) +
  geom_line(aes(color=label), size=2)+
  geom_point(aes(color=label), size=3) +
  labs(y = "Number of DEGs", x = "|log2FC| threshold") + 
  facet_grid(series ~ type,scales = "free") 


df_fc <- df[df$logFCthr==0.5,] # & df$series!="SA1035",]
df_fc <- df_fc[df_fc$type != "all",]

#### Looking at holiday reversal

unnormalized_df_fc <- df_fc
### Normalize !!!


### calculate total number in-cis per sample pair, so I can normalize
calc_total_incis <- function (file_T, file_U) {
  data_T <- read.csv(file_T)
  data_U <- read.csv(file_U)
  return (sum(abs(data_T$median_cn - data_U$median_cn) >= 1))
}

pt4_1 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA609/SA609X4XB03084_cnv.csv",   ### TO CHANGE THIS
                         file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA609/SA609X4XB03080_cnv.csv")

pt4_2 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA609/SA609X5XB03226_cnv.csv",   ## TO CHANGE
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA609/SA609X5XB03223_cnv.csv")

pt4_3 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA609/SA609X6XB03387_cnv.csv",   ## TO CHANGE
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA609/SA609X6XB03447_cnv.csv")

pt4_4 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA609/SA609X7XB03573_cnv.csv",   ## TO CHANGE
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA609/SA609X7XB03554_cnv.csv")


pt5_1 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA535/SA535X6XB03101_cnv.csv",  
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA535/SA535X6XB03099_cnv.csv")

pt5_2 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA535/SA535X7XB03304_cnv.csv",   
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA535/SA535X7XB03448_cnv.csv")

pt5_3 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA535/SA535X8XB03431_cnv.csv",   
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA535/SA535X8XB03664_cnv.csv")

pt5_5 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA535/SA535X10XB03696_cnv.csv",   
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA535/SA535X9XB03776_cnv.csv")


pt6_1 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA1035/SA1035X5XB03015_cnv.csv",  
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA1035/SA1035X5XB03021_cnv.csv")

pt6_2 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA1035/SA1035X6XB03211_cnv.csv",  
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA1035/SA1035X6XB03216_cnv.csv")

pt6_3 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA1035/SA1035X7XB03338_cnv.csv",  
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA1035/SA1035X7XB03502_cnv.csv")

pt6_4 <- calc_total_incis (file_T="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA1035/SA1035X8XB03425_cnv.csv",  
                           file_U="../../../drug_resistance_local/extranalysis/scripts/segment_cn/SA1035/SA1035X8XB03631_cnv.csv")

df_fc[df_fc$series=="Pt4" & df_fc$label=="1Rx","total_incis"] <- pt4_1
df_fc[df_fc$series=="Pt4" & df_fc$label=="2Rx","total_incis"] <- pt4_2
df_fc[df_fc$series=="Pt4" & df_fc$label=="3Rx","total_incis"] <- pt4_3
df_fc[df_fc$series=="Pt4" & df_fc$label=="4Rx","total_incis"] <- pt4_4

df_fc[df_fc$series=="Pt5" & df_fc$label=="1Rx","total_incis"] <- pt5_1
df_fc[df_fc$series=="Pt5" & df_fc$label=="2Rx","total_incis"] <- pt5_2
df_fc[df_fc$series=="Pt5" & df_fc$label=="3Rx","total_incis"] <- pt5_3
df_fc[df_fc$series=="Pt5" & df_fc$label=="5Rx","total_incis"] <- pt5_5

df_fc[df_fc$series=="Pt6" & df_fc$label=="1Rx","total_incis"] <- pt6_1
df_fc[df_fc$series=="Pt6" & df_fc$label=="2Rx","total_incis"] <- pt6_2  #### TO CHANGE pt6_2 is very small compared with others
df_fc[df_fc$series=="Pt6" & df_fc$label=="3Rx","total_incis"] <- pt6_3
df_fc[df_fc$series=="Pt6" & df_fc$label=="4Rx","total_incis"] <- pt6_4


df_fc$normalized <- df_fc$number/df_fc$total_incis



ggplot(data=df_fc, aes(x=label, y=number, fill=type)) +
  geom_bar(stat="identity") +  
  labs(y = "# DEGs", x = "Sample") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(~ series, scales="free_x", space="free")


#ggsave(paste0("deg_genes_clone", clone,".pdf"), width=7.5, heigh=3.5, useDingbats=FALSE)

ggplot(data=df_fc, aes(x=label, y=normalized, fill=type)) +
  geom_bar(stat="identity") +  
  labs(y = "Normalized # DEGs", x = "Sample") 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(~ series, scales="free_x")


#ggsave("../../figures/summary_figure/deg_genes_norm.pdf", width=7.5, heigh=2, useDingbats=FALSE)


ggplot(data=df_fc, aes(x=label, y=total_incis)) +
  geom_bar(stat="identity") +  
  labs(y = "Total in-cis", x = "Sample") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(~ series, scales="free_x")

ggsave("../../figures/summary_figure/total_incis.pdf", width=7, heigh=1.8, useDingbats=FALSE)

##########################3
### Normalize by the total number of in-cis
## I don't have the data for it


# not used
find_nonDE_genes <- function(file1, file2) {
  hol <- read.csv(file1)
  tr <- read.csv(file2)
  
  tr <- tr[tr$FDR <= 0.01 & abs(tr$logFC)>=0.5 & tr$gene_symbol %in% pgenes$gene_symbol,]
  hol1 <- hol[hol$FDR > 0.01 & abs(hol$logFC)<0.5 & hol$gene_symbol %in% pgenes$gene_symbol,]
  reversed <- intersect(tr$gene_symbol, hol1$gene_symbol)
  return(reversed)
}

# not used
find_DE_genes <- function(file1, file2) {
  hol <- read.csv(file1)
  tr <- read.csv(file2)
  
  tr <- tr[tr$FDR <= 0.01 & abs(tr$logFC)>=0.5 & tr$gene_symbol %in% pgenes$gene_symbol,]
  hol2 <- hol[hol$FDR <= 0.01 & abs(hol$logFC)>=0.5 & hol$gene_symbol %in% pgenes$gene_symbol,]
  fixed <- intersect(tr$gene_symbol, hol2$gene_symbol)
  return(fixed)
}

getDE_in_pathways <- function(de, fdr=0.01, logFC=0.5) {
  return(de[de$FDR <= fdr & abs(de$logFC)>=logFC & de$gene_symbol %in% pgenes$gene_symbol,])
}

getnonDE_in_pathways <- function(de, fdr=0.1) {
  return(de[de$FDR > fdr & de$gene_symbol %in% pgenes$gene_symbol,])
}

getDE_all <- function(de, fdr=0.01, logFC=0.5) {
  return(de[de$FDR <= fdr & abs(de$logFC)>=logFC,])
}

getnonDE_all <- function(de, fdr=0.01) {
  return(de[de$FDR > fdr,])
}


############################
############################
# function to prepare the cis and trans

get_cnv2 <- function(file) {
  cnv <- read.csv(file=file, header=TRUE)
  # add the symbol gene name
  symbs <- read.csv("Symbol_ensembl.csv")
  row.names(symbs) <- symbs$Ensembl
  cnv$Symbol <- symbs[rownames(cnv),]$Symbol
  return(cnv)
}


#####################
#####################
### Fig 6b  reversed holday genes
reversed_genes_logfc <- 0.5


#Pt4_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA609_cnv_mat.csv.gz")
#Pt5_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA535_cnv_mat.csv.gz")
#Pt6_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA1035_cnv_mat.csv.gz")


Pt4_cnv_whole <- get_cnv2(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/mapped_wholedata_SA609.csv.gz")
Pt5_cnv_whole <- get_cnv2(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/mapped_wholedata_SA535.csv.gz")
Pt6_cnv_whole <- get_cnv2(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/mapped_wholedata_SA1035.csv.gz")



find_reversed_holiday_genes <- function(fileTvsU, fileTvsH, fileHvsU, pt, cloneT, cloneU, filename=NULL) {
  deTvsU <- read.csv(fileTvsU)  ## DE
  deTvsH <- read.csv(fileTvsH)  ## DE
  deHvsU <- read.csv(fileHvsU)  ## not DE 
  
  
  ## get DE of Rx vs RxH
  deTvsH <- getDE_all(deTvsH, fdr=0.01, logFC=reversed_genes_logfc)
  
  # also use deTvsU
  #deTvsU <- getDE(deTvsU)
  
  # or include all genes from TvsU
  deTvsU <- getDE_all(deTvsU, fdr=1, logFC=0)

  #nondeHvsU <- getnonDE(deHvsU)
  #reversed <- intersect(intersect(deTvsU$gene_symbol, deTvsH$gene_symbol), nondeHvsU$gene_symbol)
  
  common <- intersect(deTvsU$gene_symbol, deTvsH$gene_symbol)
  deTvsH <- deTvsH[deTvsH$gene_symbol %in% common,]
  deTvsU <- deTvsU[deTvsU$gene_symbol %in% common,]
  m <- merge(deTvsH, deTvsU, by = "gene_symbol")
  
  # adding the new cis/trans
  pt$ct <- ifelse(pt[[cloneT]]==pt[[cloneU]], "trans", "cis")
  m <- merge(m,pt, by.x="gene_symbol", by.y="gene_ymbol")
  m <- m[!is.na(m$ct),]
  #m[is.na(m$ct),]$ct <- "unkn"
  
  #return(m)
  same_direction_cis <- sum (m$logFC.x * m$logFC.y > 0 & m$ct=="cis")
  same_direction_trans <- sum (m$logFC.x * m$logFC.y > 0 & m$ct=="trans")
  same_direction_unkn <- sum (m$logFC.x * m$logFC.y > 0 & m$ct=="unkn")
  opposite_direction_cis <- sum (m$logFC.x * m$logFC.y < 0 & m$ct=="cis")
  opposite_direction_trans <- sum (m$logFC.x * m$logFC.y < 0 & m$ct=="trans")
  opposite_direction_unkn <- sum (m$logFC.x * m$logFC.y < 0 & m$ct=="unkn")
  # to add cis and trans
  if (!is.null(filename)) {
    library(dplyr)
    m <- select(m,gene_symbol,ensembl_gene_id.x,logFC.x,FDR.x,logFC.y,FDR.y,ct)
    colnames(m) <- c("gene_symbol","ensembl_gene_id","TvsH.logFC","TvsH.FDR","TvsU.logFC","TvsU.FDR","cistrans")
    m$direction <- ifelse (m$TvsH.logFC * m$TvsU.logFC > 0, "Towards UnRx", "Away from UnRx")
    write.csv(m, file=filename, quote=FALSE, row.names=FALSE)
  }
  
  return(list(same_cis=same_direction_cis, same_trans=same_direction_trans, same_unkn=same_direction_unkn, 
              opp_cis=opposite_direction_cis, opp_trans=opposite_direction_trans, opp_unkn=opposite_direction_unkn))
}




if (clone=="aware") {
  
  
  
  
  
  reversed_X5_609 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_2_SA609_UTT_R_UUU_H_logfc_results.csv",
                                             fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_20_SA609_UTT_R_UTU_B_logfc_results.csv",
                                             fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_22_SA609_UTU_B_UUU_H_logfc_results.csv",
                                             pt=Pt4_cnv_whole, cloneT="cn_median_R", cloneU="cn_median_H",
                                             filename="manuscript_files/SA609_X5_reversed_holiday_genes.csv")
  
  
  reversed_X6_609 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_3_SA609_UTTT_R_UUUU_H_logfc_results.csv",
                                             fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_25_SA609_UTTT_R_UTTU_R_logfc_results.csv",
                                             fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_23_SA609_UTTU_R_UUUU_H_logfc_results.csv",
                                             pt=Pt4_cnv_whole, cloneT="cn_median_R", cloneU="cn_median_H",
                                             filename="manuscript_files/SA609_X6_reversed_holiday_genes.csv")
  
  reversed_X7_609 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_4_SA609_UTTTT_R_UUUUU_H_logfc_results.csv",
                                             fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_26_SA609_UTTTT_R_UTTTU_R_logfc_results.csv",
                                             fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_24_SA609_UTTTU_R_UUUUU_H_logfc_results.csv",
                                             pt=Pt4_cnv_whole, cloneT="cn_median_R", cloneU="cn_median_H",
                                             filename="manuscript_files/SA609_X7_reversed_holiday_genes.csv")
  

  ##### SA535 
  reversed_X7_535 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_5_SA535_UUTT_A_E_J_UUUU_G_logfc_results.csv",
                                                 fileTvsH="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_51_SA535_UUTT_A_E_J_UUTU_A_E_J_logfc_results.csv",
                                                 fileHvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_52_SA535_UUTU_A_E_J_UUUU_G_logfc_results.csv",
                                                 pt=Pt5_cnv_whole, cloneT="cn_median_A", cloneU="cn_median_G")   ### TODO: should be AEJ??
  
  reversed_X8_535 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_1_SA535_UUTTT_A_UUUUU_G_logfc_results.csv",
                                                 fileTvsH="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_15_SA535_UUTTT_A_UUTTU_D_logfc_results.csv",
                                                 fileHvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_16_SA535_UUTTU_D_UUUUU_G_logfc_results.csv",
                                                 pt=Pt5_cnv_whole, cloneT="cn_median_A", cloneU="cn_median_G",
                                                 filename="manuscript_files/SA535_X8_reversed_holiday_genes.csv")
  

  reversed_X10_535 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_3_SA535_UUTTTTT_A_UUUUUU_G_logfc_results.csv",
                                                  fileTvsH="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_35_SA535_UUTTTTT_A_UUTTTTU_E_logfc_results.csv",
                                                  fileHvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_36_SA535_UUTTTTU_E_UUUUUU_G_logfc_results.csv",
                                                  pt=Pt5_cnv_whole, cloneT="cn_median_A", cloneU="cn_median_G")
  
  
  
  reversed_X6_1035 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_3_SA1035_UTT_G_UUU_G_logfc_results.csv",
                                                  fileTvsH="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_31_SA1035_UTT_G_UTU_B_logfc_results.csv",
                                                  fileHvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_32_SA1035_UTU_B_UUU_G_logfc_results.csv",
                                                  pt=Pt6_cnv_whole, cloneT="cn_median_G", cloneU="cn_median_G")
  
  reversed_X7_1035 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_1_SA1035_UTTT_H_UUUU_E_logfc_results.csv",
                                                  fileTvsH="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_11_SA1035_UTTT_H_UTTU_G_logfc_results.csv",
                                                  fileHvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_12_SA1035_UTTU_G_UUUU_E_logfc_results.csv",
                                                  pt=Pt6_cnv_whole, cloneT="cn_median_H", cloneU="cn_median_E")
  
  reversed_X8_1035 <- find_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_2_SA1035_UTTTT_H_UUUUU_E_logfc_results.csv",
                                                  fileTvsH="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_21_SA1035_UTTTT_H_UTTTU_G_logfc_results.csv",
                                                  fileHvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_22_SA1035_UTTTU_G_UUUUU_E_logfc_results.csv",
                                                  pt=Pt6_cnv_whole, cloneT="cn_median_H", cloneU="cn_median_E",
                                                  filename="manuscript_files/SA1035_X8_reversed_holiday_genes.csv")
   
  
}
  
df <- NULL
df <- rbind(df, data.frame(series="Pt4", title="X5 Rx:A RxH:B UnRx:H", direction="Towards UnRx cis", number=reversed_X5_609$same_cis))
df <- rbind(df, data.frame(series="Pt4", title="X5 Rx:A RxH:B UnRx:H", direction="Towards UnRx trans", number=reversed_X5_609$same_trans))
df <- rbind(df, data.frame(series="Pt4", title="X5 Rx:A RxH:B UnRx:H", direction="Towards UnRx unknown", number=reversed_X5_609$same_unkn))
df <- rbind(df, data.frame(series="Pt4", title="X5 Rx:A RxH:B UnRx:H", direction="Away from UnRx cis", number=reversed_X5_609$opp_cis))
df <- rbind(df, data.frame(series="Pt4", title="X5 Rx:A RxH:B UnRx:H", direction="Away from UnRx trans", number=reversed_X5_609$opp_trans))
df <- rbind(df, data.frame(series="Pt4", title="X5 Rx:A RxH:B UnRx:H", direction="Away from UnRx unknown", number=reversed_X5_609$opp_unkn))

df <- rbind(df, data.frame(series="Pt4", title="X6 Rx:A RxH:A UnRx:H", direction="Towards UnRx cis",number=reversed_X6_609$same_cis))
df <- rbind(df, data.frame(series="Pt4", title="X6 Rx:A RxH:A UnRx:H", direction="Towards UnRx trans",number=reversed_X6_609$same_trans))
df <- rbind(df, data.frame(series="Pt4", title="X6 Rx:A RxH:A UnRx:H", direction="Towards UnRx unknown",number=reversed_X6_609$same_unkn))
df <- rbind(df, data.frame(series="Pt4", title="X6 Rx:A RxH:A UnRx:H", direction="Away from UnRx cis",number=reversed_X6_609$opp_cis))
df <- rbind(df, data.frame(series="Pt4", title="X6 Rx:A RxH:A UnRx:H", direction="Away from UnRx trans",number=reversed_X6_609$opp_trans))
df <- rbind(df, data.frame(series="Pt4", title="X6 Rx:A RxH:A UnRx:H", direction="Away from UnRx unknown",number=reversed_X6_609$opp_unkn))

df <- rbind(df, data.frame(series="Pt4", title="X7 Rx:A RxH:A UnRx:H", direction="Towards UnRx cis", number=reversed_X7_609$same_cis))
df <- rbind(df, data.frame(series="Pt4", title="X7 Rx:A RxH:A UnRx:H", direction="Towards UnRx trans", number=reversed_X7_609$same_trans))
df <- rbind(df, data.frame(series="Pt4", title="X7 Rx:A RxH:A UnRx:H", direction="Towards UnRx unknown", number=reversed_X7_609$same_unkn))
df <- rbind(df, data.frame(series="Pt4", title="X7 Rx:A RxH:A UnRx:H", direction="Away from UnRx cis", number=reversed_X7_609$opp_cis))
df <- rbind(df, data.frame(series="Pt4", title="X7 Rx:A RxH:A UnRx:H", direction="Away from UnRx trans", number=reversed_X7_609$opp_trans))
df <- rbind(df, data.frame(series="Pt4", title="X7 Rx:A RxH:A UnRx:H", direction="Away from UnRx unknown", number=reversed_X7_609$opp_unkn))

df <- rbind(df, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ UnRx:G", direction="Towards UnRx cis", number=reversed_X7_535$same_cis))
df <- rbind(df, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ UnRx:G", direction="Towards UnRx trans", number=reversed_X7_535$same_trans))
df <- rbind(df, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ UnRx:G", direction="Towards UnRx unknown", number=reversed_X7_535$same_unkn))
df <- rbind(df, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ UnRx:G", direction="Away from UnRx cis", number=reversed_X7_535$opp_cis))
df <- rbind(df, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ UnRx:G", direction="Away from UnRx trans", number=reversed_X7_535$opp_trans))
df <- rbind(df, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ UnRx:G", direction="Away from UnRx unknown", number=reversed_X7_535$opp_unkn))

df <- rbind(df, data.frame(series="Pt5", title="X8 Rx:A RxH:D UnRx:G", direction="Towards UnRx cis", number=reversed_X8_535$same_cis))
df <- rbind(df, data.frame(series="Pt5", title="X8 Rx:A RxH:D UnRx:G", direction="Towards UnRx trans", number=reversed_X8_535$same_trans))
df <- rbind(df, data.frame(series="Pt5", title="X8 Rx:A RxH:D UnRx:G", direction="Towards UnRx unknown", number=reversed_X8_535$same_unkn))
df <- rbind(df, data.frame(series="Pt5", title="X8 Rx:A RxH:D UnRx:G", direction="Away from UnRx cis", number=reversed_X8_535$opp_cis))
df <- rbind(df, data.frame(series="Pt5", title="X8 Rx:A RxH:D UnRx:G", direction="Away from UnRx trans", number=reversed_X8_535$opp_trans))
df <- rbind(df, data.frame(series="Pt5", title="X8 Rx:A RxH:D UnRx:G", direction="Away from UnRx unknown", number=reversed_X8_535$opp_unkn))

#df <- rbind(df, data.frame(series="SA535", title="SA535 X9", number=length(reversed_X9_535)))
df <- rbind(df, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E UnRx:G", direction="Towards UnRx cis", number=reversed_X10_535$same_cis))
df <- rbind(df, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E UnRx:G", direction="Towards UnRx trans", number=reversed_X10_535$same_trans))
df <- rbind(df, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E UnRx:G", direction="Towards UnRx unknown", number=reversed_X10_535$same_unkn))
df <- rbind(df, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E UnRx:G", direction="Away from UnRx cis", number=reversed_X10_535$opp_cis))
df <- rbind(df, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E UnRx:G", direction="Away from UnRx trans", number=reversed_X10_535$opp_trans))
df <- rbind(df, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E UnRx:G", direction="Away from UnRx unknown", number=reversed_X10_535$opp_unkn))

df <- rbind(df, data.frame(series="Pt6", title="X6 Rx:G RxH:B UnRx:G", direction="Towards UnRx cis", number=reversed_X6_1035$same_cis))
df <- rbind(df, data.frame(series="Pt6", title="X6 Rx:G RxH:B UnRx:G", direction="Towards UnRx trans", number=reversed_X6_1035$same_trans))
df <- rbind(df, data.frame(series="Pt6", title="X6 Rx:G RxH:B UnRx:G", direction="Towards UnRx unknown", number=reversed_X6_1035$same_unkn))
df <- rbind(df, data.frame(series="Pt6", title="X6 Rx:G RxH:B UnRx:G", direction="Away from UnRx cis", number=reversed_X6_1035$opp_cis))
df <- rbind(df, data.frame(series="Pt6", title="X6 Rx:G RxH:B UnRx:G", direction="Away from UnRx trans", number=reversed_X6_1035$opp_trans))
df <- rbind(df, data.frame(series="Pt6", title="X6 Rx:G RxH:B UnRx:G", direction="Away from UnRx unknown", number=reversed_X6_1035$opp_unkn))

df <- rbind(df, data.frame(series="Pt6", title="X7 Rx:H RxH:G UnRx:E", direction="Towards UnRx cis", number=reversed_X7_1035$same_cis))
df <- rbind(df, data.frame(series="Pt6", title="X7 Rx:H RxH:G UnRx:E", direction="Towards UnRx trans", number=reversed_X7_1035$same_trans))
df <- rbind(df, data.frame(series="Pt6", title="X7 Rx:H RxH:G UnRx:E", direction="Towards UnRx unknown", number=reversed_X7_1035$same_unkn))
df <- rbind(df, data.frame(series="Pt6", title="X7 Rx:H RxH:G UnRx:E", direction="Away from UnRx cis", number=reversed_X7_1035$opp_cis))
df <- rbind(df, data.frame(series="Pt6", title="X7 Rx:H RxH:G UnRx:E", direction="Away from UnRx trans", number=reversed_X7_1035$opp_trans))
df <- rbind(df, data.frame(series="Pt6", title="X7 Rx:H RxH:G UnRx:E", direction="Away from UnRx unknown", number=reversed_X7_1035$opp_unkn))

df <- rbind(df, data.frame(series="Pt6", title="X8 Rx:H RxH:G UnRx:E", direction="Towards UnRx cis", number=reversed_X8_1035$same_cis))
df <- rbind(df, data.frame(series="Pt6", title="X8 Rx:H RxH:G UnRx:E", direction="Towards UnRx trans", number=reversed_X8_1035$same_trans))
df <- rbind(df, data.frame(series="Pt6", title="X8 Rx:H RxH:G UnRx:E", direction="Towards UnRx unknown", number=reversed_X8_1035$same_unkn))
df <- rbind(df, data.frame(series="Pt6", title="X8 Rx:H RxH:G UnRx:E", direction="Away from UnRx cis", number=reversed_X8_1035$opp_cis))
df <- rbind(df, data.frame(series="Pt6", title="X8 Rx:H RxH:G UnRx:E", direction="Away from UnRx trans", number=reversed_X8_1035$opp_trans))
df <- rbind(df, data.frame(series="Pt6", title="X8 Rx:H RxH:G UnRx:E", direction="Away from UnRx unknown", number=reversed_X8_1035$opp_unkn))

df$series <- factor(df$series, levels=c("Pt4","Pt5","Pt6"))

source("ggplot_utils.R")

ggplot(data=df, aes(x=title, y=number, fill=direction)) +
  geom_bar(stat="identity") + 
  labs(y = "No. of genes", x = "Sample") +
  thesis_theme + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("darkred","coral2","burlywood", "deepskyblue4","deepskyblue","grey")) +
  theme(legend.position = "none") +    # no legend
  theme(axis.title.x=element_blank(),
    axis.text.x=element_blank()) + 
  facet_grid(~ series, scales="free_x") 

ggsave(paste0("manuscript_files/reversed_genes_clone", clone, "_logfc", reversed_genes_logfc, ".pdf"), width=3, heigh=2, useDingbats=FALSE)


####################################################
#### ACTIVATED AND REPRESSED GENES --- Fig 6a

activated_repressed_logfc <- 0.5

find_non_transient_genes <- function(fileTvsU, fileTvsH, fileHvsU, pt, cloneT, cloneU, filename=NULL) {
  deTvsU <- read.csv(fileTvsU)  ## DE
  deTvsH <- read.csv(fileTvsH)  ## DE
  deHvsU <- read.csv(fileHvsU)  ## not DE 
  
  
  ## get NON DE of Rx vs RxH
  deTvsH <- getnonDE_all(deTvsH, fdr=0.1)
  
  # get DE between Rx and UnRx
  #deTvsU <- getDE(deTvsU)
  
  # or include all genes from TvsU
  deTvsU <- getDE_all(deTvsU, fdr=0.01, logFC=activated_repressed_logfc)
  
  #nondeHvsU <- getnonDE(deHvsU)
  #reversed <- intersect(intersect(deTvsU$gene_symbol, deTvsH$gene_symbol), nondeHvsU$gene_symbol)
  
  common <- intersect(deTvsU$gene_symbol, deTvsH$gene_symbol)
  deTvsH <- deTvsH[deTvsH$gene_symbol %in% common,]
  deTvsU <- deTvsU[deTvsU$gene_symbol %in% common,]
  m <- merge(deTvsH, deTvsU, by = "gene_symbol")
  
  # adding the new cis/trans
  pt$ct <- ifelse(pt[[cloneT]]==pt[[cloneU]], "trans", "cis")
  m <- merge(m,pt, by.x="gene_symbol", by.y="gene_symbol")
  m <- m[!is.na(m$ct),]
  #if(sum(is.na(m$ct)) > 0) {
  #  m[is.na(m$ct),]$ct <- "unkn"
  #}

  #m <- m[!is.na(m$ct),]
  
  #return(m)
  activated_cis <- sum (m$logFC.y > 0 & m$ct=="cis")
  activated_trans <- sum (m$logFC.y > 0 & m$ct=="trans")
  activated_unkn <- sum (m$logFC.y > 0 & m$ct=="unkn")
  repressed_cis <- sum (m$logFC.y < 0 & m$ct=="cis")
  repressed_trans <- sum (m$logFC.y < 0 & m$ct=="trans")
  repressed_unkn <- sum (m$logFC.y < 0 & m$ct=="unkn")
  
  if (!is.null(filename)) {
    library(dplyr)
    m <- select(m,gene_symbol,ensembl_gene_id.x,logFC.x,FDR.x,logFC.y,FDR.y)
    colnames(m) <- c("gene_symbol","ensembl_gene_id","TvsH.logFC","TvsH.FDR","TvsU.logFC","TvsU.FDR")
    m$direction <- ifelse (m$TvsU.logFC > 0, "Activated", "Repressed")
    write.csv(m, file=filename, quote=FALSE, row.names=FALSE)
  }
  
  return(list(act_cis=activated_cis, act_trans=activated_trans, act_unkn=activated_unkn, 
              rep_cis=repressed_cis, rep_trans=repressed_trans, rep_unkn=repressed_unkn))
}



if (clone=="aware") {
  
  
  
  
  
  nontr_X5_609 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_2_SA609_UTT_R_UUU_H_logfc_results.csv",
                                           fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_20_SA609_UTT_R_UTU_B_logfc_results.csv",
                                           fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_22_SA609_UTU_B_UUU_H_logfc_results.csv",
                                           pt=Pt4_cnv_whole, cloneT="cn_median_R", cloneU="cn_median_H", 
                                                 filename="manuscript_files/SA609_X5_non_transient_genes.csv")
  
  ### No genes reversed at thr fdr 0.1. Just VIM reversed at thr fdr 0.05, but line plot doesn't look consistent
  
  nontr_X6_609 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_3_SA609_UTTT_R_UUUU_H_logfc_results.csv",
                                                 fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_25_SA609_UTTT_R_UTTU_R_logfc_results.csv",
                                                 fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_23_SA609_UTTU_R_UUUU_H_logfc_results.csv",
                                           pt=Pt4_cnv_whole, cloneT="cn_median_R", cloneU="cn_median_H", 
                                                 filename="manuscript_files/SA609_X6_non_transient_genes.csv")
  
  nontr_X7_609 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_4_SA609_UTTTT_R_UUUUU_H_logfc_results.csv",
                                                 fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_26_SA609_UTTTT_R_UTTTU_R_logfc_results.csv",
                                                 fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_24_SA609_UTTTU_R_UUUUU_H_logfc_results.csv",
                                           pt=Pt4_cnv_whole, cloneT="cn_median_R", cloneU="cn_median_H", 
                                                 filename="manuscript_files/SA609_X7_non_transient_genes.csv")
  
  
  ##### SA535 
  nontr_X7_535 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_5_SA535_UUTT_A_E_J_UUUU_G_logfc_results.csv",
                                                 fileTvsH="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_51_SA535_UUTT_A_E_J_UUTU_A_E_J_logfc_results.csv",
                                                 fileHvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_52_SA535_UUTU_A_E_J_UUUU_G_logfc_results.csv",
                                           pt=Pt5_cnv_whole, cloneT="cn_median_A", cloneU="cn_median_G", 
                                              filename="manuscript_files/SA535_X7_non_transient_genes.csv")
  
  nontr_X8_535 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_1_SA535_UUTTT_A_UUUUU_G_logfc_results.csv",
                                                 fileTvsH="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_15_SA535_UUTTT_A_UUTTU_D_logfc_results.csv",
                                                 fileHvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_16_SA535_UUTTU_D_UUUUU_G_logfc_results.csv",
                                           pt=Pt5_cnv_whole, cloneT="cn_median_A", cloneU="cn_median_G", 
                                                 filename="manuscript_files/SA535_X8_non_transient_genes.csv")
  
  #source("plot_genes_v6.R")
  #plot_genes_for_series(dataset="SA535", drug="cisplatin", version="v7", 
  #                      genes=reversed_X8_535, 
  #                      pathway="NOCL_HOLIDAY_REVERSAL_X8", normalization="sctransform")
  
  
  
  nontr_X10_535 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_3_SA535_UUTTTTT_A_UUUUUU_G_logfc_results.csv",
                                                  fileTvsH="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_35_SA535_UUTTTTT_A_UUTTTTU_E_logfc_results.csv",
                                                  fileHvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_36_SA535_UUTTTTU_E_UUUUUU_G_logfc_results.csv",
                                            pt=Pt5_cnv_whole, cloneT="cn_median_A", cloneU="cn_median_G", 
                                               filename="manuscript_files/SA535_X10_non_transient_genes.csv")
  
  
  
  nontr_X6_1035 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_3_SA1035_UTT_G_UUU_G_logfc_results.csv",
                                                  fileTvsH="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_31_SA1035_UTT_G_UTU_B_logfc_results.csv",
                                                  fileHvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_32_SA1035_UTU_B_UUU_G_logfc_results.csv",
                                            pt=Pt6_cnv_whole, cloneT="cn_median_G", cloneU="cn_median_G", 
                                               filename="manuscript_files/SA1035_X6_non_transient_genes.csv")
  
  nontr_X7_1035 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_1_SA1035_UTTT_H_UUUU_E_logfc_results.csv",
                                                  fileTvsH="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_11_SA1035_UTTT_H_UTTU_G_logfc_results.csv",
                                                  fileHvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_12_SA1035_UTTU_G_UUUU_E_logfc_results.csv",
                                            pt=Pt6_cnv_whole, cloneT="cn_median_H", cloneU="cn_median_E", 
                                               filename="manuscript_files/SA1035_X7_non_transient_genes.csv")
  
  nontr_X8_1035 <- find_non_transient_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_2_SA1035_UTTTT_H_UUUUU_E_logfc_results.csv",
                                                  fileTvsH="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_21_SA1035_UTTTT_H_UTTTU_G_logfc_results.csv",
                                                  fileHvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_22_SA1035_UTTTU_G_UUUUU_E_logfc_results.csv",
                                            pt=Pt6_cnv_whole, cloneT="cn_median_H", cloneU="cn_median_E", 
                                                  filename="manuscript_files/SA1035_X8_non_transient_genes.csv")
  
}

df2 <- NULL
df2 <- rbind(df2, data.frame(series="Pt4", title="X5 Rx:A RxH:B Un:H", direction="Activated cis", number=nontr_X5_609$act_cis))
df2 <- rbind(df2, data.frame(series="Pt4", title="X5 Rx:A RxH:B Un:H", direction="Activated trans", number=nontr_X5_609$act_trans))
df2 <- rbind(df2, data.frame(series="Pt4", title="X5 Rx:A RxH:B Un:H", direction="Activated unknown", number=nontr_X5_609$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt4", title="X5 Rx:A RxH:B Un:H", direction="Repressed cis", number=nontr_X5_609$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt4", title="X5 Rx:A RxH:B Un:H", direction="Repressed trans", number=nontr_X5_609$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt4", title="X5 Rx:A RxH:B Un:H", direction="Repressed unknown", number=nontr_X5_609$rep_unkn))

df2 <- rbind(df2, data.frame(series="Pt4", title="X6 Rx:A RxH:A Un:H", direction="Activated cis",number=nontr_X6_609$act_cis))
df2 <- rbind(df2, data.frame(series="Pt4", title="X6 Rx:A RxH:A Un:H", direction="Activated trans",number=nontr_X6_609$act_trans))
df2 <- rbind(df2, data.frame(series="Pt4", title="X6 Rx:A RxH:A Un:H", direction="Activated unknown",number=nontr_X6_609$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt4", title="X6 Rx:A RxH:A Un:H", direction="Repressed cis",number=nontr_X6_609$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt4", title="X6 Rx:A RxH:A Un:H", direction="Repressed trans",number=nontr_X6_609$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt4", title="X6 Rx:A RxH:A Un:H", direction="Repressed unknown",number=nontr_X6_609$rep_unkn))

df2 <- rbind(df2, data.frame(series="Pt4", title="X7 Rx:A RxH:A Un:H", direction="Activated cis", number=nontr_X7_609$act_cis))
df2 <- rbind(df2, data.frame(series="Pt4", title="X7 Rx:A RxH:A Un:H", direction="Activated trans", number=nontr_X7_609$act_trans))
df2 <- rbind(df2, data.frame(series="Pt4", title="X7 Rx:A RxH:A Un:H", direction="Activated unknown", number=nontr_X7_609$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt4", title="X7 Rx:A RxH:A Un:H", direction="Repressed cis", number=nontr_X7_609$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt4", title="X7 Rx:A RxH:A Un:H", direction="Repressed trans", number=nontr_X7_609$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt4", title="X7 Rx:A RxH:A Un:H", direction="Repressed unknown", number=nontr_X7_609$rep_unkn))

df2 <- rbind(df2, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ Un:G", direction="Activated cis", number=nontr_X7_535$act_cis))
df2 <- rbind(df2, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ Un:G", direction="Activated trans", number=nontr_X7_535$act_trans))
df2 <- rbind(df2, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ Un:G", direction="Activated unknown", number=nontr_X7_535$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ Un:G", direction="Repressed cis", number=nontr_X7_535$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ Un:G", direction="Repressed trans", number=nontr_X7_535$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt5", title="X7 Rx:AEJ RxH:AEJ Un:G", direction="Repressed unknown", number=nontr_X7_535$rep_unkn))

df2 <- rbind(df2, data.frame(series="Pt5", title="X8 Rx:A RxH:D Un:G", direction="Activated cis", number=nontr_X8_535$act_cis))
df2 <- rbind(df2, data.frame(series="Pt5", title="X8 Rx:A RxH:D Un:G", direction="Activated trans", number=nontr_X8_535$act_trans))
df2 <- rbind(df2, data.frame(series="Pt5", title="X8 Rx:A RxH:D Un:G", direction="Activated unknown", number=nontr_X8_535$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt5", title="X8 Rx:A RxH:D Un:G", direction="Repressed cis", number=nontr_X8_535$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt5", title="X8 Rx:A RxH:D Un:G", direction="Repressed trans", number=nontr_X8_535$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt5", title="X8 Rx:A RxH:D Un:G", direction="Repressed unknown", number=nontr_X8_535$rep_unkn))
#df <- rbind(df, data.frame(series="SA535", title="SA535 X9", number=length(reversed_X9_535)))
df2 <- rbind(df2, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E Un:G", direction="Activated cis", number=nontr_X10_535$act_cis))
df2 <- rbind(df2, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E Un:G", direction="Activated trans", number=nontr_X10_535$act_trans))
df2 <- rbind(df2, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E Un:G", direction="Activated unknown", number=nontr_X10_535$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E Un:G", direction="Repressed cis", number=nontr_X10_535$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E Un:G", direction="Repressed trans", number=nontr_X10_535$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt5", title="X9+1 Rx:A RxH:E Un:G", direction="Repressed unknown", number=nontr_X10_535$rep_unkn))

df2 <- rbind(df2, data.frame(series="Pt6", title="X6 Rx:G RxH:B Un:G", direction="Activated cis", number=nontr_X6_1035$act_cis))
df2 <- rbind(df2, data.frame(series="Pt6", title="X6 Rx:G RxH:B Un:G", direction="Activated trans", number=nontr_X6_1035$act_trans))
df2 <- rbind(df2, data.frame(series="Pt6", title="X6 Rx:G RxH:B Un:G", direction="Activated unknown", number=nontr_X6_1035$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt6", title="X6 Rx:G RxH:B Un:G", direction="Repressed cis", number=nontr_X6_1035$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt6", title="X6 Rx:G RxH:B Un:G", direction="Repressed trans", number=nontr_X6_1035$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt6", title="X6 Rx:G RxH:B Un:G", direction="Repressed unknown", number=nontr_X6_1035$rep_unkn))

df2 <- rbind(df2, data.frame(series="Pt6", title="X7 Rx:H RxH:G Un:E", direction="Activated cis", number=nontr_X7_1035$act_cis))
df2 <- rbind(df2, data.frame(series="Pt6", title="X7 Rx:H RxH:G Un:E", direction="Activated trans", number=nontr_X7_1035$act_trans))
df2 <- rbind(df2, data.frame(series="Pt6", title="X7 Rx:H RxH:G Un:E", direction="Activated unknown", number=nontr_X7_1035$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt6", title="X7 Rx:H RxH:G Un:E", direction="Repressed cis", number=nontr_X7_1035$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt6", title="X7 Rx:H RxH:G Un:E", direction="Repressed trans", number=nontr_X7_1035$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt6", title="X7 Rx:H RxH:G Un:E", direction="Repressed unknown", number=nontr_X7_1035$rep_unkn))

df2 <- rbind(df2, data.frame(series="Pt6", title="X8 Rx:H RxH:G Un:E", direction="Activated cis", number=nontr_X8_1035$act_cis))
df2 <- rbind(df2, data.frame(series="Pt6", title="X8 Rx:H RxH:G Un:E", direction="Activated trans", number=nontr_X8_1035$act_trans))
df2 <- rbind(df2, data.frame(series="Pt6", title="X8 Rx:H RxH:G Un:E", direction="Activated unknown", number=nontr_X8_1035$act_unkn))
df2 <- rbind(df2, data.frame(series="Pt6", title="X8 Rx:H RxH:G Un:E", direction="Repressed cis", number=nontr_X8_1035$rep_cis))
df2 <- rbind(df2, data.frame(series="Pt6", title="X8 Rx:H RxH:G Un:E", direction="Repressed trans", number=nontr_X8_1035$rep_trans))
df2 <- rbind(df2, data.frame(series="Pt6", title="X8 Rx:H RxH:G Un:E", direction="Repressed unknown", number=nontr_X8_1035$rep_unkn))

df2$series <- factor(df2$series, levels=c("Pt4","Pt5","Pt6"))

### FIg 6a

ggplot(data=df2, aes(x=title, y=number, fill=direction)) +
  geom_bar(stat="identity") + 
  labs(y = "No. of genes", x = "Sample") +
  thesis_theme + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("darkred","coral2","burlywood","deepskyblue4","deepskyblue","grey")) +
  theme(legend.position = "none") +    # no legend
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + 
  facet_grid(~ series, scales="free_x") 

ggsave(paste0("manuscript_files/activated_or_repressed_genes_clone", clone, "_logfc", activated_repressed_logfc, ".pdf"), width=3, heigh=2, useDingbats=FALSE)





find_info_reversed_holiday_genes <- function(fileTvsU, fileTvsH, fileHvsU) {
  ### same as previous function, but returns more info
  deTvsU <- read.csv(fileTvsU)  ## DE
  deTvsH <- read.csv(fileTvsH)  ## DE
  deHvsU <- read.csv(fileHvsU)  ## not DE 
  
  
  ## get DE of Rx vs RxH
  deTvsH <- getDE(deTvsH)
  
  # also use deTvsU
  #deTvsU <- getDE(deTvsU)
  
  # or include all genes from TvsU
  deTvsU <- getDE(deTvsU, fdr=1, logFC=0)
  
  #nondeHvsU <- getnonDE(deHvsU)
  #reversed <- intersect(intersect(deTvsU$gene_symbol, deTvsH$gene_symbol), nondeHvsU$gene_symbol)
  
  common <- intersect(deTvsU$gene_symbol, deTvsH$gene_symbol)
  deTvsH <- deTvsH[deTvsH$gene_symbol %in% common,]
  deTvsU <- deTvsU[deTvsU$gene_symbol %in% common,]
  m <- merge(deTvsH, deTvsU, by = "gene_symbol")
  return(m)
}


plot_degree_rev <- function (info, title) {
  library(ggrepel)
  info$col <- info$logFC.x * info$logFC.y > 0
  ggplot(info, aes(x = logFC.x, y = logFC.y)) + 
    #geom_point(size = 5, color = "#0099f9") + 
    geom_point(size = 5, aes(color = col)) +
    geom_text_repel(label = info$gene_symbol,  size=3.5) + 
    labs(
      x = "logFC Rx vs. RxH",
      y = "logFC Rx vs. UnRx",
      title = title
    )
}



## get DE of Rx vs RxH
deTvsH <- getDE_all(deTvsH, fdr=0.01, logFC=reversed_genes_logfc)

# also use deTvsU
#deTvsU <- getDE(deTvsU)

# or include all genes from TvsU
deTvsU <- getDE_all(deTvsU, fdr=1, logFC=0)

#nondeHvsU <- getnonDE(deHvsU)
#reversed <- intersect(intersect(deTvsU$gene_symbol, deTvsH$gene_symbol), nondeHvsU$gene_symbol)

common <- intersect(deTvsU$gene_symbol, deTvsH$gene_symbol)
deTvsH <- deTvsH[deTvsH$gene_symbol %in% common,]
deTvsU <- deTvsU[deTvsU$gene_symbol %in% common,]
m <- merge(deTvsH, deTvsU, by = "gene_symbol")


#return(m)
same_direction_cis <- sum (m$logFC.x * m$logFC.y > 0 & m$cnchange.y !=0)
same_direction_trans <- sum (m$logFC.x * m$logFC.y > 0 & m$cnchange.y ==0)
opposite_direction_cis <- sum (m$logFC.x * m$logFC.y < 0 & m$cnchange.y !=0)
opposite_direction_trans <- sum (m$logFC.x * m$logFC.y < 0 & m$cnchange.y ==0)





### we have to compare vs UnRx, not holiday vs others
plot_reversed_holiday_genes <- function(fileTvsU, fileTvsH, fileHvsU, title, pt, cloneT, cloneU, filename) {
  ### same as previous function, but returns more info
  deTvsU <- read.csv(fileTvsU)  ## DE
  deTvsH <- read.csv(fileTvsH)  ## DE
  deHvsU <- read.csv(fileHvsU)  ## not DE 
  
  
  
  ## get DE of Rx vs RxH
  deTvsH <- getDE_all(deTvsH, fdr=0.01, logFC=reversed_genes_logfc)
  
  # also use deTvsU
  #deTvsU <- getDE(deTvsU)
  
  # or include all genes from TvsU
  deTvsU <- getDE_all(deTvsU, fdr=1, logFC=0)
  
  #nondeHvsU <- getnonDE(deHvsU)
  #reversed <- intersect(intersect(deTvsU$gene_symbol, deTvsH$gene_symbol), nondeHvsU$gene_symbol)
  
  common <- intersect(deTvsU$gene_symbol, deTvsH$gene_symbol)
  deTvsH <- deTvsH[deTvsH$gene_symbol %in% common,]
  deTvsU <- deTvsU[deTvsU$gene_symbol %in% common,]
  m <- merge(deTvsH, deTvsU, by = "gene_symbol")  
  
  # adding the new cis/trans
  pt$ct <- ifelse(pt[[cloneT]]==pt[[cloneU]], "trans", "cis")
  m <- merge(m,pt, by.x="gene_symbol", by.y="Symbol")
  m[is.na(m$ct),]$ct <- "unkn"

  info <- m
  library(ggrepel)
  info$col <- info$logFC.x * info$logFC.y > 0
  info$label <- 0
  info$color <- 0
  
  info[info$logFC.x * info$logFC.y < 0 & info$ct=="cis",c("label","color")] <- data.frame(label="Away from UnRx cis", color="darkred")
  info[info$logFC.x * info$logFC.y < 0 & info$ct=="trans",c("label","color")] <- data.frame(label="Away from UnRx trans", color="coral2")  
  info[info$logFC.x * info$logFC.y > 0 & info$ct=="cis",c("label","color")] <- data.frame(label="Towards UnRx cis", color="deepskyblue4")
  info[info$logFC.x * info$logFC.y > 0 & info$ct=="trans",c("label","color")] <- data.frame(label="Towards UnRx trans", color="deepskyblue")
  
  print(table(info$label))
  mycol <- as.character(info$color)
  names(mycol) <- as.character(info$label)
  
  write.csv(info, file=paste0(filename, ".csv"), quote=FALSE, row.names=FALSE)
  
  ggplot(info, aes(y = logFC.x, x = logFC.y)) + 
    #geom_point(size = 5, color = "#0099f9") + 
    geom_point(size = 1, aes(color = label), alpha=0.5) +
    scale_color_manual(values=mycol) +
    #scale_color_manual(values = c(breaks ="darkred","coral2","deepskyblue4","deepskyblue",
    #                              values= "darkred","coral2","deepskyblue4","deepskyblue")) +
    thesis_theme + 
    scale_x_continuous(limits = c(-3, 3)) +
    scale_y_continuous(limits = c(-3, 3)) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    geom_text_repel(label = info$gene_symbol,  size=2, max.overlaps = 50, segment.alpha = 0.2) + 
    theme(legend.position = "none") +
    labs(
      x = "logFC Rx vs. Un",
      y = "logFC Rx vs. RxH",
      #y = "logFC RxH vs. Rx",
      #x = "logFC RxH vs. UnRx",
      title = title
    )
  ggsave(paste0(filename, ".pdf"), width=2.5, heigh=2.2, useDingbats=FALSE)
}



### we have to compare vs UnRx, not holiday vs others
plot_activated_repressed_holiday_genes <- function(fileTvsU, fileTvsH, fileHvsU, title, pt, cloneT, cloneU, filename) {
  
  deTvsU <- read.csv(fileTvsU)  ## DE
  deTvsH <- read.csv(fileTvsH)  ## DE
  deHvsU <- read.csv(fileHvsU)  ## not DE 
  
  
  ## get NON DE of Rx vs RxH
  deTvsH <- getnonDE_all(deTvsH, fdr=0.1)
  
  # get DE between Rx and UnRx
  #deTvsU <- getDE(deTvsU)
  
  # or include all genes from TvsU
  deTvsU <- getDE_all(deTvsU, fdr=0.01, logFC=activated_repressed_logfc)
  
  #nondeHvsU <- getnonDE(deHvsU)
  #reversed <- intersect(intersect(deTvsU$gene_symbol, deTvsH$gene_symbol), nondeHvsU$gene_symbol)
  
  common <- intersect(deTvsU$gene_symbol, deTvsH$gene_symbol)
  deTvsH <- deTvsH[deTvsH$gene_symbol %in% common,]
  deTvsU <- deTvsU[deTvsU$gene_symbol %in% common,]
  m <- merge(deTvsH, deTvsU, by = "gene_symbol")
  
  # adding the new cis/trans
  pt$ct <- ifelse(pt[[cloneT]]==pt[[cloneU]], "trans", "cis")
  m <- merge(m,pt, by.x="gene_symbol", by.y="Symbol")
  m[is.na(m$ct),]$ct <- "unkn"
  #m <- m[!is.na(m$ct),]
  
  
  info <- m
  library(ggrepel)
  info$col <- info$logFC.x * info$logFC.y > 0
  info$label <- 0
  info$color <- 0
  
  
  info[info$logFC.y > 0 & info$ct=="cis",c("label","color")] <- data.frame(label="Activated cis", color="darkred")
  info[info$logFC.y > 0 & info$ct=="trans",c("label","color")] <- data.frame(label="Activated trans", color="coral2")  
  info[info$logFC.y < 0 & info$ct=="cis",c("label","color")] <- data.frame(label="Repressed cis", color="deepskyblue4")
  info[info$logFC.y < 0 & info$ct=="trans",c("label","color")] <- data.frame(label="Repressed trans", color="deepskyblue")
  
  print(table(info$label))
  mycol <- as.character(info$color)
  names(mycol) <- as.character(info$label)
  
  write.csv(info, file=paste0(filename, ".csv"), quote=FALSE, row.names=FALSE)
  
  ggplot(info, aes(x = logFC.x, y = logFC.y)) + 
    #geom_point(size = 5, color = "#0099f9") + 
    geom_point(size = 1, aes(color = label), alpha=0.5) +
    scale_color_manual(values=mycol) +
    #scale_color_manual(values = c(breaks ="darkred","coral2","deepskyblue4","deepskyblue",
    #                              values= "darkred","coral2","deepskyblue4","deepskyblue")) +
    thesis_theme + 
    scale_x_continuous(limits = c(-0.17, 0.17)) +
    scale_y_continuous(limits = c(-3, 3)) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    geom_text_repel(label = info$gene_symbol,  size=2, max.overlaps = 50, segment.alpha = 0.2) + 
    theme(legend.position = "none") +
    labs(
      y = "logFC Rx vs. Un",
      x = "logFC Rx vs. RxH"
      #y = "logFC RxH vs. Rx",
      #x = "logFC RxH vs. UnRx",
      #title = title
    )
  ggsave(paste0(filename, ".pdf"), width=1.2, heigh=2.2, useDingbats=FALSE)
}


#### focus on holiday rather than treated



if (clone=="aware") {
  
  
  plot_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_2_SA609_UTT_R_UUU_H_logfc_results.csv",
                              fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_20_SA609_UTT_R_UTU_B_logfc_results.csv",
                              fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_22_SA609_UTU_B_UUU_H_logfc_results.csv",
                              title="Pt4 X5 Rx:A RxH:B Un:H", pt=Pt4_cnv, cloneT="A", cloneU="H",
                              filename=paste0("manuscript_files/transient_pt4_x5_logfc",reversed_genes_logfc))
  
  plot_activated_repressed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_2_SA609_UTT_R_UUU_H_logfc_results.csv",
                              fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_20_SA609_UTT_R_UTU_B_logfc_results.csv",
                              fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_22_SA609_UTU_B_UUU_H_logfc_results.csv",
                              title="Pt4 X5 Rx:A RxH:B Un:H", pt=Pt4_cnv, cloneT="A", cloneU="H",
                              filename=paste0("manuscript_files/act_repr_pt4_x5_logfc",reversed_genes_logfc))
  
  
  plot_activated_repressed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_3_SA609_UTTT_R_UUUU_H_logfc_results.csv",
                                         fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_25_SA609_UTTT_R_UTTU_R_logfc_results.csv",
                                         fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_23_SA609_UTTU_R_UUUU_H_logfc_results.csv",
                                         title="Pt4 X5 Rx:A RxH:B Un:H", pt=Pt4_cnv, cloneT="A", cloneU="H",
                                         filename=paste0("manuscript_files/act_repr_pt4_x6_logfc",reversed_genes_logfc))
  
  plot_activated_repressed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_4_SA609_UTTTT_R_UUUUU_H_logfc_results.csv",
                                         fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_26_SA609_UTTTT_R_UTTTU_R_logfc_results.csv",
                                         fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_24_SA609_UTTTU_R_UUUUU_H_logfc_results.csv",
                                         title="Pt4 X5 Rx:A RxH:B Un:H", pt=Pt4_cnv, cloneT="A", cloneU="H",
                                         filename=paste0("manuscript_files/act_repr_pt4_x7_logfc",reversed_genes_logfc))
  
  plot_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_3_SA609_UTTT_R_UUUU_H_logfc_results.csv",
                                         fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_25_SA609_UTTT_R_UTTU_R_logfc_results.csv",
                                         fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_23_SA609_UTTU_R_UUUU_H_logfc_results.csv",
                                         title="Pt4 X5 Rx:A RxH:B Un:H", pt=Pt4_cnv, cloneT="A", cloneU="H",
                                         filename=paste0("manuscript_files/transient_pt4_x6_logfc",reversed_genes_logfc))
  
  plot_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_4_SA609_UTTTT_R_UUUUU_H_logfc_results.csv",
                                         fileTvsH="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_26_SA609_UTTTT_R_UTTTU_R_logfc_results.csv",
                                         fileHvsU="../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_24_SA609_UTTTU_R_UUUUU_H_logfc_results.csv",
                                         title="Pt4 X5 Rx:A RxH:B Un:H", pt=Pt4_cnv, cloneT="A", cloneU="H",
                                         filename=paste0("manuscript_files/transient_pt4_x7_logfc",reversed_genes_logfc))
  
  
  
  
  ### SUPP
  plot_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_1_SA535_UUTTT_A_UUUUU_G_logfc_results.csv",
                              fileTvsH="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_15_SA535_UUTTT_A_UUTTU_D_logfc_results.csv",
                              fileHvsU="../../../drug_resistance_local/extranalysis/results/SA535-v7/comps2/scrande_SA535_16_SA535_UUTTU_D_UUUUU_G_logfc_results.csv",
                              title="Pt5 X8 Rx:A RxH:D Un:G", pt=Pt5_cnv, cloneT="A", cloneU="G",
                              filename=paste0("manuscript_files/transient_pt5_x8_logfc",reversed_genes_logfc))
  
  plot_reversed_holiday_genes(fileTvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_2_SA1035_UTTTT_H_UUUUU_E_logfc_results2.csv",
                              fileTvsH="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_21_SA1035_UTTTT_H_UTTTU_G_logfc_results2.csv",
                              fileHvsU="../../../drug_resistance_local/extranalysis/results/SA1035-v6/comps2/scrande_SA1035_22_SA1035_UTTTU_G_UUUUU_E_logfc_results2.csv",
                              title= "Pt6 X8 Rx:H RxH:G Un:E", pt=Pt6_cnv, cloneT="H", cloneU="E",
                              filename=paste0("manuscript_files/transient_pt6_x8_logfc",reversed_genes_logfc))
  
  
  
  
  ### In Supplementary
  
  
  plot_reversed_holiday_genes(fileTvsU="../results/SA609-v6/comps2/scrande_SA609_3_SA609_UTTT_R_UUUU_H_logfc_results.csv",
                              fileTvsH="../results/SA609-v6/comps2/scrande_SA609_25_SA609_UTTT_R_UTTU_R_logfc_results.csv",
                              fileHvsU="../results/SA609-v6/comps2/scrande_SA609_23_SA609_UTTU_R_UUUU_H_logfc_results.csv",
                              title="Pt4 X6 Rx:A RxH:A UnRx:H")
  
  plot_reversed_holiday_genes(fileTvsU="../results/SA609-v6/comps2/scrande_SA609_4_SA609_UTTTT_R_UUUUU_H_logfc_results.csv",
                              fileTvsH="../results/SA609-v6/comps2/scrande_SA609_26_SA609_UTTTT_R_UTTTU_R_logfc_results.csv",
                              fileHvsU="../results/SA609-v6/comps2/scrande_SA609_24_SA609_UTTTU_R_UUUUU_H_logfc_results.csv",
                              title="Pt4 X7 Rx:A RxH:A UnRx:H",
                              filename="../../figures/summary_figure/magnitude_pt5_x7.pdf")  
  
  plot_reversed_holiday_genes(fileTvsU="../results/SA535-v7/comps2/scrande_SA535_5_SA535_UUTT_A_E_J_UUUU_G_logfc_results.csv",
                              fileTvsH="../results/SA535-v7/comps2/scrande_SA535_51_SA535_UUTT_A_E_J_UUTU_A_E_J_logfc_results.csv",
                              fileHvsU="../results/SA535-v7/comps2/scrande_SA535_52_SA535_UUTU_A_E_J_UUUU_G_logfc_results.csv",
                              title="Pt5 X7 Rx:AEJ RxH:AEJ UnRx:G")
  
  plot_reversed_holiday_genes(fileTvsU="../results/SA535-v7/comps2/scrande_SA535_3_SA535_UUTTTTT_A_UUUUUU_G_logfc_results.csv",
                              fileTvsH="../results/SA535-v7/comps2/scrande_SA535_35_SA535_UUTTTTT_A_UUTTTTU_E_logfc_results.csv",
                              fileHvsU="../results/SA535-v7/comps2/scrande_SA535_36_SA535_UUTTTTU_E_UUUUUU_G_logfc_results.csv",
                              title="Pt5 X10 Rx:A RxH:D UnRx:G")  
  
  plot_reversed_holiday_genes(fileTvsU="../results/SA1035-v6/comps2/scrande_SA1035_3_SA1035_UTT_G_UUU_G_logfc_results.csv",
                              fileTvsH="../results/SA1035-v6/comps2/scrande_SA1035_31_SA1035_UTT_G_UTU_B_logfc_results.csv",
                              fileHvsU="../results/SA1035-v6/comps2/scrande_SA1035_32_SA1035_UTU_B_UUU_G_logfc_results.csv",
                              title= "Pt6 X6 Rx:G RxH:B UnRx:G")
  
  plot_reversed_holiday_genes(fileTvsU="../results/SA1035-v6/comps2/scrande_SA1035_1_SA1035_UTTT_H_UUUU_E_logfc_results.csv",
                              fileTvsH="../results/SA1035-v6/comps2/scrande_SA1035_11_SA1035_UTTT_H_UTTU_G_logfc_results.csv",
                              fileHvsU="../results/SA1035-v6/comps2/scrande_SA1035_12_SA1035_UTTU_G_UUUU_E_logfc_results.csv",
                              title= "Pt6 X7 Rx:H RxH:G UnRx:E")
  
}


##### SA535 



#### Get the common genes between edgeR and slingshot

#traj_pt4 <- read.csv("../results/SA609-v6/patternTest_diffEndTest_SA609_l1_Rx_lg2_RxH_top50.csv")
traj_trans_pt4 <- read.csv("SA609_trajectory_genes_modules/SA609_transient_genes.csv")
traj_activ_pt4 <- read.csv("SA609_trajectory_genes_modules/SA609_activated_genes.csv")
traj_repr_pt4 <- read.csv("SA609_trajectory_genes_modules/SA609_repressed_genes.csv")
edger_trans_pt4_x5 <- read.csv("../results/SA609-v6/SA609_X5_reversed_holiday_genes.csv")
edger_trans_pt4_x6 <- read.csv("../results/SA609-v6/SA609_X6_reversed_holiday_genes.csv")
edger_trans_pt4_x7 <- read.csv("../results/SA609-v6/SA609_X7_reversed_holiday_genes.csv")

edger_actrep_pt4_x5 <- read.csv("../results/SA609-v6/SA609_X5_non_transient_genes.csv")
edger_actrep_pt4_x6 <- read.csv("../results/SA609-v6/SA609_X6_non_transient_genes.csv")
edger_actrep_pt4_x7 <- read.csv("../results/SA609-v6/SA609_X7_non_transient_genes.csv")


pt4_traj_trans_genes <- traj_trans_pt4$gene_symbol
pt4_traj_activ_genes <- traj_activ_pt4$gene_symbol
pt4_traj_repr_genes <- traj_repr_pt4$gene_symbol

pt4_edger_trans_genes <- unique(rbind(edger_trans_pt4_x5,edger_trans_pt4_x6,edger_trans_pt4_x7)$gene_symbol)
pt4_edger_actrep <- rbind(edger_actrep_pt4_x5,edger_actrep_pt4_x6,edger_actrep_pt4_x7)
pt4_edger_activ_genes <- unique(pt4_edger_actrep[pt4_edger_actrep$direction=="Activated",]$gene_symbol)
pt4_edger_repr_genes <- unique(pt4_edger_actrep[pt4_edger_actrep$direction=="Repressed",]$gene_symbol)





### Example of confusion matrix with ggplot

traj_class <- factor(rep(c("Transient", "Activated","Repressed","Total"), each=4))
traj_class <- factor(traj_class, levels=c("Transient", "Activated","Repressed","Total"))
de_class <- factor(rep(c("Transient", "Activated","Repressed", "Total"), times=4))
de_class <- factor(de_class, levels=c("Transient", "Activated","Repressed","Total"))
Y      <- c(length(intersect(pt4_traj_trans_genes,pt4_edger_trans_genes)),
            length(intersect(pt4_traj_trans_genes,pt4_edger_activ_genes)),
            length(intersect(pt4_traj_trans_genes,pt4_edger_repr_genes)),
            length(pt4_traj_trans_genes),
            length(intersect(pt4_traj_activ_genes,pt4_edger_trans_genes)),
            length(intersect(pt4_traj_activ_genes,pt4_edger_activ_genes)),
            length(intersect(pt4_traj_activ_genes,pt4_edger_repr_genes)),
            length(pt4_traj_activ_genes),
            length(intersect(pt4_traj_repr_genes,pt4_edger_trans_genes)),
            length(intersect(pt4_traj_repr_genes,pt4_edger_activ_genes)),
            length(intersect(pt4_traj_repr_genes,pt4_edger_repr_genes)),
            length(pt4_traj_repr_genes),
            length(pt4_edger_trans_genes),
            length(pt4_edger_activ_genes),
            length(pt4_edger_repr_genes),
            0
)
           

traj_class <- factor(rep(c("Transient", "Activated","Repressed"), each=3))
traj_class <- factor(traj_class, levels=c("Transient", "Activated","Repressed"))
de_class <- factor(rep(c("Transient", "Activated","Repressed"), times=3))
de_class <- factor(de_class, levels=c("Transient", "Activated","Repressed"))
Y      <- c(length(intersect(pt4_traj_trans_genes,pt4_edger_trans_genes)),
            length(intersect(pt4_traj_trans_genes,pt4_edger_activ_genes)),
            length(intersect(pt4_traj_trans_genes,pt4_edger_repr_genes)),
            length(intersect(pt4_traj_activ_genes,pt4_edger_trans_genes)),
            length(intersect(pt4_traj_activ_genes,pt4_edger_activ_genes)),
            length(intersect(pt4_traj_activ_genes,pt4_edger_repr_genes)),
            length(intersect(pt4_traj_repr_genes,pt4_edger_trans_genes)),
            length(intersect(pt4_traj_repr_genes,pt4_edger_activ_genes)),
            length(intersect(pt4_traj_repr_genes,pt4_edger_repr_genes))
)
 
df <- data.frame(traj_class, de_class, Y)

library(ggplot2)
ggplot(data =  df, mapping = aes(x = traj_class, y = de_class)) +
  geom_tile(aes(fill = Y), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Y)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none")
