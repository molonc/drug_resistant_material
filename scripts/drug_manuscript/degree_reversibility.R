plot_reversed_holiday_genes <- function(fileTvsU, fileTvsH, fileHvsU, title, output_dir, datatag, plttitle) {
  ### same as previous function, but returns more info
  deTvsU <- read.csv(fileTvsU)  ## DE
  deTvsH <- read.csv(fileTvsH)  ## DE
  deHvsU <- read.csv(fileHvsU)  ## not DE 
  
  dim(deTvsH)
  dim(deHvsU)
  ## TO COME BACK and use directly HvsT
  ## get DE of Rx vs RxH
  deTvsH <- getDE(deTvsH)
  
  deHvsT <- deTvsH
  deHvsT$logFC <- -deHvsT$logFC
  
  # also use deTvsU
  #deTvsU <- getDE(deTvsU)
  
  # or include all genes from TvsU
  deHvsU <- getDE(deHvsU, fdr=1, logFC=0.5)
  
  #nondeHvsU <- getnonDE(deHvsU)
  #reversed <- intersect(intersect(deTvsU$gene_symbol, deTvsH$gene_symbol), nondeHvsU$gene_symbol)
  
  common <- intersect(deHvsU$gene_symbol, deHvsT$gene_symbol)
  deHvsT <- deHvsT[deHvsT$gene_symbol %in% common,]
  deHvsU <- deHvsU[deHvsU$gene_symbol %in% common,]
  dim(deHvsT)
  head(deHvsT)
  # From Hoa
  deHvsT <- deHvsT %>%
    select(gene_symbol,logFC)%>%
    mutate(logFC=abs(logFC))%>%
    rename(absLogFC_H_vs_T=logFC)
  
  deHvsU <- deHvsU %>%
    select(gene_symbol,logFC)%>%
    mutate(logFC=abs(logFC))%>%
    rename(absLogFC_H_vs_U=logFC)
  stat <- deHvsT %>% inner_join(deHvsU, by='gene_symbol')
  # dim(stat)  
  # View(head(stat))
  # colnames(stat)
  ht_thrs <- quantile(stat$absLogFC_H_vs_T,0.75)
  hu_thrs <- quantile(stat$absLogFC_H_vs_U,0.75)
  stat_labels <- stat %>%
    filter(absLogFC_H_vs_T>ht_thrs & absLogFC_H_vs_U>hu_thrs)%>%
    mutate(sum_stat=absLogFC_H_vs_T+absLogFC_H_vs_U)%>%
    dplyr::arrange(desc(sum_stat))
  
  stat_labels <- stat_labels[1:30,]
  dim(stat_labels)
  p <- ggplot(stat, aes(x = absLogFC_H_vs_T, y = absLogFC_H_vs_U)) + 
      geom_point(size = 3, alpha=0.6, aes(color = absLogFC_H_vs_U))+
      scale_color_continuous(low = "blue", high = "red") +
      theme_classic() + 
      geom_abline(intercept =0 , slope = 1, linetype = "dashed")+ xlim(0.5,3)+ylim(0.5,3)+
      geom_text_repel(data = stat_labels, aes(label = gene_symbol), size = 3.7,   # was 3.3
                    min.segment.length = 0, max.overlaps = Inf,
                    color='black', segment.alpha = 0.2)+
      labs(title=plttitle)+
      theme(
        legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size=11)
      )
  
  # p
  png(paste0(save_dir,datatag, "_degree_reverse.png"), height = 2*400, width=2*450,res = 2*72)
  print(p)
  dev.off()
  # m <- merge(deHvsT, deHvsU, by = "gene_symbol")
  # 
  # info <- m
  # library(ggrepel)
  # info$col <- info$logFC.x * info$logFC.y > 0
  # ggplot(info, aes(x = logFC.x, y = logFC.y)) + 
  #   #geom_point(size = 5, color = "#0099f9") + 
  #   geom_point(size = 2, alpha=0.6, aes(color = col)) +
  #   geom_text_repel(label = info$gene_symbol,  size=3.5) + 
  #   labs(
  #     x = "logFC RxH vs. Rx",
  #     y = "logFC RxH vs. UnRx",
  #     title = title
  #   )
}

get_hallmark_genes <- function(){
  ### get the pathway genes
  gmt <- data.frame(read.csv("h.all.v7.0.symbols.gmt",sep='\t',header = FALSE))
  gmt[,2] <- NULL
  
  library(reshape2)
  
  # Specify id.vars: the variables to keep but not split apart on
  pgenes <- melt(gmt, id.vars=c("V1"))
  pgenes$variable <- NULL
  colnames(pgenes) <- c("pathway","gene_symbol")
}
getDE <- function(de, fdr=0.01, logFC=0.5) {
  # return(de[de$FDR <= fdr & abs(de$logFC)>=logFC & de$gene_symbol %in% pgenes$gene_symbol,])
  return(de[de$FDR <= fdr & abs(de$logFC)>=logFC,])
}

getnonDE <- function(de, fdr=0.1) {
  # return(de[de$FDR > fdr & de$gene_symbol %in% pgenes$gene_symbol,])
  return(de[de$FDR > fdr,])
}
#### focus on holiday rather than treated
pw_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/pathway/DE-nocl/'
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/pathway/'
datatag <- 'SA609'
plot_reversed_holiday_genes(fileTvsU=paste0(pw_dir,"scran-nocl-SA609-2_SA609_UTT_UUU-SA609-SA609X5XB03230-SA609X5XB03223_logfc_results.csv"),
                            fileTvsH=paste0(pw_dir,"scran-nocl-SA609-12_SA609_UTT_UTU-SA609-SA609X5XB03230-SA609X5XB03231_logfc_results.csv"),
                            fileHvsU=paste0(pw_dir,"scran-nocl-SA609-22_SA609_UTU_UUU-SA609-SA609X5XB03231-SA609X5XB03223_logfc_results.csv"),
                            title="Pt4 X5", save_dir, datatag, 'Pt4')
# plot_reversed_holiday_genes(fileTvsU=paste0(pw_dir,"../results/SA609-v6/comps2/scran-nocl-SA609-2_SA609_UTT_UUU-SA609-SA609X5XB03230-SA609X5XB03223_logfc_results.csv",
#                                             fileTvsH="../results/SA609-v6/comps2/scran-nocl-SA609-12_SA609_UTT_UTU-SA609-SA609X5XB03230-SA609X5XB03231_logfc_results.csv",
#                                             fileHvsU="../results/SA609-v6/comps2/scran-nocl-SA609-22_SA609_UTU_UUU-SA609-SA609X5XB03231-SA609X5XB03223_logfc_results.csv",
#                                             title="Pt4 X5")
datatag <- 'SA535'
# plot_reversed_holiday_genes(fileTvsU="../results/SA535-v7/comps2/scran-nocl-SA535-3_SA535_UUTTT_UUUUU-SA535-SA535X8XB03431-SA535X8XB03664_logfc_results.csv",
#                             fileTvsH="../results/SA535-v7/comps2/scran-nocl-SA535-13_SA535_UUTTT_UUTTU-SA535-SA535X8XB03431-SA535X8XB03434_logfc_results.csv",
#                             fileHvsU="../results/SA535-v7/comps2/scran-nocl-SA535-23_SA535_UUTTU_UUUUU-SA535-SA535X8XB03434-SA535X8XB03664_logfc_results.csv",
#                             title="Pt5 X8")
plot_reversed_holiday_genes(fileTvsU=paste0(pw_dir,"scran-nocl-SA535-3_SA535_UUTTT_UUUUU-SA535-SA535X8XB03431-SA535X8XB03664_logfc_results.csv"),
                            fileTvsH=paste0(pw_dir,"scran-nocl-SA535-13_SA535_UUTTT_UUTTU-SA535-SA535X8XB03431-SA535X8XB03434_logfc_results.csv"),
                            fileHvsU=paste0(pw_dir,"scran-nocl-SA535-23_SA535_UUTTU_UUUUU-SA535-SA535X8XB03434-SA535X8XB03664_logfc_results.csv"),
                            title="Pt5 X8", save_dir, datatag, 'Pt5')
datatag <- 'SA1035'
plot_reversed_holiday_genes(fileTvsU=paste0(pw_dir,"scran-nocl-SA1035-22_SA1035_UTU_UUU-SA1035-SA1035X6XB03209-SA1035X6XB03216_logfc_results.csv"),
                            fileTvsH=paste0(pw_dir,"scran-nocl-SA1035-23_SA1035_UTTU_UUUU-SA1035-SA1035X7XB03340-SA1035X7XB03502_logfc_results.csv"),
                            fileHvsU=paste0(pw_dir,"scran-nocl-SA1035-24_SA1035_UTTTU_UUUUU-SA1035-SA1035X8XB03420-SA1035X8XB03631_logfc_results.csv"),
                            title= "Pt6 X8", save_dir, datatag, 'Pt6')
# plot_reversed_holiday_genes(fileTvsU="../results/SA1035-v6/comps2/scran-nocl-SA1035-22_SA1035_UTU_UUU-SA1035-SA1035X6XB03209-SA1035X6XB03216_logfc_results.csv",
#                             fileTvsH="../results/SA1035-v6/comps2/scran-nocl-SA1035-23_SA1035_UTTU_UUUU-SA1035-SA1035X7XB03340-SA1035X7XB03502_logfc_results.csv",
#                             fileHvsU="../results/SA1035-v6/comps2/scran-nocl-SA1035-24_SA1035_UTTTU_UUUUU-SA1035-SA1035X8XB03420-SA1035X8XB03631_logfc_results.csv",
#                             title= "Pt6 X8")