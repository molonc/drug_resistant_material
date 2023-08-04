calc_mode <- function(x) {
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}

get_library_id <- function(cell_ids, cores_use=2) {
  
  labels <- mclapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  }, mc.cores = cores_use)
  return(as.character(labels))
  
}
get_gene_coordinates <- function(segments){ #, min_pc_overlap=0.1
  # segments$chr <- paste0('chr',chrs_df$chr)
  # library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  g <- genes(txdb, single.strand.genes.only=FALSE)
  
  segments_gr <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
  overlaps <- findOverlaps(g, segments_gr, ignore.strand = TRUE)
  overlapping_ranges <- pintersect(g[queryHits(overlaps)], segments_gr[subjectHits(overlaps)], ignore.strand = TRUE, drop.nohit.ranges = TRUE)
  
  percentage_overlap <- width(overlapping_ranges) / width(segments_gr[subjectHits(overlaps)])
  # wd_overlap <- width(overlapping_ranges)
  # percentage_overlap <- width(overlapping_ranges) / width(g[queryHits(overlaps)])
  
  percentage_overlap_sum <- sapply(percentage_overlap, sum)
  # wd_overlap_sum <- sapply(wd_overlap, sum)
  gene_cn <- tibble(entrezgene = names(g)[queryHits(overlaps)], percent_overlap=percentage_overlap_sum) #wd_overlap=wd_overlap_sum
  # print(dim(gene_cn))
  # length(unique(gene_cn$entrezgene))
  # t <- gene_cn[duplicated(gene_cn$entrezgene),]
  # dim(t)
  # head(t)
  # t1 <- gene_cn[gene_cn$entrezgene=='10000',]
  gene_cn <- gene_cn %>%
    cbind(mcols(segments_gr[subjectHits(overlaps)]))
  # head(gene_cn)
  # dim(gene_cn)
  # length(unique(gene_cn$entrezgene))
  ext_rows <- gene_cn %>% # Remove the cases where a given gene that map to several regions.
    dplyr::group_by(entrezgene, cluster) %>%
    dplyr::summarise(percent_overlap=max(percent_overlap))%>% 
    dplyr::mutate(desc=paste0(entrezgene, cluster, percent_overlap))%>%
    dplyr::ungroup()%>%
    dplyr::pull(desc)
  # length(ext_rows)
  gene_cn <- gene_cn %>%
    dplyr::mutate(desc=paste0(entrezgene, cluster, percent_overlap))%>%
    dplyr::filter(desc %in% ext_rows)%>%
    dplyr::select(-desc)%>%
    dplyr::mutate(percent_overlap=round(percent_overlap, 4))
  
  
  # dim(gene_cn)
  gene_cn$ensembl_gene_id <- mapIds(org.Hs.eg.db,
                                    keys = gene_cn$entrezgene,
                                    column="ENSEMBL",
                                    keytype="ENTREZID",
                                    multiVals="first") #"CharacterList"
  # ens_mapped <- mapIds(org.Hs.eg.db,
  #             keys = gene_cn$entrezgene,
  #             column="ENSEMBL",
  #             keytype="ENTREZID",
  #             multiVals="list")
  
  ## Other way to do it
  # chrs <- c(as.character(1:22), "X")
  # t <- annotables::grch38 %>%
  #   dplyr::select(ensembl_gene_id = ensgene, entrezgene=entrez)  %>%
  #   dplyr::filter(entrezgene %in% gene_cn$entrezgene) #  & chr %in% chrs
  # #   inner_join(deg_df) %>%
  
  
  # ens_mapped <- unlist(ens_mapped)
  # ens_mapped_out <- data.frame(entrezgene=names(ens_mapped), ensembl_gene_id=ens_mapped)
  # ens_mapped_out <- ens_mapped_out %>%
  #   dplyr::filter(!is.na(ensembl_gene_id))
  # gene_cn <-  gene_cn %>% inner_join(ens_mapped_out, by='entrezgene')
  # print(dim(gene_cn))
  
  # Filter for complete entries
  gene_cn <- gene_cn %>%
    as.data.frame %>%
    drop_na()
  
  # gene_cn <- gene_cn %>%
  #   dplyr::filter(percent_overlap>min_pc_overlap)
  # print(dim(gene_cn))
  # gene_cn <- annotables::grch37 %>% 
  #   dplyr::select(ensembl_gene_id = ensgene, symbol) %>% 
  #   inner_join(gene_cn)
  # print(dim(gene_cn))
  # gene_cn$percent_overlap <- round(gene_cn$percent_overlap, 4)
  return(gene_cn)
}

get_mapping_genes <- function(segments){
  if(!grepl('chr',segments$chr[1])){
    segments <- segments %>%
      dplyr::mutate(chr=paste0("chr", chr))
  }
  unassigned_clone_lbs <- paste0('clone_',c('None','Unassigned', 'Un','un'))
  segments <- segments %>%
    # dplyr::filter(clone_id =='clone_D') %>%
    dplyr::filter(!clone_id %in% unassigned_clone_lbs) %>%
    dplyr::rename(cluster=clone_id)#%>%
  # dplyr::select(chr,start,end,cluster,copy_number,pct_pure)
  gene_cn <- get_gene_coordinates(segments)
  return(gene_cn)
}
