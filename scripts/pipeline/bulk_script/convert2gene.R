

## How to use this script
## 
## First install dependent packages: 
## BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
## BiocManager::install("dplyr")
## BiocManager::install("tximport")
## BiocManager::install("annotables")
## BiocManager::install("argparse")

## How to run script from command line:
## Option 1: 
## Do not specify --txt2gene_fn parameter, using annotable package reference
## Rscript yourdir/convert2gene.R --abundance_fn yourdir/convert2geneLevel/abundance.tsv --output_fn yourdir/convert2geneLevel/abundance2gene_using_annotables_package.csv.gz --is_mouse TRUE 

## Option 2: 
## Specify --txt2gene_fn parameter, using our build package from file
## Rscript yourdir/convert2gene.R --abundance_fn yourdir/convert2geneLevel/abundance.tsv --txt2gene_fn yourdir/convert2geneLevel/reference_tx2gene/MouseEnstx2gene.csv.gz --output_fn yourdir/convert2geneLevel/abundance2gene_using_Elena_ref.csv.gz --is_mouse TRUE 



suppressMessages({
  library("dplyr")
  library("tximport")
  library("annotables")
  library("argparse")
})
parser <- ArgumentParser(description = "Convert transcripts into gene level")
parser$add_argument('--abundance_fn', metavar='FILE', type='character',
                    help="Path to abundance file from kallisso")
parser$add_argument('--txt2gene_fn', metavar='FILE', type='character',default = NULL,
                    help="Path to reference gene set, if NULL using default setting, load a reference gene set from annotables package")
parser$add_argument('--is_mouse', type = 'character', default = 'TRUE',# mouse or human reference?
                    help="TRUE or FALSE, if not include, using default setting, as TRUE")
parser$add_argument('--output_fn', metavar='FILE', type = 'character',default = '',
                    help="Output file name, if not include, using default setting")

args <- parser$parse_args()



load_txt2genes <- function(txt2gene_fn){
  # BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
  # library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  # BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
  library("TxDb.Mmusculus.UCSC.mm10.ensGene")
  txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
  
  ## Return 142446 transcripts, and 40107 empty gene id
  # k <- keys(txdb, keytype = "TXNAME")
  # tx2gene <- select(txdb, k, "GENEID", "TXNAME")
  
  ## Return 102339 transcripts, and with available gene id
  k <- keys(txdb, keytype = "GENEID")
  tx2gene <- select(txdb, keys = k,  columns = "TXNAME", keytype = "GENEID")
  
  
  print(dim(tx2gene))
  print(head(tx2gene))
  print(sum(is.na(tx2gene$GENEID)))
  tx2gene <- tx2gene %>%
    dplyr::select(TXNAME, GENEID)
  data.table::fwrite(tx2gene, txt2gene_fn)
  return(tx2gene)
}

load_tx2gene_from_file <- function(load_from_file=TRUE, txt2gene_fn=NULL, mouse_reference=TRUE){
  if(load_from_file & !is.null(txt2gene_fn)){
    ## Loading from the our in-house built tx2gene file, Elena: version 37, Hoa: version 38
    # tx2gene <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene.csv', header=T) %>% as.data.frame()
    tx2gene <- data.table::fread(txt2gene_fn, header=T)
    if('V1' %in% colnames(tx2gene)){
      tx2gene$V1 <- NULL  
    }
    
    print(head(tx2gene))
    print(dim(tx2gene))
    tx2gene <- tx2gene %>%
      dplyr::select(TXNAME, GENEID)
    return(tx2gene)  
  }else{ #Loading from the latest version of tx2gene from annotables reference package
    # library("annotables")
    
    if(mouse_reference){
      print('Using mouse reference gene set grcm38_tx2gene')
      tx2gene <- grcm38_tx2gene
      # sum(is.na(tx2gene$enstxp))
      print(dim(tx2gene))
      colnames(tx2gene) <- c("TXNAME","GENEID")
    }else{ # 
      print('Using human reference gene set grch38_tx2gene')
      tx2gene <- grch38_tx2gene
      colnames(tx2gene) <- c("TXNAME","GENEID")
      print(dim(tx2gene))
    }
    
    return(tx2gene)
  }
  
  
}

load_kallisto_data <- function(bulk_fns, output_fn, tx2gene=NULL, mouse_reference=TRUE){
  for(f in bulk_fns){
    if(!file.exists(f)){
      print(f)
      stop('File do not exist')
      
    }
  }
  if(is.null(tx2gene)){ ## Using build version from package annotables in this case 
    tx2gene <- load_tx2gene(load_from_file=FALSE, txt2gene_fn=NULL)
  }
  # head(tx2gene)
  print(dim(tx2gene))
  txi <- tximport(bulk_fns, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=T) #, ignoreAfterBar = TRUE
  print(dim(txi$counts))
  
  # Get genes symbols here
  chrs <- c(as.character(1:22), "X") # don't take into account chr Y
  if(mouse_reference){
    print('Using mouse reference gene set to convert ensemble gene id to gene symbol')
    annots <- annotables::grcm38 %>%
      dplyr::select(ensembl_gene_id = ensgene, gene_symbol=symbol,chr) %>%
      dplyr::filter(chr %in% chrs)  
  }else{ # 
    print('Using human reference gene set to convert ensemble gene id to gene symbol')
    annots <- annotables::grch38 %>% # or annotables::grch37 version
      dplyr::select(ensembl_gene_id = ensgene, gene_symbol=symbol,chr) %>%
      dplyr::filter(chr %in% chrs)  
  }
  
  annots <- annots[!duplicated(annots$ensembl_gene_id),]
  # colnames(annots)
  
  counts_values <- as.data.frame(txi$counts)
  gn <- sapply(strsplit(rownames(counts_values),'\\.'),function(x){ # remove version for name, to convert to gene symbol
    return(x[1])
  })
  
  counts_values$ensembl_gene_id <- as.character(gn)
  counts_values <- counts_values %>% left_join(annots, by='ensembl_gene_id')
  
  counts_values <- counts_values %>%
    dplyr::rename(bulk_exp=V1)
  print(head(counts_values))
  print(dim(counts_values))
  data.table::fwrite(counts_values, output_fn)
  print('Saved output as: ')
  print(output_fn)
  # return(TRUE)
  
}


## For human reference 
# txt2gene_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene.csv'
# txt2gene_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/tx2gene_updated_06_Dec_2021.csv.gz'



main <- function(abundance_fn, output_fn='', is_mouse=TRUE, txt2gene_fn=NULL){
  if(output_fn==''){
    output_fn <- paste0(dirname(abundance_fn),'/','abundance2gene.csv.gz')
  }
  ## Parameter setting
  # input_dir <- '/home/htran/storage/images_dataset/merfish_data/bulk/convert2geneLevel/'
  ## Step 1: First building a tx2gene file, and save into folder
  # txt2gene_fn <- paste0(input_dir,'MouseEnstx2gene.csv.gz')
  # tx2gene <- load_txt2genes(txt2gene_fn) 
  # head(tx2gene)
  
  ## After that, for future run, we can just load this file, skip step 1 above
  ## Step2: Loading file 
  if(is_mouse==TRUE){
    mouse_reference <- TRUE
  }else{
    mouse_reference <- FALSE
  }
  if(is.null(txt2gene_fn)){
    tx2gene <- load_tx2gene_from_file(load_from_file=FALSE, txt2gene_fn=NULL, mouse_reference)  
  }else{
    tx2gene <- load_tx2gene_from_file(load_from_file=TRUE, txt2gene_fn, mouse_reference)
  }
  
  # print(sum(is.na(tx2gene$GENEID)))
  ## convert transcript id to gene id
  counts_genes <- load_kallisto_data(abundance_fn, output_fn, tx2gene, mouse_reference)
  print('Complete!!!')
}


# abundance_fn <- paste0(input_dir, 'abundance.tsv') 
# output_fn <- paste0(save_dir, 'abundance2gene.csv.gz')
if(args$is_mouse=='TRUE'){
  is_mouse=TRUE
}else{
  is_mouse=FALSE
}
main(args$abundance_fn, args$output_fn, is_mouse, args$txt2gene_fn)






