#' This is the entry point for the pipeline
#' 
#' 

# MA: this is similar to fix_sce
## Get rid of stuff that breaks rowRanges
rd_to_rename <- c("start", "end", "strand")
unbreak_rowranges <- function(sce) {
  for(r in rd_to_rename) {
    if(r %in% names(rowData(sce))) {
      rowData(sce)[[r]] <- NULL
    }
  }
  sce
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(argparse)
})

parser <- ArgumentParser(description = "Integrate SCEs")
parser$add_argument('--sce_human', metavar='FILE', type='character', 
                    help="Path to human SCE")

parser$add_argument('--sce_mouse', metavar='FILE', type='character', 
                    help="Path to mouse SCE")

parser$add_argument('--sample', type='character',
                    help="Sample name")

parser$add_argument('--output_figure', metavar='FILE', type = 'character',
                    help="Output figure file")

parser$add_argument('--output_csv', metavar='FILE', type = 'character',
                    help="Output dataframe with summary of # mouse cells")
parser$add_argument('--output_sce', type = 'character', metavar = 'FILE',
                    help="Output path for SCE with mouse cells removed")
print("0")
args <- parser$parse_args()

sce_human <- readRDS(args$sce_human)
# MA: validate this input -- to move this to a separate script or rule
# 1. Make sure it doesn't contain col names that break rowranges 
sce_human <- unbreak_rowranges(sce_human)

# 2. Remove the columns id, series, material from here ()
colrm <- intersect(names(colData(sce_human)), c("id", "series", "material"))
colData(sce_human) <- colData(sce_human)[,!(names(colData(sce_human)) %in% colrm)]

# 3. Sometimes the rowData(sce)$ID has NA values. Remove those rows.
sce_human <- sce_human[which(!is.na(rowData(sce_human)$ID)),]

sce_human_orig <- sce_human

## Check if sce_mouse is missing and write output files if so
if(file.info(args$sce_mouse)$size==0) {
    # Just write SCE human
    #sce_human <- unbreak_rowranges(sce_human)
    
    sce_human$is_mouse <- FALSE
    
    ggplot()
    ggsave(args$output_figure, width = 6, height = 5)
    
    tibble() %>% 
      write_csv(args$output_csv)
    saveRDS(sce_human, file = args$output_sce)
    
    quit()
} else {

  sce_mouse <- readRDS(args$sce_mouse)
  sce_mouse <- unbreak_rowranges(sce_mouse)
  
  # Used to check this
  ## if (is.character(sce_mouse) & sce_mouse == "empty")) 
  
  if(nrow(sce_human) > 40000) {
    stop("Human SCE has > 40k genes: is this aligned to mouse hybrid?")
  }
  if(nrow(sce_mouse) > 40000) {
    stop("Mouse SCE has > 40k genes: is this aligned to mouse hybrid?")
  }
  
  
  
  sce_human <- calculateQCMetrics(sce_human)
  sce_mouse <- calculateQCMetrics(sce_mouse)
  
  get_metric_df <- function(sce) {
    as.data.frame(colData(sce)) %>% 
      as_tibble() %>% 
      dplyr::select(Barcode, total_counts)
  }
  
  df_hm <- inner_join(
    get_metric_df(sce_human),
    get_metric_df(sce_mouse),
    by = "Barcode",
    suffix = c("_human", "_mouse")
  )
  
  df_hm <- dplyr::mutate(df_hm, is_mouse = total_counts_mouse > total_counts_human)
  
  plt <- ggplot(df_hm, aes(x = total_counts_human, y = total_counts_mouse, colour = is_mouse)) + 
    geom_point() +
    geom_abline()
  
  ggsave(args$output_figure, plt, width = 6, height = 5)
  
  mouse_bcs <- dplyr::filter(df_hm, is_mouse) %>% .$Barcode
  
  # MA: removing reading again so I keep the corrections at the top
  # sce_human <- readRDS(args$sce_human)
  sce_human <- sce_human_orig
  
  sce_human$is_mouse <- sce_human$Barcode %in% mouse_bcs
  
  tbl <- table(sce_human$is_mouse)
  
  as.data.frame(tbl) %>% 
    dplyr::rename(is_mouse = Var1, n_cells = Freq) %>% 
    dplyr::mutate(sample = args$sample) %>% 
    dplyr::select(sample, everything()) %>% 
    write_csv(args$output_csv)
  
  sce_human <- sce_human[, !sce_human$is_mouse]
  
  # Write outputs
  saveRDS(sce_human, file = args$output_sce)
  
  cat("Completed.\n")
}  
