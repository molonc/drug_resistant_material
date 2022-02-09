
# make_clone_palette <- function(levels) {
#   if (length(levels) <= 8) {
#     pal <- RColorBrewer::brewer.pal(max(length(levels), 3), "Dark2")
#   } else {
#     pal <- colorRampPalette(RColorBrewer::brewer.pal(max(length(levels), 3), "Dark2"))(length(levels))
#   }
#   names(pal) <- levels
#   pal <- pal[levels]
#   return(pal)
# }

cn <- ascn
tree = tree_cg
clusters = clones_df
normalize_ploidy = FALSE
normalize_tree = FALSE
branch_length = 1
spacer_cols=15
plottree = TRUE
plotcol = "state_phase"
reorderclusters = FALSE
pctcells = 0.05
library_mapping = NULL
clone_pal = NULL
sample_label_idx = 1
fillna = TRUE
frequencycutoff = 2
maxf = NULL
plotfrequency = FALSE
show_legend = TRUE
show_library_label = TRUE
plotHeatmap <- function(cn,
                        tree = NULL,
                        clusters = NULL,
                        normalize_ploidy = FALSE,
                        normalize_tree = FALSE,
                        branch_length = 1,
                        spacer_cols=20,
                        plottree = TRUE,
                        plotcol = "state",
                        reorderclusters = FALSE,
                        pctcells = 0.05,
                        library_mapping = NULL,
                        clone_pal = NULL,
                        sample_label_idx = 1,
                        fillna = TRUE,
                        frequencycutoff = 2,
                        maxf = NULL,
                        plotfrequency = FALSE,
                        show_legend = TRUE,
                        show_library_label = TRUE,
                        ...){
  
  if (is.hscn(cn) | is.ascn(cn)){
    CNbins <- cn$data
  } else{
    CNbins <- cn
  }
  
  if (!plotcol %in% c("state", "state_BAF", "state_phase", "state_AS", "state_min")){
    stop(paste0("Column name - ", plotcol, " not available for plotting, please use one of state, state_BAF, state_phase, state_AS or state_min"))
  }
  
  if (!plotcol %in% names(CNbins)){
    stop(paste0("Column name - ", plotcol, " not in CNbins data frame..."))
  }
  
  if (plotcol == "state"){
    colvals <- schnapps:::cn_colours
    legendname <- "Copy Number"
  }
  
  if (plotcol == "state_BAF"){
    colvals <- schnapps:::cn_colours_bafstate
    legendname <- "Allelic Imbalance"
  }
  
  if (plotcol == "state_AS"){
    colvals <- schnapps:::cn_colours_loh
    legendname <- "Allele Specific Copy Number"
  }
  
  if (plotcol == "state_min"){
    colvals <- schnapps:::cn_colours_minorallele
    legendname <- "Minor Allele Copy Number"
  }
  
  if (plotcol == "state_phase"){
    colvals <- schnapps:::cn_colours_phase
    legendname <- "Allelic Imbalance"
  }
  
  ncells <- length(unique(CNbins$cell_id))
  
  if (is.null(clusters)){
    ordered_cell_ids <- paste0(unique(CNbins$cell_id))
  } else{
    ordered_cell_ids <- paste0(clusters$cell_id)
  }
  
  if (is.null(tree) & is.null(clusters)){
    message("No tree or cluster information provided, clustering using HDBSCAN")
    clustering_results <- schnapps::umap_clustering(CNbins, minPts = max(round(pctcells * ncells), 2), field = "copy")
    tree <- clustering_results$tree
    tree_ggplot <- schnapps:::make_tree_ggplot(tree, as.data.frame(clustering_results$clusters), clone_pal = clone_pal)
    tree_plot_dat <- tree_ggplot$data
    message("Creating tree...")
    tree_hm <- schnapps:::make_corrupt_tree_heatmap(tree_ggplot)
    ordered_cell_ids <- schnapps:::get_ordered_cell_ids(tree_plot_dat)
    
    clusters <- clustering_results$clustering %>%
      dplyr::select(cell_id, clone_id)
  }
  
  if (!is.null(clusters)){
    if (!"clone_id" %in% names(clusters)){
      stop("No clone_id columns in clusters dataframe, you might need to rename your clusters")
    }
    if (reorderclusters == TRUE){
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
      ordered_cell_ids <- paste0(clusters$cell_id)
    }
  }
  
  if (plottree == TRUE){
    if (normalize_tree == T){
      tree <- schnapps:::format_tree(tree, branch_length)
    }
    
    tree_ggplot <- schnapps:::make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal)
    tree_plot_dat <- tree_ggplot$data
    
    message("Creating tree...")
    tree_hm <- schnapps:::make_corrupt_tree_heatmap(tree_ggplot)
    ordered_cell_ids <- schnapps:::get_ordered_cell_ids(tree_plot_dat)
  }
  
  message("Creating copy number heatmap...")
  copynumber <- schnapps:::createCNmatrix(CNbins, field = plotcol, fillna = fillna)
  if (normalize_ploidy == T){
    message("Normalizing ploidy for each cell to 2")
    copynumber <- schnapps:::normalize_cell_ploidy(copynumber)
  }
  copynumber <- schnapps:::format_copynumber(copynumber, ordered_cell_ids, spacer_cols = spacer_cols)
  clones_formatted <- schnapps:::format_clones(as.data.frame(clusters), ordered_cell_ids)
  if (!is.null(clone_pal)){
    clones_idx <- dplyr::distinct(clones_formatted, clone_id, clone_label)
    clone_pal <- clone_pal[clones_idx$clone_id]
    names(clone_pal) <- clones_idx$clone_label
  }
  copynumber_hm <- schnapps:::make_copynumber_heatmap(copynumber,
                                           clones_formatted,
                                           colvals = colvals,
                                           legendname = legendname,
                                           library_mapping = library_mapping,
                                           clone_pal = clone_pal,
                                           sample_label_idx = sample_label_idx,
                                           cutoff = frequencycutoff,
                                           maxf = maxf,
                                           plotcol = plotcol,
                                           plotfrequency = plotfrequency,
                                           show_legend = show_legend,
                                           show_library_label = show_library_label,
                                           ...)
  copynumber_hm <- schnapps:::make_copynumber_heatmap(copynumber,
                                                      clones_formatted,
                                                      colvals = colvals,
                                                      legendname = legendname,
                                                      library_mapping = library_mapping,
                                                      clone_pal = clone_pal,
                                                      sample_label_idx = sample_label_idx,
                                                      cutoff = frequencycutoff,
                                                      maxf = maxf,
                                                      plotcol = plotcol,
                                                      plotfrequency = plotfrequency,
                                                      show_legend = show_legend,
                                                      show_library_label = show_library_label)
  
  copynumber_hm <- ComplexHeatmap::Heatmap(name = legendname, 
                                           as.matrix(copynumber), col = colvals, na_col = "white", 
                                           show_row_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, 
                                           show_column_names = FALSE, bottom_annotation = schnapps:::make_bottom_annot(copynumber), 
                                           left_annotation = make_left_annot(copynumber, clones=clones_formatted, 
                                                                             library_mapping = library_mapping, clone_pal = clone_pal, 
                                                                             idx = sample_label_idx, show_legend = show_legend, 
                                                                             show_library_label = show_library_label), heatmap_legend_param = list(nrow = 4), 
                                           top_annotation = make_top_annotation_gain(copynumber, 
                                                                                     cutoff = cutoff, maxf = maxf, plotfrequency = plotfrequency, 
                                                                                     plotcol = plotcol), use_raster = TRUE, raster_quality = 5)
  
  if (plottree == TRUE){
    h <- tree_hm + copynumber_hm
  } else {
    h <- copynumber_hm
  }
  
  return(h)
}


make_left_annot <- function (copynumber, clones, library_mapping = NULL, show_library_label = TRUE, 
          clone_pal = NULL, idx = 1, show_legend = TRUE) 
{
  annot_colours <- list()
  library_labels <- schnapps:::get_library_labels(rownames(copynumber), 
                                       idx = idx)
  if (!is.null(library_mapping)) {
    library_labels <- unlist(library_mapping[library_labels])
    if (!all(library_labels %in% names(library_mapping)) == 
        FALSE) {
      warning("Not all library ids present in library to name mapping, using library IDs as annotations...")
      library_labels <- schnapps:::get_library_labels(rownames(copynumber))
    }
  }
  library_levels <- gtools::mixedsort(unique(library_labels))
  annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
  annot_colours$Sample <- annot_colours$Sample[!is.na(annot_colours$Sample)]
  library_legend_rows <- 10
  if (!is.null(clones)) {
    clone_levels <- unique(clones$clone_label)
    clone_level_none <- clone_levels[grepl("None", clone_levels)]
    clone_levels <- gtools::mixedsort(clone_levels[!grepl("None", 
                                                          clone_levels)])
    if (is.null(clone_pal)) {
      clone_pal <- schnapps:::make_clone_palette(clone_levels)
    }
    if (length(clone_level_none > 0)) {
      clone_pal[[clone_level_none]] <- clone_none_black
    }
    annot_colours$Clone <- clone_pal
    clone_label_generator <- function(index) {
      clone_label_pos <- get_clone_label_pos(clones)
      y_pos <- 1 - unlist(clone_label_pos)/nrow(clones)
      grid::grid.text(names(clone_label_pos), 0.5, y_pos, 
                      just = c("centre", "centre"))
    }
    clone_legend_rows <- 10
    if (length(clone_levels) > 10) {
      clone_legend_rows <- round(sqrt(length(clone_levels) * 
                                        4))
    }
    if (show_library_label) {
      left_annot <- ComplexHeatmap::HeatmapAnnotation(Clone = clones$clone_label, 
                                                      clone_label = clone_label_generator, Sample = library_labels, 
                                                      col = annot_colours, show_annotation_name = c(TRUE, 
                                                                                                    FALSE, TRUE), which = "row", annotation_width = grid::unit(rep(0.4, 
                                                                                                                                                                   3), "cm"), annotation_legend_param = list(Clone = list(nrow = clone_legend_rows), 
                                                                                                                                                                                                             Sample = list(nrow = library_legend_rows)), 
                                                      show_legend = show_legend)
    }
    else {
      left_annot <- ComplexHeatmap::HeatmapAnnotation(Clone = clones$clone_label, 
                                                      clone_label = clone_label_generator, col = annot_colours, 
                                                      show_annotation_name = c(TRUE, FALSE), which = "row", 
                                                      annotation_width = grid::unit(rep(0.4, 2), "cm"), 
                                                      annotation_legend_param = list(Clone = list(nrow = clone_legend_rows)), 
                                                      show_legend = show_legend)
    }
  }
  else {
    left_annot <- ComplexHeatmap::HeatmapAnnotation(Sample = library_labels, 
                                                    col = annot_colours, which = "row", simple_anno_size = grid::unit(0.4, 
                                                                                                                      "cm"), annotation_legend_param = list(Sample = list(nrow = library_legend_rows), 
                                                                                                                      ), show_legend = show_legend)
  }
  return(left_annot)
}

