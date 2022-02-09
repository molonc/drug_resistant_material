# pkgs <- c("factoextra",  "NbClust","Rfast") 
# install.packages(pkgs)
# BiocManager::install("Rfast")
# library(factoextra)
# library(NbClust)
# library(Rfast)
suppressPackageStartupMessages({
    library(tidyverse)
    library(dplyr)
    library(ggplot2)
    library(SingleCellExperiment)
    library(WGCNA)
})

source('/home/htran/Projects/farhia_project/rnaseq/method_testing/network_utils.R')
# Convert from 1 vector to encoding feature 0 1 0, 0 0 1,...
# Input: meta_cells_df should contain a column 'treatment_status'
encode_metadata <- function(meta_cells_df, save_dir, datatag){
  # BiocManager::install("vtreat")
  # library(vtreat)
  tz <- vtreat::designTreatmentsZ(meta_cells_df, c("treatment_status"))
  new_df <- vtreat::prepare(tz, meta_cells_df)
  # View(head(new_df))
  # sum(new_df$treatment_status_lev_x_X6_UUT==1)  
  print(dim(new_df))
  print(colnames(new_df))
  names(new_df) <- gsub('treatment_status_lev_x_','',names(new_df))
  new_df <- new_df %>%
            dplyr::select(-treatment_status_catP)
  new_df$cell_id <- meta_cells_df$cell_id
  write.csv(new_df,paste0(save_dir, datatag,'_meta_cells.csv'), quote = F, row.names = F)
  rownames(new_df) <- new_df$cell_id
  new_df$cell_id <- NULL
  write.csv(meta_cells_df,paste0(save_dir, datatag,'_meta_cells_treated_df.csv'), quote = F, row.names = F)
  return(new_df)
}
remove_outliers <- function(datExpr, save_dir, datatag, tag){
  # rn <- rownames(datExpr)
  # datExpr <- apply(datExpr, 2, as.numeric)
  # rownames(datExpr) <- rn
  print(dim(datExpr))
  # Run this to check if there are gene outliers
  gsg = WGCNA::goodSamplesGenes(datExpr, verbose = 3)
  print(gsg$allOK)
  if(gsg$allOK){
    print('There is no outlier cells and outlier genes')
  }
  
  # If the last statement returns TRUE, all genes have passed the cuts. 
  # If not, we remove the offending genes and samples from the data with the following:
  outlier_genes <- NULL
  outlier_cells <- NULL
  if(!gsg$allOK)
  {
    if(sum(!gsg$goodGenes)>0){
      outlier_genes <- names(datExpr)[!gsg$goodGenes]
      print(paste0("Number of outlier genes: ", length(outlier_genes)))
      printFlush(paste("Removing genes:", paste(outlier_genes, collapse= ", ")));
    }
      
    if(sum(!gsg$goodSamples)>0){
      outlier_cells <- rownames(datExpr)[!gsg$goodSamples]
      print(paste0("Number of outlier cells: ", length(outlier_cells)))
      printFlush(paste("Removing samples:", paste(outlier_cells, collapse=", ")))
    }
        
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    print(dim(datExpr))
  }
  
  data.table::fwrite(x = datExpr, paste0(save_dir, datatag,'_',tag,'_datExpr.csv.gz'),row.names = T)
  # datExpr <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_datExpr.csv.gz'))
  # datExpr <- data.frame(datExpr)
  # rownames(datExpr) <- as.character(datExpr$V1)
  # datExpr$V1 <- NULL
  
  return(list(outlier_genes=outlier_genes, outlier_cells=outlier_cells,
              datExpr=datExpr))
  
}
select_soft_threshold_power <- function(datExpr, save_dir, option=1){
  # library(WGCNA)
  # library(flashClust)
  WGCNA::enableWGCNAThreads()
  if(option==1){
    powers = c(c(1:10), seq(from=12, to=20, by=2)) #choosing a set of soft-thresholding powers
    sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function
    
    png(paste0(save_dir, "soft_threshold_power.png"), height =2*400, width= 2*900,res = 2*72)
      # sizeGrWindow(9,5)
      par(mfrow= c(1,2))
      cex1=0.9
      p <- plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
                xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", 
                type= "n", main= paste("Scale independence"))
      p <- p + text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
      p <- p + abline(h=0.90, col="red")
      p1 <- plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
      p1 <- p1 + text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
      # summary_plt <- cowplot::plot_grid(p, p1, ncol = 2, align = 'h')
      # print(summary_plt)
    dev.off()
    
  }else{
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    png(paste0(save_dir, "soft_threshold_power.png"), height =2*400, width= 2*900,res = 2*72)
      p <- plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
                type='n', main = paste('Scale independence'))
      p <- p + text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                    labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
      print(p)
    dev.off()
    
  }
  
  
}

extract_genes_correlation_network <- function(results_mtx, results, datExpr, save_dir, datatag, softPowerThres=8, cluster_use='blue'){
  # Recalculate topological 
  # overlapTOM = WGCNA::TOMsimilarityFromExpr(datExpr, power = softPowerThres)
  # Read in the annotation file
  probeToGene = data.frame(results$gene_cluster$ens_gene, results$gene_cluster$gene)

  # Select module
  # module = "brown";# Select module probes
  probes = names(datExpr)
  inModule = (results$dynamicColors==cluster_use)  #moduleColors
  modProbes = probes[inModule]
  # Select the corresponding Topological Overlap
  adjMat = results_mtx$TOM[inModule, inModule]
  dimnames(adjMat) = list(modProbes, modProbes)
  # Export the network into an edge list file VisANT can read, original code
  # vis = WGCNA::exportNetworkToVisANT(modTOM, file = paste(save_dir, "VisANTInput-", cluster_use, ".txt", sep=""),
  #                             weighted = TRUE, threshold = 0,
  #                             probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol))

  adjMat = as.matrix(adjMat)
  adjMat[is.na(adjMat)] = 0
  nRow = nrow(adjMat)
  WGCNA::checkAdjMat(adjMat, min = -1, max = 1)
  probes = dimnames(adjMat)[[1]]
  if (!is.null(probeToGene)) {
    probes2genes = match(probes, probeToGene[, 1])
    if (sum(is.na(probes2genes)) > 0) 
      stop("Error translating probe names to gene names: some probe names could not be translated.")
    probes = probeToGene[probes2genes, 2]
  }
  rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
  colMat = matrix(c(1:nRow), nRow, nRow)
  adjDst = as.dist(adjMat)
  dstRows = as.dist(rowMat)
  dstCols = as.dist(colMat)
  maxNConnections <- 100  # keep the top 100 biggest weight connection between 2 genes
  threshold <- as.numeric(quantile(abs(adjDst),probs=0.5))  # threshold, keep biggest weight connection  
  if (is.null(maxNConnections)) 
    maxNConnections = length(adjDst)
  ranks = rank(-abs(adjDst), na.last = TRUE, ties.method = "first")
  edges = abs(adjDst) > threshold & ranks <= maxNConnections
  nEdges = sum(edges)
  print(paste0('Number of nEdges: ', nEdges))
  weighted = TRUE
  edges_df = data.frame(from = probes[dstRows[edges]], to = probes[dstCols[edges]], 
                          direction = rep(0, nEdges), #method = rep("M0039", nEdges), 
                          weight = if (weighted) 
                            adjDst[edges]
                          else rep(1, nEdges))
  print(dim(edges_df))
  write.csv(edges_df, file = paste0(save_dir, "edges_genes_viz_", cluster_use, ".csv"), quote = FALSE, row.names = FALSE)
  
  return(list(edges_df=edges_df, node_df=results$gene_cluster))
  
}


# node_df <- gene_cluster
# source('/home/htran/Projects/farhia_project/rscript/method_testing/network_utils.R')
construct_graph_v2 <- function(edges_df, node_df, save_dir, datatag, cluster_use){
  
  colnames(node_df)[grepl('logFC',colnames(node_df))] <- 'logFC'
  colnames(node_df)[which(colnames(node_df)=='gene')] <- 'gene_symb'
  # View(head(node_df))
  # node_df$logfc
  rownames(node_df) <- node_df$gene_symb
  g = read_ltm_tree(edges_df)
  
  print('Load String DB...')
  string_net_out_fn <- paste0(save_dir,datatag,'_string_net_stat.rds')
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  string_net_stat <- load_stringDB(string_net_out_fn, input_dir, node_df, save_dir, datatag, 300, F)
  string_network <- string_net_stat$string_network
  protein_infos <- string_net_stat$protein_info
  rownames(protein_infos) <- protein_infos$protein_external_id
  string_connection <- paste0(protein_infos[string_network$protein1,'preferred_name'],'_',
                              protein_infos[string_network$protein2,'preferred_name'])


  is_stringdb <- c('NOT_INCLUDED_STRINGDB','INCLUDED_STRINGDB') #
  names(is_stringdb) <- c(FALSE, TRUE)
  edges_df$is_stringdb <- is_stringdb[as.factor(paste0(edges_df$from,'_', edges_df$to) %in% string_connection)]
  print(summary(edges_df$is_stringdb))
  edges_df$is_stringdb <- factor(edges_df$is_stringdb, levels=c('NOT_INCLUDED_STRINGDB','INCLUDED_STRINGDB'))
  
  
  
  igraph::E(g)$weight <- as.numeric(edges_df$weight) 
  igraph::E(g)$width <- as.numeric(edges_df$weight) 
  # igraph::E(g)$width <- as.numeric(edges_df$width)
  igraph::E(g)$is_stringdb <- as.character(edges_df$is_stringdb)
  
  # igraph::V(g)$size <- node_df[igraph::V(g)$name,'logFC']
  # linetypes <- c(2,1)
  # # edge_attr(g)
  # edge_colrs <- c("grey","black")
  # names(edge_colrs) <- c('NOT_INCLUDED_STRINGDB','INCLUDED_STRINGDB')
  # edge_colrs <- edge_colrs[as.character(unique(edges_df$is_stringdb))]
  
  edge_colrs <- c("NOT_INCLUDED_STRINGDB" = "grey","INCLUDED_STRINGDB" = "black")
  linetypes <- c("NOT_INCLUDED_STRINGDB" = 2,"INCLUDED_STRINGDB" = 1)
  #igraph::E(g)$lty <- linetypes[as.factor(E(g)$is_stringdb)]
  # summary(as.factor(E(g)$is_string))
  # E(g)$color <- edge_colrs[as.factor(E(g)$is_stringdb)]
  edge_layer <- ggraph::geom_edge_link(aes_(color = ~is_stringdb, linetype= ~is_stringdb), show.legend =T, alpha=.7) #
                
  
  
  igraph::V(g)$color <-  node_df[igraph::V(g)$name,'logFC']
  igraph::V(g)$size <- 20
  
  # summary(as.factor(E(g)$is_string))
  # E(g)$color <- edge_colrs[as.factor(E(g)$is_string)]
  
  # E(g)$width <- E(g)$width / 3
  # g <- igraph::set_graph_attr(g, "layout", igraph::layout_with_kk(g))
  width_plot <- 120
  if(width_plot < 80){
    width_plot <- 500
    height_plot <- 500
  }else{
    width_plot <- width_plot * 8
    height_plot <- 500
  }
  
  p <- ggraph::ggraph(g, layout='kk', circular = FALSE) +
    edge_layer +
    ggraph::scale_edge_color_manual(values=edge_colrs) + #name=names(edge_colrs),
    ggraph::scale_edge_linetype_manual(values=linetypes) + 
    ggraph::geom_node_point(aes_(color=~as.numeric(as.character(color))), size=6) + #
    scale_colour_gradient2(name = "logFC", low = "green",
                           mid = "blue", high = "red") 
  
  # TO DO: if label is meta gene or appear in stringdb network or pathway --> show label
  p <- p + ggraph::geom_node_text(aes_(label=~name), repel=TRUE, size=7)
  
  p <- p + ggtitle(paste0(datatag," \n Gene cluster: ",cluster_use)) + 
    theme(
      plot.title=element_text(color="black", size=15, hjust = 0.5, face='bold')
    )
  # p
  png(paste0(save_dir, datatag,'_',cluster_use,'.png'), height = 2*height_plot, width=2*width_plot,res = 2*72)
  print(p)
  dev.off()
  return(p)

}
#build an adjacency "correlation" matrix
compute_adjacency_correlation_mtx <- function(datExpr, save_dir, softPowerThres=18, option=1, save_data=T){
  if(option==1){
    WGCNA::enableWGCNAThreads()
    adjacency = WGCNA::adjacency(datExpr, power = softPowerThres, type = "signed") #specify network type
    print(head(adjacency))
    
    # Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
    #translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
    TOM = WGCNA::TOMsimilarity(adjacency, TOMType="signed") # specify network type
    diss_mtx = 1-TOM
    results <- list(diss_mtx=diss_mtx, TOM=TOM, adjacency=adjacency)
    save(results, file= paste0(save_dir, "results_mtx.RData"))
    return(results)
    
  }else{
    s = abs(WGCNA::bicor(datExpr))
    #calculation of adjacency matrix
    beta = 3  # need to pre-define beta value
    TOM = s^beta
    #dissimilarity measure
    diss_mtx = 1-TOM
    results <- list(diss_mtx=diss_mtx, TOM=TOM, adjacency=s)
    if(save_data){
      save(results, file= paste0(save_dir, "results_mtx.RData"))
    }
    
    return(results)
  }  
  
}
#Generate a clustered gene tree
generate_cluster_gene_tree_v1 <- function(meta_genes, datExpr, diss_mtx, save_dir, softPowerThres, minClusterSz=10, save_data=T){
  # Rfast::Dist(hclust_mtx, method = "euclidean")
  geneTree = flashClust::flashClust(as.dist(diss_mtx), method="average")
  # geneTree = hclust(as.dist(diss_mtx), method="average")
  
  plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
  #This sets the minimum number of genes to cluster into a module
  dynamicMods = dynamicTreeCut::cutreeDynamic(dendro=geneTree, distM=diss_mtx, deepSplit=2, #cutHeight=0.97,
                              pamRespectsDendro=FALSE, minClusterSize=minClusterSz)
  dynamicColors= labels2colors(dynamicMods)
  
  
  
  png(paste0(save_dir, "clustering_WGCNA_v1.png"), height =2*400, width= 2*900,res = 2*72)
  plotDendroAndColors(geneTree, dynamicColors, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, main='')
  dev.off()
  
  # softPowerThres <- 12
  MEList= moduleEigengenes(datExpr, colors=dynamicColors, softPower=softPowerThres)
  MEs= MEList$eigengenes
  MEDiss= 1-cor(MEs)
  METree= flashClust::flashClust(as.dist(MEDiss), method="average")
  
  png(paste0(save_dir, "eigengenes_clustering.png"), height =2*400, width= 2*400,res = 2*72)
  # sizeGrWindow(7, 6)
    plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
    # MEDissThres = 0.35
    # abline(h=MEDissThres, col = "red")
  dev.off()
  # Plot the cut line into the dendrogram
  # 
  # Call an automatic merging function
  # merge = WGCNA::mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  # mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  # mergedMEs = merge$newMEs
  
  # png(paste0(save_dir, "clustering_WGCNA_v1_merged.png"), height =2*400, width= 2*900,res = 2*72)
  # plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),
  # dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)#
  # dev.off()
  
  results <- list(dynamicMods=dynamicMods, dynamicColors=dynamicColors, 
                  MEList=MEList, MEs=MEs, MEDiss=MEDiss, METree=METree)
  
  
  gene_cluster <- data.frame(ens_gene=names(datExpr), cluster_id=results$dynamicMods, 
                             cluster=results$dynamicColors, stringsAsFactors=F)
  gene_cluster <- gene_cluster %>% left_join(meta_genes, by=c("ens_gene"))
  
  write.csv(gene_cluster, paste0(save_dir,'genes_clusters.csv'), quote=F, row.names = F)
  results[['gene_cluster']] <- gene_cluster
  print(table(results$dynamicMods, results$dynamicColors))
  
  if(save_data){
    save(results, file= paste0(save_dir, "Network_v1.RData"))
  }
  #plot the result with phytools package
  # par(mar=c(2,2,2,2))
  # phytools::plot.phylo.to.map(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
  # tiplabels(frame = 'circle',col='black', text=rep('',length(unique(dynamicMods))), 
  #           bg = levels(as.factor(dynamicColors)))
  return(results)
  
}

generate_cluster_gene_tree_v2 <- function(datExpr, diss_mtx, save_dir, minClusterSz=30){
  library(ape)
  #create gene tree by average linkage hierarchical clustering 
  geneTree = hclust(as.dist(diss_mtx), method = 'average')
  
  #module identification using dynamic tree cut algorithm
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = diss_mtx, deepSplit = 4, pamRespectsDendro = FALSE,
                              minClusterSize = minClusterSz)
  #assign module colours
  dynamicColors = labels2colors(dynamicMods)
  
  #plot the dendrogram and corresponding colour bars underneath
  png(paste0(save_dir, "_clustering_WGCNA_v2.png"), height =2*450, width= 2*900,res = 2*72)
  plotDendroAndColors(geneTree, dynamicColors, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, main='')
  dev.off()
  #calculate eigengenes
  MEs = moduleEigengenes(datExpr, colors = dynamicColors, excludeGrey = FALSE)$eigengenes
  
  #calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  
  #cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = 'average')
  
  results <- list(dynamicMods=dynamicMods, dynamicColors=dynamicColors, 
                  MEs=MEs, MEDiss=MEDiss, METree=METree)
  save(results, file= paste0(save_dir, "Network_v2.RData"))
  
  
  #plot the result with phytools package
  # par(mar=c(2,2,2,2))
  # plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
  # tiplabels(frame = 'circle',col='black', text=rep('',length(unique(dynamicMods))), bg = levels(as.factor(dynamicColors)))
  # 
  return(results)
}  

correlate_genes_with_treatment_effect <- function(results, datExpr, meta_cells, save_dir){
  #Define number of genes and samples
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  #Recalculate MEs with color labels
  MEs0 = WGCNA::moduleEigengenes(datExpr, results$gene_cluster$cluster)$eigengenes  #results$dynamicColors
  MEs = WGCNA::orderMEs(MEs0)
  moduleTraitCor = WGCNA::cor(MEs, meta_cells, method = "spearman")
  moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
  
  
  
  #Print correlation heatmap between modules and traits
  textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep= "")
  dim(textMatrix)= dim(moduleTraitCor)
  par(mar= c(6, 8.5, 3, 3))
  
  
  #display the corelation values with a heatmap plot
  #INCLUE THE NEXT LINE TO SAVE TO FILE
  #pdf(file="heatmap.pdf")
  png(paste0(save_dir, "hm_genes_treatment_correlation.png"), height =2*400, width= 2*500,res = 2*72)
  labeledHeatmap(Matrix= moduleTraitCor,
                 xLabels= names(meta_cells),
                 yLabels= names(MEs),
                 ySymbols= names(MEs),
                 colorLabels= FALSE,
                 colors= blueWhiteRed(50),
                 textMatrix= textMatrix,
                 setStdMargins= TRUE,
                 cex.text= 0.9,
                 zlim= c(-1,1),
                 main= paste("Gene clusters treatment effect correlation"))
  
  #INCLUE THE NEXT LINE TO SAVE TO FILE
  dev.off()
  
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  # meta_cells <- meta_cells %>%
  #               dplyr::select(-X4_UT)
  geneTraitSignificance = as.data.frame(cor(datExpr, meta_cells, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  
  names(geneTraitSignificance) = paste("GS.", names(meta_cells), sep="");
  names(GSPvalue) = paste("p.GS.", names(meta_cells), sep="");
  
  avg_gs <- round(rowMeans(abs(geneTraitSignificance)),4)
  
  
  names(datExpr)
  moduleColors <- results$gene_cluster$cluster
  for(module in unique(moduleColors)){
    # names(datExpr)[moduleColors==module]
    column = match(module, modNames);
    moduleGenes = moduleColors==module;
    
    png(paste0(save_dir, module,"_gs_module_membership.png"), height =2*400, width= 2*500,res = 2*72)
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "cluster"),
                       ylab = "Gene significance",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    
    dev.off()
  }
  
  
  probeToGene = data.frame(ens_gene_id=results$gene_cluster$ens_gene, 
                           gene_symbol=results$gene_cluster$gene, 
                           gene_type=results$gene_cluster$Gene_Type)
  rownames(probeToGene) <- probeToGene$ens_gene_id
  # Select module
  # module = "brown";# Select module probes
  probes = names(datExpr)
  # Create the starting data frame
  geneInfo0 = data.frame(ens_gene = probes,
                         gene = probeToGene[probes,'gene_symbol'],
                         gene_type = probeToGene[probes,'gene_type'],
                         cluster = moduleColors,
                         avg_gs=avg_gs,
                         geneTraitSignificance,
                         GSPvalue)#
  dim(geneInfo0)
  # dim(geneInfo0)
  # colnames(geneInfo0)
  # # Order modules by their significance for weight
  # modOrder = order(-abs(cor(MEs, meta_cells, use = "p")));
  # # Add module membership information in the chosen order
  # for (mod in 1:ncol(geneModuleMembership))
  # {
  #   oldNames = names(geneInfo0)
  #   geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, mod], 
  #                          MMPvalue[, mod]);
  #   names(geneInfo0) = c(oldNames, paste("MM.", modNames[mod], sep=""),
  #                        paste("p.MM.", modNames[mod], sep=""))
  # }
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  # get only positive linear trend genes
  gene_type_obs <- c('In_cis_Increase_UpRegulated','In_cis_Decrease_DownRegulated','In_trans_DownRegulated','In_trans_UpRegulated')
  geneInfo0 <- geneInfo0 %>%
    dplyr::filter(gene_type %in% gene_type_obs)
  
  geneOrder = order(geneInfo0$cluster, -abs(geneInfo0$avg_gs));
  geneInfo = geneInfo0[geneOrder, ]
  # View(geneInfo[1:20,])
  write.csv(geneInfo, file = paste0(save_dir, "gene_significant_GS_MM.csv"), row.names = F, quote = F)
  
  # Get top 10 genes in each cluster
  rs <- round(rowMeans(abs(moduleTraitCor[,2:ncol(moduleTraitCor)])),4)
  top_genes_ls <- list()
  for(c in unique(moduleColors)){
    gene_tmp_cis <- geneInfo %>%
      dplyr::filter(cluster==c & grepl('In_cis*',gene_type))
    gene_tmp_cis = gene_tmp_cis[order(gene_tmp_cis$cluster, -abs(gene_tmp_cis$avg_gs)), ]
    gene_tmp_cis <- gene_tmp_cis[1:5,]
    gene_tmp_cis$cluster_signif <- as.numeric(rs[paste0("ME",c)])
    if(nrow(gene_tmp_cis)>0){
      top_genes_ls[[paste0(c,'_cis')]] <- gene_tmp_cis
    }
    
    gene_tmp_trans <- geneInfo %>%
      dplyr::filter(cluster==c & grepl('In_trans*',gene_type))
    gene_tmp_trans = gene_tmp_trans[order(gene_tmp_trans$cluster, -abs(gene_tmp_trans$avg_gs)), ]
    gene_tmp_trans <- gene_tmp_trans[1:5,]
    gene_tmp_trans$cluster_signif <- as.numeric(rs[paste0("ME",c)])
    if(nrow(gene_tmp_trans)>0){
      top_genes_ls[[paste0(c,'_trans')]] <- gene_tmp_trans
    }
    
  }
  
  top_genes_cluster_df <- do.call(rbind, top_genes_ls)
  dim(top_genes_cluster_df)
  # View(top_genes_cluster_df)
  write.csv(top_genes_cluster_df, file = paste0(save_dir, "top5_gene_significant_GS_MM.csv"), row.names = F, quote = F)
  
}


# viz_gene <- function(original_mtx, MEs){
#   #Calculate module membership
#   MM = abs(bicor(original_mtx, MEs))
#   
#   #plot individual module of interest (MOI)
#   MOI = 3 #T cell differentiation co-expression module
#   plot(-log10(GS.lscore[modules==MOI,1]), MM[modules==MOI,MOI], pch=20,
#        cex=(GS.lscore[modules==MOI,4]/max(GS.lscore[,4],na.rm=TRUE))*4,
#        xlab='p-value (-log10) lymphocyte score', ylab='membership to module 3')
#   abline(v=-log10(0.05), lty=2, lwd=2)
# }

# compute_mean_exp <- function(norm_data, meta_cells_df){
#   exp_mean_df <- norm_data %>% 
#     # convert to long format
#     pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
#     # join with sample info table
#     full_join(meta_cells_df, by = ("cell_id")) #%>% 
#   # filter to retain only genes of interest
#   # filter(gene %in% observed_genes) %>% 
#   # for each gene
#   exp_mean_df$exp <- as.numeric(exp_mean_df$exp)
#   exp_mean_df <- exp_mean_df %>%
#     dplyr::group_by(gene) %>% 
#     # scale the cts column
#     dplyr::mutate(exp_scaled = (exp - mean(exp, na.rm=T))/sd(exp,na.rm=T)) 
#   
#   print(summary(exp_mean_df$exp_scaled))
#   # data.table::fwrite(x = exp_mean_df, paste0(save_dir, datatag,'_',tag,'_exp_mean_df.csv.gz'))
#   # data.table::fwrite(x = norm_data, paste0(save_dir, datatag,'_',tag,'_norm_data.csv.gz'))
#   # data.table::fwrite(x = meta_cells_df, paste0(save_dir, datatag,'_',tag,'_meta_cells.csv.gz'))
#   # 
#   gexp_cls <- exp_mean_df %>% 
#     # for each gene, summary by treatment_status
#     group_by(gene, treatment_status) %>%
#     # calculate the mean (scaled) cts
#     summarise(mean_exp_scaled = mean(exp_scaled, na.rm=T),
#               nrep = n()) %>% 
#     ungroup()
#   
#   return(gexp_cls)
# }

compute_mean_exp_v2 <- function(norm_data, meta_cells_df){
  cells_use <- colnames(norm_data)
  cells_use <- cells_use[cells_use!='gene']
  # length(cells_use)
  downsample_ratio <- 0.7
  gexp_cls_ls <- list()
  for(i in seq(10)){
    cells_to_sample <- sample(cells_use, floor(ncol(norm_data)*downsample_ratio), replace = FALSE)
    # print(cells_to_sample[1:2])
    exp_mean_df <- norm_data[,c(cells_to_sample,'gene')] %>% 
      # convert to long format
      pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
      # join with sample info table
      left_join(meta_cells_df, by = ("cell_id")) #%>% 
    # filter to retain only genes of interest
    # filter(gene %in% observed_genes) %>% 
    # for each gene
    exp_mean_df$exp <- as.numeric(exp_mean_df$exp)
    
    
    # print(summary(exp_mean_df$exp))
    gexp_cls <- exp_mean_df %>% 
      # for each gene, summary by treatment_status
      group_by(gene, treatment_status) %>%
      # calculate the mean (scaled) cts
      summarise(mean_exp_scaled_i = mean(exp, na.rm=T),
                nrep = n()) %>% 
      ungroup()
    
    gexp_cls_ls[[as.character(i)]] <- gexp_cls
  }
  
  s <- do.call(rbind, gexp_cls_ls)
  s <- as.data.frame(s)
  gexp_cls <- s %>%
    dplyr::group_by(gene, treatment_status) %>% 
    dplyr::summarise(mean_exp_scaled = mean(mean_exp_scaled_i, na.rm=T)) 
  
  
  # data.table::fwrite(x = exp_mean_df, paste0(save_dir, datatag,'_',tag,'_exp_mean_df.csv.gz'))
  # data.table::fwrite(x = norm_data, paste0(save_dir, datatag,'_',tag,'_norm_data.csv.gz'))
  # data.table::fwrite(x = meta_cells_df, paste0(save_dir, datatag,'_',tag,'_meta_cells.csv.gz'))
  # 
  return(gexp_cls)
}

compute_mean_zscore_exp_v3 <- function(norm_data, meta_cells_df){
  cells_use <- colnames(norm_data)
  cells_use <- cells_use[cells_use!='gene']
  # length(cells_use)
  downsample_ratio <- 0.7
  gexp_cls_ls <- list()
  for(i in seq(10)){
    cells_to_sample <- sample(cells_use, floor(ncol(norm_data)*downsample_ratio), replace = FALSE)
    # print(cells_to_sample[1:2])
    exp_mean_df <- norm_data[,c(cells_to_sample,'gene')] %>% 
      # convert to long format
      pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
      # join with sample info table
      left_join(meta_cells_df, by = ("cell_id")) #%>% 
    # filter to retain only genes of interest
    # filter(gene %in% observed_genes) %>% 
    # for each gene
    exp_mean_df$exp <- as.numeric(exp_mean_df$exp)
    
    exp_mean_df <- exp_mean_df %>%
      dplyr::group_by(gene) %>% 
      # scale the cts column
      dplyr::mutate(exp_scaled = (exp - mean(exp, na.rm=T))/sd(exp,na.rm=T)) 
    # print(summary(exp_mean_df$exp))
    gexp_cls <- exp_mean_df %>% 
      # for each gene, summary by treatment_status
      group_by(gene, treatment_status) %>%
      # calculate the mean (scaled) cts
      summarise(mean_exp_scaled_i = mean(exp_scaled, na.rm=T),
                nrep = n()) %>% 
      ungroup()
    
    gexp_cls_ls[[as.character(i)]] <- gexp_cls
  }
  
  s <- do.call(rbind, gexp_cls_ls)
  s <- as.data.frame(s)
  gexp_cls <- s %>%
    dplyr::group_by(gene, treatment_status) %>% 
    dplyr::summarise(mean_exp_scaled = mean(mean_exp_scaled_i, na.rm=T)) 
  
  
  # data.table::fwrite(x = exp_mean_df, paste0(save_dir, datatag,'_',tag,'_exp_mean_df.csv.gz'))
  # data.table::fwrite(x = norm_data, paste0(save_dir, datatag,'_',tag,'_norm_data.csv.gz'))
  # data.table::fwrite(x = meta_cells_df, paste0(save_dir, datatag,'_',tag,'_meta_cells.csv.gz'))
  # 
  return(gexp_cls)
}

compute_mean_exp <- function(norm_data, meta_cells_df){
  exp_mean_df <- norm_data %>% 
    # convert to long format
    pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
    # join with sample info table
    full_join(meta_cells_df, by = ("cell_id")) #%>% 
  # filter to retain only genes of interest
  # filter(gene %in% observed_genes) %>% 
  # for each gene
  exp_mean_df$exp <- as.numeric(exp_mean_df$exp)
  exp_mean_df <- exp_mean_df %>%
    dplyr::group_by(gene) %>% 
    # scale the cts column
    dplyr::mutate(exp_scaled = (exp - mean(exp, na.rm=T))/sd(exp,na.rm=T)) 
  
  print(summary(exp_mean_df$exp_scaled))
  # data.table::fwrite(x = exp_mean_df, paste0(save_dir, datatag,'_',tag,'_exp_mean_df.csv.gz'))
  # data.table::fwrite(x = norm_data, paste0(save_dir, datatag,'_',tag,'_norm_data.csv.gz'))
  # data.table::fwrite(x = meta_cells_df, paste0(save_dir, datatag,'_',tag,'_meta_cells.csv.gz'))
  # 
  gexp_cls <- exp_mean_df %>% 
    # for each gene, summary by treatment_status
    group_by(gene, treatment_status) %>%
    # calculate the mean (scaled) cts
    summarise(mean_exp_scaled = mean(exp_scaled, na.rm=T),
              nrep = n()) %>% 
    ungroup()
  
  return(gexp_cls)
}

plot_genes_exp <- function(gene_cluster, meta_genes, 
                           norm_data, norm_data_untreated,meta_cells_untreated,
                           meta_cells_df, obs_clones, 
                           save_dir, datatag, clone_aware=TRUE, plttitle=''){
  tag <- paste(obs_clones, collapse = '_') # tag <- 'H_A_K_I_J' SA1035
  if(datatag=='SA1035'){
    tag <- 'H_A_K_I_J'
  }
  gexp_cls <- compute_mean_zscore_exp_v3(norm_data, meta_cells_df)  #compute_mean_exp
  print(summary(gexp_cls$mean_exp_scaled))
  gexp_cls <- gexp_cls %>% inner_join(gene_cluster, by = "gene")
  data.table::fwrite(x = gexp_cls, paste0(save_dir, datatag,'_',tag,'_gexp_cls.csv'))
  # saveRDS(norm_data, paste0(save_dir, datatag,'_',tag,'_norm_data.rds'))
  # saveRDS(gexp_cls, file=paste0(save_dir, datatag,'_',tag,'_gexp_cls.rds'))
  print("Plot gene clusters:")
  pcm_treated <- plot_gene_clusters(gexp_cls, tag, datatag, obs_clones, save_dir, clone_aware,
                     xl="", yl="Mean z-score gene", trim_data=F, plttitle)
  png(paste0(save_dir,datatag,tag,"_treated_gene_regression.png"), height =2*300, width= 2*300+250,res = 2*72)
  print(pcm_treated)
  dev.off()
  
  # saveRDS(pcm_treated, paste0(save_dir,datatag, "pcm_treated_plt.rds"))
  
  tag <- paste(obs_clones_untreated, collapse = '_') # tag <- 'H_A_K_I_J' SA1035
  gexp_cls_untreated <- compute_mean_zscore_exp_v3(norm_data_untreated, meta_cells_untreated)

  print(summary(gexp_cls_untreated$mean_exp_scaled))
  gexp_cls_untreated <- gexp_cls_untreated %>% inner_join(gene_cluster, by = "gene")
  # data.table::fwrite(x = gexp_cls, paste0(save_dir, datatag,'_',tag,'_gexp_cls.csv'))
  # saveRDS(gexp_cls, file=paste0(save_dir, datatag,'_',tag,'_gexp_cls.rds'))
  print("Plot gene clusters:")
  pcm_untreated <- plot_gene_clusters(gexp_cls_untreated, tag, datatag, obs_clones, save_dir, clone_aware,
                                    xl="", yl="Mean z-score gene", trim_data=F, plttitle)

  # pcm_untreated
  p <- cowplot::plot_grid(pcm_treated, pcm_untreated, ncol=1, align = 'v') #label=c('Drug Cycle Treatment','Untreated')
  png(paste0(save_dir,datatag, "total_gene_cluster_exp_z-score.png"), height =2*570, width= 2*100*length(unique(gene_cluster$cluster))+250,res = 2*72)
  print(p)
  dev.off()
  
}

plot_genes_exp_v3 <- function(gene_cluster, meta_genes, norm_data,
                           meta_cells_df, obs_clones, 
                           save_dir, datatag, clone_aware=TRUE){
  tag <- paste(obs_clones, collapse = '_') # tag <- 'H_A_K_I_J' SA1035
  if(datatag=='SA1035'){
    tag <- 'H_A_K_I_J'
  }
  gexp_cls <- compute_mean_exp(norm_data, meta_cells_df)
  print(summary(gexp_cls$mean_exp_scaled))
  gexp_cls <- gexp_cls %>% inner_join(gene_cluster, by = c("gene"="ens_gene"))
  data.table::fwrite(x = gexp_cls, paste0(save_dir, datatag,'_',tag,'_gexp_cls.csv'))
  # saveRDS(gexp_cls, file=paste0(save_dir, datatag,'_',tag,'_gexp_cls.rds'))
  print("Plot gene clusters:")
  pcm_treated <- plot_gene_clusters(gexp_cls, tag, datatag, obs_clones, save_dir, clone_aware,
                                    xl="Time point - treatment status", yl="Mean z-score value", trim_data=F)
  return(pcm_treated)
  
}
plot_actual_genes_exp <- function(gene_cluster, meta_genes, norm_data, norm_data_untreated,meta_cells_untreated,
                           meta_cells_df, obs_clones, 
                           save_dir, datatag, clone_aware=TRUE){
  tag <- paste(obs_clones, collapse = '_') # tag <- 'H_A_K_I_J' SA1035
  if(datatag=='SA1035'){
    tag <- 'H_A_K_I_J'
  }
  # gexp_cls <- compute_mean_exp(norm_data, meta_cells_df)
  gexp_cls <- compute_mean_exp_v2(norm_data, meta_cells_df)
  print(summary(gexp_cls$mean_exp_scaled))
  gexp_cls <- gexp_cls %>% inner_join(gene_cluster, by = "gene")
  data.table::fwrite(x = gexp_cls, paste0(save_dir, datatag,'_',tag,'_gexp_cls.csv'))
  # saveRDS(gexp_cls, file=paste0(save_dir, datatag,'_',tag,'_gexp_cls.rds'))
  print("Plot gene clusters:")
  # pcm_treated <- plot_gene_clusters(gexp_cls, tag, datatag, obs_clones, save_dir, clone_aware, 
  #                    xl="Time point - treatment status", yl="Mean z-score value", trim_data=F)
  pcm_treated <- plot_gene_clusters(gexp_cls, tag, datatag, obs_clones, save_dir, clone_aware, 
                                    xl="Time point - treatment status", yl="Mean normalized exp value", trim_data=F)
  
  tag <- paste(obs_clones_untreated, collapse = '_') # tag <- 'H_A_K_I_J' SA1035
  # gexp_cls_untreated <- compute_mean_exp(norm_data_untreated, meta_cells_untreated)
  gexp_cls_untreated <- compute_mean_exp_v2(norm_data_untreated, meta_cells_untreated)
  
  print(summary(gexp_cls_untreated$mean_exp_scaled))
  gexp_cls_untreated <- gexp_cls_untreated %>% inner_join(gene_cluster, by = "gene")
  # data.table::fwrite(x = gexp_cls, paste0(save_dir, datatag,'_',tag,'_gexp_cls.csv'))
  # saveRDS(gexp_cls, file=paste0(save_dir, datatag,'_',tag,'_gexp_cls.rds'))
  print("Plot gene clusters:")
  # pcm_untreated <- plot_gene_clusters(gexp_cls_untreated, tag, datatag, obs_clones, save_dir, clone_aware, 
  #                                   xl="Time point - treatment status", yl="Mean z-score value", trim_data=F)
  pcm_untreated <- plot_gene_clusters(gexp_cls_untreated, tag, datatag, obs_clones, save_dir, clone_aware, 
                                      xl="Time point - treatment status", yl="Mean normalized exp value", trim_data=F)
  # pcm_untreated
  p <- cowplot::plot_grid(pcm_treated, pcm_untreated, ncol=1, align = 'v') #label=c('Drug Cycle Treatment','Untreated')
  png(paste0(save_dir,datatag, "total_gene_cluster_actual_normexp.png"), height =2*570, width= 2*100*length(unique(gene_cluster$cluster))+250,res = 2*72)
  print(p)
  dev.off()
  
}
gene_clustering <- function(meta_genes, norm_data, meta_cells_df, obs_clones, 
                   save_dir, datatag, clone_aware=TRUE){
  tag <- paste(obs_clones, collapse = '_') # tag <- 'H_A_K_I_J' SA1035
  
  exp_mean_df <- norm_data %>% 
    # convert to long format
    pivot_longer(!gene, names_to = "cell_id", values_to = "exp")  %>% 
    # join with sample info table
    full_join(meta_cells_df, by = ("cell_id")) %>% 
    # filter to retain only genes of interest
    # filter(gene %in% observed_genes) %>% 
    # for each gene
    group_by(gene) %>% 
    # scale the cts column
    mutate(exp_scaled = (exp - mean(exp, na.rm=T))/sd(exp)) 
  
  print(summary(exp_mean_df$exp_scaled))
  # exp_mean_df$exp_scaled <- ifelse(exp_mean_df$exp_scaled>2,2,
  #                                  ifelse(exp_mean_df$exp_scaled<-2,-2,exp_mean_df$exp_scaled))
  print(summary(exp_mean_df$exp_scaled))
  data.table::fwrite(x = exp_mean_df, paste0(save_dir, datatag,'_',tag,'_exp_mean_df.csv.gz'))
  data.table::fwrite(x = norm_data, paste0(save_dir, datatag,'_',tag,'_norm_data.csv.gz'))
  data.table::fwrite(x = meta_cells_df, paste0(save_dir, datatag,'_',tag,'_meta_cells.csv.gz'))
  # 
  gexp_cls <- exp_mean_df %>% 
    # for each gene, summary by treatment_status
    group_by(gene, treatment_status) %>%
    # calculate the mean (scaled) cts
    summarise(mean_exp_scaled = mean(exp_scaled, na.rm=T),
              nrep = n()) %>% 
    ungroup()
  
  
  dim(gexp_cls)
  length(unique(gexp_cls$gene))
  print(summary(gexp_cls$mean_exp_scaled))
  
  # Create a matrix
  print("Compute distance matrix:")
  hclust_mtx <- norm_data %>%
    select(-gene) %>%
    as.matrix()
  
  # assign rownames
  rownames(hclust_mtx) <- norm_data$gene
  dim(hclust_mtx)
  # hclust_mtx <- hclust_mtx[observed_genes, ]
  hclust_mtx <- hclust_mtx %>% 
    # transpose the matrix so genes are as columns
    t() %>% 
    # apply scalling to each column of the matrix (genes)
    scale() %>% 
    # transpose back so genes are as rows again
    t()
  
  
  print("Gene Regression:")
  gene_dist <- Rfast::Dist(hclust_mtx, method = "euclidean")
  data.table::fwrite(x = gene_dist, paste0(save_dir, datatag,'_',tag,'_gene_dist.csv.gz'))
  
  fitclusters <- NbClust::NbClust(hclust_mtx, diss=gene_dist, distance = NULL, 
                                  min.nc = 6, max.nc = 10, method = "kmeans", index = "silhouette") #index = "all"
  
  
  saveRDS(fitclusters, file=paste0(save_dir, datatag,'_',tag,'_fit_clusters.rds'))
  
  
  # gene_hclust <- hclust(gene_dist, method = "complete")
  # The default `plot()` function can be used to produce a simple dendrogram
  # plot(gene_hclust, labels = FALSE)
  # abline(h = 110, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
  # nbclusters <- 10
  # gene_cluster <- cutree(gene_hclust, k = nbclusters) %>% 
  #   # turn the named vector into a tibble
  #   enframe() 
  
  
  gene_cluster <- fitclusters$Best.partition %>% tibble::enframe()
  gene_cluster <- as.data.frame(gene_cluster)
  colnames(gene_cluster) <- c('gene', 'cluster')

  print(summary(as.factor(gene_cluster$cluster)))
  gene_cluster <- gene_cluster %>% left_join(meta_genes, by=c("gene"))
  data.table::fwrite(x = gene_cluster, paste0(save_dir, datatag,'_',tag,'_gene_cluster.csv'))
  # data.table::fwrite(x = gene_cluster, paste0(save_dir, datatag,'_',tag,'_gene_cluster_classified.csv'))
  
  gexp_cls <- gexp_cls %>% inner_join(gene_cluster, by = "gene")
  data.table::fwrite(x = gexp_cls, paste0(save_dir, datatag,'_',tag,'_gexp_cls.csv'))
  # saveRDS(gexp_cls, file=paste0(save_dir, datatag,'_',tag,'_gexp_cls.rds'))
  print("Plot gene clusters:")
  plot_gene_clusters(gexp_cls, tag, datatag, obs_clones, save_dir, clone_aware, 
                     xl="Time point - treatment status", yl="Mean z-score value", trim_data=F)
  
  print(dim(gexp_cls))
  print(summary(gexp_cls$mean_exp_scaled))
  print("Gene Regression:")
  # gene_regression_SA1035(exp_mean_df, gene_cluster, save_dir, paste0(datatag,'_',tag))
  # exp_mean_df <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_exp_mean_df.csv.gz'))
  # # View(head(exp_mean_df))
  # gene_cluster <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_gene_cluster.csv')) #
  # meta_cells_df <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_meta_cells.csv.gz'))
  # meta_cells_df <- as.data.frame(meta_cells_df)
  # gene_cluster <- as.data.frame(gene_cluster)
  # View(head(meta_cells_df))
  # length(intersect(unique(exp_mean_df$gene), gene_cluster$gene))
  
  # gene_cluster <- gene_regression_SA1035(exp_mean_df, gene_cluster, save_dir, paste0(datatag,'_',tag))
  gene_cluster <- gene_regression_SA535(exp_mean_df, gene_cluster, save_dir, paste0(datatag,'_',tag))
  table(gene_cluster$classify_gene, gene_cluster$cluster)
  # obs_clusters <- gene_cluster[gene_cluster$classify_gene %in% c('mono_incrc','mono_decrc'),'cluster']
  # print("Plot gene important clusters:")
  # if(!is.null(obs_clusters)){
  #   plot_gene_each_cluster(gexp_cls, datatag, obs_clones, as.character(obs_clusters), save_dir,
  #                          xl="Time point - treatment status", yl="Mean z-score value")
  # }
  
  
 
  # head(trans_cts_cluster)
  
}

plot_heatmap <- function(gene_cluster, norm_data, meta_cells_df){
  unique(gene_cluster$cluster)
  gene_cluster$cluster <- paste0('Cluster ',gene_cluster$cluster)
  rownames(gene_cluster) <- gene_cluster$ens_gene
  rownames(meta_cells_df) <- meta_cells_df$cell_id
  gene_cluster$ens_gene[1:3]
  norm_data <- norm_data %>%
    select(-gene) %>%
    as.matrix()
  row_cls <- gene_cluster[rownames(norm_data),'cluster']
  col_cls <- meta_cells_df[colnames(norm_data),'treatment_status']
  p <- ComplexHeatmap::Heatmap(norm_data, na_col = "white",
                               show_column_names=F,
                               show_row_names = F,
                               cluster_rows=F,cluster_columns=F,
                               name = "Genes Clusters", 
                               # row_order = sort(rownames(test)),
                               row_split= row_cls,
                               row_title_rot = 0,
                               # row_gap = unit(2, "mm"),
                               column_split = col_cls, 
                               column_title = "Resistant Genes Regression across time",
                               # column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 10),
                               row_names_gp = grid::gpar(fontsize = 10),
                               show_heatmap_legend = T,
                               # top_annotation=top_anno,
                               # left_annotation = left_anno,
                               # cell_fun = cell_func,
                               row_dend_reorder=F)
  # p
  
  
}

gene_regression_SA1035 <- function(exp_mean_df, gene_cluster, save_dir, datatag){
  # gene_cluster <- as.data.frame(gene_cluster)
  # exp_mean_df <- as.data.frame(exp_mean_df)
  epsilon <- 0.05 # error rate between different batches
  rownames(gene_cluster) <- gene_cluster$gene
  gene_cluster$classify_gene <- ''
  gene_cluster$slope <- NA
  for(cl in unique(gene_cluster$cluster)){
    print(paste0("Observe clone: ", cl))
    genes_use <- gene_cluster[gene_cluster$cluster==cl,]$gene
    tmp <- exp_mean_df[exp_mean_df$gene %in% genes_use,]
    # tmp <- exp_mean_df %>%
    #               dplyr::filter(gene %in% as.character(genes_use))
    tmp$treatment_status <- as.factor(tmp$treatment_status)
    dim(tmp)
    rg <- lm(exp_scaled ~ treatment_status, data = tmp)  
    saveRDS(rg, paste0(save_dir, datatag,'_cluster_',cl,'_lm.rds'))
    gene_cluster[genes_use,'slope'] <- rg$coefficients[2] 
    print(summary(rg))
    coeffs <- coef(rg)
    # Get means in each treatment
    ts_med <- tapply(tmp$exp_scaled, tmp$treatment_status, median)
    # print(ts_means)
    c2 <- coeffs[2] > 0
    c3 <- coeffs[3] > 0
    c21 <- ts_med[2] - ts_med[1] >= -epsilon
    c32 <- ts_med[3] - ts_med[2] >= -epsilon
    c31 <- ts_med[3] - ts_med[1] >= -epsilon
    # c43 <- ts_med[4] - ts_med[3] >= -epsilon
    
    c12 <- ts_med[1] - ts_med[2] >= -epsilon
    c23 <- ts_med[2] - ts_med[3] >= -epsilon
    c13 <- ts_med[1] - ts_med[3] >= -epsilon
    # c34 <- ts_med[3] - ts_med[4] >= -epsilon
    if(c21 & c31){
      if(c32){
        gene_type <- 'mono_incrc'
      } else{
        gene_type <- 'incrc'
      }
      print(gene_type)
    } else if(!c2 & !c3 & c12 & c13){
      if(c23){
        gene_type <- 'mono_decrc'
      }else{
        gene_type <- 'decrc'
      }
      print(gene_type)
      
    } else{
      gene_type <- 'undefined'
    }
    
    gene_cluster[genes_use,'classify_gene'] <- gene_type
    print(paste0("Cluster ",cl," and gene type is: ",gene_type))
    print(" ")
  }
  
  write.csv(gene_cluster, paste0(save_dir, datatag,'_gene_clusters_classified.csv'), quote = F, row.names = F)
  return(gene_cluster)
}

gene_regression_SA535 <- function(exp_mean_df, gene_cluster, save_dir, datatag){
  # gene_cluster <- as.data.frame(gene_cluster)
  # exp_mean_df <- as.data.frame(exp_mean_df)
  epsilon <- 0.05 # error rate between different batches
  rownames(gene_cluster) <- gene_cluster$gene
  gene_cluster$classify_gene <- ''
  gene_cluster$slope <- NA
  for(cl in unique(gene_cluster$cluster)){
    print(paste0("Observe clone: ", cl))
    genes_use <- gene_cluster[gene_cluster$cluster==cl,]$gene
    tmp <- exp_mean_df[exp_mean_df$gene %in% genes_use,]
    # tmp <- exp_mean_df %>%
    #               dplyr::filter(gene %in% as.character(genes_use))
    tmp$treatment_status <- as.factor(tmp$treatment_status)
    dim(tmp)
    rg <- lm(exp_scaled ~ treatment_status, data = tmp)  
    saveRDS(rg, paste0(save_dir, datatag,'_cluster_',cl,'_lm.rds'))
    gene_cluster[genes_use,'slope'] <- rg$coefficients[2] 
    print(summary(rg))
    coeffs <- coef(rg)
    # Get means in each treatment
    ts_med <- tapply(tmp$exp_scaled, tmp$treatment_status, median)
    # print(ts_means)
    c2 <- coeffs[2] > 0
    c3 <- coeffs[3] > 0
    c21 <- ts_med[2] - ts_med[1] >= -epsilon
    c32 <- ts_med[3] - ts_med[2] >= -epsilon
    c31 <- ts_med[3] - ts_med[1] >= -epsilon
    # c43 <- ts_med[4] - ts_med[3] >= -epsilon
    
    c12 <- ts_med[1] - ts_med[2] >= -epsilon
    c23 <- ts_med[2] - ts_med[3] >= -epsilon
    c13 <- ts_med[1] - ts_med[3] >= -epsilon
    # c34 <- ts_med[3] - ts_med[4] >= -epsilon
    if(c21 & c31){
      if(c2 & c3 & c32){
        gene_type <- 'mono_incrc'
      } else{
        gene_type <- 'incrc'
      }
      print(gene_type)
    } else if(!c2 & !c3 & c12 & c13){
      if(c23){
        gene_type <- 'mono_decrc'
      }else{
        gene_type <- 'decrc'
      }
      print(gene_type)
      
    } else{
      gene_type <- 'undefined'
    }
    
    
    gene_cluster[genes_use,'classify_gene'] <- gene_type
    print(paste0("Cluster ",cl," and gene type is: ",gene_type))
    print(" ")
  }
  
  
  
  write.csv(gene_cluster, paste0(save_dir, datatag,'_gene_clusters_classified.csv'), quote = F, row.names = F)
  return(gene_cluster)
}


gene_regression <- function(exp_mean_df, gene_cluster, save_dir, datatag){
  gene_cluster <- as.data.frame(gene_cluster)
  exp_mean_df <- as.data.frame(exp_mean_df)
  epsilon <- 0.1 # error rate between different batches
  rownames(gene_cluster) <- gene_cluster$gene
  gene_cluster$classify_gene <- ''
  gene_cluster$slope <- NA
  for(cl in unique(gene_cluster$cluster)){
    print(paste0("Observe clone: ", cl))
    genes_use <- gene_cluster[gene_cluster$cluster==cl,]$gene
    tmp <- exp_mean_df[exp_mean_df$gene %in% genes_use,]
    # tmp <- exp_mean_df %>%
    #               dplyr::filter(gene %in% as.character(genes_use))
    tmp$treatment_status <- as.factor(tmp$treatment_status)
    dim(tmp)
    rg <- lm(exp_scaled ~ treatment_status, data = tmp)  
    saveRDS(rg, paste0(save_dir, datatag,'_cluster_',cl,'_lm.rds'))
    gene_cluster[genes_use,'slope'] <- rg$coefficients[2] 
    print(summary(rg))
    coeffs <- coef(rg)
    # Get means in each treatment
    ts_means <- tapply(tmp$exp_scaled, tmp$treatment_status, mean)
    # print(ts_means)
    c2 <- coeffs[2] > 0
    c3 <- coeffs[3] > 0
    c4 <- coeffs[4] > 0
    
    c32 <- ts_means[3] - ts_means[2] >= -epsilon
    c43 <- ts_means[4] - ts_means[3] >= -epsilon
    
    c21 <- ts_means[2] - ts_means[1] >= -epsilon
    c23 <- ts_means[2] - ts_means[3] >= -epsilon
    c34 <- ts_means[3] - ts_means[4] >= -epsilon
    if(c2 & c4){
      if(c32 & c43){
        gene_type <- 'mono_incrc'
      } else{
        gene_type <- 'incrc'
        
      }
    } else if(!c21 & !c4){
      if(c23 & c34){
        gene_type <- 'mono_decrc'
      } else{
        gene_type <- 'decrc'
      }
    } else if(c21 & !c3 & !c4 & c23 & c34){
      gene_type <- 'decrc'
    } else{
      gene_type <- 'unknown'
    }
    
    gene_cluster[genes_use,'classify_gene'] <- gene_type
    print(paste0("Cluster ",cl," and gene type is: ",gene_type))
    print(" ")
  }
  
  
  
  write.csv(gene_cluster, paste0(save_dir, datatag,'_gene_clusters_classified.csv'), quote = F, row.names = F)
  return(gene_cluster)
}

plot_gene_clusters <- function(gexp_cls, tag, datatag, obs_clones, save_dir, clone_aware,
                               xl="Time point - treatment status", yl="Mean z-score value", trim_data=TRUE, plttitle=''){
  print(summary(gexp_cls$mean_exp_scaled))
  # gexp_cls$gene_type <- ifelse(grepl('In_cis',gexp_cls$Gene_Type),'In_cis','In_trans')
  gexp_cls$mean_exp_scaled <- ifelse(as.numeric(gexp_cls$mean_exp_scaled)>4,4, gexp_cls$mean_exp_scaled)
  # if(trim_data){
  #   gexp_cls$mean_exp_scaled <- ifelse(as.numeric(gexp_cls$mean_exp_scaled)>-1.5,-1.5, gexp_cls$mean_exp_scaled)
  # }
  nbclusters <- length(unique(gexp_cls$cluster))
  if(clone_aware){
    if(is.null(tag)){
      tag <- paste0(' in clone ',paste(obs_clones, collapse = '_'))
    }
    tag1 <- paste0('_clone_aware_',paste(obs_clones, collapse = '_'))
  }else{
    tag <- ''
    tag1 <- '_without_clone_aware'
  }
  # plttitle <- paste0(datatag, ": DE genes across time cls ",tag)
  
  # lowpass.spline <- smooth.spline(gexp_cls$treatment_status, gexp_cls$mean_exp_scaled, spar = 0.6) ## Control spar for amount of smoothing
  
  gexp_cls$cluster <- paste0('Cluster ',gexp_cls$cluster) # nb of genes in each cluster here 
  pcm <- gexp_cls %>% 
    ggplot(aes(treatment_status, mean_exp_scaled)) +
    geom_boxplot() 
  edge_colrs <- c("In_cis"="black", "In_trans"="grey51")
  pcm <- gexp_cls %>% 
    ggplot(aes(treatment_status, mean_exp_scaled)) +
    # geom_line(aes(group = gene, color=gene_type), alpha = 0.7) +  #geom_line  #in case incis intrans
    geom_line(aes(group = gene), color='black', alpha = 0.7) +
    scale_color_manual(values=edge_colrs) + 
    # geom_boxplot() +
    # geom_line(predict(lowpass.spline, gexp_cls$treatment_status), col = "red", lwd = 2) + 
    geom_line(stat = "summary", fun = "median", colour = "brown", size = 1,
              aes(group = 1)) +
    facet_grid(cols = vars(cluster))
  
  pcm <- pcm + labs(x=xl, y=yl, title=plttitle)
  pcm <- pcm + theme(plot.title = element_text(color="black", size=13, hjust = 0.5),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.y = element_text(color="black", size=10),
                     axis.text.x = element_text(color="black", size=8, angle=90),
                     axis.title = element_text(color="black", size=12)
  )#egend.text = element_text(color="black", size=6)
  # png(paste0(save_dir,datatag,tag, tag1 ,"_gene_regression.png"), height =2*300, width= 2*100*nbclusters+250,res = 2*72)
  #   print(pcm)
  # dev.off()
  
  return(pcm)
  
}

plot_genes_network <- function(){
  foldChange <- gene_cluster1$logFC_UUTTTT_S_T_UUUUUU_J_Q
  names(foldChange) <- gene_cluster1$gene
  geneSets <- list(brown=gene_cluster1$gene)
  cnetplot_v2(geneSets, datatag, tag, save_dir, foldChange)
}
plot_gene_each_cluster <- function(gexp_cls, datatag, obs_clones, obs_clusters, save_dir,
                                   xl="Time point - treatment status", yl="Mean z-score value"){
  tag <- paste0(' in clone ',paste(obs_clones, collapse = ' '))
  # if(clone_aware){
  #   tag <- paste0(' in clone ',paste(obs_clones, collapse = ' '))
  #   tag1 <- paste0('_clone_aware_',paste(obs_clones, collapse = '_'))
  # }else{
  #   tag <- ''
  #   tag1 <- '_without_clone_aware'
  # }
  plttitle <- paste0("Resistant genes across treatment time points ",tag)
  # lowpass.spline <- smooth.spline(gexp_cls$treatment_status, gexp_cls$mean_exp_scaled, spar = 0.6) ## Control spar for amount of smoothing
  plot_ls <- list()
  for(cl in obs_clusters){
    pltsubtitle <- paste0("Genes Cluster ",cl)
    gexp_cls_tmp <- gexp_cls %>%
      dplyr::filter(cluster==cl)
    pcm <- gexp_cls_tmp %>% 
      ggplot(aes(treatment_status, mean_exp_scaled)) +
      geom_line(aes(group = gene, color=gene), alpha = 0.3) +  #geom_line  #
      # geom_smooth(method = "loess") +
      # geom_line(predict(lowpass.spline, gexp_cls$treatment_status), col = "red", lwd = 2) + 
      geom_line(stat = "summary", fun = "median", colour = "brown", size = 1,
                aes(group = 1)) #+
    # facet_grid(cols = vars(cluster))
    
    pcm <- pcm + labs(x=xl, y=yl, title=pltsubtitle)
    pcm <- pcm + theme(plot.title = element_text(color="black", size=13, hjust = 0.5),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                       axis.text.y = element_text(color="black", size=10),
                       axis.text.x = element_text(color="black", size=8, angle=90),
                       axis.title = element_text(color="black", size=12),
                       legend.text = element_text(color="black", size=5)
    )#
    plot_ls[[cl]] <- pcm
  }
  main_plot <- cowplot::plot_grid(plotlist = plot_ls,
                                  ncol = 1,
                                  align = 'v'
  )
  
  
  png(paste0(save_dir,datatag,tag ,"_gene_regression.png"), height = 2*600, width=2*1200,res = 2*72)
  print(main_plot)
  dev.off()
  
}


get_top_up_genes <- function(de_genes, ntop=NULL, minLogFC=0.25, FDR_cutoff=0.01, pValueThrs=0.5){
  de_genes <- de_genes %>%
    dplyr::filter(logFC>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
  de_genes <- de_genes[order(de_genes$logFC,decreasing = T),] 
  if(ntop > nrow(de_genes)){
    ntop = nrow(de_genes)
  }
  if(!is.null(ntop)){
    return(de_genes[1:ntop,])  
  }else{  # select all resistant genes with logFC > 0
    return(de_genes)
  }
  
}

downsample_cells <- function(meta_cells_df,downsample_ratio=0.5,
                             thres_small_clone=500){
  # meta_cells_df <- as.data.frame(meta_cells_df)
  ts <- unique(meta_cells_df$treatmentSt)
  thres_minority <- 20
  ext_cells <- c()
  for(t in ts){
    tmp <- meta_cells_df %>%
      dplyr::filter(treatmentSt==t)
    
    if(!is.null(tmp) & nrow(tmp)<thres_small_clone & nrow(tmp)>thres_minority){   #if clone contain less than 15 cells, remove this clone from consideration
      cells_to_sample <- tmp$cell_id
    } else if(!is.null(tmp) & nrow(tmp) >= thres_small_clone){
      cells_to_sample <- sample(tmp$cell_id, floor(nrow(tmp)*downsample_ratio),replace = FALSE)
    } else{
      cells_to_sample <- NULL
      
      if(!is.null(tmp)){
        print(paste0('DEBUG: double check this step, treatment st: ',t, ' nb cells: ',nrow(tmp)))
      } else{
        print(paste0('DEBUG: double check this step, treatment st: ',t))
      }
    }
    print(paste0('Extract ',length(cells_to_sample),' from treatment st: ', t))
    
    if(length(cells_to_sample)>0){
      ext_cells <- c(ext_cells, cells_to_sample)
    }
    
  }
  print(length(ext_cells))
  return(meta_cells_df[ext_cells,])
}


get_treatment_status <- function(status) {
  labels <- sapply(strsplit(status, "-"), function(x) {
    return(x[3])
  })
  return(as.character(labels))
}

set_theme <- function(p, plttitle, subtitle, rv_yaxis=F){
  p <- p + labs(title = plttitle, subtitle = subtitle)
  if(rv_yaxis){
    p <- p + ylab("")
  }
  p <- p + theme(plot.title = element_text(color="black", size=12, hjust = 0.5), #, face="bold"
                 plot.subtitle =  element_text(color="black", size=9, hjust = 0.5), 
                 # legend.text = element_blank(),
                 legend.position="none")
  return(p)
}


plot_3series_genes_clustering <- function(){
  library(ggplot2)
  
  # legend.position = c(0.8, 0.2)
  
  
  datatag <- 'SA609'
  gene_type <- 'increase'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/gene_regression_',gene_type,'/')
  p609_increase <- readRDS(paste0(save_dir, datatag, 'pcm_treated_plt.rds'))
  p609_increase <- set_theme(p609_increase, 'SA609 cisplatin: up-regulated DE genes', 'in Rx: UTTTT-R vs UnRx: UUUUU-H')
  
  gene_type <- 'decrease'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/gene_regression_',gene_type,'/')
  p609_decrease <- readRDS(paste0(save_dir, datatag, 'pcm_treated_plt.rds'))
  p609_decrease <- set_theme(p609_decrease, 'SA609 cisplatin: down-regulated DE genes', 'in Rx: UTTTT-R vs UnRx: UUUUU-H', T)
  
  datatag <- 'SA1035'
  gene_type <- 'increase'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/gene_regression_',gene_type,'/')
  p1035_increase <- readRDS(paste0(save_dir, datatag, 'pcm_treated_plt.rds'))
  p1035_increase <- set_theme(p1035_increase, 'SA1035 cisplatin: up-regulated DE genes', 'in Rx: UTTTT-H vs UnRx: UUUUU-E')
  
  
  gene_type <- 'decrease'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/gene_regression_',gene_type,'/')
  p1035_decrease <- readRDS(paste0(save_dir, datatag, 'pcm_treated_plt.rds'))
  p1035_decrease <- set_theme(p1035_decrease, 'SA1035 cisplatin: down-regulated DE genes', 'in Rx: UTTTT-H vs UnRx: UUUUU-E', T)
  
  datatag <- 'SA535'
  gene_type <- 'increase'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/gene_regression_',gene_type,'_cis/')
  p535_increase_cis <- readRDS(paste0(save_dir, datatag, 'pcm_treated_plt.rds'))
  p535_increase_cis <- set_theme(p535_increase_cis, 'SA535 cisplatin: up-regulated DE genes', 'in Rx: UUTTT-S_T vs UnRx: UUUUU-J')
  
  gene_type <- 'decrease'
  datatag <- 'SA535'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/gene_regression_',gene_type,'_cis/')
  p535_decrease_cis <- readRDS(paste0(save_dir, datatag, 'pcm_treated_plt.rds'))
  # p535_decrease_cis <- set_theme(p535_decrease_cis, 'SA535 cisplatin: down-regulated DEG in Rx vs UnRx', 'Rx: UUTTT-S_T vs UnRx: UUUUU-J')
  # p535_decrease_cis <- p535_decrease_cis + labs(title = 'SA535 cisplatin: down-regulated DE genes', subtitle = 'in Rx: UUTTT-S_T vs UnRx: UUUUU-J')
  p <- p535_decrease_cis
  
  p535_decrease_cis$data$gene_type <- paste0(p535_decrease_cis$data$gene_type,'blahblah')
  # length(p535_decrease_cis$layers[1])
  p535_decrease_cis <- set_theme(p535_decrease_cis, 'SA535 cisplatin: down-regulated DE genes', 'in Rx: UUTTT-S_T vs UnRx: UUUUU-J')
  
  
  p <- p + theme(plot.title = element_text(color="black", size=12, hjust = 0.5),#, face="bold"
                                                 plot.subtitle =  element_text(color="black", size=9, hjust = 0.5),
                                                 legend.title  = element_text(color="black", size=15, hjust = 0.5),
                                                 legend.text = element_text(color="black", size=15, hjust = 0.5)) +
    guides(fill = guide_legend(title="Gene Type"))
  lg <- cowplot::get_legend(p)
  pc <- ggdraw() + draw_plot(lg)
  combine_plt <- plot_grid(p535_decrease_cis, lg, rel_widths = c(2.1,1))
  combine_plt
  
  gene_type <- 'increase'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/gene_regression_',gene_type,'_cx/')
  p535_increase_cx <- readRDS(paste0(save_dir, datatag, 'pcm_treated_plt.rds'))
  p535_increase_cx <- set_theme(p535_increase_cx, 'SA535 CX5461: up-regulated DE genes', 'in Rx: UXXXX-U vs UnRx: UUUUU-J')
  
  gene_type <- 'decrease'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/gene_regression_',gene_type,'_cx/')
  p535_decrease_cx <- readRDS(paste0(save_dir, datatag, 'pcm_treated_plt.rds'))
  p535_decrease_cx <- set_theme(p535_decrease_cx, 'SA535 CX5461: down-regulated DE genes', 'in Rx: UXXXX-U vs UnRx: UUUUU-J',T)
  
  plotls <- list(p609_increase, p609_decrease, p1035_increase, p1035_decrease, 
                 p535_increase_cis, combine_plt, p535_increase_cx, p535_decrease_cx)
  
  lbs <- c('A','B','C','D','E','F','G','H')
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/gene_regression_increase_cis/'
  p <- cowplot::plot_grid(plotlist = plotls, ncol=2, align = 'h', labels = lbs, rel_widths = c(1.02, 1)) 
  png(paste0(save_dir, "genes_clustering.png"), height =2*1000, width= 2*800,res = 2*72)
  print(p)
  dev.off()
  
}


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
