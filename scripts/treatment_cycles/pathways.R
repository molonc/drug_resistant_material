## May 18, 2021
# modifying to incorporate many samples et




library(reshape2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(cowplot)

set <- "HALLMARK"
#set <- "KEGG"
#sthr <- "0.0001"
#sthr <- "0.001"
sthr <- 0.01
#sthr <- 0.05
#scol <- "PValue"
scol <- "FDR"
ver <- 'v6'

# The KEGG results don't make much sense
#set <- "KEGG"

what <- "all_pathways_by_clone_scran"

#minfc <- "0.0"
minfc <- "0.25"    ## IN the paper Fig 5 we used log2 fold change 0.25
#minfc <- "0.5"
pfdr <- 0.05
#pfdr <- 0.5
#pfdr <- 0.01

fbody <- paste0(set,"-",scol,"-",sthr,"-MINFC-",minfc,"-PFDR-",pfdr)
pt4_name <- "-SA609-v6-"
pt5_name <- "-SA535_cisplatin-v7-"
pt6_name <- "-SA1035-v6-"

## First keep the consistent, than other abundant clones
if (what == "all_pathways_by_clone_scran") {
  short_names <- c(
    paste0("scrande_SA609_11"),
    paste0("scrande_SA609_1"),
    paste0("scrande_SA609_22"),   # holiday B
    paste0("scrande_SA609_2"),
    paste0("scrande_SA609_31"),   # holiday R
    paste0("scrande_SA609_3"),
    paste0("scrande_SA609_23"),   # holiday
    paste0("scrande_SA609_4"),
    paste0("scrande_SA609_24"),   # holiday
    
    #paste0("scrande_SA535_4"),
    #paste0("scrande_SA535_5"),
    #paste0("scrande_SA535_52"),   # holiday
    paste0("scrande_SA535_1"),
    paste0("scrande_SA535_16"),   # holiday
    paste0("scrande_SA535_3"),
    paste0("scrande_SA535_78"),   # holiday A
    paste0("scrande_SA535_36"),   # holiday E
    
    #paste0("scran-nocl-SA1035-1"),
    paste0("scrande_SA1035_41"),
    paste0("scrande_SA1035_3"),
    paste0("scrande_SA1035_32"),   # holiday
    paste0("scrande_SA1035_51"),   # 
    paste0("scrande_SA1035_12"),   # holiday
    paste0("scrande_SA1035_52"),   # 
    paste0("scrande_SA1035_22"),   # holiday
    paste0("scrande_SA1035_1"),
    paste0("scrande_SA1035_36"),   # holiday  
    paste0("scrande_SA1035_2"),
    paste0("scrande_SA1035_34")    # holiday
  )
  
  files <- paste0(short_names[grep("SA609", short_names)], pt4_name, fbody)
  files <- c(files, paste0(short_names[grep("SA535", short_names)], pt5_name, fbody))
  files <- c(files, paste0(short_names[grep("SA1035", short_names)], pt6_name, fbody))
  
  labels1 <- c("X4 1 Rx:B / UnRx:H",
               "X4 2 Rx:A / UnRx:H",
               "X5 1 RxH:B / UnRx:H",
              "X5 2 Rx:A / UnRx:H",
              "X5 3 RxH:A / UnRx:H",
              "X6 Rx:A / UnRx:H",
              "X6 RxH:A / UnRx:H",
              "X7 Rx:A / UnRx:H",
              "X7 RxH:A / UnRx:H",

              "X8 Rx:A / UnRx:G",
              "X8 RxH:D / UnRx:G",
              "X10 Rx:A / UnRx:G",
              "X10 RxH:A / UnRx:G",
              "X10 RxH:E / UnRx:G",

              # Feb 11, 2024: replaced clone D with B
              # "X5 Rx:D / UnRx:A",
              "X5 Rx:B / UnRx:A",
              "X6 Rx:G / UnRx:G",
              "X6 RxH:B / UnRx:G",
              "X7 Rx:G / UnRx:E",
              "X7 RxH:G / UnRx:E",
              "X8 Rx:G / UnRx:E",
              "X8 RxH:G / UnRx:E",
              "X7 Rx:H / UnRx:E",
              "X7 RxH:H / UnRx:E",
              "X8 Rx:H / UnRx:E",
              "X8 RxH:H / UnRx:E"
  )
  
  groups <- c("Pt4\nX4", "Pt4\nX4", "Pt4\nX5", "Pt4\nX5", "Pt4\nX5", "Pt4\nX6", "Pt4\nX6", "Pt4\nX7", "Pt4\nX7", 
              #"Pt5\nX6", "Pt5\nX7", "Pt5\nX7", 
              "Pt5\nX8", "Pt5\nX8", "Pt5\nX10", "Pt5\nX10", "Pt5\nX10",
              "Pt6\nX5", "Pt6\nX6",  "Pt6\nX6", "Pt6\nX7G", "Pt6\nX7G", "Pt6\nX8G", "Pt6\nX8G", "Pt6\nX7H", "Pt6\nX7H",  "Pt6\nX8H", "Pt6\nX8H"
  )
  
  passages <- c("X4", "X4", "X5", "X5", "X5", "X6", "X6", "X7", "X7", 
              #"Pt5\nX6", "Pt5\nX7", "Pt5\nX7", 
              "X8", "X8", "X10", "X10", "X10",
              "X5", "X6",  "X6", "X7", "X7", "X8", "X8", "X7", "X7",  "X8", "X8"
  )
  series <- c(rep("Pt4", 9), rep("Pt5", 5), rep("Pt6", 11))
  
  filename <- paste0("pathways_minfc",minfc,"_per_clone_scran_pfdr",pfdr,"_defdr", sthr,"_minpath2.pdf")
}

input_data <- data.frame(file_header=short_names, labels_detailed=labels1, passage=passages, series=series)

write.csv(file="input_files.csv", input_data, quote = FALSE)

nseries <- length(files) ###5  ## 3 series

#fgroups <- files

m2 <- NULL
data <- list()

for (i in 1:length(files)) {
  dir <- "../../materials/gsea_results/"
  
  file <- paste0(dir, files[i], "-pathway.csv")
  data[[i]] <- read.csv(file, header=TRUE)
  colnames(data[[i]])[3] <- labels1[i]   ## labels1
}

m2 <- merge(x = data[[1]][,c(1,3)], y = data[[2]][,c(1,3)], by = "Term", all = TRUE) 
for (i in 3:length(files)) {
  m2 <- merge(x = m2, y = data[[i]][,c(1,3)], by = "Term", all = TRUE)  
}


for (i in 1:nrow(m2)) {
  m2[i,"Term2"] <- gsub('HALLMARK_', '', m2$Term[i])
}


m2[is.na(m2)] <- 0

m2$Term2 <- as.factor(m2$Term2)
m2$Term2 <- factor(m2$Term2, levels=m2$Term2[order(rowSums(m2[,c(2:(ncol(m2)-1))]!=0))])
#m2$Term2 <- factor(m2$Term2, levels=m2$Term2[order(rowSums(m2[,c(2,3,4,5)]!=0))])
##m2$Term2 <- factor(m2$Term2, levels=m2$Term2[order(sum(m2[,2]!=0)*400 + sum(m2[,3]!=0)*300)])



#### NOv 24, 2021: selecting only the pathways that appear more than twice

all_pathways_m2 <- m2
m2 <- m2[rowSums(m2[,c(2:(ncol(m2)-1))]!=0)>2,]


m3 <- melt(m2)
m3$variable <- as.factor(m3$variable)
m3$group <- rep(groups, each=length(unique(m3$Term2)))
#m3$group <- rep(groups, each=nrow(m3)/ngroups)
for (i in 1:nrow(m3)) {
  m3[i,"variable2"] <- gsub('\\.\\S', '', m3$variable[i])
}
m3$variable2 <- factor(m3$variable2, levels = labels1)

## making the pathway order manual


pathway_order <- c( "IL2_STAT5_SIGNALING",
                    "COAGULATION",
                    "APICAL_SURFACE",
                    "MYOGENESIS",
                    "UV_RESPONSE_DN",
                    "UNFOLDED_PROTEIN_RESPONSE",
                    "NOTCH_SIGNALING",
                    "BILE_ACID_METABOLISM",
                    "IL6_JAK_STAT3_SIGNALING",
                    "COMPLEMENT",
                    "ALLOGRAFT_REJECTION",
                    "GLYCOLYSIS",
                    "INFLAMMATORY_RESPONSE",
                    "REACTIVE_OXYGEN_SPECIES_PATHWAY",
                    "FATTY_ACID_METABOLISM",
                    "MTORC1_SIGNALING",
                    "MYC_TARGETS_V2",
                    "MYC_TARGETS_V1",
                    "OXIDATIVE_PHOSPHORYLATION",
                    "INTERFERON_ALPHA_RESPONSE",
                    "INTERFERON_GAMMA_RESPONSE",
                    "ESTROGEN_RESPONSE_LATE",
                    "P53_PATHWAY",
                    "UV_RESPONSE_UP",
                    "XENOBIOTIC_METABOLISM",
                    "CHOLESTEROL_HOMEOSTASIS",
                    "ESTROGEN_RESPONSE_EARLY",
                    "KRAS_SIGNALING_DN",
                    "MITOTIC_SPINDLE",
                    "PEROXISOME",
                    "KRAS_SIGNALING_UP",
                    "G2M_CHECKPOINT",
                    "E2F_TARGETS",
                    "TGF_BETA_SIGNALING",
                    "SPERMATOGENESIS",
                    "HYPOXIA",
                    "TNFA_SIGNALING_VIA_NFKB",
                    "EPITHELIAL_MESENCHYMAL_TRANSITION",
                    "APOPTOSIS"
)

               
m3$group <- factor(m3$group, levels=unique(groups))
m3$Term2 <- factor(m3$Term2, levels = pathway_order)

## Adding all genes
### We may not need this any more
#for (var in unique(m3$variable2)) {
#  m3 <- rbind(m3,data.frame(Term="",
#                            Term2="All differentially expressed genes",
#                            variable="",
#                            value=0,
#                            group="Resistant vs. sensitive clones",
#                            variable2=var))
#}

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")), space="Lab")
m3[m3$value==0,]$value <- NA
# rescale the values so they go over the same range positive and negative
minv <- min(m3$value, na.rm = T)
maxv <- min(5,max(m3$value, na.rm = T))
newlim <- max(abs(minv),abs(maxv))
m3$value <- rescale(m3$value, to = c(-newlim, newlim))

source("ggplot_utils.R")
plot1 <- ggplot(data = m3, aes(variable2, Term2, fill = value))+
  geom_tile(color = "white")+
  #scale_fill_gradient2(low = "blue", high = "red", mid = "grey95", 
  #                     midpoint = 0, space = "Lab", 
  #                     name="NES") +
  scale_fill_gradientn(colours = myPalette(100), na.value = "grey95", name="NES") + 
  #theme_minimal()+ 
  thesis_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  #coord_fixed() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  #scale_x_discrete(position = "top") +
  #facet_wrap(. ~ group, scales="free_x")
  facet_grid(. ~ group, scales="free_x", space="free") +
  theme(panel.spacing = unit(0.5, "mm", data = NULL))
#ggsave(paste0("pathways", set, "-", scol, "-", sthr ,"-", ver, ".pdf"), width=6.5, heigh=7, useDingbats=FALSE)
#print(plot1)

#legend1 <- get_legend(plot1)
#as_ggplot(legend1)
#ggsave("legend1.pdf")

print(plot1)

ggsave(paste0("manuscript_files/", filename), width=7.1, heigh=4.7, useDingbats=FALSE)

### Save the data
m4 <- m3[,c("variable2","Term2","group","value")]
colnames(m4) <- c("comparison", "pathway", "group", "NES")
m4$group <- gsub('\n', ' ', m4$group)

# This is Supplementary Table 8
write.csv(file=paste0("manuscript_files/pathways.csv"), m4, quote = FALSE)




#### Now add the cis/trans


get_cnv <- function(file) {
  cnv <- read.csv(file=file, header=TRUE, row.names=1)
  # add the symbol gene name
  symbs <- read.csv("Symbol_ensembl.csv")
  row.names(symbs) <- symbs$Ensembl
  cnv$Symbol <- symbs[rownames(cnv),]$Symbol
  return(cnv)
}

Pt4_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA609_cnv_mat.csv.gz")
Pt5_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA535_cnv_mat.csv.gz")
Pt6_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA1035_cnv_mat.csv.gz")

pathway_genes <- NULL

prepare_cis_trans <- function(cnv_df, series, data_num, clone1, clone2, type="Rx") {
  
  #degenes <- read.csv(de_file)
  
  #degenes <- degenes[degenes$FDR <= 0.01,]
  
  #####
  pathway <- data[[data_num]]
  cistrans <- cnv_df
  dataset <- labels1[data_num]

  for (term in pathway$Term) {
    cistrans <- cnv_df
    allgenes <- pathway[pathway$Term==term,]$genes
    mygenes <- unlist(strsplit(allgenes,";"))
    #mygenes <- intersect(mygenes,degenes$gene_symbol)  ## I don't think this is needed
    cistrans <- cistrans[cistrans$Symbol %in% mygenes,]
    # select only the columns I need
    cistrans <- na.omit(cistrans[,c(clone1,clone2,"Symbol")])
    cistrans$genetype <- ifelse(cistrans[,1] == cistrans[,2],"trans","cis")
    cistrans$series <- series
    cistrans$pathway <- term
    
    prop_cis <- sum(cistrans$genetype=="cis")/nrow(cistrans)
    prop_trans <- sum(cistrans$genetype=="trans")/nrow(cistrans)
    
    genes_cis <- paste(mygenes[cistrans$genetype=="cis"],collapse = ";")
    genes_trans <- paste(mygenes[cistrans$genetype=="trans"],collapse = ";")
    #nes <- pathway[pathway$Term==term,]$nes
    
    pathway_genes <<- rbind(pathway_genes, data.frame(cistrans[,c("series","pathway","Symbol","genetype")], 
                                                      clone_Rx_or_RxH=clone1, clone_UnRx=clone2, 
                                                      comparison=paste0(type, " vs. UnRx")))
    print(head(pathway_genes))

    ## all pathways together
    ct <- rbind(ct, data.frame(term=term, series=series, whichpath="all", dataset=dataset, Proportion=prop_cis, type=paste0(type, " in cis")))  #, genes=genes_cis))
    ct <- rbind(ct, data.frame(term=term, series=series, whichpath="all", dataset=dataset, Proportion=prop_trans, type=paste0(type, " in trans")))  #, genes=genes_trans))

    ## add the "up pathways
    if (pathway[pathway$Term==term,dataset] > 0) {
      ct <- rbind(ct, data.frame(term=term, series=series, whichpath="up", dataset=paste0(dataset, " up"), Proportion=prop_cis, type=paste0(type, " in cis up")))  #, genes=genes_cis))
      ct <- rbind(ct, data.frame(term=term, series=series, whichpath="up", dataset=paste0(dataset, " up"), Proportion=prop_trans, type=paste0(type, " in trans up")))  #, genes=genes_trans))
    }
    if (pathway[pathway$Term==term,dataset] < 0) {
      ct <- rbind(ct, data.frame(term=term, series=series, whichpath="down", dataset=paste0(dataset, " down"), Proportion=prop_cis, type=paste0(type, " in cis down")))  #, genes=genes_cis))
      ct <- rbind(ct, data.frame(term=term, series=series, whichpath="down", dataset=paste0(dataset, " down"), Proportion=prop_trans, type=paste0(type, " in trans down")))  #, genes=genes_trans))
    }
    
  }
  return(ct)
}

ct <- NULL
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X4", data_num=1, clone1="B", clone2="H", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X4", data_num=2, clone1="A", clone2="H", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X5", data_num=3, clone1="B", clone2="H", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X5", data_num=4, clone1="A", clone2="H", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X5", data_num=5, clone1="A", clone2="H", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X6", data_num=6, clone1="A", clone2="H", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X6", data_num=7, clone1="A", clone2="H", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X7", data_num=8, clone1="A", clone2="H", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt4_cnv, series="Pt4 X7", data_num=9, clone1="A", clone2="H", type="RxH")

#ct <- prepare_cis_trans(cnv_df=Pt5_cnv, series="Pt5", data_num=8, clone1="A", clone2="G")
#ct <- prepare_cis_trans(cnv_df=Pt5_cnv, series="Pt5", data_num=9, clone1="A", clone2="G")
#ct <- prepare_cis_trans(cnv_df=Pt5_cnv, series="Pt5", data_num=10, clone1="A", clone2="G")
ct <- prepare_cis_trans(cnv_df=Pt5_cnv, series="Pt5 X8", data_num=10, clone1="A", clone2="G", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt5_cnv, series="Pt5 X8", data_num=11, clone1="D", clone2="G", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt5_cnv, series="Pt5 X9+1", data_num=12, clone1="A", clone2="G", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt5_cnv, series="Pt5 X9+1", data_num=13, clone1="A", clone2="G", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt5_cnv, series="Pt5 X9+1", data_num=14, clone1="E", clone2="G", type="RxH")

ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X5", data_num=15, clone1="D", clone2="A", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X6", data_num=16, clone1="G", clone2="G", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X6", data_num=17, clone1="B", clone2="G", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X7G", data_num=18, clone1="G", clone2="E", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X7G", data_num=19, clone1="G", clone2="E", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X7G+1", data_num=20, clone1="G", clone2="E", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X7G+1", data_num=21, clone1="G", clone2="E", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X7H", data_num=22, clone1="H", clone2="E", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X7H", data_num=23, clone1="H", clone2="E", type="RxH")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X8H", data_num=24, clone1="H", clone2="E", type="Rx")
ct <- prepare_cis_trans(cnv_df=Pt6_cnv, series="Pt6 X8H", data_num=25, clone1="H", clone2="E", type="RxH")


# correcting the series labels so it makes more sense
pathway_genes[pathway_genes$series=="Pt5 X9+1",]$series <- "Pt5 X10"
pathway_genes[pathway_genes$series=="Pt6 X7G+1",]$series <- "Pt6 X8G"
### This is Supplementary Table 7
write.csv(file=paste0("manuscript_files/pathway_genes.csv"), pathway_genes, quote = FALSE)



### plot data

## plot all pathways

ct_all <- ct[ct$whichpath=="all",]
ct_up_down <- ct[ct$whichpath!="all",]



##### SUpplementary Figure 12
###### Plot only the last time point, for a supplementary figure that includes cis/trans proportions per pathway

small_ptable <- m3[m3$group %in% c("Pt4\nX7","Pt5\nX10","Pt6\nX8G","Pt6\nX8H"),]
uniq <- small_ptable[!is.na(small_ptable$value),]
small_ptable <- small_ptable[small_ptable$Term2 %in% unique(uniq$Term2),]

ggplot(data = small_ptable, aes(variable2, Term2, fill = value))+
  geom_tile(color = "white")+
  #scale_fill_gradient2(low = "blue", high = "red", mid = "grey95", 
  #                     midpoint = 0, space = "Lab", 
  #                     name="NES") +
  scale_fill_gradientn(colours = myPalette(100), na.value = "grey95", name="NES") + 
  #theme_minimal()+ 
  thesis_theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  #coord_fixed() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  #scale_x_discrete(position = "top") +
  #facet_wrap(. ~ group, scales="free_x")
  facet_grid(. ~ group, scales="free_x", space="free") +
  theme(panel.spacing = unit(0.5, "mm", data = NULL))

ggsave(paste0("manuscript_files/pathways_late.pdf"), width=4.4, heigh=5.5, useDingbats=FALSE)


## Plot the cis and trans proportions of genes for each pathway above
small_ct <-ct[ct$series %in% c("Pt4 X7","Pt5 X9+1","Pt6 X7G+1","Pt6 X8H"),]
small_ct <- small_ct[small_ct$whichpath=="all",]
small_ct[small_ct$type=="Rx in cis",]$type <- "in cis"
small_ct[small_ct$type=="Rx in trans",]$type <- "in trans"
small_ct[small_ct$type=="RxH in cis",]$type <- "in cis"
small_ct[small_ct$type=="RxH in trans",]$type <- "in trans"

for (i in 1:nrow(small_ct)) {
  small_ct[i,"Term2"] <- gsub('HALLMARK_', '', small_ct$term[i])
}

small_ct <- small_ct[small_ct$Term2 %in% small_ptable$Term2,]
small_ct$Term2 <- factor(small_ct$Term2, levels = pathway_order)
small_ct$dataset <- factor(small_ct$dataset, levels = unique(small_ct$dataset))

library(stringr)
ct_summary <- read.csv("revision_35_summary.csv")
p4_avg <- ct_summary[ct_summary$series=="Pt4",]$avg_pct_genes
p4_sd <- ct_summary[ct_summary$series=="Pt4",]$sd_pct_genes
p5_avg <- ct_summary[ct_summary$series=="Pt5",]$avg_pct_genes
p5_sd <- ct_summary[ct_summary$series=="Pt5",]$sd_pct_genes
p6_avg <- ct_summary[ct_summary$series=="Pt6",]$avg_pct_genes
p6_sd <- ct_summary[ct_summary$series=="Pt6",]$sd_pct_genes

ct_summary <- ct_summary[ct_summary$comp_type=="Rx expand vs UnRx" & ct_summary$gt=="trans gene",]

## Not using dashlines any more
#small_ct[str_detect(small_ct$series,"Pt4"),]$dashline1 <- (p4_avg - p4_sd*3)/100
#small_ct[str_detect(small_ct$series,"Pt4"),]$dashline2 <- (p4_avg + p4_sd*3)/100
#small_ct[str_detect(small_ct$series,"Pt5"),]$dashline1 <- (p5_avg - p5_sd*3)/100
#small_ct[str_detect(small_ct$series,"Pt5"),]$dashline2 <- (p5_avg + p5_sd*3)/100
#small_ct[str_detect(small_ct$series,"Pt6"),]$dashline1 <- (p6_avg - p6_sd*3)/100
#small_ct[str_detect(small_ct$series,"Pt6"),]$dashline2 <- (p6_avg + p6_sd*3)/100


small_ct[small_ct$type=="in cis",]$type <- "cis"
small_ct[small_ct$type=="in trans",]$type <- "trans"

small_ct$type <- factor(small_ct$type, levels=c("trans", "cis"))

ggplot(data=small_ct, aes(x=Proportion, y=Term2, fill=type)) +
  geom_bar(stat="identity") +
  #coord_flip() +
  # geom_text(
  #   aes(x = term, y = prop, label = prop, group = dataset), 
  #   hjust = -0.5, size = 2,
  #   position = position_dodge(width = 1),
  #   inherit.aes = TRUE
  # ) + 
  #geom_tile(color="white") +
  #theme_minimal() + 
  thesis_theme + 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1))+
  #geom_text(size = 3, position = position_stack(hjust = 0.5)) +
  #geom_text(aes(label = prop, x = pos), size = 3) + 
  #scale_color_manual("Gene type", breaks = c(in-cis, in-trans),
  #                   values=c("red3", "royalblue3")) +
  #scale_fill_manual("Gene type", values = c("in cis" = "brown", "in trans" = "cornflowerblue")) +
  scale_fill_manual("Gene type", values = c("cis" = "chocolate", "trans" = "blue2")) +
  #coord_fixed() +
  theme(axis.title.y=element_blank()) +
  #geom_vline(aes(xintercept=dashline1), linetype="dashed", color = "orange") +
  #geom_vline(aes(xintercept=dashline2), linetype="dashed", color = "orange") +
  facet_grid(. ~ dataset, scales="free_x", space="free") +
  theme(panel.spacing = unit(0.5, "mm", data = NULL)) +
  scale_x_continuous(breaks = c(0,1))
  #theme(axis.text.y = element_blank(),
  #      axis.ticks.y = element_blank(),
  #      axis.title.y = element_blank() )


mean(small_ct[small_ct$series=="Pt4 X7" & small_ct$type=="cis",]$Proportion)
mean(small_ct[small_ct$series=="Pt5 X9+1" & small_ct$type=="cis",]$Proportion)
mean(small_ct[small_ct$series %in% c("Pt6 X7G+1","Pt6 X8H") & small_ct$type=="cis",]$Proportion)

# Find Pt4 pathway with max cis
small_ct[small_ct$type=="cis" & small_ct$Proportion==max(small_ct[small_ct$series=="Pt4 X7" & small_ct$type=="cis",]$Proportion),]
# Find Pt4 pathway with min cis
small_ct[small_ct$type=="cis" & small_ct$Proportion==min(small_ct[small_ct$series=="Pt4 X7" & small_ct$type=="cis",]$Proportion),]

ggsave(paste0("manuscript_files/pathways_cis_trans.pdf"), width=5.7, heigh=4.7, useDingbats=FALSE)



############# Making the boxplots for Figure 5 panel b


ggplot(ct_all, aes(x = dataset, y = Proportion)) + 
  geom_boxplot(aes(fill=type, color=type)) +
  scale_fill_manual(name= "type", values = c("darkorange4", "darkolivegreen","darkorange2","darkolivegreen3"))+
  scale_color_manual(name= "type", values = c("darkorange4", "darkolivegreen","darkorange2","darkolivegreen3"))+
  #scale_color_manual(name = "type", values = c("black", "black", "azure4", "azure4"))+
  #geom_dotplot(
  #  aes(fill = type, color = type), trim = FALSE,
  #  binaxis='y', stackdir='center', dotsize = 0.3,
  #  position = position_dodge(0.8)
  #) +
  labs(x="Sample", y="Proportion") +
  thesis_theme +
  #theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90)) +
  background_grid() + 
  #scale_fill_manual(values = c("Rx Down" = "blue2", "Rx Up" = "darkred", "RxH Down" = "cornflowerblue", "RxH Up" = "coral1")) +
  #scale_fill_manual(values = c("Rx Down" = "blue2", "Rx Up" = "darkred", "RxH Down" = "cornflowerblue", "RxH Up" = "coral1")) +
  facet_grid(~ series, scales="free_x", space="free")


# Mean Pt4
mean(ct_all[ct_all$series %in% c("Pt4 X4", "Pt4 X5", "Pt4 X6", "Pt4 X7") & ct_all$type %in% c("Rx in cis", "RxH in cis"),]$Proportion)
# 33%

mean(ct_all[ct_all$series %in% c("Pt5 X8", "Pt5 X9+1") & ct_all$type %in% c("Rx in cis", "RxH in cis"),]$Proportion)
# 19%

mean(ct_all[ct_all$series %in% c("Pt6 X5", "Pt6 X6", "Pt6 X7G", "Pt6 X7G+1", "Pt6 X7H", "Pt6 X8H") & ct_all$type %in% c("Rx in cis", "RxH in cis"),]$Proportion)
# 3.7%
     
ggsave(paste0("manuscript_files/pathways_minfc",minfc,"_boxplots.pdf"), width=5.8, heigh=2.1, useDingbats=FALSE)

####################################
## Make a plot with the percentages of fixed/reverted/new holiday vs treated
### Figure 5d

calc_perc_reverted <- function (df, series, title, tr, hol) {
  valsT <- get(tr,all_pathways_m2)
  valsH <- get(hol,all_pathways_m2)
  total <- sum(valsT | valsH)
  reversed <- sum(valsT * valsH < 0)   # How many changed color
  overlap <- sum(valsT * valsH > 0)    # How many stayed the same color
  transient <- sum(valsT != 0) - sum(valsT & valsH)
  new <- sum(valsH != 0) - sum(valsT & valsH)
  #df <- rbind(df, data.frame(title=title, feature="total", value=total))
  #df <- rbind(df, data.frame(title=title, feature="overlap", value=overlap))
  #df <- rbind(df, data.frame(title=title, feature="reverted", value=reverted))
  #df <- rbind(df, data.frame(title=title, feature="new", value=new))
  df <- rbind(df, data.frame(series=series, title=title, feature="% Rx only", value=transient/total*100))
  df <- rbind(df, data.frame(series=series, title=title, feature="% RxH only", value=new/total*100))
  df <- rbind(df, data.frame(series=series, title=title, feature="% both same direction", value=overlap/total*100))
  df <- rbind(df, data.frame(series=series, title=title, feature="% both reverse direction", value=reversed/total*100))
  #df <- rbind(df, data.frame(title=title, total=total, overlap=overlap, reverted=reverted,new=new,
  #                           poverlap=overlap/total, preverted=reverted/total, pnew=new/total))
  return(df)
}


clone <- "aware"

if (clone=="aware") {
  
  df <- NULL
  df <- calc_perc_reverted (df, series="Pt4", title="Pt4 1 X5 HAB", tr="X5 2 Rx:A / UnRx:H", hol="X5 1 RxH:B / UnRx:H")
  df <- calc_perc_reverted (df, series="Pt4", title="Pt4 2 X5 HAA", tr="X5 2 Rx:A / UnRx:H", hol="X5 3 RxH:A / UnRx:H")
  df <- calc_perc_reverted (df, series="Pt4", title="Pt4 3 X6 HAA", tr="X6 Rx:A / UnRx:H", hol="X6 RxH:A / UnRx:H")
  df <- calc_perc_reverted (df, series="Pt4", title="Pt4 4 X 7 HAA", tr="X7 Rx:A / UnRx:H", hol="X7 RxH:A / UnRx:H")
  
  #df <- calc_perc_reverted (df, series="Pt5", title="Pt5 X7", tr="X7 Rx:AEJ vs. UnRx:G", hol="X7 RxH:AEJ vs. UnRx:G")
  df <- calc_perc_reverted (df, series="Pt5", title="Pt5 1 X8 GAD", tr="X8 Rx:A / UnRx:G", hol="X8 RxH:D / UnRx:G")
  #df <- calc_perc_reverted (df, series="SA535", title="SA535 X9", tr="SA535 Rx X9 vs. UnRx X9", hol="SA535 RxH X9 vs. UnRx X9")
  df <- calc_perc_reverted (df, series="Pt5", title="Pt5 2 X10 GAA", tr="X10 Rx:A / UnRx:G", hol="X10 RxH:A / UnRx:G")
  df <- calc_perc_reverted (df, series="Pt5", title="Pt5 3 X10 GAE", tr="X10 Rx:A / UnRx:G", hol="X10 RxH:E / UnRx:G")
  
  df <- calc_perc_reverted (df, series="Pt6", title="Pt6 1 X6 GGB", tr="X6 Rx:G / UnRx:G", hol="X6 RxH:B / UnRx:G")
  df <- calc_perc_reverted (df, series="Pt6", title="Pt6 2 X7 EGG", tr="X7 Rx:G / UnRx:E", hol="X7 RxH:G / UnRx:E")
  df <- calc_perc_reverted (df, series="Pt6", title="Pt6 3 X8 EGG", tr="X8 Rx:G / UnRx:E", hol="X8 RxH:G / UnRx:E")
  df <- calc_perc_reverted (df, series="Pt6", title="Pt6 4 X7 EHH", tr="X7 Rx:H / UnRx:E", hol="X7 RxH:H / UnRx:E")
  df <- calc_perc_reverted (df, series="Pt6", title="Pt6 5 X8 EHH", tr="X8 Rx:H / UnRx:E", hol="X8 RxH:H / UnRx:E")
  
  df$feature <- factor(df$feature, levels=c("% Rx only", "% RxH only", "% both same direction", "% both reverse direction"))
  df$series <- factor(df$series, levels=c("Pt4", "Pt5", "Pt6"))
}

library("wesanderson")

ggplot(data=df, aes(x=title, y=value, fill=feature)) +
  geom_bar(stat="identity") + 
  labs(x="Sample", y="% pathways") +
  scale_fill_manual(values = rev(wes_palette("FantasticFox1"))) +
  thesis_theme + 
  background_grid() + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(~ series, scales="free_x", space="free")

ggsave(paste0("manuscript_files/pathways_minfc",minfc,"_reversed_pathways_clone", clone, ".pdf"), width=6, heigh=2, useDingbats=FALSE)



###############################

## Now figuring out the number of pathways that evolved 
## Figure 5e


calc_num_pathways <- function (num, series, title, comp, type="Rx") {
  # type can be "Rx" or "RxH"
  vals <- get(comp,all_pathways_m2)
  num <- rbind (num, data.frame(series=series, title=paste0(title, " ", type), feature=paste0(type, " Up"), value=sum(vals > 0)))
  num <- rbind (num, data.frame(series=series, title=paste0(title, " ", type), feature=paste0(type, " Down"), value=sum(vals < 0)))
  return(num)
}




if (clone=="aware") {
  num <- NULL
  num <- calc_num_pathways (num, series="Pt4 X4", title="Pt4 1 X4 HB", comp="X4 1 Rx:B / UnRx:H", type="Rx")
  num <- calc_num_pathways (num, series="Pt4 X4", title="Pt4 2 X4 HA", comp="X4 2 Rx:A / UnRx:H", type="Rx")
  num <- calc_num_pathways (num, series="Pt4 X5", title="Pt4 1 X5 HB", comp="X5 1 RxH:B / UnRx:H", type="RxH")
  num <- calc_num_pathways (num, series="Pt4 X5", title="Pt4 2 X5 HA", comp="X5 2 Rx:A / UnRx:H", type="Rx")
  num <- calc_num_pathways (num, series="Pt4 X5", title="Pt4 3 X5 HA", comp="X5 3 RxH:A / UnRx:H", type="RxH")
  num <- calc_num_pathways (num, series="Pt4 X6", title="Pt4 X6 HA", comp="X6 Rx:A / UnRx:H", type="Rx")
  num <- calc_num_pathways (num, series="Pt4 X6", title="Pt4 X6 HA", comp="X6 RxH:A / UnRx:H", type="RxH")
  num <- calc_num_pathways (num, series="Pt4 X7", title="Pt4 X7 HA", comp="X7 Rx:A / UnRx:H", type="Rx")
  num <- calc_num_pathways (num, series="Pt4 X7", title="Pt4 X7 HA", comp="X7 RxH:A / UnRx:H", type="RxH")
  
  #num <- calc_num_pathways (num, series="Pt5", title="Pt5 X6", tr="X6 Rx:AG vs. UnRx:G", hol=NULL)
  #num <- calc_num_pathways (num, series="Pt5", title="Pt5 X7", tr="X7 Rx:AEJ vs. UnRx:G", hol="X7 RxH:AEJ vs. UnRx:G")
  num <- calc_num_pathways (num, series="Pt5 X8", title="Pt5 X8 GA", comp="X8 Rx:A / UnRx:G", type="Rx")
  num <- calc_num_pathways (num, series="Pt5 X8", title="Pt5 X8 GD", comp="X8 RxH:D / UnRx:G", type="RxH")
  num <- calc_num_pathways (num, series="Pt5 X9+1", title="Pt5 X9+1 GA", comp="X10 Rx:A / UnRx:G", type="Rx")
  num <- calc_num_pathways (num, series="Pt5 X9+1", title="Pt5 X9+1 GA", comp="X10 RxH:A / UnRx:G", type="RxH")
  num <- calc_num_pathways (num, series="Pt5 X9+1", title="Pt5 X9+1 GE", comp="X10 RxH:E / UnRx:G", type="RxH")
  
  num <- calc_num_pathways (num, series="Pt6 X5", title="Pt6 X5 AB", comp="X5 Rx:B / UnRx:A", type="Rx")
  num <- calc_num_pathways (num, series="Pt6 X6", title="Pt6 X6 GG", comp="X6 Rx:G / UnRx:G", type="Rx")
  num <- calc_num_pathways (num, series="Pt6 X6", title="Pt6 X6 GzB", comp="X6 RxH:B / UnRx:G", type="RxH")
  num <- calc_num_pathways (num, series="Pt6 X7G", title="Pt6 X7 EG", comp="X7 Rx:G / UnRx:E", type="Rx")
  num <- calc_num_pathways (num, series="Pt6 X7G", title="Pt6 X7 EG", comp="X7 RxH:G / UnRx:E", type="RxH")
  num <- calc_num_pathways (num, series="Pt6 X8G", title="Pt6 X8 EG", comp="X8 Rx:G / UnRx:E", type="Rx")
  num <- calc_num_pathways (num, series="Pt6 X8G", title="Pt6 X8 EG", comp="X8 RxH:G / UnRx:E", type='RxH')
  num <- calc_num_pathways (num, series="Pt6 X7H", title="Pt6 X7 EH", comp="X7 Rx:H / UnRx:E", type="Rx")
  num <- calc_num_pathways (num, series="Pt6 X7H", title="Pt6 X7 EH", comp="X7 RxH:H / UnRx:E", type="RxH")
  num <- calc_num_pathways (num, series="Pt6 X8H", title="Pt6 X8 EH", comp="X8 Rx:H / UnRx:E", type="Rx")
  num <- calc_num_pathways (num, series="Pt6 X8H", title="Pt6 X8 EH", comp="X8 RxH:H / UnRx:E", type="RxH")
  
  num$series <- factor(num$series, levels=unique(num$series))
  num$feature <- factor(num$feature, levels=c("Rx Up", "Rx Down", "RxH Up", "RxH Down"))
}



ggplot(data=num, aes(x=title, y=value, fill=feature)) +
  geom_bar(stat="identity") + 
  labs(x="Sample", y="# pathways") +
  thesis_theme + 
  background_grid() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("Rx Down" = "blue2", "Rx Up" = "darkred", "RxH Down" = "cornflowerblue", "RxH Up" = "coral1")) +
  facet_grid(~ series, scales="free_x", space="free")

ggsave(paste0("manuscript_files/pathways_minfc", minfc, "_number_pathways_clone", clone, ".pdf"), width=5.75, height=2.1, useDingbats=FALSE)





