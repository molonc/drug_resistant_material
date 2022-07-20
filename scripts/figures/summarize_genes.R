# Summarize how many genes changed over treatment and how many reversed

library(reshape2)
library(dplyr)
library(ggplot2)

stepFC <- 0.25
FDRthr <- 0.01
maxFC <- 1.5

source("ggplot_utils.R")

clone <- "aware"
#clone <- "unaware"


### get the pathway genes
gmt <- data.frame(read.csv("h.all.v7.0.symbols.gmt",sep='\t',header = FALSE))
gmt[,2] <- NULL

# Specify id.vars: the variables to keep but not split apart on
pgenes <- melt(gmt, id.vars=c("V1"))
pgenes$variable <- NULL
colnames(pgenes) <- c("pathway","gene_symbol")


#######################
#######################
### I need from here on



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
logfc_threshold <- 0.5


#Pt4_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA609_cnv_mat.csv.gz")
#Pt5_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA535_cnv_mat.csv.gz")
#Pt6_cnv <- get_cnv(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/SA1035_cnv_mat.csv.gz")

if (!exists("Pt4_cnv_whole")) {
  Pt4_cnv_whole <- get_cnv2(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/mapped_wholedata_SA609.csv.gz")
}
if (!exists("Pt5_cnv_whole")) {
  Pt5_cnv_whole <- get_cnv2(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/mapped_wholedata_SA535.csv.gz")
}
if (!exists("Pt6_cnv_whole")) {
  Pt6_cnv_whole <- get_cnv2(file="../../materials/dlp_cnv/Fig5_Mirela_cnv/mapped_wholedata_SA1035.csv.gz")
}


#### Use the meta file
#meta <- read.csv("../../materials/comparisons/comparisons_transient.csv")
meta <- read.csv("../../materials/comparisons/comparisons_transient_all.csv")
dir <- "../../materials/comparisons/"


colorcode <-read.csv("../../materials/umap_figs/figs/colorcode_total_v3.csv.gz")

find_reversed_holiday_genes <- function(patient, passage, group=NA, type="transient") {
  ## type can be transient or nontransient
  if (patient=="Pt4") {
    pt <- Pt4_cnv_whole
  } else if (patient=="Pt5") {
    pt <- Pt5_cnv_whole
  } else if (patient=="Pt6") {
    pt <- Pt6_cnv_whole
  }
  thismeta <- meta[meta$series==patient & meta$passage==passage,]
  if (!is.na(group)) {
    thismeta <- thismeta[thismeta$group==group,]
  }
  fileTvsU <- paste0(dir, thismeta[thismeta$comp=="TvsU","filename"])
  fileTvsH <- paste0(dir, thismeta[thismeta$comp=="TvsH","filename"])
  fileHvsU <- paste0(dir, thismeta[thismeta$comp=="HvsU","filename"])
  cloneT <- paste0("cn_median_",thismeta[thismeta$comp=="TvsU","clone1"])
  cloneU <- paste0("cn_median_",thismeta[thismeta$comp=="TvsU","clone2"])
  
  
  deTvsU <- read.csv(fileTvsU)  ## DE
  deTvsH <- read.csv(fileTvsH)  ## DE
  deHvsU <- read.csv(fileHvsU)  ## not DE 
  
  
  if (type=="transient") {
    print ("Transient")
    ## get DE of Rx vs RxH
    deTvsH <- getDE_all(deTvsH, fdr=0.01, logFC=logfc_threshold)
    # also use deTvsU
    #deTvsU <- getDE(deTvsU)
    
    # or include all genes from TvsU
    deTvsU <- getDE_all(deTvsU, fdr=1, logFC=0)
    
  } else {
    print ("Nontransient")
    ## get NON DE of Rx vs RxH
    deTvsH <- getnonDE_all(deTvsH, fdr=0.1)
    
    # get DE between Rx and UnRx
    #deTvsU <- getDE(deTvsU)
    
    # or include all genes from TvsU
    deTvsU <- getDE_all(deTvsU, fdr=0.01, logFC=logfc_threshold)
    
  }

  


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
  m <- select(m,gene_symbol,ensembl_gene_id.x,logFC.x,FDR.x,logFC.y,FDR.y,ct)
  colnames(m) <- c("gene_symbol","ensembl_gene_id","TvsH.logFC","TvsH.FDR","TvsU.logFC","TvsU.FDR","cistrans")

  if (type=="transient") {
    m$direction <- ifelse (m$TvsH.logFC * m$TvsU.logFC > 0, "Towards UnRx", "Away from UnRx")
  } else {
    m$direction <- ifelse(m$TvsU.logFC >0, "Activated", "Repressed")
  }
  m$patient <- patient
  m$passage <- passage
  m$group <- group
  
  return(m)
  #m[is.na(m$ct),]$ct <- "unkn"
  
  #return(m)
  
  #return(list(same_cis=same_direction_cis, same_trans=same_direction_trans, same_unkn=same_direction_unkn, 
  #            opp_cis=opposite_direction_cis, opp_trans=opposite_direction_trans, opp_unkn=opposite_direction_unkn))
}





get_bar_data <- function(data, series, type="transient") {
  df <- NULL
  print(series)
  mp <- meta[meta$series==series,] 
  mp[mp$series=="Pt4" & mp$clone1=="R","clone1"] <-"A"
  mp[mp$series=="Pt4" & mp$clone2=="R","clone2"] <-"A"
  passages <- unique(data[data$patient==series,"passage"])
  for (pa in passages) {
    print(pa)
    groups <- unique(data[data$patient==series & data$passage==pa,"group"])
    for (group in groups) {
      title <- paste0(ifelse(pa=="X10", "X9+1", pa), " UnRx:",mp[mp$passage==pa & mp$group==group,"clone2"][1],
                      " Rx:",mp[mp$passage==pa & mp$group==group,"clone1"][1], 
                      " RxH:",mp[mp$passage==pa & mp$group==group,"clone1"][3])
      
      un <- mp[mp$passage==pa & mp$group==group & mp$comp=="TvsU","clone2"]
      rx <- mp[mp$passage==pa & mp$group==group & mp$comp=="TvsU","clone1"]
      rxh <- mp[mp$passage==pa & mp$group==group & mp$comp=="TvsH","clone2"]
      if (type=="transient") {
        df <- rbind(df, data.frame(series=series, passage=pa, title=title, direction="Towards UnRx cis", group=group,
                                   UnRx=un, Rx=rx, RxH=rxh,
                                   number=nrow(data[data$patient==series & data$passage==pa & data$group==group & data$direction=="Towards UnRx" & data$cistrans=="cis",])))
        df <- rbind(df, data.frame(series=series, passage=pa, title=title, direction="Towards UnRx trans", group=group,
                                   UnRx=un, Rx=rx, RxH=rxh,
                                   number=nrow(data[data$patient==series & data$passage==pa & data$group==group & data$direction=="Towards UnRx" & data$cistrans=="trans",])))
        df <- rbind(df, data.frame(series=series, passage=pa, title=title, direction="Away from UnRx cis", group=group,
                                   UnRx=un, Rx=rx, RxH=rxh,
                                   number=nrow(data[data$patient==series & data$passage==pa & data$group==group & data$direction=="Away from UnRx" & data$cistrans=="cis",])))
        df <- rbind(df, data.frame(series=series, passage=pa, title=title, direction="Away from UnRx trans", group=group,
                                   UnRx=un, Rx=rx, RxH=rxh,
                                   number=nrow(data[data$patient==series & data$passage==pa & data$group==group & data$direction=="Away from UnRx" & data$cistrans=="trans",])))
      
      } else {
        df <- rbind(df, data.frame(series=series, passage=pa, title=title, direction="Activated cis", group=group,
                                   UnRx=un, Rx=rx, RxH=rxh,
                                   number=nrow(data[data$patient==series & data$passage==pa & data$group==group & data$direction=="Activated" & data$cistrans=="cis",])))
        df <- rbind(df, data.frame(series=series, passage=pa, title=title, direction="Activated trans", group=group,
                                   UnRx=un, Rx=rx, RxH=rxh,
                                   number=nrow(data[data$patient==series & data$passage==pa & data$group==group & data$direction=="Activated" & data$cistrans=="trans",])))
        df <- rbind(df, data.frame(series=series, passage=pa, title=title, direction="Repressed cis", group=group,
                                   UnRx=un, Rx=rx, RxH=rxh,
                                   number=nrow(data[data$patient==series & data$passage==pa & data$group==group & data$direction=="Repressed" & data$cistrans=="cis",])))
        df <- rbind(df, data.frame(series=series, passage=pa, title=title, direction="Repressed trans", group=group,
                                   UnRx=un, Rx=rx, RxH=rxh,
                                   number=nrow(data[data$patient==series & data$passage==pa & data$group==group & data$direction=="Repressed" & data$cistrans=="trans",])))
      }
    }
  }
  
  return(df)
}



make_gene_plot <- function(info, type="transient", title) {
  
  library(ggrepel)
  info$label <- 0
  info$color <- 0
  
  if (type=="transient") {
    info[info$TvsH.logFC * info$TvsU.logFC < 0 & info$cistrans=="cis",c("label","color")] <- data.frame(label="Away from UnRx cis", color="darkred")
    info[info$TvsH.logFC * info$TvsU.logFC < 0 & info$cistrans=="trans",c("label","color")] <- data.frame(label="Away from UnRx trans", color="coral2")  
    info[info$TvsH.logFC * info$TvsU.logFC > 0 & info$cistrans=="cis",c("label","color")] <- data.frame(label="Towards UnRx cis", color="deepskyblue4")
    info[info$TvsH.logFC * info$TvsU.logFC > 0 & info$cistrans=="trans",c("label","color")] <- data.frame(label="Towards UnRx trans", color="deepskyblue")
  } else {
    info[info$TvsU.logFC > 0 & info$cistrans=="cis",c("label","color")] <- data.frame(label="Activated cis", color="darkred")
    info[info$TvsU.logFC > 0 & info$cistrans=="trans",c("label","color")] <- data.frame(label="Activated trans", color="coral2")  
    info[info$TvsU.logFC < 0 & info$cistrans=="cis",c("label","color")] <- data.frame(label="Repressed cis", color="deepskyblue4")
    info[info$TvsU.logFC < 0 & info$cistrans=="trans",c("label","color")] <- data.frame(label="Repressed trans", color="deepskyblue")
  }
  
  mycol <- as.character(info$color)
  names(mycol) <- as.character(info$label)  
  
  
  p <- ggplot(info, aes(x = TvsH.logFC, y = TvsU.logFC)) + 
    #geom_point(size = 5, color = "#0099f9") + 
    geom_point(size = 1, aes(color = label), alpha=0.5) +
    scale_color_manual(values=mycol)
    #scale_color_manual(values = c(breaks ="darkred","coral2","deepskyblue4","deepskyblue",
    #                              values= "darkred","coral2","deepskyblue4","deepskyblue")) +
    #thesis_theme + 
  if (type=="transient") {
    p <- p +
      scale_x_continuous(limits = c(-3, 3)) 
  } else {
    p <- p +
      scale_x_continuous(limits = c(-0.8, 0.8)) 
  }

  p <- p +
    scale_y_continuous(limits = c(-3, 3)) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    geom_text_repel(label = info$gene_symbol,  size=1.2, max.overlaps = 50, segment.alpha = 0.2) + 
    theme_bw() + 
    theme(legend.position = "none") +
    labs(
      x = "logFC Rx vs. RxH",
      y = "logFC Rx vs. UnRx",
      title = title
    )
  print(p)
  if (type=="transient") {
    ggsave(paste0("manuscript_files/transient_", title, ".pdf"), width=2.5, heigh=2.2, useDingbats=FALSE)
  } else {
    ggsave(paste0("manuscript_files/nontransient_", title, ".pdf"), width=1.3, heigh=2.2, useDingbats=FALSE)
    #ggsave(paste0("manuscript_files/nontransient_", title, ".pdf"), width=5.3, heigh=5.2, useDingbats=FALSE)
  }
}


make_plot <- function (type="transient") {
  
  data <- NULL
  
  data <- find_reversed_holiday_genes(patient="Pt4", passage="X5", group="1", type=type)
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt4", passage="X5", group="2", type=type))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt4", passage="X6", group="1", type=type))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt4", passage="X7", group="1", type=type))
  ##### SA535  
  #data <- rbind (data, find_reversed_holiday_genes(patient="Pt5", passage="X7", group="1"))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt5", passage="X8", group="1", type=type))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt5", passage="X10", group="1", type=type))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt5", passage="X10", group="2", type=type))
  # SA1035
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt6", passage="X6", group="1", type=type))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt6", passage="X7", group="1", type=type))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt6", passage="X7", group="2", type=type))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt6", passage="X8", group="1", type=type))
  data <- rbind (data, find_reversed_holiday_genes(patient="Pt6", passage="X8", group="2", type=type))
  
  print(head(data))
  
  df <- NULL
  df <- rbind(df, get_bar_data(data,"Pt4",type))
  df <- rbind(df, get_bar_data(data,"Pt5",type))
  df <- rbind(df, get_bar_data(data,"Pt6",type))
  

  
  df$series <- factor(df$series, levels=c("Pt4","Pt5","Pt6"))
  

  write.csv(file=paste0("manuscript_files/data_file_", type,".csv"), data, quote = FALSE)
  
  
  ggplot(data=df, aes(x=title, y=number, fill=direction)) +
    geom_bar(stat="identity") + 
    labs(y = "No. of genes", x = "Sample") +
    #thesis_theme + 
    scale_fill_manual(values = c("darkred","coral2","deepskyblue4","deepskyblue")) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "none") +    # no legend
    #theme(axis.title.x=element_blank(),
    #  axis.text.x=element_blank()) + 
    facet_grid(~ series, scales="free_x", space="free") 
  
  ggsave(paste0("manuscript_files/", type, "_genes_clone", clone, "_logfc", logfc_threshold, ".pdf"), width=2.95, heigh=3.2, useDingbats=FALSE)

  #######################
  # Now make the annotation heatmap  
  if(type=="transient") {
    
    hm <- df[df$direction=="Towards UnRx cis",c("series","passage","UnRx","Rx","RxH")]
    levels(hm$series) <- c("SA609","SA535","SA1035")
    
    hm$id <- paste0("C",1:nrow(hm))
    rownames(hm) <- paste0("C",1:nrow(hm))
    
    hm <- melt(hm, id=c("series","passage","id"))
    colnames(hm) <- c("series","passage","id","treatment","clone")
    hm[hm$series=="SA609" & hm$clone=="R","clone"] <- "A"
    hm$treatment <- factor(hm$treatment, levels=levels(hm$treatment))
    hm$id <- factor(hm$id, levels=paste0("C",1:nrow(hm)))
    
    # add hm$colour
    hm <- merge(x=hm,y=colorcode,by.x=c("series","clone"), by.y=c("datatag","clone_id"), all.x=TRUE)
    hm$group <- paste0(hm$series,"_",hm$clone)
    
    col <- as.character(hm$colour)
    names(col) <- as.character(hm$group)  
    
    
    ggplot(hm, aes(x=id, y=treatment, fill=group)) +
      geom_tile(color="white", lwd=1, linetype=1) +
      geom_text(aes(label=clone), color="white", size=3) +
      scale_fill_manual(values=col) +
      theme_minimal() +
      theme(legend.position = "none") +    # no legend
      #coord_fixed() + 
      facet_grid(~ series, scales="free_x", space="free") 
    
    
    ggsave(paste0("manuscript_files/three-way-comparison-pretty-colors.pdf"), width=3.1, heigh=1.3, useDingbats=FALSE)
  }
  
  for (pt in unique(data$patient)) {
    print(pt)
    for (pass in unique(data[data$patient==pt,"passage"])) {
      print(pass)
      for (gr in unique(data[data$patient==pt & data$passage==pass,"group"])) {
        print(gr)
        mydata <- data[data$patient==pt & data$passage==pass & data$group==gr,]
        print(dim(mydata))
        if (type=="transient") {
          make_gene_plot(mydata, type="transient", title=paste0(pt," ", pass, " ", gr))
        } else {
          make_gene_plot(mydata, type="nontransient", title=paste0(pt," ", pass, " ", gr))
        }
      }
    }
  }
  
}


make_plot(type="transient")
make_plot(type="nontransient")




### selecting some genes
df <- read.csv("manuscript_files/data_file_transient.csv")
df <- df[df$patient == "Pt4",]
sel <- df[abs(df$TvsU.logFC) > 1 & abs(df$TvsH.logFC) > 1,]
sel <- sel[sel$direction=="Towards UnRx",]
sel$gene_symbol


d1 <- read.csv("manuscript_files/data_file_nontransient.csv")
d1$type <- "Treatment induced or repressed"
d2 <- read.csv("manuscript_files/data_file_transient.csv")
d2$type <- "Treatment holiday diverged"
d <- rbind(d1,d2)
write.csv(file=paste0("manuscript_files/dynamic_genes.csv"), d, quote = FALSE)

