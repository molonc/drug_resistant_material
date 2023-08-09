### Apr 27, 2023
## correlation between cn change and logFC of DE comparisons

#source("ggplot_utils.R")
library(wesanderson)
library(ggplot2)
library(readr)
library(tidyverse)
#library(feather)


## Mar 16: to come back to the clones that are considered. Based on clones to ignore in clonealign? Or based on the dlp frequency?


meta <- read.csv("../../materials/comparisons/comparisons_final_Hoa.csv")
dir <- "../../materials/comparisons/results2/"

compute_correlation <- function(file, name, type) {
  
  dedata <- read.csv(paste0(dir,file))
  #dedata <- read.csv("../../../drug_resistance_local/extranalysis/results/SA609-v6/comps2/scrande_SA609_4_SA609_UTTTT_R_UUUUU_H_logfc_results2.csv")
  dedata <- dedata[dedata$FDR<=0.01 & abs(dedata$logFC)>0.25,]
  # keeping only whole-integer cnchanges
  ## dedata <- dedata[dedata$cnchange%%1==0,]
  # in the case of ""scrande_SA1035_32_SA1035_UTU_B_UUU_G_logfc_results2.csv", cnchange is continuous
  dedata$cnchange <- round(dedata$cnchange)
  # setting limits to 3 up and down
  if (length(dedata[dedata$cnchange>3,]$cnchange) > 0) {
    dedata[dedata$cnchange>3,]$cnchange <- 3
  }
  if (length(dedata[dedata$cnchange< -3,]$cnchange) > 0) {
    dedata[dedata$cnchange< -3,]$cnchange <- -3
  }
  mymax <- max(dedata$cnchange)
  mymin <- min(dedata$cnchange)
  
  dedata$cnchange <- as.factor(dedata$cnchange)
  title <- sub("vs.","/",name)
  ggplot(dedata, aes(x=cnchange, y=logFC)) + geom_boxplot(alpha=0.7) +
    stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      y = "log2 fold change",
      x = "CN change",
      title=title
    )
  ggsave(paste0("manuscript_files/correlation_", name,".pdf"), width=2.5, heigh=2, useDingbats=FALSE)
  
    
  bygene_cor <- cor(dedata$logFC,as.numeric(dedata$cnchange))
  
  ### calculating correlation with means
  meancor_data <- NULL
  for (i in mymin:mymax) {
    meancor_data <- rbind(meancor_data, data.frame(mean_logFC=mean(dedata[dedata$cnchange==i,]$logFC), cnchange=i))
  }
  meancor_data <- na.omit(meancor_data)
  mean_cor <- cor(meancor_data$mean_logFC, meancor_data$cnchange)
  
  
  mediancor_data <- NULL
  for (i in mymin:mymax) {
    mediancor_data <- rbind(mediancor_data, data.frame(mean_logFC=median(dedata[dedata$cnchange==i,]$logFC), cnchange=i))
  }
  mediancor_data <- na.omit(mediancor_data)
  median_cor <- cor(mediancor_data$mean_logFC, mediancor_data$cnchange)
  
  return(data.frame(name=name, type=type, bygene_cor=bygene_cor, mean_cor=mean_cor, median_cor=median_cor, filename=file))

}

correlations <- NULL

for (i in 1:nrow(meta)) {
  filename <- meta[i,"result_fn"]
  filename <- sub("results.csv","results2.csv", filename)
  filename <- sub("B-C","B_C", filename)
  if(filename %in% c("scrande_SA1035_3_SA1035_UTT_G_UUU_G_logfc_results2.csv",      #only CN 0
                     "scrande_SA1035_12_SA1035_UTTU_G_UUUU_E_logfc_results2.csv",   # only CN 0 and only 6 of CN-1 and CN1
                     "scrande_SA1035_41_SA1035_UT_D_UU_A_logfc_results2.csv"))      # onlu CN 0 and 2 CN-1
  {    
    next
  }
  name <- paste0(meta[i,"series"], " ", meta[i,"labels_detailed"])
  if (meta[i,"comp_type"]=="untreated") {
    type <- "UnRx / UnRx"
  } else if(meta[i,"labels_short"]=="Rx") {
    type <- "Rx / UnRx"
  } else if(meta[i,"labels_short"]=="RxH") {
    type <- "RxH / UnRx"
  }
  
  print(paste0(name, " ", type, " ", filename))
  
  correlations <- rbind(correlations, compute_correlation(file=filename, name=name, type=type))
  
}




### Now make the summary plots
ggplot(correlations, aes(x=type, y=mean_cor)) + geom_boxplot(alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(0.7,1)+
  labs(
    x = "Comparison type",
    y = "Pearson correlation (means)"
  )
ggsave(paste0("manuscript_files/correlation_summary_mean.pdf"), width=3, heigh=2, useDingbats=FALSE)



ggplot(correlations, aes(x=type, y=median_cor)) + geom_boxplot(alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(0.7,1)+
  labs(
    x = "Comparison type",
    y = "Pearson correlation (medians)"
  )
ggsave(paste0("manuscript_files/correlation_summary_median.pdf"), width=3, heigh=2, useDingbats=FALSE)



