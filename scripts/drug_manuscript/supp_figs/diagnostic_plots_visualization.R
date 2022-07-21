
p609 <- p609 + theme(legend.position = 'bottom')
p_total <- cowplot::plot_grid(p609, p1035, p535_cis, p535_cx, nrow = 2, align='hv')
p_total <- cowplot::plot_grid(p609, p1035, p535_cis, p535_cx, ncol = 1, align='v')
png(paste0(save_dir,"edgeR_evaluation_3series_vertical_08_Feb.png"), height = 2*1100, width=2*700,res = 2*72)
print(p_total)
dev.off()

my_font <- "Helvetica"

lg_pos <- "bottom"
# lg_pos <- "none"
thesis_theme <- ggplot2::theme(
  text = element_text(size = 8, hjust = 0.5, family=my_font),
  axis.title.x = element_text(size=8, hjust = 0.5, family=my_font),
  axis.title.y = element_text(size=8, hjust = 0.5, family=my_font),
  # axis.text.x = element_text(color="black",size=7, hjust = 0.5, family=my_font, angle = 90),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size=7, hjust = 0.5, family=my_font),
  plot.title = element_text(color="black",size=9, face="bold", hjust=0, family=my_font),
  plot.subtitle = element_text(size=8, hjust=0, family=my_font),
  legend.title=element_text(size=7, hjust = 0.5, family=my_font),
  legend.text=element_text(size=7, hjust = 0.5, family=my_font),
  strip.text.x = element_text(size=9, family=my_font),
  strip.text.y = element_text(size=9, family=my_font),
  legend.spacing.x = unit(0.1, 'mm'),
  legend.spacing.y = unit(0.1, 'mm'),
  legend.key.height=unit(1,"line"),
  legend.position = lg_pos,
  panel.grid.major = element_blank(), panel.grid.minor = element_blank()
)
library(ggExtra)
library(ggplot2)

res_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/normalization_evaluation/'
datatag <- 'SA609'
c_609 <- readRDS(paste0(res_dir, 'SA609/edgeR_corr_eval_SA609:_UTTTT-A_vs_UUUUU-H.rds'))

# ylb <- 'SCTransform (Rx X7:cloneA - UnRx X7:cloneH)'
ylb <- 'SCTransform normalized exp: Res(X7:A) - Sen(X7:H)'
# xlb <- 'edgeR log2FC (Rx X7:cloneA vs. UnRx X7:cloneH)'
xlb <- 'edgeR log2FC: Res(X7:A) vs. Sen(X7:H)'
cr <- cor.test(c_609$data$exp_mean, c_609$data$logFC, method = "spearman")
plttitle <- paste0(datatag,': edgeR - SCTransform')
subtitle <- paste0('Correlation:',round(cr$estimate,3),', p-value:',round(cr$p.value,3))
c_609_v2 <- c_609 +
  labs(y=ylb, x=xlb, title=plttitle, subtitle = subtitle) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.3) + 
  thesis_theme + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.3) 
  
# c_609_v2 <- c_609 + theme(legend.position="none") + geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.3) + 
  
c_609_v2 <- ggMarginal(c_609_v2, type="boxplot", size=10, color="#636363") 


get_statistic <- function(exp, datatag, tag){
  desc <- paste0(datatag,':',tag,': ')
  desc <- paste0(desc,' mean:',round(mean(exp),2),' sd:',round(sd(exp),2),' median: ', round(median(exp),2))
  print(desc)
}
get_statistic(c_609$data$exp_mean, datatag, 'SCTransform')
get_statistic(c_609$data$logFC, datatag,'edgeR')



c_1035 <- readRDS(paste0(res_dir, 'SA1035/edgeR_corr_eval_SA1035:_UTTTT-H_vs_UUUUU-E.rds'))
datatag <- 'SA1035'
get_statistic(c_1035$data$exp_mean, datatag, 'SCTransform')
get_statistic(c_1035$data$logFC, datatag,'edgeR')


ylb <- 'SCTransform normalized exp: Res(X8:H) - Sen(X8:E)'
xlb <- 'edgeR log2FC: Res(X8:H) vs. Sen(X8:E)'
cr <- cor.test(c_1035$data$exp_mean, c_1035$data$logFC, method = "spearman")
plttitle <- paste0(datatag,': edgeR - SCTransform')
subtitle <- paste0('Correlation:',round(cr$estimate,3),', p-value:',round(cr$p.value,3))

c_1035_v2 <- c_1035 +
  labs(y=ylb, x=xlb, title=plttitle, subtitle = subtitle) +
  # labs(y='SCTransform (Rx X8:cloneH - UnRx X8:cloneE)', x='edgeR log2FC (Rx X8:cloneH vs UnRx X8:cloneE)') +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.3) + 
  thesis_theme + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.3) 

c_1035_v2 <- ggMarginal(c_1035_v2, type="boxplot", size=10, color="#636363")


c_535_cis <- readRDS(paste0(res_dir, 'SA535_cisplatin/edgeR_corr_eval_SA535_cisplatin:_UUTTTT-T_vs_UUUUUU-J.rds'))
# c_535_cis_v2 <- c_535_cis + theme(legend.position="bottom") + 
#   geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.3) + 
#   geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.3) 
datatag <- 'SA535:Cisplatin'
get_statistic(c_535_cis$data$exp_mean, datatag, 'SCTransform')
get_statistic(c_535_cis$data$logFC, datatag,'edgeR')
ylb <- 'SCTransform normalized exp: Res(X10:T) - Sen(X9:J)'
xlb <- 'edgeR log2FC: Res(X10:T) vs. Sen(X9:J)'
cr <- cor.test(c_535_cis$data$exp_mean, c_535_cis$data$logFC, method = "spearman")
plttitle <- paste0(datatag,': edgeR - SCTransform')
subtitle <- paste0('Correlation:',round(cr$estimate,3),', p-value:',round(cr$p.value,3))

c_535_cis_v2 <- c_535_cis +
  labs(y=ylb, x=xlb, title=plttitle, subtitle = subtitle) +
  # labs(y='SCTransform (Rx X10:cloneT - UnRx X9:cloneJ)', x='edgeR log2FC (Rx X10:cloneT vs UnRx X9:cloneJ)',
  #      title="SA535:Cisplatin: edgeR - SCTransform \n Correlation: 0.78") +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.3) + 
  thesis_theme + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.3) 

c_535_cis_v2 <- ggMarginal(c_535_cis_v2, type="boxplot", size=10, color="#636363")



c_535_cx <- readRDS(paste0(res_dir, 'SA535_CX5461/edgeR_corr_eval_SA535_CX5461:_UXXXX-U_vs_UUUUU-J.rds'))
datatag <- 'SA535:CX5461'
get_statistic(c_535_cx$data$exp_mean, datatag, 'SCTransform')
get_statistic(c_535_cx$data$logFC, datatag,'edgeR')

ylb <- 'SCTransform normalized exp: Res(X8:U) - Sen(X9:J)'
xlb <- 'edgeR log2FC: Res(X8:U) vs. Sen(X9:J)'
cr <- cor.test(c_535_cx$data$exp_mean, c_535_cx$data$logFC, method = "spearman",exact=FALSE)
plttitle <- paste0(datatag,': edgeR - SCTransform')
subtitle <- paste0('Correlation:',round(cr$estimate,3),', p-value:',round(cr$p.value,3))

c_535_cx_v2 <- c_535_cx +
  labs(y=ylb, x=xlb, title=plttitle, subtitle = subtitle) +
  # labs(y='SCTransform (Rx X8:cloneU - UnRx X9:cloneJ)', x='edgeR log2FC (Rx X8:cloneU vs UnRx X9:cloneJ)',
  #      title="SA535:CX5461: edgeR - SCTransform \n Correlation: 0.79") +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.3) + 
  thesis_theme + 
  geom_vline(xintercept=0, linetype="dashed", color = "red", size=0.3) 

c_535_cx_v2 <- ggMarginal(c_535_cx_v2, type="boxplot", size=10, color="#636363")

# c_609_v2 <- c_609 + theme(legend.position = 'none')
# c_535_cis_v2 <- c_535_cis + theme(legend.position = 'none')
# p_corr <- cowplot::plot_grid(c_609_v2, c_1035, c_535_cis_v2, c_535_cx, nrow = 2, rel_widths = c(1, 1.35, 1, 1.35),
#                              labels = c('a','b','c','d'))
p_corr <- cowplot::plot_grid(c_609_v2, c_1035_v2, c_535_cis_v2, c_535_cx_v2, 
                             nrow = 2, rel_heights = c(1, 1.18, 1, 1.18),
                             labels = c('a','b','c','d'))

# png(paste0(res_dir,"diagnostic_plot_edgeR_SCTransform_correlation_09_Feb.png"), height = 2*850, width=2*800,res = 2*72)
# print(p_corr)
# dev.off()

ggsave(paste0(res_dir,"diagnostic_plot_edgeR_SCTransform_correlation_26_Feb.pdf"),
       plot = p_corr,
       height = 7.5,
       width = 6.5,
       useDingbats=F)

ggsave(paste0(res_dir,"diagnostic_plot_edgeR_SCTransform_correlation_26_Feb.png"),
       plot = p_corr,
       height = 7.5,
       width = 6.5,
       # useDingbats=F,
       type = "cairo-png",
       dpi=200
)


rc_609_v2 <- rc_609 + theme(legend.position = 'none')
rc_535_cis_v2 <- rc_535_cis + theme(legend.position = 'none')
p_rcorr <- cowplot::plot_grid(rc_609_v2, rc_1035, rc_535_cis_v2, rc_535_cx, nrow = 2, rel_widths = c(1, 1.35, 1, 1.35),
                              labels = c('a','b','c','d'))

png(paste0(res_dir,"diagnostic_plot_edgeR_raw_data_correlation_08_Feb.png"), height = 2*600, width=2*800,res = 2*72)
print(p_rcorr)
dev.off()









