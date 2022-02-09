for(gene in unique(genedf$name)) {
  print(gene)
  ## untreated
  genelm <- genedf[genedf$name==gene & genedf$condition %in% c("UUU","UUUU","UUUUU"),]
  linearMod <- lm(mean ~ time, data=genelm)
  genedf[genedf$name==gene & genedf$condition %in% c("UUU","UUUU","UUUUU"),"slope"] <- linearMod$coefficients[2]
  ## treated
  genelm <- genedf[genedf$name==gene & genedf$condition %in% c("UTT","UTTT","UTTTT"),]  
  linearMod <- lm(mean ~ time, data=genelm)
  genedf[genedf$name==gene & genedf$condition %in% c("UTT","UTTT","UTTTT"),"slope"] <- linearMod$coefficients[2]  
  mon <- (genelm[2,"mean"]-genelm[1,"mean"]) * (genelm[3,"mean"]-genelm[2,"mean"])
  if (mon > 0) {
    m <- TRUE
  } else {
    m <- FALSE
  }
  genedf[genedf$name==gene,"monotonic"] <- m
}


# This above plots the increasing in T / decreasing in U
genedf <- genedf[genedf$monotonic==TRUE,]
genedf_ch <- genedf[!is.na(genedf$slope),]
genes_t_up <- genedf_ch[genedf_ch$slope>0 & genedf_ch$treatment=="Rx","name"]
genes_u_down <- genedf_ch[genedf_ch$slope <= 0 & genedf_ch$treatment=="Un","name"]
common_genes_up <- intersect(genes_t_up, genes_u_down)
mygenes <- NULL
for (gene in common_genes_up) {
  slopeT <- genedf[genedf$name==gene & genedf$condition=="UTT","slope"]
  slopeU <- genedf[genedf$name==gene & genedf$condition=="UUU","slope"]
  mygenes <- rbind(mygenes, data.frame(name=gene, slopeT=slopeT, slopeU=slopeU))
}
write.csv(mygenes, file=paste0(outdir,"/", dataset, "_monotonically_increasing_genes.csv"), quote=FALSE, row.names=FALSE)
labels <- mygenes[mygenes$slopeT>=0.25 | mygenes$slopeU < -0.05,]
png(filename=paste0(outdir,"/", dataset, "_monotonically_increasing_genes.png"), width = 1600, height = 1200)   #, units="in", res=300, pointsize=4)
ggplot(mygenes, aes(x=slopeT, y=slopeU)) + 
  geom_point(size=3,color="black", alpha=0.5) +
  geom_text_repel(data=labels, aes(label = name), size=5) + 
  labs(title=paste0(dataset, " genes mon. increasing in T, decreasing in U"),
       y = "Linear slope untreated",
       x = "Linear slope treated") +
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(size=20),  
    axis.text.y = element_text(size=20),
    plot.title = element_text(size=30, hjust=0))  
dev.off()  