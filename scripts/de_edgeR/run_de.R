## New runs on 1 Dec 2020, version 6
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp30", series="SA535", version="v6", title="SA535", sample1 = "SA535X6XB03101", clone1="R", clonelabel1="X6 UUT R", sample2="SA535X6XB03101", clone2="Q", clonelabel2="X6 UUT Q"))
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp31", series="SA535", version="v6", title="SA535", sample1 = "SA535X6XB03101", clone1="R", clonelabel1="X6 UUT R", sample2="SA535X6XB03101", clone2="J", clonelabel2="X6 UUT J"))



## comparisons for the pathway analysis
## the following are comparisons that help us understand how the pathways change with more treatments
# comp10 UT vs U
# comp11 UTT vs UT
# comp12 UTTT vs UTT
# comp13 UTTTT vs UTTT

## SA609
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp10", series="SA609", version="v6", title="SA609", sample1 = "SA609X4XB003083", clonelabel1="X4 UT", sample2="SA609X3XB01584", clonelabel2="X3 U"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp11", series="SA609", version="v6", title="SA609", sample1 = "SA609X5XB03230", clonelabel1="X5 UTT", sample2="SA609X4XB003083", clonelabel2="X4 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp12", series="SA609", version="v6", title="SA609", sample1 = "SA609X6XB03404", clonelabel1="X6 UTTT", sample2="SA609X5XB03230", clonelabel2="X5 UTT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp13", series="SA609", version="v6", title="SA609", sample1 = "SA609X7XB03505", clonelabel1="X7 UTTTT", sample2="SA609X6XB03404", clonelabel2="X6 UTTT"))

## to add some of the reproducible sammples
# comp15 should be replicate of comp10

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp15", series="SA609", version="v6", title="SA609", sample1 = "SA609X4XB03084", clonelabel1="X4 UT", sample2="SA609X3XB01584", clonelabel2="X3 U"))
# comp16 should be replicate of comp11
# NOTE: comp16 and comp17 don't really makke sense because the compared sample are from different branches
# I could compare UTTTT vs UT from the replicate branch,
# 3573 (UTTTT replicate) could be compared with 3083 (UT), this would be a replicate of comparing 3505 with 3083
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp16", series="SA609", version="v6", title="SA609", sample1 = "SA609X5XB03230", clonelabel1="X5 UTT", sample2="SA609X4XB03084", clonelabel2="X4 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp17", series="SA609", version="v6", title="SA609", sample1 = "SA609X7XB03573", clonelabel1="X7 UTTTT", sample2="SA609X6XB03404", clonelabel2="X6 UTTT"))



## SA535cisplatin
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp10", series="SA535", version="v6", title="SA535_cisplatin", sample1 = "SA535X6XB03101", clonelabel1="X6 UUT", sample2="SA535X5XB02895", clonelabel2="X5 UU"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp11", series="SA535", version="v6", title="SA535_cisplatin", sample1 = "SA535X7XB03304", clonelabel1="X7 UUTT", sample2="SA535X6XB03101", clonelabel2="X6 UUT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp12", series="SA535", version="v6", title="SA535_cisplatin", sample1 = "SA535X8XB03431", clonelabel1="X8 UUTTT", sample2="SA535X7XB03304", clonelabel2="X7 UUTT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp13", series="SA535", version="v6", title="SA535_cisplatin", sample1 = "SA535X9XB03617", clonelabel1="X9 UUTTTT", sample2="SA535X8XB03431", clonelabel2="X8 UUTTT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp14", series="SA535", version="v6", title="SA535_cisplatin", sample1 = "SA535X10XB03696", clonelabel1="X10 UUTTTTT", sample2="SA535X9XB03617", clonelabel2="X9 UUTTTT"))

## SA535 CX
# comp10 and comp11 are replicates, 2891 and 2894 are derived from the same 2498
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp10", series="SA535", version="v6", title="SA535_CX5461", sample1 = "SA535X5XB02891", clonelabel1="X5 UX", sample2="SA535X4XB02498", clonelabel2="X4 U"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp11", series="SA535", version="v6", title="SA535_CX5461", sample1 = "SA535X5XB02894", clonelabel1="X5 UX", sample2="SA535X4XB02498", clonelabel2="X4 U"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp12", series="SA535", version="v6", title="SA535_CX5461", sample1 = "SA535X6XB03200", clonelabel1="X6 UXX", sample2="SA535X5XB02891", clonelabel2="X5 UX"))
# comp13 doesn't make sense because they are not direct relatives
#rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp13", series="SA535", version="v6", title="SA535_CX5461", sample1 = "SA535X6XB03200", clonelabel1="X6 UXX", sample2="SA535X5XB02894", clonelabel2="X5 UX"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp14", series="SA535", version="v6", title="SA535_CX5461", sample1 = "SA535X8XB03548", clonelabel1="X8 UXXXX", sample2="SA535X6XB03200", clonelabel2="X6 UXX"))

## SA1035
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp10", series="SA1035", version="v6", title="SA1035", sample1 = "SA1035X5XB03015", clonelabel1="X5 UT", sample2="SA1035X4XB02879", clonelabel2="X4 U"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp11", series="SA1035", version="v6", title="SA1035", sample1 = "SA1035X6XB03211", clonelabel1="X6 UTT", sample2="SA1035X5XB03015", clonelabel2="X5 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp12", series="SA1035", version="v6", title="SA1035", sample1 = "SA1035X7XB03338", clonelabel1="X7 UTTT", sample2="SA1035X6XB03211", clonelabel2="X6 UTT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp13", series="SA1035", version="v6", title="SA1035", sample1 = "SA1035X8XB03425", clonelabel1="X8 UTTTT", sample2="SA1035X7XB03338", clonelabel2="X7 UTTT"))


# New runs on 5 Oct 2020, version 5
# comp1 resistant at UTTT vs sensitive at UUUU
# comp2 Resistant clones in the last treated timepoint vs. sensitive clone in the last untreated timepoint
# comp3 whole sample last treated timepoint vs all sample last untreated timepoint
# comp4 Resistant clones in the last treated time point vs the same resistant clones in the first treated time point
# comp5 Resistant clones in the last treated time point UTTTT vs the same resistant clones in UTT

# comp8 whole sample last treated time point vs first treated timepoint (UTTTT vs UT)
# comp9 whole sample last treated time point vs U (UTTTT vs U)



## look at the number of de genes
#de535 <- read.csv("../results/SA535-v5/de/comp2-SA535X9XB03617-A_D-SA535X8XB03664-G_logfc_results.csv")
#de609 <- read.csv("../results/SA609-v5/de/comp2-SA609X7XB03505-R-SA609X7XB03554-C_H_logfc_results.csv")
#de1035 <- read.csv("../results/SA1035-v5/de/comp2-SA1035X8XB03425-G_H-SA1035X8XB03631-E_logfc_results.csv")

#dim(de609[de609$FDR <= 0.01,])
#dim(de535[de535$FDR <= 0.01,])
#dim(de1035[de1035$FDR <= 0.01,])

## SA609
#rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp1", series="SA609", version="v5", title="SA609", sample1 = "SA609X6XB03404", clone1="R", clonelabel1="X6 UTTT R", sample2="SA609X6XB03447", clone2="C_H", clonelabel2="X6 UUUU C,H"))
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp2", series="SA609", version="v5", title="SA609", sample1 = "SA609X7XB03505", clone1="R", clonelabel1="X7 UTTTT R", sample2="SA609X7XB03554", clone2="C_H", clonelabel2="X7 UUUUU C,H"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp3", series="SA609", version="v5", title="SA609", sample1 = "SA609X7XB03505", clonelabel1="X7 UTTTT", sample2="SA609X7XB03554", clonelabel2="X7 UUUUU"))
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp4", series="SA609", version="v5", title="SA609", sample1 = "SA609X7XB03505", clone1="R", clonelabel1="X7 UTTTT R", sample2="SA609X4XB003083", clone2="R", clonelabel2="X4 UT R"))
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp5", series="SA609", version="v5", title="SA609", sample1 = "SA609X7XB03505", clone1="R", clonelabel1="X7 UTTTT R", sample2="SA609X5XB03230", clone2="R", clonelabel2="X5 UTT R"))

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp8", series="SA609", version="v5", title="SA609", sample1 = "SA609X7XB03505", clonelabel1="X7 UTTTT", sample2="SA609X4XB003083", clonelabel2="X4 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp9", series="SA609", version="v5", title="SA609", sample1 = "SA609X7XB03505", clonelabel1="X7 UTTTT", sample2="SA609X3XB01584", clonelabel2="X3 U"))

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp10", series="SA609", version="v5", title="SA609", sample1 = "SA609X4XB003083", clonelabel1="X4 UT", sample2="SA609X3XB01584", clonelabel2="X3 U"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp11", series="SA609", version="v5", title="SA609", sample1 = "SA609X5XB03230", clonelabel1="X5 UTT", sample2="SA609X4XB003083", clonelabel2="X4 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp12", series="SA609", version="v5", title="SA609", sample1 = "SA609X6XB03404", clonelabel1="X6 UTTT", sample2="SA609X5XB03230", clonelabel2="X5 UTT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp13", series="SA609", version="v5", title="SA609", sample1 = "SA609X7XB03505", clonelabel1="X7 UTTTT", sample2="SA609X6XB03404", clonelabel2="X6 UTTT"))

# untreated
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp20", series="SA609", version="v5", title="SA609", sample1 = "SA609X7XB03554", clonelabel1="X7 UUUUU", sample2="SA609X3XB01584", clonelabel2="X3 U"))


## SA535
# clonealign results are not that great for comp1, D_E clone
# rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp1", series="SA535", version="v5", title="SA535", sample1 = "SA535X8XB03431", clone1="A_D", clonelabel1="X8 UTTT A,D", sample2="SA535X8XB03664", clone2="G", clonelabel2="X8 UUUU G"))
# comp2 This has to be redone when I have UUUUU, until then compare with UUUU
#rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp1", series="SA535", version="v5", title="SA535", sample1 = "SA535X8XB03431", clone1="A_D", clonelabel1="X8 UTTT A,D", sample2="SA535X8XB03664", clone2="G", clonelabel2="X8 UUUU G"))

rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp2", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03425", clone1="G_H", clonelabel1="X8 UTTTT G,H", sample2="SA1035X8XB03631", clone2="E", clonelabel2="X8 UUUUU E"))

rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp2", series="SA535", version="v5", title="SA535", sample1 = "SA535X9XB03617", clone1="A_D", clonelabel1="X9 UTTTT A,D", sample2="SA535X8XB03664", clone2="G", clonelabel2="X8 UUUU G"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp3", series="SA535", version="v5", title="SA535", sample1 = "SA535X9XB03617", clonelabel1="X9 UTTTT", sample2="SA535X8XB03664", clonelabel2="X8 UUUU"))
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp4", series="SA535", version="v5", title="SA535", sample1 = "SA535X9XB03617", clone1="A_D", clonelabel1="X9 UTTTT A,D", sample2="SA535X6XB03101", clone2="A_D", clonelabel2="X6 UT A,D"))
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp5", series="SA535", version="v5", title="SA535", sample1 = "SA535X9XB03617", clone1="A_D", clonelabel1="X9 UTTTT A,D", sample2="SA535X7XB03304", clone2="A_D", clonelabel2="X7 UTT A,D"))

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp8", series="SA535", version="v5", title="SA535", sample1 = "SA535X9XB03617", clonelabel1="X9 UTTTT", sample2="SA535X6XB03101", clonelabel2="X6 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp9", series="SA535", version="v5", title="SA535", sample1 = "SA535X9XB03617", clonelabel1="X9 UTTTT", sample2="SA535X5XB02895", clonelabel2="X5 U"))

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp10", series="SA535", version="v5", title="SA535", sample1 = "SA535X6XB03101", clonelabel1="X6 UT", sample2="SA535X5XB02895", clonelabel2="X5 U"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp11", series="SA535", version="v5", title="SA535", sample1 = "SA535X7XB03304", clonelabel1="X7 UTT", sample2="SA535X6XB03101", clonelabel2="X6 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp12", series="SA535", version="v5", title="SA535", sample1 = "SA535X8XB03431", clonelabel1="X8 UTTT", sample2="SA535X7XB03304", clonelabel2="X7 UTT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp13", series="SA535", version="v5", title="SA535", sample1 = "SA535X9XB03617", clonelabel1="X9 UTTTT", sample2="SA535X8XB03431", clonelabel2="X8 UTTT"))

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp20", series="SA535", version="v5", title="SA535", sample1 = "SA535X8XB03664", clonelabel1="X8 UUUU", sample2="SA535X5XB02895", clonelabel2="X5 U"))


## SA1035

rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp2", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03425", clone1="G_H", clonelabel1="X8 UTTTT G,H", sample2="SA1035X8XB03631", clone2="E", clonelabel2="X8 UUUUU E"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp3", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03425", clonelabel1="X8 UTTTT", sample2="SA1035X8XB03631", clonelabel2="X8 UUUUU"))
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp4", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03425", clone1="G", clonelabel1="X8 UTTTT G", sample2="SA1035X5XB03015", clone2="G", clonelabel2="X5 UT G"))
rmarkdown::render("differential-expression-by-clone.Rmd", params = list(cno="comp5", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03425", clone1="G", clonelabel1="X8 UTTTT G", sample2="SA1035X6XB03211", clone2="G", clonelabel2="X6 UTT G"))

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp8", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03425", clonelabel1="X8 UTTTT", sample2="SA1035X5XB03015", clonelabel2="X5 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp9", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03425", clonelabel1="X8 UTTTT", sample2="SA1035X4XB02879", clonelabel2="X4 U"))

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp10", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X5XB03015", clonelabel1="X5 UT", sample2="SA1035X4XB02879", clonelabel2="X4 U"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp11", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X6XB03211", clonelabel1="X6 UTT", sample2="SA1035X5XB03015", clonelabel2="X5 UT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp12", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X7XB03338", clonelabel1="X7 UTTT", sample2="SA1035X6XB03211", clonelabel2="X6 UTT"))
rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp13", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03425", clonelabel1="X8 UTTTT", sample2="SA1035X7XB03338", clonelabel2="X7 UTTT"))

rmarkdown::render("differential-expression-by-sample.Rmd", params = list(cno="comp20", series="SA1035", version="v5", title="SA1035", sample1 = "SA1035X8XB03631", clonelabel1="X8 UUUUU", sample2="SA1035X4XB02879", clonelabel2="X4 U"))



## combine pngs

library(png)
library(grid)
library(ggplot2)
library(gridExtra)

get_plots <- function(files) {
  plots <- lapply(ll <- files,function(x){
    print(x)
    img <- as.raster(readPNG(x))
    rasterGrob(img, interpolate = FALSE)
  })
  return(plots)
}

###############
files <- Sys.glob("../results/monotonically_changed_genes_SA1035/expression_INCREASING_*.png")
plots <- get_plots(files)
ggsave("../results/monotonically_changed_genes_SA1035/SA1035_expression_INCREASING.pdf", marrangeGrob(grobs=plots, nrow=3, ncol=2))
###############
files <- Sys.glob("../results/monotonically_changed_genes_SA1035/expression_DECREASING_*.png")
plots <- get_plots(files)
ggsave("../results/monotonically_changed_genes_SA1035/SA1035_expression_DECREASING.pdf", marrangeGrob(grobs=plots, nrow=3, ncol=2))


