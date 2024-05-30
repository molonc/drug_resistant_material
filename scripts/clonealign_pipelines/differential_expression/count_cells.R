library(SingleCellExperiment)

dir1 <- "/cellassign/fitness-scrna/mirela/data/scrna_human_normed_annotated/"
dir2 <- "/cellassign/fitness-scrna/mirela/data/scrna_human"
dir3 <- "/cellassign/fitness-scrna/mirela/results/tenxqc_normed_yesdoublets/outputs/preprocess/sce_annotated"

count_cells <- function(dir, pattern) {
  files <- list.files(path=dir,pattern = pattern,full.names = TRUE)
  total <- 0
  for (file in files) {
    print(file)
    sce <- readRDS(file)
    n <- ncol(sce)
    print(n)
    total <- total + n
  }
  print(paste0("Number of files ", length(files)))
  return(total)
}

c <- count_cells(dir2, "*.rdata")
print(paste0("Total in ", dir2, " = ", c))

c <- count_cells(dir1, "*.rdata")
print(paste0("Total in ", dir1, " = ", c))

c <- count_cells(dir3, "*.rds")
print(paste0("Total in ", dir3, " = ", c))
