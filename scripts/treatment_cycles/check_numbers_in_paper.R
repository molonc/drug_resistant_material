cells <- read.csv("SuppTable1_number_of_cells_per_sample.csv", header=TRUE)
print(paste0("Total filtered: ", sum(cells$Counts_no_doublets)))

print(paste0("Pt1 initial: ", sum(cells[cells$Patient=="Pt1",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt1",]$Counts_no_doublets)))
print(paste0("Pt2 initial: ", sum(cells[cells$Patient=="Pt2",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt2",]$Counts_no_doublets)))
print(paste0("Pt3 initial: ", sum(cells[cells$Patient=="Pt3",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt3",]$Counts_no_doublets)))

print(paste0("Pt4 UnRx initial: ", sum(cells[cells$Patient=="Pt4" & cells$Drug=="Untreated",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt4" & cells$Drug=="Untreated",]$Counts_no_doublets)))
print(paste0("Pt4 Rx initial: ", sum(cells[cells$Patient=="Pt4" & cells$Drug=="Cisplatin",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt4" & cells$Drug=="Cisplatin",]$Counts_no_doublets)))
print(paste0("Pt4 RxH initial: ", sum(cells[cells$Patient=="Pt4" & cells$Drug=="Cisplatin-holiday",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt4" & cells$Drug=="Cisplatin-holiday",]$Counts_no_doublets)))

print(paste0("Pt5 UnRx initial: ", sum(cells[cells$Patient=="Pt5" & cells$Drug=="Untreated",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt5" & cells$Drug=="Untreated",]$Counts_no_doublets)))
print(paste0("Pt5 Rx initial: ", sum(cells[cells$Patient=="Pt5" & cells$Drug=="Cisplatin",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt5" & cells$Drug=="Cisplatin",]$Counts_no_doublets)))
print(paste0("Pt5 RxH initial: ", sum(cells[cells$Patient=="Pt5" & cells$Drug=="Cisplatin-holiday",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt5" & cells$Drug=="Cisplatin-holiday",]$Counts_no_doublets)))

print(paste0("Pt6 UnRx initial: ", sum(cells[cells$Patient=="Pt6" & cells$Drug=="Untreated",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt6" & cells$Drug=="Untreated",]$Counts_no_doublets)))
print(paste0("Pt6 Rx initial: ", sum(cells[cells$Patient=="Pt6" & cells$Drug=="Cisplatin",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt6" & cells$Drug=="Cisplatin",]$Counts_no_doublets)))
print(paste0("Pt6 RxH initial: ", sum(cells[cells$Patient=="Pt6" & cells$Drug=="Cisplatin-holiday",]$Initial_counts), " filtered: ", sum(cells[cells$Patient=="Pt6" & cells$Drug=="Cisplatin-holiday",]$Counts_no_doublets)))

print(paste0("Total UnRx initial: ", sum(cells[cells$Drug=="Untreated",]$Initial_counts), " filtered: ", sum(cells[cells$Drug=="Untreated",]$Counts_no_doublets)))
print(paste0("Total Rx initial: ", sum(cells[cells$Drug=="Cisplatin",]$Initial_counts), " filtered: ", sum(cells[cells$Drug=="Cisplatin",]$Counts_no_doublets)))
print(paste0("Total RxH initial: ", sum(cells[cells$Drug=="Cisplatin-holiday",]$Initial_counts), " filtered: ", sum(cells[cells$Drug=="Cisplatin-holiday",]$Counts_no_doublets)))



#############
## Figure 4
data <- read.csv("manuscript_files/dynamic_genes.csv", header=TRUE)

data <- data[data$type=="Genes that remain fixed after drug withdrawal",]


# to measure the last time point, group 1 for Pt5
table(data[data$patient=="Pt5" & data$passage=="X10" & data$group=="1",]$direction)
## 43%

# to measure the last time point, group 1 for Pt6
table(data[data$patient=="Pt6" & data$passage=="X8" & data$group=="1",]$direction)
## 34%
