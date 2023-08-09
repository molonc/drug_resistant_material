
data <- read.csv("all_tumour_growths.csv")

data[data$Rx=="Un treated",]$Rx <- "UnRx"
data[data$Rx=="Cisplatin",]$Rx <- "Rx"
data[data$Rx=="Holiday",]$Rx <- "RxH"

get_growth <- function(tgro) {
  tgro <- as.numeric(tgro[!is.na(tgro)])
  tgro2 <- tgro
  
  tgro <- tgro[-1]
  tgro2 <- tgro2[-length(tgro2)]
  return((tgro-tgro2)/3)
}

fdata <- NULL

for (i in 1:nrow(data)) {
  data[i,"Patient"]
  tgro <- data[i,6:length(data[i,])]  ## remove the first 2 points so we include only from the treatment day, 
  ## but to take this info from a file, usually it is the third point
  vals <- get_growth(tgro)
  fdata <- rbind(fdata, data.frame(Patient=data[i,"Patient"],Cycle=data[i,"Cycle"],Rx=data[i,"Rx"],val=vals))
}

fdata$Cycle <- as.factor(fdata$Cycle)
fdata$Rx <- factor(fdata$Rx, levels=c("UnRx","Rx","RxH"))

Pt4data = fdata[fdata[,"Patient"]=="Pt4",]
levels(Pt4data$Cycle)[match("1",levels(Pt4data$Cycle))] <- "X4"
levels(Pt4data$Cycle)[match("2",levels(Pt4data$Cycle))] <- "X5"
levels(Pt4data$Cycle)[match("3",levels(Pt4data$Cycle))] <- "X6"
levels(Pt4data$Cycle)[match("4",levels(Pt4data$Cycle))] <- "X7"

library(ggplot2)
source("ggplot_utils.R")
# Basic box plot
ggplot(Pt4data, aes(x=Cycle, y=val, fill=Rx)) + 
  geom_boxplot(outlier.shape=NA, position = position_dodge(0.8,preserve = "single")) + 
  labs(x="Cycle", y=Growth~rate~per~day~(mm^3)) +
  scale_y_continuous(limits = c(-100, 250)) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  scale_fill_manual(values=c("#6d6e71", "#3d4d9b", "#f9de4f")) + 
  thesis_theme + 
  theme(axis.title.x=element_blank())
ggsave(paste0("tumour_growth_rate_Pt4.pdf"), width=6, heigh=2, useDingbats=FALSE)


Pt6data = fdata[fdata[,"Patient"]=="Pt6",]
levels(Pt6data$Cycle)[match("1",levels(Pt6data$Cycle))] <- "X5"
levels(Pt6data$Cycle)[match("2",levels(Pt6data$Cycle))] <- "X6"
levels(Pt6data$Cycle)[match("3",levels(Pt6data$Cycle))] <- "X7"
levels(Pt6data$Cycle)[match("4",levels(Pt6data$Cycle))] <- "X8"

ggplot(Pt6data, aes(x=Cycle, y=val, fill=Rx)) + 
  geom_boxplot(outlier.shape=NA, position = position_dodge(0.8,preserve = "single")) + 
  labs(x="Cycle", y=Growth~rate~per~day~(mm^3)) +
  scale_y_continuous(limits = c(-100, 250)) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  scale_fill_manual(values=c("#6d6e71", "#3d4d9b", "#f9de4f")) + 
  thesis_theme + 
  theme(axis.title.x=element_blank())

ggsave(paste0("tumour_growth_rate_Pt6.pdf"), width=6, heigh=2, useDingbats=FALSE)

