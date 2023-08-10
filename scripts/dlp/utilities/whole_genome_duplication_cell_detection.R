
library(dplyr)
metrics_fn <- '/home/htran/storage/raw_DLP/drug_resistance_DLP/SA501/A95621B/annotation/A95621B_metrics.csv.gz'
df <- data.table::fread(metrics_fn)
dim(df)
summary(as.factor(df$state_mode))
table(df$mean_copy, df$state_mode)
dim(df)
df <- df %>%
  dplyr::filter(quality>=0.75)
sum(df$state_mode==4)
19/529
library(ggplot2)
p <- ggplot(df, aes(x=mean_copy, y=state_mode)) + 
  geom_point(size=1)
p
