# library(dplyr)
# input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/comparisons/'
# df1 <- data.table::fread(paste0(input_dir, 'comparisons_drug_res_v4.csv')) %>% as.data.frame()
# 
# df2 <- data.table::fread(paste0(input_dir, 'comparisons_final.csv'), header=T) %>% as.data.frame()
# 
# dim(df2)
# View(head(df2))
# View(head(df1))


# '/home/htran/storage/images_dataset/merfish_data/evaluation_metrics_results/'
# '/home/htran/storage/images_dataset/merfish_data/evaluation_metrics_results/XP2059/case1_zSlices647'

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(argparse)
  # library(stringr)
})

parser <- ArgumentParser(description = "Plotting evaluation metrics")
parser$add_argument('--eval_criteria', type = 'character', default='1', help="Input dir") #'Across rounds' or 
parser$add_argument('--series', type = 'character', help="Input XP series ID, XPXXXX",default='XP')
parser$add_argument('--input_dir', type = 'character', help="Input dir")
parser$add_argument('--output_dir', type = 'character', help="Output dir",default=NULL)

args <- parser$parse_args()

if(args$eval_criteria=='1'){
  xplt <- 'Across rounds evaluation'  
}else(){
  xplt <- 'z-slices evaluation'  
}

# input_dir <-'/home/htran/storage/images_dataset/merfish_data/evaluation_metrics_results/XP2059/case1_zSlices647/counting_spots_25_r1/'
# FOV <- '025'
# z <- '05'
prefix <- 'merFISH_0'
input_dir <- args$input_dir
# Read all csv files histogram and combine them into 1 file for plotting
fns <- list.files(input_dir, pattern = '*_hist.csv')
fns

counts <- tibble::tibble()
for(f in fns){
  df <- data.table::fread(paste0(input_dir,f)) %>% as.data.frame()
  # f <- gsub(paste0('_',FOV,'_',z,'_hist.csv'),'',f)
  f <- gsub('_hist.csv','',f)
  f <- gsub(prefix,'R',f)
  df$desc <- f
  counts <- dplyr::bind_rows(counts, df)
}
# dim(counts)
# head(counts)

# Counting signals from 250 to 255 intensity values
stat <- counts %>%
  dplyr::group_by(desc) %>%
  dplyr::summarise(Count=sum(Count))
# stat

# Counting signals=255 intensity values only
# stat_255 <- counts %>%
#   dplyr::filter(Value==255) %>%
#   dplyr::group_by(desc) %>%
#   dplyr::summarise(Count=sum(Count))
# stat_255

datatag <- args$series
# datatag <- 'XP2059'
p <- ggplot(stat, aes(x=desc, y=Count)) +
  geom_bar(stat="identity", fill="steelblue", width=0.4)+
  theme_minimal() +
  theme(axis.text.x = element_text(color="black", size=12, angle=90),
        axis.title.x = element_text(color="black", size=13)) + 
  labs(title=paste0('Counting pixels spots ', datatag), x=xplt, y='#pixels-spots')
# p  
if(is.null(args$output_dir)){
  output_dir <- args$input_dir
  # output_dir <- input_dir
}else{
  output_dir <- args$output_dir
  dir.create(output_dir)
}
png(paste0(output_dir,datatag, "_spots_count.png"),height = 2*400, width=2*650,res = 2*72)
print(p)
dev.off()
