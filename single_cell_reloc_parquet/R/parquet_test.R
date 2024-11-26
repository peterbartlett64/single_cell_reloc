require('arrow', 'ggplot2')

setwd("E:/ALL_FINAL/Combined_by_perc")
tf = "LSM12_selected.parquet"
df <- read_parquet(tf)
lib