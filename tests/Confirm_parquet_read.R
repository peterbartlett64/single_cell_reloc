confirmation <- function(){
  directory = "D:/Microfluidics/RESULTS_ALL/Most_final_collected/Combined_by_perc/Col_with_Abund/Merged"
  dir(directory)
  setwd(directory)
  df = read_parquet("Final_wAbund.parquet")
  head(df)
  return(TRUE)
}
