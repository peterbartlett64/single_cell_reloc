# require(ggplot2, ggExtra, arrow, ggpointdensity, dplyr, ggstatsplot, ggsci, dplyr, hrbrthemes, ggsci, scales, stringr)
library(cowplot)
library(tidylog)
library(ggrepel)
library(forcats)
library(CGPfunctions)
# library(networkD3)
library(ggalluvial)
# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(alluvial)
library(dtplyr)
library(arrow)

#Quick header so that files are not getting overwritten with successive tests
if (exists("it_counter")){
  it_counter <- it_counter + 1
} else {
  it_counter <- 0
}


alluvial_ReLocScore <- function(df, protein, mod_name = 'ReLocScore'){
  date <- Sys.Date() #always call the date in function creation to stop overwriting
  used <- read_parquet(df, as_data_frame = F) %>% 
    filter(Protein == protein) %>%
    collect()
  curr <- var_ranks(used, 'Loc_Score') %>% 
    # sankey_var_it()
    
  alluvium_plot <- ggplot(data = curr %>%
                            distinct() %>%
                            mutate(Protein = fct_reorder(Protein, Percentage_reloc_less)),
                          mapping = aes(x = Frames_post_treatment,
                                        stratum = hexile_rank,
                                        alluvium = Protein,
                                        fill = hexile_rank,
                                        label = hexile_rank))+
    scale_fill_npg()+
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = 'darkgray', width = 1/10)+
    geom_stratum(width = .05)+
    theme_sankey()+
    theme(legend.position = "bottom")+
    ggtitle('Global movment of protein through hexile ranks')
  
  ggsave(sprintf("%s_%s_Hex_rank_Alluvium_proteins.png", date, mod_name), plot = alluvium_plot, width = 12, height = 6)
  return(alluvium_plot)
  
}

df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>% 
  filter(Protein %in% c('EXO1')) %>%
  collect()



alluvial_temp <- alluvial_ReLocScore("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", 'EXO1')
