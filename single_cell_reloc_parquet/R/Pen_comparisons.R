# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-02-01
#
# Script Name:
#
# Script Description:
#
#
# SETUP ------------------------------------
library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(ggsci)
library(dplyr)
library(hrbrthemes)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(scales)
library(stringr)
library(cowplot)
library(ggstatsplot)
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
library(gt)
library(data.table)
library(purrr)
library(binr)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

#Read in the dataframe and correct the concatonated protein name

name = "All_pos_pt_t_less_percentages_melt"
fix = "-perc_yet"

#, The below reads in the file to be used based on the dictionary of used values and then get the file.
#. Because there was updated penetrance values, Testing with those.
read_correct <- function(name, fix){
  #Read in the data. There were corrections made so that the all this selecting stuff does not have to done in the future
  penetrance_values_95 <- read_parquet(sprintf("D:/ALL_FINAL/95th_percentile/%s.parquet", name), as_data_frame = T) %>%
    mutate(Protein = stringr::str_replace(Protein, fix, "")) %>%
    mutate(Protein = as_factor(Protein))
  # penetrance_values_95$From  <- as_factor("95")

  penetrance_values_99 <- read_parquet(sprintf("D:/ALL_FINAL/99th_percentile/%s.parquet", name), as_data_frame = T)  %>%
    mutate(Protein = stringr::str_replace(Protein, fix, "")) %>%
    mutate(Protein = as_factor(Protein))
  # penetrance_values_99$From  <- as_factor("99")

  #Read in the table of values to use.
  selector <- read_parquet("D:/ALL_FINAL/percentage_to_use_lib.parquet",  as_data_frame = T) %>%
    mutate(Protein = as_factor(Protein))

  selector_95  <- selector %>% filter(Selected_series == "95")
  selector_99  <- selector %>% filter(Selected_series == "99")

  #Join the tables together
  penetrance_values_95 <- left_join(selector_95, penetrance_values_95, by = "Protein")
  penetrance_values_99 <- left_join(selector_99, penetrance_values_99, by = "Protein")

  penetrance_values <- penetrance_values_95 %>%
    bind_rows(penetrance_values_99)
  return(penetrance_values)
}

space_rank = function(value, space, bins){
  # return((ceiling(value/(space/bins))))
  return((bins + 1) - (ceiling(value/(space/bins))))
}

hexile_rank <- function(df, var){
  df <- df %>%
    group_by(Frames_post_treatment) %>%
    mutate(hexile_rank = ntile(Percentage_reloc_less, 6)) %>%
    mutate(bin_rank_space = space_rank(Percentage_reloc_less, 100, 6)) %>%
    mutate(bin_rank_space = ceiling(Percentage_reloc_less / 10)) %>% #This is to place each protein the possible space
    drop_na()
  return(df)
}

#Rank the pene values
pen_ranks <- function(df, space, bins){
  df <- df %>%
    group_by(Frames_post_treatment) %>%
    mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>%
    mutate(cumm_rank = cume_dist(Percentage_reloc_less)) %>%
    mutate(bin_rank_space = space_rank(Percentage_reloc_less, space, bins)) %>%
    mutate(hexile_rank = space_rank(cumm_rank, space = 1, bins = 6)) %>%
    mutate(bin_rank_space = as.integer(bin_rank_space)) %>%
    ungroup()
    # mutate(bin_rank_space = as_factor(bin_rank_space))
    # mutate(bin_rank_space = fct_reorder(bin_rank_space, integer()))


    # mutate(bin_rank_space = ceiling(Percentage_reloc_less / (space / bins))) #%>% #This is to place each protein the possible space
    # drop_na()
  return(df)
}



penetrance_values <- read_parquet("D:/ALL_FINAL/Combined_by_perc/new_percs.parquet") %>% #read_correct("All_pos_pt_t_less_percentages_melt", "-perc_yet") %>%
  rename(Percentage_reloc_less = 'Yet_perc')%>%
  filter(!(Protein %in% c("LCD1", "RAD5", 'DDC2d0224r1', 'Contam.'))) %>%
  pen_ranks(space = 100, bins = 5)



  # mutate(bin_rank_space = ((5 + 1) - (ceiling(Percentage_reloc_less/(100/5)))))
  # pen_ranks(., space = 100, bins = 5)
  # group_by(Frames_post_treatment) %>%
  # mutate(hexile_rank = ntile(Percentage_reloc_less, 6)) %>%
  # mutate(bin_rank_space = space_rank(Percentage_reloc_less, 100, 6)) %>%
  # mutate(bin_rank_space = ceiling(Percentage_reloc_less / 10)) %>% #This is to place each protein the possible space
  # drop_na()

subset_controller <- function(df, min_val = NULL){
  min_val <- ifelse(is.null(min_val), df %>%
                      group_by(Protein) %>%
                      summarise(max = max(Frames_post_treatment)) %>%
                      filter(max >= 30) %>%
                      ungroup() %>%
                      summarise(min = min(max)), min_val)
  subbed <- filter(penetrance_values, Frames_post_treatment %in% c(0, min_val))
  return(subbed)
}

#. I don't think that this is being used
# new <- subset_controller(penetrance_values) %>%
#   mutate(Time_post_treatment = Frames_post_treatment *7.5) %>%
#   arrange(Time_post_treatment) %>%
#   mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>%
#   mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>%
#   mutate(Time_post_treatment = as_factor(Time_post_treatment)) %>%
#   mutate(hexile_rank = as_factor(hexile_rank)) %>%
#   #Put the bin_rank_space in ascending order regardless of when they show up in the data
#   mutate(bin_rank_space = as_factor(bin_rank_space))
  #Put these in ascending order regardless of when they show up in the data

  # mutate(bin_rank_space = fct_relevel(bin_rank_space, sort(integer())))



# mutate(Time_post_treatment = as_factor(Time_post_treatment))
# new$Frames_post_treatment <- as.character(new$Frames_post_treatment)


# newggslopegraph(dataframe = new,
#   Times = Frames_post_treatment,
#   Measurement = Percentage_reloc_less,
#   Grouping = Protein)

# is_alluvia_form(as.data.frame(new), axes = 3:6, silent = TRUE)
# ggplot(as.data.frame(new),
#        aes(y = Percentage_reloc_less, axis1 = Frames_post_treatment, axis2 = hexile_rank)) +
#   geom_alluvium(aes(fill = Admit), width = 1/12) +
#   geom_stratum(width = 1/12, fill = "black", color = "grey") +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
#   scale_fill_brewer(type = "qual", palette = "Set1") +
#   ggtitle("UC Berkeley admissions and rejections, by sex and department")
direction_sankey <- function(start, end){
  return(paste0(start, "->", end))
}

three_bin <- function(val){
  #Don't forget that this is a rank, so high values are low in rankings
  if(val >= 5){
    return("Low")
  }else if(val <= 4 & val>= 3){
    return("Mid")
  }else if(val <= 2){
    return("High")
  }
  else{
    return(NA)

  }
}



sankey_it <- function(df){
  sankeyed_data <- df%>%
    mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>%
    group_by(Protein) %>%
    arrange(Frames_post_treatment, .by_group = T) %>%
    # group_by(Protein) %>%
    # ungroup() %>%
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>%
    mutate(three_bin = map(hexile_rank, three_bin)) %>%
    mutate(nxt_three = three_bin[row_number()+1]) %>%

    #* Generate the next hexile rank
    mutate(nxt_hex = hexile_rank[row_number()+1]) %>%
    #* Generate the next bin rank
    mutate(nxt_decSp = bin_rank_space[row_number()+1]) %>%

    #Get the difference between the prior frame present (there can be spaced for subseting) and the current
    mutate(diff_decSp = bin_rank_space - nxt_decSp) %>%
    mutate(diff_decSp = replace_na(diff_decSp, 0)) %>% #This is to deal with the final value where there is nothing to compare against
    mutate(diff_hex = hexile_rank - nxt_hex) %>%
    mutate(diff_hex = replace_na(diff_hex, 0)) %>%

    mutate(diff_hex_dup = diff_hex) %>%

    mutate(hexile_cpy = hexile_rank) %>%
    mutate(hexile_rank = as_factor(hexile_rank))%>%

    mutate(hexile_rank = recode_factor(hexile_rank, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
    mutate(bin_rank_space = as_factor(bin_rank_space)) %>%
    #Put these in ascending order regardless of when they show up in the data
    # mutate(bin_rank_space = fct_relevel(bin_rank_space, sort(integer()))) %>%
    # mutate(hexile_rank = fct_reorder(hexile_rank, hexile_cpy, .desc = T)) %>%
    # mutate(nxt_hex = fct_reorder(nxt_hex, Percentage_reloc_less, .desc = T)) %>%

    mutate(nxt_x = Frames_post_treatment[row_number()+1]) %>%
    mutate(diff_hex = as_factor(diff_hex)) %>%
    mutate(diff_hex = recode_factor(diff_hex, '0' = "No Change", '1' = "Up 1", '2' = "Up 2", '3' = "Up 3", '4' = "Up 4", '5' = "Up 5", '-1' = "Down 1", '-2' = "Down 2", '-3' = "Down 3", '-4' = "Down 4", '-5' = "Down 5")) %>%

    mutate(nxt_hex = as_factor(nxt_hex)) %>%
    mutate(nxt_hex  = recode_factor(nxt_hex, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
    mutate(Changes_full = map2(hexile_rank, nxt_hex, direction_sankey)) %>%
    mutate(Changes_full = unlist(Changes_full)) %>%
    mutate(Changes_full = as_factor(Changes_full)) %>%

    mutate(Changes_three = map2(three_bin, nxt_three, direction_sankey)) %>%
    mutate(Changes_three = unlist(Changes_three)) %>%
    mutate(Changes_three = as_factor(Changes_three)) %>%
    ungroup()
    # mutate(diff_hex = fct_reorder(diff_hex, diff_hex_dup, .desc = T)) %>%
    # drop_na() #This is to get rid of the last row for each protein. As written this can be used for several time points instead of just the first and last

  return(sankeyed_data)
}


sankey_var_it <- function(df, var){
  sankeyed_data <- df%>%
    mutate(var = as.integer(var)) %>%
    group_by(Protein) %>%
    arrange(Frames_post_treatment, .by_group = T) %>%
    # group_by(Protein) %>%
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>%
    mutate(Time_post_treatment = as_factor(Time_post_treatment)) %>%
    mutate(nxt_hex = hexile_rank[row_number()+1]) %>%
    mutate(nxt_decSp = bin_rank_space[row_number()+1]) %>%
    mutate(diff_decSp = bin_rank_space - nxt_decSp) %>%
    mutate(changes = map2_chr(direction_sankey, bin_rank_space, nxt_decSp)) %>%
    mutate(diff_decSp = replace_na(diff_decSp, 0)) %>%

    mutate(diff_hex = hexile_rank - nxt_hex) %>%
    mutate(diff_hex = replace_na(diff_hex, 0)) %>%
    # mutate(diff_hex_dup = diff_hex) %>%

    # mutate(hexile_cpy = hexile_rank) %>%
    mutate(hexile_rank = as_factor(hexile_rank))%>%
    mutate(hexile_rank = recode_factor(hexile_rank, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
    # mutate(bin_rank_space = as_factor(bin_rank_space)) %>%
    #Put these in ascending order regardless of when they show up in the data
    # mutate(bin_rank_space = fct_relevel(bin_rank_space, sort(integer()))) %>%
    mutate(nxt_hex = as_factor(nxt_hex)) %>%
    # mutate(hexile_rank = fct_reorder(hexile_rank, integer(), .desc = T)) %>%
    # mutate(nxt_hex = fct_reorder(nxt_hex, var, .desc = T)) %>%
    mutate(nxt_x = Frames_post_treatment[row_number()+1]) %>%
    mutate(diff_hex = as_factor(diff_hex)) %>%
    mutate(diff_hex = recode_factor(diff_hex, '0' = "No Change", '1' = "Up 1", '2' = "Up 2", '3' = "Up 3", '4' = "Up 4", '5' = "Up 5", '-1' = "Down 1", '-2' = "Down 2", '-3' = "Down 3", '-4' = "Down 4", '-5' = "Down 5")) %>%
    # mutate(diff_hex = fct_reorder(diff_hex, integer(), .desc = T)) %>%
    # mutate(nxt_hex  = recode_factor(hexile_rank, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
    # mutate(changes = map2(hexile_rank, nxt_hex, direction_sankey)) %>%
    # mutate(changes = unlist(changes)) %>%
    # mutate(changes = as_factor(changes)) %>%
    ungroup()
  # drop_na() #This is to get rid of the last row for each protein. As written this can be used for several time points instead of just the first and last

  return(sankeyed_data)
}

sankeyed_data <- sankey_it(subset_controller(penetrance_values, 32) %>%
                             mutate(Time_post_treatment = Frames_post_treatment *7.5) #. This would normally be present. If it is not, then make it again
                           )

head(sankeyed_data) #just a quick look

#try something new


#, Funciton to generage requested Alluvium Flow with variable name
protein_alluvium_hexile <- function(df, mod_name = ""){
  date <- Sys.Date()
  alluvium_plot <- ggplot(data = df %>%
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
  ggsave(sprintf("%s_%s_Hex_rank_Alluvium_proteins.eps", date, mod_name), plot = alluvium_plot, width = 12, height = 6)
  return(alluvium_plot)
}

protein_alluvium <- function(df, var_name){
  date <- Sys.Date()
  alluvium_plot <- ggplot(data = df %>%
                            distinct() %>%
                            mutate(Protein = fct_reorder(Protein, Percentage_reloc_less)),
                          mapping = aes(x = Frames_post_treatment,
                                        stratum = .data[[var_name]],
                                        alluvium = Protein,
                                        fill = .data[[var_name]],
                                        label = .data[[var_name]]))+
    scale_fill_npg()+
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = 'darkgray', width = 1/10)+
    geom_stratum(width = .05)+
    theme_sankey()+
    theme(legend.position = "bottom")+
    ggtitle(sprintf('Global movment of protein through %s', var_name))

  ggsave(sprintf("%s_%s_Hex_rank_Alluvium_proteins.png", date, var_name), plot = alluvium_plot, width = 12, height = 6)
  ggsave(sprintf("%s_%s_Hex_rank_Alluvium_proteins.eps", date, var_name), plot = alluvium_plot, width = 12, height = 6)
  return(alluvium_plot)
}

# sankeyed_data %>%
#   map(.x = ., .y = c('hexile_rank', 'bin_rank_space'),
#        ~protein_alluvium(.x, .y)
#   )

protein_alluvium(sankeyed_data, 'hexile_rank')
protein_alluvium(sankeyed_data, 'bin_rank_space')



get_member_table <- function(df, by = 'hexile_rank'){
  df <- df %>%
    filter(!(Protein %like% 'DDC2')) %>%
    filter(!(Protein %like% 'd0210')) %>%
    # filter((as.integer(diff_hex) < -1 )| (as.integer(diff_hex) > 1)) %>%
    filter(Frames_post_treatment == '0') %>%
    drop_na() %>%
    mutate(Changes_three = fct_relevel(Changes_three, "High->High", "High->Mid", "High->Low", "Mid->High", "Mid->Mid", "Mid->Low", "Low->High", "Low->Mid", "Low->Low")) %>% 
    mutate(Changes_full = fct_reorder(Changes_full, hexile_cpy)) %>%
    as.data.table()
    # filter((diff_hex_dup > 0 )| (diff_hex_dup < 0))
    # mutate(hexile_rank = fct_relevel(hexile_rank, sort(integer())))#Remove the proteins which did not move
  table <- df[, .("Protein Members" = list(unlist(unique(Protein))), N = .N), by = by] %>%
    # arrange(Changes_three) %>%
    # arrange(Changes_full) %>% 
    gt()
    

  gtsave(data = table, filename = sprintf("%s_%s_penHet_table.html", date, by))
  gtsave(data = table, filename =  sprintf("%s_%s_penHet_table.pdf", date, by))

    # row_group_order(groups = levels(df$hexile_rank))
  return(table)
}


get_member_table(sankeyed_data, by = 'diff_hex')
get_member_table(sankeyed_data, by = 'bin_rank_space')
get_member_table(sankeyed_data, by = 'Changes_full')
get_member_table(sankeyed_data, by = 'Changes_three')



#testing again
df <- sankeyed_data %>%
  filter(!(Protein %like% 'DDC2')) %>%
  filter(!(Protein %like% 'd0210')) %>%
  # filter((as.integer(diff_hex) < -1 )| (as.integer(diff_hex) > 1)) %>%
  filter(Frames_post_treatment == '0') %>%
  mutate(changes = fct_reorder(changes, hexile_cpy)) %>%
  drop_na() %>%
  as.data.table()



sankey_shift <- ggplot(data = sankeyed_data %>% make_long(hexile_rank, diff_hex),
       mapping = aes(x = Frames_post_treatment,
                     next_x = nxt_x,
                     node = hexile_rank,
                     next_node = nxt_hex,
                     group = Protein,
                     fill = factor(hexile_rank)))+
  geom_sankey()+
  geom_sankey_label()+
  theme_sankey()
sankey_shift

sankey_sub <- sankeyed_data %>%
  filter((diff_hex_dup > 0 )| (diff_hex_dup < 0)) #Remove the proteins which did not move


get_types <- function(df){
  bet_het_n <- df %>%
    filter(diff_hex_dup < 0) %>%
    nrow()
  bet_het_members <- df %>%
    filter(diff_hex_dup < 0) %>%
    unique(.$Protein)
  #   reframe(N = nrow(.), Proteins = unique(.$Protein))
  # misc_het <- df %>%
  #   filter(diff_hex_dup > 0) %>%
  #   reframe(N = nrow(.), Proteins = unique(.$Protein))
  # homogenous <- df %>%
  #   filter(diff_hex_dup == 0) %>%
  #   reframe(N = nrow(.), Proteins = unique(.$Protein))

  return(table_new)
}

x <- get_types(sankeyed_data)
sankeyed_data %>%
  filter(diff_hex_dup < 0) %>%
  distinct(Protein) %>%
  as.list()

table_results <- gt(sankeyed_data) %>%
  tab_header(
    title = "Assigned Schemas",
    subtitle = "Movement of proteins through hexile ranks"
  ) %>%
  tab_source_note(
    source_note = "Data from the penetrance values non-shifters removed"
  ) %>%
  tab_stubhead(label = "Relocation Schema") %>%
  tab_row_group(label = 'Homogenous',
                rows = matches(filter(sankeyed_data, diff_hex_dup < 0)$Protein)) #%>%
  # tab_row_group(label = 'Bet-hedging Relocation',
  #               rows = slice(filter(sankeyed_data, diff_hex_dup < 0))) %>%
  # tab_row_group(label = 'Miscellaneous Relocation',
  #               rows = slice(filter(sankeyed_data, diff_hex_dup > 0)))
  # tab_spanner(
  #   label = "Description of Members",
  #   columns = vars(Protein),
  #   locations = cells_body(rows())
  # )




long_sankey <- sankey_sub %>%   #sankeyed_data %>%
  ungroup() %>%
  filter(Frames_post_treatment == '0')%>%
  make_long(., hexile_rank, diff_hex) %>%
  drop_na(node) %>%
  mutate(node = as_factor(node))
long_sankey$node <- factor(long_sankey$node,
                           levels=c('Sixth','Fifth','Fourth','Third',
                                    'Second','Top', "Down 5", "Down 4",
                                    "Down 3", "Down 2", "Down 1", "Up 1",
                                    "Up 2", "Up 3", "Up 4", "Up 5")) # Got rid of the No Change for subset
                           # levels=c('Sixth','Fifth','Fourth','Third','Second','First', "Down 5", "Down 4", "Down 3", "Down 2", "Down 1", "No Change", "Up 1", "Up 2", "Up 3", "Up 4", "Up 5"))

sankey_shift <- ggplot(data = long_sankey,
       mapping = aes(x = x,
                     next_x = next_x,
                     node = node,
                     next_node = next_node,
                     fill = node,
                     label = node))+
  geom_sankey()+
  # geom_sankey_label()+
  theme_sankey()+
  geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
  # scale_fill_npg()
sankey_shift

ggsave(sprintf('%s_Sankey_shift.png', date), plot = sankey_shift, width = 12, height = 6)




df <- mtcars %>%
  make_long(cyl, vs, am, gear, carb)



penetrance_values$Frames_post_treatment <- as_factor(penetrance_values$Frames_post_treatment)


#Function to compare the starting penetrance of each protein. Get an idea of the extent of priming.

pen_comparison <- function(df, time_sub){
  ggplot(data = filter(df, Time_post_treatment == time_sub) %>%
                         distinct() %>%
                         mutate(Protein = fct_reorder(Protein, Percentage_reloc_less)),
         mapping = aes(x = Protein, y = Percentage_reloc_less))+
    geom_point(mapping = aes(color = Selected_series), alpha = 0.5) + #inherit the aesthetics from the parent ggplot
    theme(legend.position = "none",
          axis.text.x=element_blank())+
    scale_fill_npg()
    # geom_text_repel(aes(label = Protein))#remove the legend
}



pen_comparison(penetrance_values, 0)
pen_comparison(penetrance_values, 1)


#Testing the

df_testing <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F)





#Create the plot of proteins which have a GO term of interest
# ggplot(data = filter(df, Protein %in% Protein %in% c("RFA1", "HTA2", "RDH54", "POL30",
#                                               "NPL4", "CHK1", "MRC1", "MSH3", "MMS21",
#                                               "RAD51", "RAD24", "RTT107", "RRD1", "SLD2",
#                                               "DOA1", "ECO1", "RPN4", "SUB2", "RAD57",
#                                               "DBF4", "PPH3", "RAD55", "CDC1", "RAD9",
#                                               "XRS2", "LRS4", "SLD3", "INO80", "RAD54",
#                                               "SAE2", "ZIP2", "ARP4", "DPB11", "SRS2",
#                                               "CDC6", "RFC2", "SLX4", "TOP3", "NEJ1",
#                                               "RAD33", "RAD52", "NAM7", "YKU80", "PSO2",
#                                               "SGS1", "MRE11", "YKU70", "MGS1", "RAD50",
#                                               "RFC3", "RFC4", "RTS1", "EXO1", "ULS1",
#                                               "REV1", "RAD53", "DDC1", "RIM1")),
#        mapping = aes(x = Frames_post_treatment, y = Loc_score, group = 'Protein'))+
#   geom_line() #inherit the aesthetics from the header