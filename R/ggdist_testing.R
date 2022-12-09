require(pacman)
pacman:: p_load(ggdist,dplyr, GGally, ggplot2, ggthemes, ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, stringr, tidyr, ggpmisc)
# pacman:: p_load(httpgd)
# here() #Set the working directory to the current folder
# hgd()
# hgd_browse()

setwd("D:/Microfluidics/RESULTS_ALL/Most_final_collected/Combined_by_perc/Col_with_Abund") #!This is temporary

file = read.csv("EXO1_test.csv")
head(file)

# penetrance <- ggplot(data = file, mapping = aes(x = Protein, y = Percentage_reloc, fill = Percentage_reloc_less)) +
#   geom_text(aes(label = Percentage_reloc_less)) +
#   geom_tile() +
#   xlab(label = "Protein") +
#   ylab(label = "Selected Series")
# penetrance

file$Reloc_yet = file$Reloc_yet + 1

df = data.frame(
 group = c("a", "b", "c"),
 value = rnorm(
   300,
   mean = c(1, 2, 3),
   sd = c(1, 1.5, 1)
 )
)

df2 <- df
df$past = 1
df2$past = 2
df3 <- rbind(df,df2)

file %>% 
  file$Frame_x = as.character(file$Frame_x) %>% 
  ggplot(aes(
    y = Frame_x,
    x = Loc_score))+
  stat_eye(aes(fill = stat(1 < x))) +
  geom_dotsinterval(side = 'bottom', scale = 1.2, height = 0.6, fill = file$Reloc_yet)+
  scale_fill_manual(values = c("#d9d9d9fe", "skyblue"))+
  scale_color_brewer(palette = "Dark2")+
  theme_clean()

df3 %>% 
  ggplot(aes(
    y = group,
    x = value))+
  stat_eye(aes(fill = stat(1 < x))) +
  geom_dotsinterval(side = 'bottom', scale = 1.2, height = 0.6, fill = df3$past)+
  scale_fill_manual(values = c("#d9d9d9fe", "skyblue"))+
  scale_color_brewer(palette = "Dark2")+
  theme_clean()

# hgd_browse()
