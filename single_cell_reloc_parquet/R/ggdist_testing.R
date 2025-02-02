install.packages(pacman)
require(pacman)
pacman:: p_load(ggdist,dplyr, GGally, ggplot2, ggthemes, ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, stringr, tidyr, ggpmisc)

#, Change the figure display
# pacman:: p_load(httpgd)
# here() #Set the working directory to the current folder
# hgd()
# hgd_browse()

#!This is temporary and will be set based on input from pipeline
setwd("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code") #* Test with a single file


#, Read in the file for testing/run
file = read.csv("ACE2.csv", sep = ',')
# file = read.parquet("EXO1_test.par")

head(file)
describe(file)

# penetrance <- ggplot(data = file, mapping = aes(x = Protein, y = Percentage_reloc, fill = Percentage_reloc_less)) +
#   geom_text(aes(label = Percentage_reloc_less)) +
#   geom_tile() +
#   xlab(label = "Protein") +
#   ylab(label = "Selected Series")
# penetrance
file$Reloc_yet = 0
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
  # file$Frame_x = as.character(file$Frame_x) %>%
  ggplot(aes(
    y = Frame_x,
    x = Loc_score))+
  stat_eye(aes(fill = stat(1 < x))) +
  geom_dotsinterval(side = 'bottom', scale = 1.2, height = 0.6, fill = file$Reloc_yet) +
  scale_fill_manual(values = c("#d9d9d9fe", "skyblue")) +
  scale_color_brewer(palette = "Dark2")+
  theme_clean()
hgd_browse()

df3 %>%
  ggplot(aes(
    y = group,
    x = value))+
  stat_eye(aes(fill = stat(1 < x))) +
  geom_dotsinterval(side = 'bottom', scale = 1.2, height = 0.6, fill = df3$past)+
  scale_fill_manual(values = c("#d9d9d9fe", "skyblue"))+
  scale_color_brewer(palette = "Dark2")+
  theme_clean()

class(df3$value)
# hgd_browse()
