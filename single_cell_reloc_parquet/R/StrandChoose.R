library(dplyr)
library(purrr)

f_choose <- function(sel_col, forward, reverse){
  if (sel_col == 'C>T'){
    return(forward)
  }else if (sel_col == 'C>G'){
    return(reverse)
  }
}

sig <- c('C>G', 'C>T', 'C>T', 'C>G')
fwd <- c('Fwd', 'Fwd', 'Fwd', 'Fwd') # This can be replaced with the actual top strand
rev <- c('Rev', 'Rev', 'Rev', 'Rev') # '' except reverse
df <- data.frame(sig, fwd, rev)


df_strand_chosen <- df %>% 
  mutate(new_column = pmap(list(sig, fwd, rev), f_choose))
