#OpenWorm Omega Turn analysis

library(dplyr)
library(ggplot2)
library(data.table)
library(rhdf5)
library(cowplot)
library(ggstatsplot)
library(tidyquant)

setwd("D:/APW Lab/Food ON OFF/trial 5/opentest2/together")
temp = intersect(list.files(pattern = "\\.hdf5$", full.names = T), list.files(pattern = "features", full.names = T))

my_data = list()
meta = ''
df = ''
worms = list()

for (i in unique(temp)){
  df = h5read(file = i, name = "features_means")
  df = as.data.frame(df)
  #add column with i
  df = dfx
  df = df %>% select(1:3, 5, 694, 661, 708) %>% filter(as.numeric(n_valid_skel) > 300)
  df$file.ID <- i
  worms[[i]] = df %>% distinct(worm_index) %>% pull(worm_index)
  my_data[[i]] = df
  
}

h5closeAll()
agg = data.frame(rbindlist(my_data)) %>% filter(as.numeric(n_valid_skel) > 400) %>%
  mutate(file = stringr::str_extract(file.ID, "(?<=/)[^_]+")) %>% 
  mutate(food = case_when(grepl('f', file) ~ 'food',
                          grepl('s', file) ~ 'starve')) %>%
  mutate(geno = case_when(grepl('a', file) ~ 'AMPK KO',
                          grepl('n', file) ~ 'WT')) %>%
  mutate(group = case_when(food == 'food' & geno == 'WT' ~ 'WT Food',
                           food == 'starve' & geno == 'WT' ~ 'WT Starve',
                           food == 'food' & geno == 'AMPK KO' ~ 'AMPK KO Food',
                           food == 'starve' & geno == 'AMPK KO' ~ 'AMPK KO Starve'))
agg$omega_turns_frequency[is.nan(agg$omega_turns_frequency)]=0
agg = agg %>% mutate(omega = omega_turns_frequency * n_frames, upsilon = upsilon_turns_frequency * n_frames)
avg.agg = agg %>% group_by(group, worm_index) %>%
  summarise(meanO = mean(omega), meanOF = mean(omega_turns_frequency), 
            meanU = mean(upsilon), meanUF = mean(upsilon_turns_frequency), meanD = mean(worm_dwelling),
            meanPC = mean(path_curvature), mean.absPC = mean(path_curvature)) %>% as.data.frame()

ggplot(avg.agg, aes(x = group, y = meanOF)) + geom_boxplot() + geom_marg
ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanO, 
                           pairwise.display = "all", p.adjust.method = "fdr", 
                           marginal.type = "density", 
                           title = "Omega Turn Frequency Between WT Worms ± Food")
