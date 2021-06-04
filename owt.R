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
my_data2= list()
meta = ''
df = ''
df2 = ''
dfx
worms = list()

for (i in unique(temp)){
  df = h5read(file = i, name = "features_means")
  df = as.data.frame(df)
  #add column with i
  dfx = df
  df = df %>% select(1:3, 5, 694, 661, 708) %>% filter(as.numeric(n_valid_skel) > 300)
  df$file.ID <- i
  worms[[i]] = df %>% distinct(worm_index) %>% pull(worm_index)
  my_data[[i]] = df
  if (i == './nf4_features.hdf5') {df.nf4 = as.data.frame(df)}
  
  df2 = h5read(file = i, name = 'features_timeseries')
  df2 = as.data.frame(df2)
  df2x = df2
  df2 = df2 %>% select(1:4,57)
  df2$file.ID = i
  my_data2[[i]] = df2
}

h5closeAll()
{ agg = data.frame(rbindlist(my_data)) %>% filter(as.numeric(n_valid_skel) > 400) %>%
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
agg$upsilon_turns_frequency[is.nan(agg$upsilon_turns_frequency)]=0
agg = agg %>% mutate(omega = omega_turns_frequency * n_frames, upsilon = upsilon_turns_frequency * n_frames) %>%
  mutate(worm_ID = paste0(worm_index,".",file)) }
avg.agg = agg %>% group_by(group, worm_ID) %>%
  summarise(meanO = mean(omega), meanOF = mean(omega_turns_frequency), 
            meanU = mean(upsilon), meanUF = mean(upsilon_turns_frequency), meanD = mean(worm_dwelling),
            meanPC = mean(path_curvature), mean.absPC = mean(path_curvature)) %>% as.data.frame()

# ggplot(avg.agg, aes(x = group, y = meanOF)) + geom_boxplot()
{
pof = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanOF, 
                           pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE, ylab = "Omega Turn Frequency (1/s)", xlab = '')
po = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanO, 
                                  pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE, ylab = "Omega Turns", xlab = '')
puf = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanUF, 
                                  pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE, ylab = "Upsilon Turn Frequency (1/s)", xlab = '')
pu = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanU, 
                                  pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE, ylab = "Upsilon Turns", xlab = '') }


agg2 = data.frame(rbindlist(my_data2)) %>%
  mutate(file = stringr::str_extract(file.ID, "(?<=/)[^_]+")) %>% 
  mutate(food = case_when(grepl('f', file) ~ 'food',
                          grepl('s', file) ~ 'starve')) %>%
  mutate(geno = case_when(grepl('a', file) ~ 'AMPK KO',
                          grepl('n', file) ~ 'WT')) %>%
  mutate(group = case_when(food == 'food' & geno == 'WT' ~ 'WT Food',
                           food == 'starve' & geno == 'WT' ~ 'WT Starve',
                           food == 'food' & geno == 'AMPK KO' ~ 'AMPK KO Food',
                           food == 'starve' & geno == 'AMPK KO' ~ 'AMPK KO Starve')) %>%
  mutate(skeleton_fit = case_when(skeleton_id == -1 ~ NaN, TRUE ~ skeleton_id)) %>%
  mutate(worm_ID = paste0(worm_index,".",file)) %>% group_by(worm_ID) %>% arrange(group, file, worm_ID, timestamp) %>% 
  mutate(difference_f = timestamp - lead(timestamp),
    difference_mm = motion_modes - lead(motion_modes)) %>%
  mutate(reversal = case_when(difference_f == 1 | 2 | 3 | 4 | 5 & difference_mm == 0 ~ 0,
                              difference_f == 1 | 2 | 3 | 4 | 5 & difference_mm == 1 | 2 ~ 1,
                              TRUE ~ NaN)) %>%
  mutate(RevNum = sum(reversal)) %>% arrange(group, file, worm_ID, timestamp) %>% na.omit() %>% select(1,10,16) %>% distinct()

ggplot(agg2, aes(x = group, y = RevNum)) + geom_boxplot()
ggplot(agg2, aes(x = timestamp, y = RevNum)) + geom_line() + facet_wrap(~ group)
pr = ggstatsplot::ggbetweenstats(data = agg2, x = group, y = RevNum, 
                                 pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE,
                                 xlab = "Group", ylab = "Reversals")

agg.agg = left_join(agg, agg2) %>% ungroup() %>% filter(!is.na(RevNum)) %>% mutate(RevFreq = 10 * RevNum / n_frames)
prf = ggstatsplot::ggbetweenstats(data = agg.agg, x = group, y = RevFreq, 
                                  pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE,
                                  xlab = "Group", ylab = "Reversal Frequency (1/s)")

pt = ggplot() + labs(title = "Movement Statistics in WT Worms ± Food")
plots = plot_grid(pof, po, puf, pu, prf, pr, labels = "AUTO", align = 'hv', ncol = 2)
full_plot = plot_grid(pt, plots, ncol = 1, rel_heights = c(0.04, 1))
full_plot
