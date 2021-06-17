#OpenWorm Omega Turn analysis for Anezka's mitoOFF experiments

library(dplyr)
library(ggplot2)
library(data.table)
library(rhdf5)
library(cowplot)
library(ggstatsplot)
library(tidyquant)

#for /f "tokens=*" %a in ('dir /b') do ren "%a" "u%a"

{ setwd("D:/APW Lab/Anezka/APW208/Anezka_redo/Tierpsy-210615/Results")
temp = intersect(list.files(pattern = "\\.hdf5$", full.names = T), list.files(pattern = "features", full.names = T))
my_data = list()
my_data2= list()
meta = ''
df = ''
df2 = ''
dfx = ''
worms = list() }

for (i in unique(temp)){
  print(i)
  df = h5read(file = i, name = "features_means")
  df = as.data.frame(df)
  #add column with i
  dfx = df
  df = df %>% select(1:3, 5, 694, 661, 708) %>% filter(as.numeric(n_valid_skel) > 0)
  df$file.ID <- i
  worms[[i]] = df %>% distinct(worm_index) %>% pull(worm_index)
  my_data[[i]] = df
  # if (i == './ALDP-210603_features.hdf5') {df.nf4 = as.data.frame(df)}
  
  df2 = h5read(file = i, name = 'features_timeseries')
  df2 = as.data.frame(df2)
  df2x = df2
  df2 = df2 %>% select(1:4,57)
  df2$file.ID = i
  my_data2[[i]] = df2
}

h5closeAll()

{ 
  agg = data.frame(rbindlist(my_data)) %>% filter(as.numeric(n_valid_skel) > 600) %>%
    mutate(file = stringr::str_extract(file.ID, "(?<=/)[^_]+")) %>% 
    mutate(worm_ID = paste0(worm_index,".",file)) %>%
    mutate(ATR = case_when(grepl('A', file) ~ 'ATR',
                            grepl('X', file) ~ 'unt')) %>%
    mutate(light = case_when(grepl('L', file) ~ 'Light',
                           grepl('N', file) ~ 'Dark')) %>%
    mutate(food = case_when(grepl('F', file) ~ 'Food',
                            grepl('D', file) ~ 'Off-Food')) %>%
    mutate(starve = case_when(grepl('S', file) ~ 'Starved',
                              grepl('P', file) ~ 'Pre')) %>%
    mutate(group = case_when(ATR == 'ATR' & light == 'Light' & food == 'Off-Food' & starve == 'Pre' ~ 'A+L+F-',
                             ATR == 'ATR' & light == 'Light' & food == 'Food' & starve == 'Pre' ~ 'A+L+F+',
                             ATR == 'ATR' & light == 'Dark' & food == 'Off-Food' & starve == 'Pre' ~ 'A+L-F-',
                             ATR == 'ATR' & light == 'Dark' & food == 'Food' & starve == 'Pre' ~ 'A+L-F+',
                             ATR == 'unt' & light == 'Dark' & food == 'Off-Food' & starve == 'Pre' ~ 'A-L-F-',
                             ATR == 'unt' & light == 'Dark' & food == 'Food' & starve == 'Pre' ~ 'A-L-F+',
                             ATR == 'unt' & light == 'Light' & food == 'Off-Food' & starve == 'Pre' ~ 'A-L+F-',
                             ATR == 'unt' & light == 'Light' & food == 'Food' & starve == 'Pre' ~ 'A-L+F+')) %>%
    filter(file != 'ALFP1-210603')
  agg$omega_turns_frequency[is.nan(agg$omega_turns_frequency)]=0
  agg$upsilon_turns_frequency[is.nan(agg$upsilon_turns_frequency)]=0
  agg = agg %>% mutate(omega = omega_turns_frequency * n_frames, upsilon = upsilon_turns_frequency * n_frames)

avg.agg = agg %>% group_by(group, worm_ID) %>%
  summarise(grouping = food,
            meanO = sum(omega), meanOF = sum(omega_turns_frequency),
            meanU = mean(upsilon), meanUF = mean(upsilon_turns_frequency), meanD = mean(worm_dwelling),
            meanPC = mean(path_curvature), mean.absPC = mean(path_curvature)) %>% as.data.frame()
}
write.csv(avg.agg, 'avg.agg600f.filterbad.csv')

{
  pof = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanOF, type = 'np', 
                                    pairwise.display = "significant", p.adjust.method = "holm", 
                                    bf.message = FALSE, ylab = "Omega Turn Frequency (1/s)", xlab = '')
  pof
  }
  
{
agg2 = data.frame(rbindlist(my_data2)) %>%
  mutate(file = stringr::str_extract(file.ID, "(?<=/)[^_]+")) %>% 
  mutate(ATR = case_when(grepl('A', file) ~ 'ATR',
                         grepl('X', file) ~ 'unt')) %>%
  mutate(light = case_when(grepl('L', file) ~ 'Light',
                           grepl('N', file) ~ 'Dark')) %>%
  mutate(food = case_when(grepl('F', file) ~ 'Food',
                          grepl('D', file) ~ 'Off-Food')) %>%
  mutate(starve = case_when(grepl('S', file) ~ 'Starved',
                            grepl('P', file) ~ 'Pre')) %>%
  mutate(group = case_when(ATR == 'ATR' & light == 'Light' & food == 'Off-Food' & starve == 'Pre' ~ 'A+L+F-',
                           ATR == 'ATR' & light == 'Light' & food == 'Food' & starve == 'Pre' ~ 'A+L+F+',
                           ATR == 'ATR' & light == 'Dark' & food == 'Off-Food' & starve == 'Pre' ~ 'A+L-F-',
                           ATR == 'ATR' & light == 'Dark' & food == 'Food' & starve == 'Pre' ~ 'A+L-F+',
                           ATR == 'unt' & light == 'Dark' & food == 'Off-Food' & starve == 'Pre' ~ 'A-L-F-',
                           ATR == 'unt' & light == 'Dark' & food == 'Food' & starve == 'Pre' ~ 'A-L-F+',
                           ATR == 'unt' & light == 'Light' & food == 'Off-Food' & starve == 'Pre' ~ 'A-L+F-',
                           ATR == 'unt' & light == 'Light' & food == 'Food' & starve == 'Pre' ~ 'A-L+F+')) %>%
  mutate(skeleton_fit = case_when(skeleton_id == -1 ~ NaN, TRUE ~ skeleton_id)) %>%
  mutate(worm_ID = paste0(worm_index,".",file)) %>% group_by(worm_ID) %>% arrange(worm_ID, timestamp) %>% na.omit() %>%
  mutate(difference_f = timestamp - lag(timestamp),
         difference_mm = motion_modes - lead(motion_modes)) %>%
  mutate(reversal = case_when(difference_f == 1 & difference_mm == 2 ~ 1,
                              TRUE ~ 0)) %>%
  mutate(RevNum = sum(reversal)) %>% arrange(group, file, worm_ID, timestamp) %>% na.omit() %>% select(c(12,14,18)) %>% distinct()

agg.agg = left_join(agg, agg2) %>% ungroup() %>% filter(!is.na(RevNum)) %>% mutate(RevFreq = (10 * RevNum) / n_frames)
}
# ggplot(agg2, aes(x = timestamp, y = reversal)) + geom_line() + facet_wrap(~ group)
prf = ggstatsplot::ggbetweenstats(data = agg.agg, x = group, y = RevFreq, 
                                  pairwise.comparisons = F, p.adjust.method = "fdr", bf.message = FALSE,
                                  xlab = "Group", ylab = "Reversal Frequency (1/s)")
prf
ggplot(agg.agg, aes(x = group, y = RevFreq)) + geom_boxplot()
# write.csv(agg.agg, 'rev.agg.agg.csv')


