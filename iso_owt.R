#OpenWorm Omega Turn analysis

library(dplyr)
library(ggplot2)
library(data.table)
library(rhdf5)
library(cowplot)
library(ggstatsplot)
library(tidyquant)

setwd("D:/APW Lab/Tierp test vids/Jun/Results")
temp = intersect(list.files(pattern = "\\.hdf5$", full.names = T), list.files(pattern = "features", full.names = T))

my_data = list()
my_data2= list()
meta = ''
df = ''
df2 = ''
dfx = ''
worms = list()

for (i in unique(temp)){
  df = h5read(file = i, name = "features_means")
  df = as.data.frame(df)
  #add column with i
  dfx = df
  df = df %>% select(1:3, 5, 694, 661, 708) %>% filter(as.numeric(n_valid_skel) > 1)
  df$file.ID <- i
  worms[[i]] = df %>% distinct(worm_index) %>% pull(worm_index)
  my_data[[i]] = df
  # if (i == './nf4_features.hdf5') {df.nf4 = as.data.frame(df)}
  
  df2 = h5read(file = i, name = 'features_timeseries')
  df2 = as.data.frame(df2)
  df2x = df2
  df2 = df2 %>% select(1:7,40,57)
  df2$file.ID = i
  my_data2[[i]] = df2
}

h5closeAll()
agg = data.frame(rbindlist(my_data)) %>%
  mutate(file = stringr::str_extract(file.ID, "(?<=/)[^_]+")) %>% 
  mutate(group = case_when(grepl('0', file) ~ '0%',
                          grepl('2', file) ~ '2%',
                          grepl('3', file) ~ '3%',
                          grepl('4', file) ~ '4%')) %>%
  mutate(worm_ID = paste0(worm_index,".",file)) %>% dplyr::filter(as.numeric(n_valid_skel) > 50)
# agg$omega_turns_frequency[is.nan(agg$omega_turns_frequency)]=0
# agg$upsilon_turns_frequency[is.nan(agg$upsilon_turns_frequency)]=0
# agg = agg %>% mutate(omega = omega_turns_frequency * n_frames, upsilon = upsilon_turns_frequency * n_frames)
# avg.agg = agg %>% group_by(group, worm_ID) %>%
#   summarise(meanO = mean(omega), meanOF = mean(omega_turns_frequency), 
#             meanU = mean(upsilon), meanUF = mean(upsilon_turns_frequency), meanD = mean(worm_dwelling),
#             meanPC = mean(path_curvature), mean.absPC = mean(path_curvature)) %>% as.data.frame()
ggplot(agg, aes(x = n_valid_skel, fill = group, after_stat(count / max(count)), color = group)) + geom_density() + 
  geom_vline(xintercept = c(100, 300, 900)) + facet_wrap(~group) + ylab("Frequency") + xlab("Frame Length of Valid Skeletons")
ggplot(agg, aes(x = n_valid_skel, group = group)) + geom_histogram(binwidth = 90) + facet_wrap(~ group)

# 
# 
# # ggplot(avg.agg, aes(x = group, y = meanOF)) + geom_boxplot()
# {
# pof = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanOF, 
#                            pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE, ylab = "Omega Turn Frequency (1/s)", xlab = '')
# po = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanO, 
#                                   pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE, ylab = "Omega Turns", xlab = '')
# puf = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanUF, 
#                                   pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE, ylab = "Upsilon Turn Frequency (1/s)", xlab = '')
# pu = ggstatsplot::ggbetweenstats(data = avg.agg, x = group, y = meanU, 
#                                   pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE, ylab = "Upsilon Turns", xlab = '') }


agg2 = data.frame(rbindlist(my_data2)) %>%
  mutate(file = stringr::str_extract(file.ID, "(?<=/)[^_]+")) %>% mutate(worm_ID = paste0(worm_index,".",file)) %>%
  mutate(group = case_when(grepl('0', file) ~ '0%',
                           grepl('2', file) ~ '2%',
                           grepl('3', file) ~ '3%',
                           grepl('4', file) ~ '4%')) %>%
  mutate(skeleton_fit = case_when(skeleton_id == -1 ~ NaN, TRUE ~ skeleton_id)) %>%
  group_by(worm_ID) %>% arrange(group, file, worm_ID, timestamp) %>% ungroup() %>%
  dplyr::filter(worm_ID %in% agg$worm_ID) %>% dplyr::filter(midbody_width < 4.5) %>%
  # mutate(difference_f = timestamp - lead(timestamp),
  #   difference_mm = motion_modes - lead(motion_modes)) %>%
  # mutate(reversal = case_when(difference_f == 1 | 2 | 3 | 4 | 5 & difference_mm == 0 ~ 0,
  #                             difference_f == 1 | 2 | 3 | 4 | 5 & difference_mm == 1 | 2 ~ 1,
  #                             TRUE ~ NaN)) %>%
  # mutate(RevNum = sum(reversal)) %>% 
  arrange(group, file, worm_ID, timestamp) %>% na.omit() %>% 
  select(2,4,7,8,12,13) %>% distinct() %>% mutate(time = timestamp / 540) %>% arrange(group, worm_ID, time) %>% ungroup()

agg2mean = agg2 %>% group_by(group, timestamp) %>%
  summarise(width = mean(midbody_width), speed.abs = mean(abs(midbody_speed)), speed = mean(midbody_speed), time = time, n = n()) %>% distinct()
agg2mean_speed = agg2mean %>% select(time, group, speed.abs)

# agg.upd = agg %>% filter(worm_ID %in% agg2)
ggplot(agg2, aes(x = midbody_width, group = group, after_stat(count / max(count)), color = group, fill = group)) + geom_density() + facet_wrap(~ group) + 
  geom_vline(xintercept = 4.5) + xlab("Worm midbody width (px)") + ylab("Frequency")

ggplot(agg2, aes(x = time, y = abs(midbody_speed), group = group, color = group)) + geom_point(aes(group = group), size = 0.3, alpha = 0.3) + 
  geom_smooth(se = F) +
  facet_wrap(~ group) + #+ geom_line(data = agg2mean_speed, aes(x = timestamp, y = speed.abs, color = group), size = 0.5, alpha = 1) + geom_smooth(se = F)
  xlab("Time (min)") + ylab("|Worm midbody speed| (px / s)")


ggplot(agg2mean, aes(x = timestamp, group = group, color = group, y = speed.abs)) + geom_line() + geom_smooth() + facet_wrap(~ group)
ggplot(agg2mean, aes(x = width, group = group, after_stat(count / max(count)), color = group, fill = group)) + geom_density() + facet_wrap(~ group) + 
  geom_vline(xintercept = 4.5) + xlab("Worm midbody width (px)") + ylab("Frequency")
ggplot(agg2mean, aes(x = time, y = n, group = group, color = group)) + geom_point() + geom_line() + facet_wrap(~ group) + 
  xlab("Time (min)") + ylab("Number of unique worms per frame")
# write.csv(agg2, 'agg2.csv')

agg2sum = agg2 %>% mutate(time = timestamp / 540) %>% group_by(group, timestamp) %>%
  summarise(speed.abs = mean(abs(midbody_speed)), speed = mean(midbody_speed), time = time, n = n()) %>% distinct()

agg.agg = left_join(agg, agg2, by = c("worm_ID", "group")) %>% ungroup() #%>% filter(!is.na(RevNum)) %>% mutate(RevFreq = 10 * RevNum / n_frames)

ggplot(agg2sum, aes(x = time, y = n, group = group, color = group)) + geom_point() + geom_line() + facet_wrap(~ group) + 
  xlab("Time (min)") + ylab("Number of unique worms per frame")
ggplot(agg2sum, aes(x = time, y = speed.abs, group = group, color = group)) + geom_point() + facet_wrap(~ group) + 
  xlab("Time (min)") + ylab("|Speed| (px / s)")



  ggplot(agg2, aes(x = timestamp, y = RevNum)) + geom_line() + facet_wrap(~ group)
pr = ggstatsplot::ggbetweenstats(data = agg2, x = group, y = RevNum, 
                                 pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE,
                                 xlab = "Group", ylab = "Reversals")

agg.agg = left_join(agg, agg2, by = c("worm_ID", "group")) %>% ungroup() %>% filter(!is.na(RevNum)) %>% mutate(RevFreq = 10 * RevNum / n_frames)
prf = ggstatsplot::ggbetweenstats(data = agg.agg, x = group, y = RevFreq, 
                                  pairwise.display = "all", p.adjust.method = "fdr", bf.message = FALSE,
                                  xlab = "Group", ylab = "Reversal Frequency (1/s)")

pt = ggplot() + labs(title = "Movement Statistics in WT Worms ± Food")
plots = plot_grid(pof, po, puf, pu, prf, pr, labels = "AUTO", align = 'hv', ncol = 2)
full_plot = plot_grid(pt, plots, ncol = 1, rel_heights = c(0.04, 1))
full_plot
