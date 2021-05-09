library(dplyr)
library(ggplot2)
library(data.table)
library(rhdf5)
library(cowplot)
library(ggstatsplot)
library(tidyquant)

#set wd
#setwd("C:/Users/shona/OneDrive/Desktop/n2 vs 269/n2 vs apw269/final")
setwd("D:/APW Lab/Food ON OFF/trial 5/Together")
temp = intersect(list.files(pattern = "\\.hdf5$", full.names = T), list.files(pattern = "featuresN", full.names = T))

my_data = list()
meta = ''
df = ''
worms = list()

for (i in unique(temp)){
  meta <- h5read(file = i, name = "trajectories_data")
  meta <- as.data.frame(meta)
  df = h5read(file = i, name = "timeseries_data")
  df = as.data.frame(df)
  #add column with i
  
  meta = meta %>% filter(worm_label == 1)
  df = df %>% filter(worm_index %in% meta$worm_index_joined)
  df = df %>% mutate_at(vars(worm_index), funs(meta$worm_index_manual[match(., meta$worm_index_joined)]))
  df = df %>% select(1:4)
  df$file.ID <- i
  worms[[i]] = df %>% distinct(worm_index) %>% pull(worm_index)
  my_data[[i]] = df
  
}

h5closeAll()
agg = data.frame(rbindlist(my_data)) %>% 
  mutate(file = stringr::str_extract(file.ID, "(?<=/)[^_]+")) %>% 
  mutate(food = case_when(grepl('f', file) ~ 'food',
                          grepl('s', file) ~ 'starve')) %>%
  mutate(geno = case_when(grepl('a', file) ~ 'AMPK KO',
                          grepl('n', file) ~ 'WT')) %>%
  mutate(group = case_when(food == 'food' & geno == 'WT' ~ 'WT Food',
                           food == 'starve' & geno == 'WT' ~ 'WT Starve',
                           food == 'food' & geno == 'AMPK KO' ~ 'AMPK KO Food',
                           food == 'starve' & geno == 'AMPK KO' ~ 'AMPK KO Starve'))
avg.agg = agg %>% filter(timestamp > 200) %>% group_by(group, worm_index) %>% na.omit() %>% 
  summarise(mean = mean(abs(speed))) %>% as.data.frame()
ggstatsplot::ggwithinstats(data = avg.agg, x = group, y = mean, 
                            pairwise.display = "all", p.adjust.method = "fdr")

# write.csv(avg.agg, 'avg.agg.csv')




##Filter for movement, add time instead of frame, correct for fps diff, label genos or trials
n22.data.ff = n22.data.f %>% select(1:4) %>% filter(!between(timestamp, 0, 48)) %>% mutate(time = timestamp / 48) %>% mutate(trial = "WT 1") %>% mutate(geno = 'WT')
n23.data.ff = n23.data.f %>% select(1:4) %>% filter(!between(timestamp, 0,1050)) %>% mutate(time = timestamp / 45) %>% mutate(trial = "WT 2")%>% mutate(geno = 'WT')
n24.data.ff = n24.data.f %>% select(1:4) %>% filter(!between(timestamp, 0,350)) %>% mutate(time = timestamp / 52) %>% mutate(trial = "WT 3")%>% mutate(geno = 'WT')

#Concatenate trials together for plotting
full.data2 = ''
full.data2 = rbindlist(list(n22.data.ff, n23.data.ff, n24.data.ff))
fd2 = full.data2 %>% na.omit() %>% group_by(trial, worm_index, time) %>% summarise(., avg = mean(abs(speed)))

ggplot(fd2, aes(x = time, y = avg, color = as.factor(trial))) + ylim(NA,10) + geom_line() + 
  geom_smooth(color='black') #+ facet_wrap(~ trial)

#Plot one trial
ggplot(n2.data.ff, aes(x = timestamp, y = speed, color = as.factor(worm_index), group = as.factor(worm_index))) + geom_line() 

#plot all trials, avg worms, color by trial
ggplot(full.data, aes(x = time, y = speed, group = trial, color = trial)) + 
  #geom_point(alpha = 0.1, na.rm = TRUE) + 
  geom_smooth(fullrange=FALSE, na.rm = TRUE) + facet_wrap(~ trial)

#Presentation figs, abs(speed) and separate by worm and by trial
p1 = ggplot(full.data, aes(time, y = abs(speed), group = trial, color = trial)) + geom_smooth(fullrange=FALSE, na.rm = TRUE) + 
  facet_wrap(~ as.factor(worm_index))
p2 = ggplot(full.data, aes(time, y = angular_velocity, group = trial, color = trial)) + 
  geom_smooth(fullrange=FALSE, na.rm = TRUE) + facet_wrap(~ as.factor(worm_index))
plot_grid(p1, p2, labels = "AUTO", align = 'hv')

# #Speed vs abs(speed)
# ggplot(full.data, aes(timestamp, abs(speed))) + geom_smooth(fullrange=FALSE, na.rm = TRUE)
# ggplot(full.data, aes(timestamp, y = angular_velocity, group = trial, color = trial)) +
#          geom_smooth(fullrange=FALSE, na.rm = TRUE)


#Compare angular velocity to speed per frame per worm to see correlation using ggstatsplot
angvel.speed.comp = ggstatsplot::ggscatterstats(
  data = full.data,
  x = speed,
  y = angular_velocity,
  type = 'parametric',
  bf.message = FALSE,
  results.subtitle = FALSE,
  marginal.type = 'density',
  title = 'Relationship between worm velocity and angular velocity across trials',
  subtitle = 'P(spearman) = 0.02, P(pearson) = 0.06'
  )
angvel.speed.comp

