library(dplyr)
library(ggplot2)
library(data.table)
library(rhdf5)
library(cowplot)
library(ggstatsplot)
library(tidyquant)

#set wd
#setwd("C:/Users/shona/OneDrive/Desktop/n2 vs 269/n2 vs apw269/final")
setwd("D:/APW Lab/Food ON OFF/trial 3/Results")
temp = intersect(list.files(pattern = "\\.hdf5$", full.names = T), list.files(pattern = "featuresN", full.names = T))

my_data = list()
meta = ''
df = ''


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
  my_data[[i]] = df
  
}

h5closeAll()


#TEST
n2.hd5 = "rb off 1_featuresN.hdf5"
n22.meta = h5read("rb off 1_featuresN.hdf5", "trajectories_data")
n22.data = h5read("rb off 1_featuresN.hdf5", "timeseries_data")
h5closeAll()
n22.meta.f = n22.meta %>% filter(worm_label == 1)
n22.data.f = n22.data %>% filter(worm_index %in% n22.meta.f$worm_index_joined)
n22.data.f %>% pull(worm_index) %>% unique()
n22.data %>% pull(worm_index) %>% unique()
my_data[[7]] %>% pull(worm_index) %>% unique()

# # 
# # #load hd5 files
# # #Here, we had 4 trials, I just deleted the first trial. Variables labeled n (for n2 WT worms) + trial#
# n2.hd5 = "n2 food 2_featuresN.hdf5"
# # n23.hd5 = "n2 3_featuresN.hdf5"
# # n24.hd5 = "n2 4_featuresN.hdf5"
# # 
# # #Extract metadata with manual curated worms and movement data for timeseries per worm
# n22.meta = h5read("n2 food 2_featuresN.hdf5", "trajectories_data")
# n22.data = h5read("n2 food 2_featuresN.hdf5", "timeseries_data")
# # n23.meta = h5read("n2 3_featuresN.hdf5", "trajectories_data")
# # n23.data = h5read("n2 3_featuresN.hdf5", "timeseries_data")
# # n24.meta = h5read("n2 4_featuresN.hdf5", "trajectories_data")
# # n24.data = h5read("n2 4_featuresN.hdf5", "timeseries_data")
# 
# #close all files just in case
# h5closeAll()
# 
# #Trim data to only manually curated worms
# n22.meta.f = n22.meta %>% filter(worm_label == 1)
# n22.data.f = n22.data %>% filter(worm_index %in% n22.meta.f$worm_index_joined)
# n23.meta.f = n23.meta %>% filter(worm_label == 1)
# n23.data.f = n23.data %>% filter(worm_index %in% n23.meta.f$worm_index_joined)
# n24.meta.f = n24.meta %>% filter(worm_label == 1)
# n24.data.f = n24.data %>% filter(worm_index %in% n24.meta.f$worm_index_joined)
# 
# #Limit actual timestamp data to only those with manual curation
# n22.data.f = n22.data.f %>% mutate_at(vars(worm_index), funs(n22.meta.f$worm_index_manual[match(., n22.meta.f$worm_index_joined)]))
# n23.data.f = n23.data.f %>% mutate_at(vars(worm_index), funs(n23.meta.f$worm_index_manual[match(., n23.meta.f$worm_index_joined)]))
# n24.data.f = n24.data.f %>% mutate_at(vars(worm_index), funs(n24.meta.f$worm_index_manual[match(., n24.meta.f$worm_index_joined)]))
# 
# #Pull individual worms
# n22.data.f %>% pull(worm_index) %>% unique()
# n23.data.f %>% pull(worm_index) %>% unique()
# n24.data.f %>% pull(worm_index) %>% unique()

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

