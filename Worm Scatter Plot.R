install.packages("viridis")  # Install
install.packages("mclust")
install.packages("xlsx")
library("viridis")    
library(tidyverse)
library(janitor)
library(readxl)
library(reshape2)
library(ggplot2)
library(ggpointdensity)
library(mclust)
library(rhdf5)
library("xlsx")

Worm.2=NULL

### Extract Data from HDF5 
A = H5Fopen(file.choose())
B = A$'trajectories_data'
C = B[, c("timestamp_raw", "worm_index_joined","coord_x","coord_y")]
Worm= na.omit(C)
names(Worm)= c("timebin", "animal_number","x","y") #rename the columns to make it easier to code
change <- as.numeric(0.4)
adjust= Worm %>% mutate(animal_number=animal_number+change)

Worm.2<- rbind(adjust, Worm.2)


# tracks.df <- read_excel("data/GPR_Data_TEMPLATE.xlsx")
# # tracks.df<-tracks.df %>% drop_na()

#####sample rate of video
tim<-10  #seconds/frame
#########
subsamp<-10 #subsample rate seconds

tracks.df<-as.data.frame(Worm.2)
anim<- unique(tracks.df$animal_number)
time<-unique(tracks.df$timebin)
tracks22.df<-NULL
tracks3.df<-NULL
behav2.df<-NULL
dY<-NULL
dX<-NULL
p<-NULL
p_sum<-NULL
p_sum2<-NULL
ptot<-NULL
str<-NULL
sx<-NULL

tracks.df <- tracks.df %>%
  #####GREG Correction 
  mutate(x=x/1002)%>%
  mutate(y=y/1002)%>%
  ##########
mutate(y = 2.748-y)

tracks2.df<-tracks.df

  for( i in anim){
    tracks22.df<- tracks2.df %>% filter(animal_number == i)
    tracks22.df<- tracks22.df[seq(1, nrow(tracks22.df), subsamp), ]
    tracks22.df<-tracks22.df%>%
      mutate('lin vel'=sqrt((x-lag(x))^2+(y-lag(y))^2))
    
    ####Corrects each track to center of assay based on resolution of camera  
    yy<-head(tracks22.df$y,1)
    xx<-head(tracks22.df$x,1)
    ###1.92, 1.374 center point of resolution in cm 
    dX<-as.numeric(1.92-xx)
    dY<-as.numeric(1.374-yy)
    tracks22.df<-tracks22.df %>% 
      mutate(X= dX+x, Y=dY+y)
    
    
    tracks22.df<-tracks22.df%>%
      #converts to mm/s
      
      
      #####
    mutate(`lin vel`=(`lin vel`/tim)*10)%>%
      
      mutate(dy= y - lag(y))%>%
      mutate(dx= x - lag(x))%>%
      mutate(q1= abs(dy/dx))%>%
      mutate(q2= abs(dx/dy))%>%
      group_by(animal_number)%>%
      mutate(angle = case_when(
        dy==0 & dx>0 ~180,
        dy == 0 & dx < 0 ~ 0,
        dx == 0 & dy > 0 ~ 90,
        dx == 0 & dy < 0 ~ 270,
        dx < 0 & dy < 0 ~ (270 + (180 / pi) * atan(q2)),
        dx > 0 & dy < 0 ~ (180 + (180 / pi) * atan(q1)),
        dx > 0 & dy > 0 ~ (90 + (180 / pi) * atan(q2)),
        dx < 0 & dy > 0 ~ ((180 / pi) * atan(q1))
      )) %>%
      mutate(dthet= (angle - lead(angle))) %>% 
      mutate(theta= case_when
             (abs(dthet) >180 ~ abs(dthet)-360,
               abs(dthet) <180 ~ dthet)) %>% 
      mutate(ang= (abs(theta)))
    
    tracks3.df<-rbind(tracks3.df,tracks22.df)
  }

# tracks3.df<-tracks3.df %>% mutate(turns=case_when(ang<=30 ~"run",
#                                                   ang>30 & ang<150 ~"curve",
#                                                   ang>=150 & ang<=180~"turn"
# ))

####Turn Probabilities 
# for(q in time){
#   behav.df<-tracks3.df %>% filter(timebin == q)
#   for(i in anim){
#     behav2.df<-behav.df %>% filter(animal_number==i)
#     behav2.df<-drop_na(behav2.df)
#     p_sum<-count(behav2.df, turns)
#     ptot<-length(behav2.df$turns)
#     str<-head(behav2.df$strain,1)
#     sx<-head(behav2.df$sex,1)
#     p_sum<-p_sum %>%
#       mutate(prob= n/ptot,timestamp=q,sex=sx,strain=str) %>% 
#       select(strain,sex,timestamp,animal_number,turns,n,prob)
#     p_sum2<-rbind(p_sum2,p_sum)   
#     
#   }
# }
# turn_prob<-p_sum2 %>% as.data.frame(arrange(sex,turns,timestamp,animal_number))
# write.xlsx(turn_prob, file = "output/turn_probabilities.xlsx", sheetName = "turn_prob",col.names = TRUE, row.names=TRUE)


th<-as.data.frame(tracks3.df$theta)
av<-as.data.frame(tracks3.df$ang)
lv<-as.data.frame(tracks3.df$`lin vel`)
th<-th %>% drop_na()
av<-av %>% drop_na()
lv<-lv %>% drop_na()
tmin<-min(abs(th))  
tmax<-max(abs(th))
amax<-max(abs(av))
amin<-min(abs(av))
lvmax<-max(abs(lv))
lvmin<-min(abs(lv))

tracks_plot<-tracks3.df

### Scatter Plot
tracks4<- tracks3.df
clust.df<-tracks4 %>%ungroup() %>%  select(timebin,`lin vel`,ang)
angvel_speed<-clust.df%>%mutate(ang=abs(ang))
angvel_speed<-drop_na(angvel_speed)
ggplot(angvel_speed, aes(x=angvel_speed$ang,y=angvel_speed$`lin vel`))+geom_point()

# tracks_plot<-tracks_plot %>% 
#   mutate(seconds=case_when(timebin==0 ~seconds,
#                            timebin==10 ~seconds+600,
#                            timebin==20 ~seconds+1200
#   ))

# tracks_plot2<-filter(tracks_plot, sex =='h')

# tracks_plot<-filter(tracks_plot, animal_number ==1)
# tracks_plot<-tracks_plot %>% drop_na()

tracks_plot%>% ggplot(aes(x=X, y=Y))+
  geom_point(size=1.25, color="light blue")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  # theme_bw()+
  theme(axis.line = element_line(colour = "black"))

# tracks3.df%>% ggplot(aes(x=seconds, y=ang, color=sex))+
#   geom_point(size=1.25)

##### overlay parameters on track

#######plot linear velocity
tracks_plot<-filter(tracks3.df, sex =='h')
tracks_plot<-filter(tracks_plot, timebin ==0)



### Plot Tracks by animal
tracks_plot%>% ggplot(aes(x=X, y=Y,group= animal_number, na.rm = FALSE))+
  geom_path(size=1.25,aes(color = animal_number))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  # theme_bw()+
   theme(axis.line = element_line(colour = "black"))+
   # scale_color_gradient(low="blue", high="orange")+
  scale_y_continuous(limits = c(0,2.748))+
  scale_x_continuous(limits = c(0,3.84))


#####plot linear velocity
tracks_plot%>% ggplot(aes(x=X, y=Y, color= abs(`lin vel`),group= animal_number, na.rm = FALSE))+
  geom_path(size=1.25)+
  scale_color_gradient(low="blue", high="orange",limits=c(lvmin,lvmax))+
   scale_y_continuous(limits = c(0,2.748))+
   scale_x_continuous(limits = c(0,3.84))


###### plot turn angle over track
tracks_plot%>% ggplot(aes(x=X, y=Y, color= ang, group= animal_number,g,na.rm = TRUE))+
  geom_path(size=1.25)+
  scale_color_gradient(low="blue", high="orange",limits=c(tmin,tmax))+
  scale_y_continuous(limits = c(0,2.748))+
  scale_x_continuous(limits = c(0,3.84))


#######################Scatter Plot #######################################


tracks4<- tracks3.df
clust.df<-tracks4 %>%ungroup() %>%  select(timebin,`lin vel`,ang)
angvel_speed<-clust.df%>%mutate(ang=abs(ang))
angvel_speed<-drop_na(angvel_speed)
ggplot(angvel_speed, aes(x=angvel_speed$ang,y=angvel_speed$`lin vel`))+geom_point()
ggplot(angvel_speed, aes(x=angvel_speed$ang,y=angvel_speed$`lin vel`, color=as.factor(timebin), fill=as.factor(timebin)))+geom_point()+facet_wrap("timebin")


#######Clustering###################################
tracks4<-tracks3.df %>% filter(sex=='h')
clust.df<-tracks4 %>%ungroup() %>%  select(ang,`lin vel`)
angvel<-drop_na(clust.df)
mod1<-Mclust(angvel)


summary(mod1)

plot(mod1)
angvel$clust<- mod1$classification

angvel$uncertainty<- mod1$uncertainty

angvel <- angvel%>% arrange(clust,ang)


ggplot(angvel, aes(x=angvel$ang,y=angvel$`lin vel`, color=as.factor(clust), fill=as.factor(clust)))+geom_point()
