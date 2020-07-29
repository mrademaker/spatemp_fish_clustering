library(dplyr)
library(ggpubr)
#### Exploring cluster properties ####
df = read.csv('/Users/markrademaker/Projects/Spatio_temporal_clustering/Data/cod_cluster_gdf.csv')

####################################
# BASED ON FEATURE RANGE (MIN MAX) #
####################################

#Total feature range
df_max_time = max(df$YearQuarter)
df_min_time = min(df$YearQuarter)
df_range_time = abs(df_max_time - df_min_time)
df_mean_time = mean(df$YearQuarter)

df_max_lat = max(df$ShootLat)
df_min_lat = min(df$ShootLat)
df_range_lat = abs(df_max_lat - df_min_lat)
df_mean_lat = mean(df$ShootLat)

df_max_lon = max(df$ShootLong)
df_min_lon = min(df$ShootLong)
df_range_lon = abs(df_max_lon - df_min_lon)
df_mean_lon = mean(df$ShootLong)

df_max_jbiom = max(df$Gadus_morhua_juvenile)
df_min_jbiom = min(df$Gadus_morhua_juvenile)
df_range_jbiom = abs(df_max_jbiom - df_min_jbiom)
df_mean_jbiom = mean(df$Gadus_morhua_juvenile)

df_max_abiom = max(df$Gadus_morhua_adult)
df_min_abiom = min(df$Gadus_morhua_adult)
df_range_abiom = abs(df_max_abiom - df_min_abiom)
df_mean_abiom = mean(df$Gadus_morhua_adult)

#Initiate empty list
list_obs=list()

#Feature range per cluster
for (group in seq(1:32)){
  print(group)
  dfx = df %>% filter(clusters == group)
  
  min_time = min(dfx$YearQuarter)
  max_time = max(dfx$YearQuarter)
  time_range = max_time - min_time
  time_prop  = time_range / df_range_time
     
  min_lat = min(dfx$ShootLat)
  max_lat = max(dfx$ShootLat)
  lat_range = max_lat - min_lat
  lat_prop = lat_range / df_range_lat
  
  min_lon = min(dfx$ShootLong)
  max_lon = max(dfx$ShootLong)
  lon_range = max_lon - min_lon
  lon_prop = lon_range / df_range_lon
  
  min_jbiom = min(dfx$Gadus_morhua_juvenile)
  max_jbiom = max(dfx$Gadus_morhua_juvenile)
  jbiom_range = max_jbiom - min_jbiom
  jbiom_prop  = jbiom_range / df_range_jbiom
  
  min_abiom = min(dfx$Gadus_morhua_adult)
  max_abiom = max(dfx$Gadus_morhua_adult)
  abiom_range = max_abiom - min_abiom
  abiom_prop = abiom_range / df_range_abiom
  
  obs = list(group,
             min_time,max_time,time_range,time_prop,
             min_lat,max_lat,lat_range,lat_prop,
             min_lon,max_lon,lon_range,lon_prop,
             min_jbiom,max_jbiom,jbiom_range,jbiom_prop,
             min_abiom,max_abiom,abiom_range,abiom_prop)
  #store as list of list
  list_obs[[group]] = obs
}
#List of lists to dataframe
df_sort=as.data.frame(do.call(rbind, lapply(list_obs, as.vector)))
colnames(df_sort)=c('cluster',
                    'min_yquarter','max_yquarter','yquarter_range','yquarter_prop',
                    'min_lat','max_lat','lat_range','lat_prop',
                    'min_lon','max_lon','lon_range','lon_prop',
                    'min_jbiom','max_jbiom','jbiom_range','jbiom_prop',
                    'min_abiom','max_abiom','abiom_range','abiom_prop'
                    )

### Boxplot of CLUSTER Properties ######
library(tidyverse)
library(ggplot2)
library(wesanderson)

#First subset columns of interest
df = df_sort[,c("cluster","yquarter_prop","lat_prop","lon_prop","jbiom_prop","abiom_prop")]
#Transform data into long format and set appropriate data type
data_long <- gather(df, variable, measurement, yquarter_prop:abiom_prop, factor_key=TRUE)
data_long$cluster = as.numeric(data_long$cluster)
data_long$measurement = as.numeric(data_long$measurement)

#Create boxplot
data_long %>% 
  ggplot(aes(x=variable,y=measurement, fill=variable)) +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 5))+
  #scale_fill_viridis_d()+
  theme_classic()+
  geom_boxplot(show.legend=FALSE) +
  geom_point(size=4,alpha=0.1,show.legend=FALSE)+
  theme(axis.ticks = element_blank())+
  scale_x_discrete(name = "\nCluster features",labels=c("yquarter_prop" = "Year Quarter",
                            "lat_prop" = "Latitude",
                            "lon_prop" = "Longitude",
                            "jbiom_prop" = "Juv. Biomass  ",
                            "abiom_prop" = "Adult Biomass"))+
  scale_y_continuous(name="Cluster feature range (max-min) / Total feature range (max-min)\n")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 18, angle = 90, hjust = .5, vjust = .5, face = "bold"))

  ggsave('/Users/markrademaker/Projects/Spatio_temporal_clustering/Boxplot cluster properties range (min max).png',dpi=300)

####################################
# HISTOGRAM BIOMASS ALL DATA       #
#################################### 
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ggpubr)

# randomly sampled dataset #
df = read.csv('/Users/markrademaker/Projects/Spatio_temporal_clustering/Data/cod_cluster_gdf.csv')
# full dataset #
#df = read.csv('/Users/markrademaker/Projects/Spatio_temporal_clustering/Data/Gadus_morhua_1991_2019_juv_adult_biomass.csv')

#First subset columns of interest
df = df[,c("clusters","YearQuarter","ShootLat","ShootLong","Gadus_morhua_juvenile","Gadus_morhua_adult")]
#if full dataset
#df = df[,c("YearQuarter","ShootLat","ShootLong","Gadus_morhua_juvenile","Gadus_morhua_adult")]

colnames(df)=c('cluster',"Year_Quarter","Latitude","Longitude","Juv_Biomass","Adult_Biomass")
#if full dataset
#colnames(df)=c("Year_Quarter","Latitude","Longitude","Juv_Biomass","Adult_Biomass")

#Create histogram adult biomass
pa = ggplot(df, aes(x=log10((Adult_Biomass+0.01)))) +
  theme_classic()+
  geom_histogram(colour="black",fill="#FF6666",alpha=0.5,bins=32)+
  theme(axis.ticks = element_blank())+
  scale_y_continuous(breaks=c(0,50,100,150,200,250,300),expand=c(0,0),name="Count")+
  scale_x_continuous(name="Log10(Adult Biomass + 0.01) (kg)")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 18, angle = 90, hjust = .5, vjust = .5, face = "bold"))
ggsave('/Users/markrademaker/Projects/Spatio_temporal_clustering/histogram properties full data adult biomass.png',dpi=300)
#Create histogram juvenile biomass
pj= ggplot(df, aes(x=log10(Juv_Biomass+0.01))) +
  theme_classic()+
  geom_histogram(colour="black",fill="#FF6666",alpha=0.5,bins=32)+
  theme(axis.ticks = element_blank())+
  scale_y_continuous(breaks=c(0,50,100,150,200,250,300),expand=c(0,0),name="Count")+
  scale_x_continuous(name="Log10(Juvenile Biomass + 0.01) (kg)")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 18, angle = 90, hjust = .5, vjust = .5, face = "bold"))
ggsave('/Users/markrademaker/Projects/Spatio_temporal_clustering/histogram properties full data juvenile biomass.png',dpi=300)

#Side by side plot
ggarrange(pa,pj)
ggsave('/Users/markrademaker/Projects/Spatio_temporal_clustering/histogram properties full data.png',dpi=300)

### Now density plot instead ####
library(ggplot2)
library(scales)

#Adult biomass
#Create histogram adult biomass
pa = ggplot(df, aes(x=Adult_Biomass)) +
  theme_classic()+
  geom_histogram(aes(y=..count..),color="black",fill="#FF6666",alpha=0.2,bins=32)+
  #geom_density(aes(y=..density..),color="darkblue",fill="lightblue",alpha=0.2)+
  scale_x_continuous(breaks=c(0,5,10,25,50,100,250,500,1000,2000,4000), trans="log1p", expand=c(0,0),name="Adult Biomass (kg)") +
  #scale_y_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),expand=c(0,0),name="Count\n") +
  scale_y_continuous(breaks=c(0,50,100,150,200,250,300),expand=c(0,0),name="Count")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 18, angle = 90, hjust = .5, vjust = .5, face = "bold"))
pa

#Create histogram juvenile biomass
pj = ggplot(df, aes(x=Juv_Biomass)) +
  theme_classic()+
  geom_histogram(aes(y=..count..),color="black",fill="#FF6666",alpha=0.2,bins=32)+
  scale_x_continuous(breaks=c(0,5,10,25,50,100,250,500,1000,2000,4000), trans="log1p", expand=c(0,0),name="Juvenile Biomass (kg)") +
  #scale_y_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),expand=c(0,0),name="Count\n") +
  scale_y_continuous(breaks=c(0,50,150,200,250,300),expand=c(0,0),name="Count\n")+
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 18, angle = 90, hjust = .5, vjust = .5, face = "bold"))
pj
#Side by side plot
ggarrange(pa,pj)
ggsave('/Users/markrademaker/Projects/Spatio_temporal_clustering/histogram properties full data kgs.png',dpi=300)
