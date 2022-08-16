#packages
library(easyr)
library(ggplot2)
library(ggpubr)
setwd("C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/simulations")
a <- read.csv("C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/simulations/model_1_k_5.0_D_0.01_r_0.0_run_1_sim-duration_306.0.csv")
b <- read.csv("C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/simulations/model_1_k_50.0_D_0.01_r_0.0_run_1_sim-duration_306.0.csv")
c <- read.csv("C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/simulations/model_1_k_100.0_D_0.01_r_0.0_run_1_sim-duration_306.0.csv") 
#data wrangling
wrangle <- function(df){
  df$dead. <- tobool(df$dead.)
  df$cell_index <- as.factor(df$cell_index)
  #for some reason some cells born dead at timepoint 0 continue to have their ages update even though they are dead, we will remove them
  df <- df[!(df$time_point=="timepoint_1_of_306" & df$dead.==TRUE),]
  df$k <- as.factor(df$k)
  return(df)
}
a <- wrangle(a)
b <- wrangle(b)
c <- wrangle(c)
#plot_1 cell trajectory. 
#
a_dead <- a[a$dead.==TRUE,]

#plot_3 cell_cycle length
#for every cell index find earliest timepoint u 
#data_temp$
age_histogram <- function(data_temp){
ggplot(data_temp,(aes(x=age,col=k)))+
  geom_freqpoly(bins=100)+
  xlim(0.8,1.5)+
  annotate("label",x=Inf,y=Inf,label=paste("D=",as.character(data_temp$D[1])),vjust=1,hjust=1)+labs(title="Cell_Cycle_Length")
}
data_temp <- 
target_area_histogram <- function(data_temp){
ggplot(data_temp,(aes(x=A0,col=k)))+
  geom_freqpoly(bins=100)+
  xlim(0,3)+
  annotate("label",x=Inf,y=Inf,label=paste("D=",as.character(data_temp$D[1])),vjust=1,hjust=1)+labs(title="Target_Area")
}
apical_area_histogram <- function(data_temp){
  ggplot(data_temp,(aes(x=apical_area,color=k)))+
    geom_density(aes(y=..density..),size=1)+
    xlim(-0.05,1)+
    annotate("label",x=Inf,y=Inf,label=paste("D=",as.character(data_temp$D[1])),vjust=1,hjust=1)+labs(title="Apical_Area")
  
}
data_temp <- rbind(a,b,c)
data_temp_dead <- data_temp[data_temp$dead.==TRUE,]
data_temp_final_timepoint <- data_temp[data_temp$time_point=="timepoint_306_of_306",]
plot_c_1 <- target_area_histogram(data_temp_final_timepoint)
plot_c_2 <- age_histogram(data_temp_dead)
plot_c_3 <- apical_area_histogram(data_temp_final_timepoint)


#now will do d=0
a <- read.csv("C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/simulations/model_1_k_5.0_D_0.0_r_0.0_run_1_sim-duration_306.0.csv")
b <- read.csv("C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/simulations/model_1_k_50.0_D_0.0_r_0.0_run_1_sim-duration_306.0.csv")
c <- read.csv("C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/simulations/model_1_k_100.0_D_0.0_r_0.0_run_1_sim-duration_306.0.csv")

a <- wrangle(a)
b <- wrangle(b)
c <- wrangle(c)
data_temp <- rbind(a,b,c)
data_temp_dead <- data_temp[data_temp$dead.==TRUE,]
data_temp_final_timepoint <- data_temp[data_temp$time_point=="timepoint_306_of_306",]
plot_b_1 <- target_area_histogram(data_temp_final_timepoint)
plot_b_2 <- age_histogram(data_temp_dead)
plot_b_3 <- apical_area_histogram(data_temp_final_timepoint)

plot_b_1
plot_c_1
plot_b_2
plot_c_2
plot_b_3
plot_c_3

multi.page <- ggarrange(plot_b_1,plot_b_2,plot_b_3,nrow=1,ncol=1)
ggexport(multi.page,filename="C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/output/figure_8_panel_b.pdf")

multi.page <- ggarrange(plot_c_1,plot_c_2,plot_c_3,nrow=1,ncol=1)
ggexport(multi.page,filename="C:/Users/Aarya/OneDrive - Imperial College London/school/year3/Final_Year_Project/Git/vertex_model/examples/output/figure_8_panel_c.pdf")

#figure 9
a <- 
b
c
d
e
