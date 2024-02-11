#Author: Xiaoxuan Li
#12012022
library(corrplot)
library(tidyverse)
library(psych)
library(xlsx)
library(car)
library(MASS)
library(rsq)
library(ggpmisc)
library(ggrepel)
library(dplyr)
library(rlist)
library(data.table)
library(lemon)
library(lmtest)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(mltools)
library(tidyr)
library(stats)
library(lidR)
library(tools) 
library(raster)
library(ggpointdensity)
library(viridis)
library(grid)
library(readxl)
library(ehaGoF)
library(Metrics)
library(rgdal)
library(caret)
library(randomForest)
library(regclass)
windowsFonts(A = windowsFont("Times New Roman"))



# ------------------------------------------------------------------------------------------------ #
#scatterplot --- biomass pred vs chm 
dir <- "E:\\GEDI\\Result\\Result\\All_b_s_p_cc_m_biom.csv"

Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
#Data <- Data[Data$Combo !=0,]
#Data$Combo <- Data$Tree
Data <- Data[Data$status == "Leaf-on",]

Data <- Data[Data$ql4 ==1,]
#Data$Combo <- Data$Combo/0.049
#Data$Combo <- 9.8441*(Data$MEAN*Data$CC)
Data$Combo <- 6.85*Data$MEAN^0.95
Data <- Data[Data$Combo >0,]
reg <- lm(formula = sqrt(Combo) ~ 
            site+Beam+sensitivity+landsat_treecover+fhd_normal+
            sqrt(RH_98)+GEDI_cc_1+GEDI_cc_2, 
          data = Data)

Data$pred <- predict(reg, Data)^2

#stats
summary(reg)
r2_1 = round(summary(reg)$r.squared,3)
RMSE_1 <- sqrt(mean((Data$Combo - Data$pred)^2))
rRMSE_1 <- 100*sqrt(mean((Data$pred - Data$Combo)^2))/mean(Data$Combo)




p1<- ggplot(Data, aes(x=Combo, y=pred))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=15,y=80,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2_1 
           ),parse=TRUE) + 
  annotate("text",x=15,y=70,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE_1,3),"Mg/ha",             
             "\n" , " %RMSE: ", round(rRMSE_1,3),"%",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("Footprint level CHM derived AGB (Leaf-on) (Mg/ha)"), 
       y=expression("Multi-variate GEDI predicted AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




dir <- "E:\\GEDI\\Result\\Result\\All_b_s_p_cc_m_biom.csv"

Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
#Data <- Data[Data$Combo !=0,]
#Data$Combo <- Data$Tree
Data <- Data[Data$status == "Leaf-off",]

Data <- Data[Data$ql4 ==1,]
#Data$Combo <- Data$Combo/0.049
#Data$Combo <- 9.8441*(Data$MEAN*Data$CC)
Data$Combo <- 6.85*Data$MEAN^0.95
Data <- Data[Data$Combo >0,]
reg <- lm(formula = sqrt(Combo) ~ 
            site+Beam+sensitivity+landsat_treecover+fhd_normal+
            sqrt(RH_98)+GEDI_cc_1+GEDI_cc_2, 
          data = Data)

Data$pred <- (predict(reg, Data)^2)
#stats
summary(reg)
r2_1 = round(summary(reg)$r.squared,3)
RMSE_1 <- sqrt(mean((Data$Combo - Data$pred)^2))
rRMSE_1 <- 100*sqrt(mean((Data$pred - Data$Combo)^2))/mean(Data$Combo)

p2<- ggplot(Data, aes(x=Combo, y=pred))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=15,y=80,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2_1 
           ),parse=TRUE) + 
  annotate("text",x=15,y=70,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE_1,3),"Mg/ha",             
             "\n" , " %RMSE: ", round(rRMSE_1,3),"%",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("Footprint level CHM derived AGB (Leaf-off) (Mg/ha)"), 
       y=expression("Multi-variate GEDI predicted AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggarrange(p1,p2)
out = "E:\\Biomass\\CSIR\\figure\\agbd_chm_mean_v2.jpg"
#out = "E:\\Biomass\\CSIR\\figure\\agbd_H_CC.jpg"
#out = "E:\\Biomass\\CSIR\\figure\\agbd_individual.jpg"
ggsave(out,height=12, width=24, dpi=600)



# ------------------------------------------------------------------------------------------------ #
#assumptions
dir <- "E:\\GEDI\\Result\\Result\\All_b_s_p_cc_m_biom.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
#Data <- Data[Data$Combo !=0,]
#Data$Combo <- Data$Tree
Data <- Data[Data$status == "Leaf-on",]

Data <- Data[Data$ql4 ==1,]
#Data$Combo <- Data$Combo/0.049
#Data$Combo <- 9.8441*(Data$MEAN*Data$CC)
Data$Combo <- 6.85*Data$MEAN^0.95
Data <- Data[Data$Combo >0,]

reg <- lm(formula = sqrt(Combo) ~ 
            sqrt(RH_98)+GEDI_cc_1+GEDI_cc_2+fhd_normal,
          data = Data)
summary(reg)
#vif(reg)^2
Data$pred <- predict(reg, Data)^2


#1/26/2023
# ------------------------------------------------------------------------------------------------ #
#test Colgan 2012 CC*H
dir <- "E:\\GEDI\\Result\\Result\\All_b_s_p_cc05_m_biom.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]

Data$Combo <- 25.61*(Data$MEAN*Data$CC0) + 3.95
Data$Combo <- 9.8441*(Data$MEAN*Data$CC0)

Data <- Data[Data$Combo >0,]


ggplot(Data, aes(x=Combo, y=AGBD_1))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("Footprint level CHM derived AGB (Leaf-on) (Mg/ha)"), 
       y=expression("On-orbit GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


reg <- lm(formula = sqrt(Combo) ~ 
            sqrt(RH_98)+GEDI_cc_1+GEDI_cc_2+fhd_normal, 
          data = Data)
summary(reg)
#vif(reg)^2
Data$pred <- predict(reg, Data)^2

ggplot(Data, aes(x=Combo, y=pred))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("Footprint level CHM derived AGB (Leaf-off) (Mg/ha)"), 
       y=expression("Multi-variate GEDI predicted AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



