#Author: Xiaoxuan Li
#01282022
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
#assumptions
dir <- "E:\\GEDI\\Result\\Result\\All_b_s_p_cc_m_biom.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- Data[Data$status == "Leaf-on",]

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




#1/28/2023 
# ------------------------------------------------------------------------------------------------ #
#scatterplot csir 25m vs. ALS MCH biomass
filedir <- "E:\\Biomass\\CSIR\\result\\CSIR_25m\\biom_CSIR_25m_combo_1.csv"
Data = read.csv(filedir)
Data <- na.omit(Data)

Data$ALS <- Data$Ind
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
MD <- mean(Data$ALS - Data$CSIR)

p1<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="black",shape=19, size=3)+
  theme_bw()+
  annotate("text",x=0,y=143,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=120,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " Bias: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  coord_cartesian(xlim = c(0, 160),ylim =  c(0, 160))+
  labs(x=expression("25 m CSIR AGB (Mg/ha)"), 
       y=expression("Individual tree based AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


Data$ALS <- Data$MCH
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
MD <- mean(Data$ALS - Data$CSIR)

p2<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="black",shape=19, size=3)+
  theme_bw()+
  annotate("text",x=0,y=143,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=120,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " Bias: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  coord_cartesian(xlim = c(0, 160),ylim =  c(0, 160))+
  labs(x=expression("25 m CSIR AGB (Mg/ha)"), 
       y=expression("Area based AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

