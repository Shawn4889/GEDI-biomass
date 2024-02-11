#Author: Xiaoxuan Li
#02082023
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
library(DAAG)
library(matrixStats)
library(FactoMineR)
library(ggcorrplot)
library(devtools)
library(factoextra)
library(purrr)
library(ggridges)
library(randomForest)
library(hexbin)
library(colorRamp2)
windowsFonts(A = windowsFont("Times New Roman"))


# ------------------------------------------------------------------------------------------------ #
#test multicollinearity
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$RH1_98 >=2.35,]
#Data <- Data[Data$solar_elevation <0,]
#Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 1.386052*Data$MEAN^1.996775

Data <- Data[Data$Combo >0,]


Data$RH1_50 <- Data$RH1_50 +100
reg <- glm(formula = Combo ~ RH1_50 + GEDI_PAVD + GEDI_PAI + GEDI_FHD + GEDI_Cover + RH1_98,
          data = Data)
vif(reg)
p = Data[, c("GEDI_PAI", "GEDI_FHD", "GEDI_Cover", "RH1_50", "RH1_75", "RH1_98")]
cor(p)


reg <- glm(formula = Combo ~ GEDI_PAVD+RH1_50+RH1_98,
           data = Data)
vif(reg)

p = Data[, c("Combo","RH1_50", "RH1_98","GEDI_PAVD")]
cor(p)
write.table(p,file='C:\\Users\\Shawn\\Desktop\\test.csv',sep = ",",col.names=TRUE)


Data$preds <- predict(reg,newdata=Data)
cor(Data$preds, Data$Combo)^2


#test importance
spit = c(train = .8, test = .2)
gs = sample(cut(seq(nrow(Data)), nrow(Data)*cumsum(c(0,spit)),labels = names(spit)))
res = split(Data, gs)
sapply(res, nrow)/nrow(Data)

rf <- randomForest(Combo ~ 
                     GEDI_Cover + RH1_98, 
                   data = res$train,ntree = 100,importance=TRUE, nodesize = 20)
importance(rf)
varImpPlot(rf)
varImp(rf)


#test correlation
reg <- lm(formula = Combo ~ GEDI_Cover+RH1_50+RH1_98,data = Data)
summary(reg)
vif(reg)
# ------------------------------------------------------------------------------------------------ #
#predicted GEDI L4A vs. updated ALS CHM biomass
#improved GEDI L4A vs. ALS CHM -- individual tree method 
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data$Combo <- Data$Tree
#Data <- Data[Data$solar_elevation <0,]
#Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data$Combo <- Data$Combo/0.049
Data <- Data[Data$Combo > 0,]
Data <- Data[Data$RHS_98 >=2.35,]

#RH98 only
train_control <- trainControl(method = "cv", number = 10, savePredictions=TRUE)

model <- train(sqrt(Combo) ~
                 sqrt(RH1_98), 
               data = Data, 
               method = "lm", 
               trControl = train_control)

print(model)
pred <- model$pred

R2s_1 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
              function(x) summary(lm(x))$r.sq))

RMSEs_1 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
                function(x) sqrt(mean(lm(x)$residuals^2))))

MDs_1 <- aggregate(x= (pred$pred)^2- (pred$obs)^2,     
                   by = list(pred$Resample),      
                   FUN = mean)


mean(R2s_1)
mean(RMSEs_1)
mean(MDs_1[,2])

#more variables were added
train_control <- trainControl(method = "cv", number = 10, savePredictions=TRUE)
model <- train(sqrt(Combo) ~GEDI_Cover+sqrt(RH1_98), 
               data = Data, method = "lm", trControl = train_control)

pred <- model$pred

R2s_2 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
              function(x) summary(lm(x))$r.sq))
RMSEs_2 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
                function(x) sqrt(mean(lm(x)$residuals^2))))
MDs_2 <- aggregate(x= (pred$pred)^2- (pred$obs)^2,     
                   by = list(pred$Resample),      
                   FUN = mean)

mean(R2s_2)
mean(RMSEs_2)
mean(MDs_2[,2])

r2 = mean(R2s_2)
RMSE <- mean(RMSEs_2)
MD <- mean(MDs_2[,2])

p1 <- ggplot(pred, aes(x=obs^2, y=pred^2))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=15,y=90,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(r2,3) 
           ),parse=TRUE) + 
  annotate("text",x=15,y=80,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Individual tree based AGB (Mg/ha)"), 
       y=expression("Predicted GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pred <- pred[(pred$obs)^2 < 80,]
pred$biases = (pred$pred)^2 - (pred$obs)^2
breakbin = round(seq(0,80,10),2)
pred$group <- cut((pred$obs)^2,breaks = breakbin,dig.lab=5)
pred %>% count(group)
aggregate(pred$biases, list(pred$group), mean)
aggregate(pred$biases, list(pred$group), mean)[,2]/
  aggregate((pred$obs)^2, list(pred$group), mean)[,2]*100

b1 <- ggplot(pred, aes(x=group, y=biases,group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 80))+
  scale_y_continuous(minor_breaks = round(seq(-80,80,20),digits = 1),
                     breaks = round(seq(-80,80,20),digits = 1))+
  theme_bw()+
  labs(x = "25 m Individual tree based AGB (Mg/ha)",
       y="Bias of predicted GEDI L4A AGB (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))



#GEDI L4A vs. ALS CHM -- area based method  H only
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
#Data <- Data[Data$RHS_98 >=2.35,]
#Data <- Data[Data$solar_elevation <0,]
#Data <- Data[Data$status == "Leaf-off",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 5.2*Data$MEAN^1.06
Data$Combo <- Data$Combo*490.625/625
Data <- Data[Data$Combo >0 & Data$Combo < 60,]
#RH98 only
train_control <- trainControl(method = "cv", number = 10, savePredictions=TRUE)

model <- train(sqrt(Combo) ~
                 sqrt(RH1_98), 
               data = Data, 
               method = "lm", 
               trControl = train_control)

pred <- model$pred

R2s_3 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
              function(x) summary(lm(x))$r.sq))

RMSEs_3 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
                function(x) sqrt(mean(lm(x)$residuals^2))))
MDs_3 <- aggregate(x= (pred$pred)^2- (pred$obs)^2,     
                   by = list(pred$Resample),      
                   FUN = mean)

mean(R2s_3)
mean(RMSEs_3)
mean(MDs_3[,2])

#more variables were added

train_control <- trainControl(method = "cv", number = 10, savePredictions=TRUE)
model <- train(sqrt(Combo) ~GEDI_Cover+sqrt(RH1_98), 
               data = Data, method = "lm", trControl = train_control)

pred <- model$pred

R2s_4 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
              function(x) summary(lm(x))$r.sq))
RMSEs_4 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
                function(x) sqrt(mean(lm(x)$residuals^2))))

MDs_4 <- aggregate(x= (pred$pred)^2- (pred$obs)^2,     
                   by = list(pred$Resample),      
                   FUN = mean)

mean(R2s_4)
mean(RMSEs_4)
mean(MDs_4[,2])

r2 = mean(R2s_4)
RMSE <- mean(RMSEs_4)
MD <- mean(MDs_4[,2])

p2 <- ggplot(pred, aes(x=obs^2, y=pred^2))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=15,y=90,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(r2,3)
           ),parse=TRUE) + 
  annotate("text",x=15,y=80,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Area-based (H only) AGB (Mg/ha)"), 
       y=expression("Predcited GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


pred <- pred[(pred$obs)^2 < 40,]
pred$biases = (pred$pred)^2 - (pred$obs)^2
breakbin = round(seq(0,40,10),2)
pred$group <- cut((pred$obs)^2,breaks = breakbin,dig.lab=5)
pred %>% count(group)/nrow(Data)*100
pred %>% count(group)
aggregate(pred$biases, list(pred$group), mean)
aggregate(pred$biases, list(pred$group), mean)[,2]/
  aggregate((pred$obs)^2, list(pred$group), mean)[,2]*100

b2 <- ggplot(pred, aes(x=group, y=biases,group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 80))+
  scale_y_continuous(minor_breaks = round(seq(-80,80,20),digits = 1),
                     breaks = round(seq(-80,80,20),digits = 1))+
  theme_bw()+
  labs(x = "25 m Area-based (H only) AGB (Mg/ha)",
       y="Bias of predicted GEDI L4A AGB (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))



#GEDI L4A vs. ALS CHM -- area based method  H*CC
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$RHS_98 >=2.35,]
#Data <- Data[Data$solar_elevation <0,]
#Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 9.623*Data$MEAN*Data$CC15 + 3.76
Data <- Data[Data$Combo >0,]
Data$Combo <- Data$Combo*490.625/625
#RH98 only
train_control <- trainControl(method = "cv", number = 10, savePredictions=TRUE)

model <- train(sqrt(Combo) ~
                 sqrt(RH1_98), 
               data = Data, 
               method = "lm", 
               trControl = train_control)

pred <- model$pred

R2s_3 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
              function(x) summary(lm(x))$r.sq))

RMSEs_3 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
                function(x) sqrt(mean(lm(x)$residuals^2))))
MDs_3 <- aggregate(x= (pred$pred)^2- (pred$obs)^2,     
                   by = list(pred$Resample),      
                   FUN = mean)

mean(R2s_3)
mean(RMSEs_3)
mean(MDs_3[,2])
#more variables were added
train_control <- trainControl(method = "cv", number = 10, savePredictions=TRUE)

model <- train(sqrt(Combo) ~GEDI_Cover+sqrt(RH1_98), 
               data = Data, method = "lm", trControl = train_control)

pred <- model$pred

R2s_4 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
              function(x) summary(lm(x))$r.sq))
RMSEs_4 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
                function(x) sqrt(mean(lm(x)$residuals^2))))

MDs_4 <- aggregate(x= (pred$pred)^2- (pred$obs)^2,     
                   by = list(pred$Resample),      
                   FUN = mean)

mean(R2s_4)
mean(RMSEs_4)
mean(MDs_4[,2])

r2 = mean(R2s_4)
RMSE <- mean(RMSEs_4)
MD <- mean(MDs_4[,2])

p3 <- ggplot(pred, aes(x=obs^2, y=pred^2))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=15,y=90,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(r2,3)
           ),parse=TRUE) + 
  annotate("text",x=15,y=80,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Area-based (H*CC) AGB (Mg/ha)"), 
       y=expression("Predicted GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


pred <- pred[(pred$obs)^2 < 80,]
pred$biases = (pred$pred)^2 - (pred$obs)^2
breakbin = round(seq(0,80,10),2)
pred$group <- cut((pred$obs)^2,breaks = breakbin,dig.lab=5)
pred %>% count(group)
aggregate(pred$biases, list(pred$group), mean)
aggregate(pred$biases, list(pred$group), mean)[,2]/
  aggregate((pred$obs)^2, list(pred$group), mean)[,2]*100

b3 <- ggplot(pred, aes(x=group, y=biases,group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 80))+
  scale_y_continuous(minor_breaks = round(seq(-80,80,20),digits = 1),
                     breaks = round(seq(-80,80,20),digits = 1))+
  theme_bw()+
  labs(x = "25 m Area-based (H*CC) AGB (Mg/ha)",
       y="Bias of predicted GEDI L4A AGB (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))



ggarrange(p1,p2,p3,ncol=3)
out = "E:\\Biomass\\CSIR\\figure\\L4A_CHM_biomass_improved.jpg"
ggsave(out,height=12, width=36, dpi=600)



ggarrange(b1,b2,b3,ncol=3)
out = "E:\\Biomass\\CSIR\\figure\\L4A_CHM_biomass_improved_biasplot.jpg"
ggsave(out,height=12, width=36, dpi=600)



ggarrange(p2,b2,ncol=2)
out = "E:\\Biomass\\CSIR\\figure\\L4A_CHM_biomass_biasplot_H_improved.jpg"
ggsave(out,height=12, width=24, dpi=600)


# ------------------------------------------------------------------------------------------------ #

#test RH98 2.35m to Biomass limitation
#test phenology condition and bias boxplot
#subsample

#GEDI L4A vs. ALS CHM -- area based method  H only
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
#Data <- Data[Data$RHS_98 <=2.35,]
#Data <- Data[Data$solar_elevation <0,]
#Data <- Data[Data$status == "Leaf-off",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 1.675*Data$MEAN^1.89
#Data$Combo <-  6.85*Data$MEAN^0.95
Data$Combo <- Data$Combo*490.625/625
Data <- Data[Data$Combo >0 & Data$Combo < 80,]


#subsample
Data <- Data %>% group_by(gr=cut(Combo, breaks= seq(0, 80, by=1))) %>% 
  arrange(as.numeric(gr))
Data <- Data %>% group_by(gr) %>% slice_sample(n=50)


#more variables were added
train_control <- trainControl(method = "cv", number = 10, savePredictions=TRUE)
model <- train(sqrt(Combo) ~GEDI_Cover+sqrt(RHS_98), 
               data = Data, method = "lm", trControl = train_control)
pred <- model$pred

#apply minimum biomass threshold
#pred <- pred[(pred$obs)^2 >= 0.566,]
#pred <- pred[(pred$obs)^2 >= 1.535,]

R2s_2 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
              function(x) summary(lm(x))$r.sq))
RMSEs_2 <- c(by(data.frame((pred$obs)^2, (pred$pred)^2), pred$Resample, 
                function(x) sqrt(mean(lm(x)$residuals^2))))
MDs_2 <- aggregate(x= (pred$pred)^2- (pred$obs)^2,     
                   by = list(pred$Resample),      
                   FUN = mean)

mean(R2s_2)
mean(RMSEs_2)
mean(MDs_2[,2])

r2 = mean(R2s_2)
RMSE <- mean(RMSEs_2)
MD <- mean(MDs_2[,2])


#check minimum biomass threshold
merge <- cbind((model$pred$pred)^2,(model$pred$obs)^2,model$trainingData)
merge <- merge[merge$RHS_98 <=2.35,]
summary(merge$`(model$pred$obs)^2`)



#scatterplot
ggplot(pred, aes(x=obs^2, y=pred^2))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=15,y=90,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(r2,3)
           ),parse=TRUE) + 
  annotate("text",x=15,y=80,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(pred))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Area-based (H only) AGB (Mg/ha)"), 
       y=expression("Predcited GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#bias boxplot

pred <- pred[(pred$obs)^2 < 80,]
pred$biases = (pred$pred)^2 - (pred$obs)^2
breakbin = round(seq(0,80,10),2)
pred$group <- cut((pred$obs)^2,breaks = breakbin,dig.lab=5)
pred %>% count(group)
aggregate(pred$biases, list(pred$group), mean)
aggregate(pred$biases, list(pred$group), mean)[,2]/
  aggregate((pred$obs)^2, list(pred$group), mean)[,2]*100

ggplot(pred, aes(x=group, y=biases,group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 80))+
  scale_y_continuous(minor_breaks = round(seq(-80,80,20),digits = 1),
                     breaks = round(seq(-80,80,20),digits = 1))+
  theme_bw()+
  labs(x = "25 m Area-based (H only) AGB (Mg/ha)",
       y="Bias of predicted GEDI L4A AGB (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))




# ------------------------------------------------------------------------------------------------ #
#scatterplot 1ha csir vs. area based and ind biomass
filedir <- "E:\\Biomass\\CSIR\\result\\CSIR\\biom_CSIR_combo.csv"
Data = read.csv(filedir)
Data <- Data[Data$CSIR >0,]
Data<-Data[!(Data$Name1=="Ven_S68_SP_5_Venetia_h.tif"),]
Data$ALS = Data$Tree
mean(Data$CSIR)
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
Mean <- mean(Data$CSIR)
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
rRMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))/mean(Data$CSIR)*100
MD <- mean(Data$ALS - Data$CSIR)
p1<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,
              color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+ 
  theme_bw()+
  geom_smooth(method="lm",se = F) +
  annotate("text",x=0,y=152,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=130,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " MSD: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  annotate("text",x=125,y=155,hjust = 0,size = 15,family= "A", label= "(a)") + 
  coord_cartesian(xlim = c(0, 160),ylim =  c(0, 160))+
  labs(x=expression("1 ha AGBD"["field_1ha"]~"(Mg/ha)"), 
       y=expression("1 ha AGBD"["ALS_tree"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


filedir <- "E:\\Biomass\\CSIR\\result\\CSIR\\biom_CSIR_combo.csv"
Data = read.csv(filedir)
Data <- Data[Data$CSIR >0,]
Data<-Data[!(Data$Name1=="Ven_S68_SP_5_Venetia_h.tif"),]
reg <- nls(CSIR~ a*H^b, start = list(a=1, b=1),data=Data)
Data$preds <- predict(reg,newdata=Data)
Data$ALS = Data$preds
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
Mean <- mean(Data$CSIR)       
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
rRMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))/mean(Data$CSIR)*100
MD <- mean(Data$ALS - Data$CSIR)

p2<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,
              color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  geom_smooth(method="lm",se = F) +
  annotate("text",x=0,y=152,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=130,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " MSD: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  annotate("text",x=125,y=155,hjust = 0,size = 15,family= "A", label= "(b)") + 
  coord_cartesian(xlim = c(0, 160),ylim =  c(0, 160))+
  labs(x=expression("1 ha AGBD"["field_1ha"]~"(Mg/ha)"), 
       y=expression("1 ha AGBD"["ALS_area_MCH"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



filedir <- "E:\\Biomass\\CSIR\\result\\CSIR\\biom_CSIR_combo.csv"
Data = read.csv(filedir)
Data <- Data[Data$CSIR >0 ,]
Data<-Data[!(Data$Name1=="Ven_S68_SP_5_Venetia_h.tif"),]
Data$HCC <- Data$H*Data$CC
reg <- lm(CSIR~ HCC, data=Data)
Data$preds <- predict(reg,newdata=Data)
Data$ALS = Data$preds
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
Mean <- mean(Data$CSIR)       
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
rRMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))/mean(Data$CSIR)*100
MD <- mean(Data$ALS - Data$CSIR)

p3<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,
              color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  geom_smooth(method="lm",se = F) +
  annotate("text",x=0,y=152,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=130,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " MSD: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  annotate("text",x=125,y=155,hjust = 0,size = 15,family= "A", label= "(c)") + 
  coord_cartesian(xlim = c(0, 160),ylim =  c(0, 160))+
  labs(x=expression("1 ha AGBD"["field_1ha"]~"(Mg/ha)"), 
       y=expression("1 ha AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggarrange(p1,p2,p3, nrow = 1,ncol = 3)
out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\CSIR_field_CHM_1ha.jpg"
ggsave(out,height=8, width=24, dpi=600)



# ------------------------------------------------------------------------------------------------ #
#scatterplot 25m csir vs. area based and ind biomass
filedir <- "E:\\Biomass\\CSIR\\result\\CSIR_25m\\biom_CSIR_25m_combo_1_h_cc.csv"
Data = read.csv(filedir)
Data <- na.omit(Data)
Data$ALS <-  Data$Ind
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np13"),]
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np15"),]
Data <- Data[Data$CSIR >0 ,]
#stats
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
rRMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))/mean(Data$CSIR)*100
MD <- mean(Data$ALS - Data$CSIR)

p1<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  annotate("text",x=0,y=122,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=110,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,3),"%",
             "\n" , " Sample size: 203")) + 
  annotate("text",x=100,y=120,hjust = 0,size = 15,family= "A", label= "(a)") + 
  coord_cartesian(xlim = c(0, 120),ylim =  c(0, 120))+
  labs(x=expression("AGBD"["field_25m"]~"(Mg/ha)"), 
       y=expression("25 m AGBD"["ALS_tree"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




filedir <- "E:\\Biomass\\CSIR\\Biomass_2018\\CSIR_25m.xlsx"
Data = read_excel(filedir)
Data <- na.omit(Data)
Data$ALS <-  0.3502*Data$MEAN05^2.5405 
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np13"),]
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np15"),]
Data<-Data[!(Data$plot=="Ven_S68_SP_5_np10"),]
Data <- Data[Data$agbd_ha >0,]

#stats
r2 = round(cor(Data$ALS, Data$agbd_ha,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ALS - Data$agbd_ha)^2))
MD <- mean(Data$ALS - Data$agbd_ha)
rRMSE <- sqrt(mean((Data$ALS - Data$agbd_ha)^2))/mean(Data$agbd_ha)*100

p2<- ggplot(Data, aes(x=agbd_ha, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  annotate("text",x=0,y=122,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=110,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,3)," Mg/ha",
             "\n" , " Sample size: 203")) + 
  annotate("text",x=100,y=120,hjust = 0,size = 15,family= "A", label= "(b)") + 
  coord_cartesian(xlim = c(0, 120),ylim =  c(0, 120))+
  labs(x=expression("AGBD"["field_25m"]~"(Mg/ha)"), 
       y=expression("25 m AGBD"["ALS_area_MCH"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



filedir <- "E:\\Biomass\\CSIR\\Biomass_2018\\CSIR_25m.xlsx"
Data = read_excel(filedir)
Data <- na.omit(Data)
Data$ALS <-  9.0665*Data$HCC
Data <- Data[Data$CC15 <1.1,]

Data<-Data[!(Data$plot=="Ven_S63_SP_5_np13"),]
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np15"),]
Data<-Data[!(Data$plot=="Ven_S68_SP_5_np10"),]
Data <- Data[Data$agbd_ha >0,]
#stats
r2 = round(cor(Data$ALS, Data$agbd_ha,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ALS - Data$agbd_ha)^2))
MD <- mean(Data$ALS - Data$agbd_ha)
rRMSE <- sqrt(mean((Data$ALS - Data$agbd_ha)^2))/mean(Data$agbd_ha)*100

p3<- ggplot(Data, aes(x=agbd_ha, y=ALS,col=CC15))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  #geom_point(aes(color=CC15),shape=19, size=3, alpha=0.5)+
  geom_pointdensity(shape=19, size=3, alpha=0.5)+
  scale_color_viridis(direction = -1)+
  theme_bw()+
  annotate("text",x=0,y=122,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=110,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,3)," %",
             "\n" , " Sample size: 203")) + 
  annotate("text",x=100,y=120,hjust = 0,size = 15,family= "A", label= "(c)") + 
  coord_cartesian(xlim = c(0, 120),ylim =  c(0, 120))+
  labs(x=expression("AGBD"["field_25m"]~"(Mg/ha)"), 
       y=expression("25 m AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"))+
  theme(legend.position = c(0.9,0.3),legend.key.height=unit(1.5,"cm"))+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\CSIR_field_25m_cc.jpg"
ggsave(out,height=10, width=10, dpi=600)


ggarrange(p1,p2,p3,nrow = 1,ncol = 3)
out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\CSIR_field_25m.jpg"
ggsave(out,height=10, width=30, dpi=600)





ggplot(Data, aes(x=agbd_ha, y=CC15,col=CC15))+ 
  geom_pointdensity(shape=19, size=3, alpha=1)+
  scale_color_viridis(direction = -1)+
  theme_bw()+
  annotate("text",x=100,y=120,hjust = 0,size = 15,family= "A", label= "(c)") + 
  coord_cartesian(xlim = c(0, 120),ylim =  c(0, 10))+
  labs(x=expression("Mean canopy height (m)"), 
       y=expression("25 m AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"))+
  theme(legend.position = c(0.9,0.3),legend.key.height=unit(1.5,"cm"))+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\CSIR_field_25m_cc2.jpg"
ggsave(out,height=10, width=10, dpi=600)



ggplot(Data, aes(x=agbd_ha, y=MEAN05,col=MEAN05))+ 
  geom_pointdensity(shape=19, size=3, alpha=1)+
  scale_color_viridis(direction = -1)+
  theme_bw()+
  coord_cartesian(xlim = c(0, 120),ylim =  c(0, 10))+
  labs(y=expression("Mean canopy height (m)"), 
       x=expression("25 m AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"))+
  theme(legend.position = c(0.9,0.3),legend.key.height=unit(1.5,"cm"))+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())






ggplot(Data, aes(x=CC15, y=MEAN05))+ 
  geom_pointdensity(shape=19, size=3, alpha=1)+
  scale_color_viridis(direction = -1)+
  theme_bw()+
  coord_cartesian(xlim = c(0, 120),ylim =  c(0, 10))+
  labs(y=expression("Mean canopy height (m)"), 
       x=expression("25 m AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"))+
  theme(legend.position = c(0.9,0.3),legend.key.height=unit(1.5,"cm"))+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plt.scatter(Data$CC15,Data$MEAN05,c=Data$agbd_ha)


png("E:\\Biomass\\CSIR\\Result_0618\\2018\\CSIR_field_25m_3d.png", width = 12, height = 12, units = 'in', res = 600)
par(mgp=c(3,1,0),mar=c(4,5,4,4), family = "A",oma=c(2,2,2,2),cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
scatter3D(Data$CC15,Data$MEAN05,Data$agbd_ha,pch = 19,cex = 2,colkey = TRUE,  
          bty = "u",col.panel ="white",phi = 0,alpha=1,theta=60,surface=TRUE,
          col.grid = "gray50", xlab="CC (%)", ylab="MCH (m)", zlab="AGBD (Mg/ha)")

dev.off()

# ------------------------------------------------------------------------------------------------ #
#als height and cover distribution
Data1 <- read.csv("E:\\Biomass\\CSIR\\ALS_25\\Dny_25.csv",header=T)
Data2 <- read.csv("E:\\Biomass\\CSIR\\ALS_25\\Ven_25.csv",header=T)
Data3 <- read.csv("E:\\Biomass\\CSIR\\ALS_25\\Wel_25.csv",header=T)

Data <- rbind(Data1,Data2,Data3)
Data<-Data1
Data<-Data2
Data<-Data3


#Data$AGB_MCH <- 1.386052*Data$MEAN^1.996775
#Data$AGB_CC <-  10.08957*Data$MEAN*Data$SUM -2.504157


Data <- Data[Data$AGB_MCH >0 & Data$AGB_MCH <200,]
Data <- Data[Data$AGB_CC >0 & Data$AGB_MCH <200,]
Data <- Data[Data$SUM >0,]
Data <- Data[Data$MEAN >0,]

Data$group_AGBD_MCH <- cut(Data$AGB_MCH,
                       breaks = c(0,5,10,15,20,25,30,35,40,45,50,55,60,70,200),
                       dig.lab=3)
Data$group_AGBD_CC <- cut(Data$AGB_CC,
                           breaks = c(0,5,10,15,20,25,30,35,40,45,50,55,60,70,200),
                           dig.lab=3)

ggplot(data=Data) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_histogram(aes(x=group_AGBD_MCH,y = stat(count) / sum(count)*100),
                 color =  "white",fill = "red",stat="count",alpha = 0.2) +
  geom_histogram(aes(x=group_AGBD_CC,y = stat(count) / sum(count)*100),
                 color =  "white",fill = "blue",stat="count",alpha = 0.2) +
  coord_cartesian(ylim =  c(0, 90))+
  xlab("ALS AGB (Mg/ha)") + ylab("Frequency (%)")+
  scale_x_discrete(labels=c("0-5","5-10","10-15","15-20","20-25",
                            "25-30","30-35","35-40","40-45","45-50",
                            "50-55","55-60","60-70",">70"))+
  scale_y_continuous(minor_breaks = round(seq(0,90,10),digits = 1),
                     breaks = round(seq(0,90,10),digits = 1))+
  scale_fill_manual(name="Legend",values=c("red","blue"),labels=c("a","b"))+
  theme(text=element_text(family="A"))+
  theme(text=element_text(size=25),axis.text=element_text(size=25))


out = "E:\\Biomass\\CSIR\\figure\\ALS_AGB_histogram wel.jpg"
ggsave(out,height=9, width=18, dpi=600)


# ------------------------------------------------------------------------------------------------ #
#three model comparison
#Local updated model
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- Data[Data$sensitivity_a1 > 0.95,]
Data$Combo <- 1.386052*Data$MEAN^1.996775
#Data$Combo <- Data$Combo*490.625/625
Data <- Data[Data$Combo >0,]

Data$RH1_50 <- Data$RH1_50 +100
#hist(Data$Combo)

Data <- Data %>% group_by(gr=cut(Combo, breaks= seq(0, 100, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

rsqd_local <- 0
bias_local <- 0
RMSE_local <- 0
bias_group_local <- list(rep(0, 11))
rbias_group_local <- list(rep(0, 11))

rsqd_orbit <- 0
bias_orbit <- 0
RMSE_orbit <- 0
bias_group_orbit <- list(rep(0, 11))
rbias_group_orbit <- list(rep(0, 11))

rsqd_dbt <- 0
bias_dbt <- 0
RMSE_dbt <- 0
bias_group_dbt <- list(rep(0, 11))
rbias_group_dbt <- list(rep(0, 11))

#subsampling
for (i in 1:100){
  
  Data <- Data %>% group_by(gr) %>% slice_sample(n = 50)
  random_sample <- createDataPartition(Data$Combo, p = 0.7, list = FALSE)
  training_dataset  <- Data[random_sample, ]
  testing_dataset <- Data[-random_sample, ]
  
  #Local updated 
  reg <- glm(formula = Combo ~ GEDI_PAVD+RH1_50+RH1_98,data = testing_dataset)
  testing_dataset$preds <- predict(reg,newdata=testing_dataset)
  #Add correction factor
  #C <- mean(reg$fitted.values^2)/mean(testing_dataset$Combo)
  #testing_dataset$preds <- C*((testing_dataset$preds)^2)
  #testing_dataset$preds <- (testing_dataset$preds)^2
  rsqd_local[i] <- round(cor(testing_dataset$Combo,testing_dataset$preds)^2,3)
  bias_local[i] <- mean(testing_dataset$preds - testing_dataset$Combo)
  RMSE_local[i] <- sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))

  testing_dataset <- testing_dataset[testing_dataset$Combo < 100 & testing_dataset$Combo < 100,]
  testing_dataset$biases = testing_dataset$preds - testing_dataset$Combo
  breakbin = round(seq(0,100,5))
  testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)
  bias_group_local <- append(bias_group_local,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]))
  rbias_group_local <- append(rbias_group_local,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]/
                                                       aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100))
  
  #On-orbit GEDI L4A
  rsqd_orbit[i] = round(cor(testing_dataset$AGBD, testing_dataset$Combo,method = "pearson")^2,3)
  RMSE_orbit[i] <- sqrt(mean((testing_dataset$AGBD - testing_dataset$Combo)^2))
  bias_orbit[i] <- mean(testing_dataset$AGBD - testing_dataset$Combo)
  
  testing_dataset$biases = testing_dataset$AGBD - testing_dataset$Combo
  bias_group_orbit <- append(bias_group_orbit,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]))
  rbias_group_orbit <- append(rbias_group_orbit,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]/
                                                       aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100))
  
  #DBT Africa GEDI L4A
  testing_dataset$AGBD_dbt <- 1.092*(-118.408 + 1.957*sqrt(testing_dataset$RH1_50) 
                                 + 9.962 * sqrt(testing_dataset$RH1_98 + 100))^2
  rsqd_dbt[i] = round(cor(testing_dataset$AGBD_dbt, testing_dataset$Combo,method = "pearson")^2,3)
  RMSE_dbt[i] <- sqrt(mean((testing_dataset$AGBD_dbt - testing_dataset$Combo)^2))
  bias_dbt[i] <- mean(testing_dataset$AGBD_dbt - testing_dataset$Combo)
  
  testing_dataset$biases = testing_dataset$AGBD_dbt - testing_dataset$Combo
  bias_group_dbt <- append(bias_group_dbt,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]))
  rbias_group_dbt <- append(rbias_group_dbt,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]/
                                                       aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100))
  
}

mean(rsqd_local)
sd(rsqd_local)
mean(bias_local)
sd(bias_local)
mean(RMSE_local)
sd(RMSE_local)
n<- nrow(testing_dataset)

p1<- ggplot(testing_dataset, aes(x=round(Combo,3), y=round(preds,3)))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=87,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_local),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(mean(RMSE_local),3)," Mg/ha",
             "\n" , " Bias: ", round(mean(bias_local),3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(c)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=expression("AGB"["ALS_area_MCH"]~"(Mg/ha)"), 
       y=expression("AGB"["SAS"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


breakbin = round(seq(0,60,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

table(testing_dataset$group)/nrow(testing_dataset)*100
df_bias_group_local <- do.call(rbind, bias_group_local)
df_bias_group_local <- df_bias_group_local[-1,]
colMeans(df_bias_group_local, na.rm = FALSE, dims = 1)
colSds(df_bias_group_local, na.rm = FALSE, dims = 1)
df_rbias_group_local <- do.call(rbind, rbias_group_local)
df_rbias_group_local <- df_rbias_group_local[-1,]
colMeans(df_rbias_group_local, na.rm = FALSE, dims = 1)
colSds(df_rbias_group_local, na.rm = FALSE, dims = 1)
testing_dataset$biases <- testing_dataset$preds - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)


b1<- ggplot(testing_dataset, aes(x=group, y=biases, group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-60, 60))+
  annotate("text",x=10,y=50,hjust = 0,size = 15,family= "A", label= "(c)") + 
  scale_y_continuous(minor_breaks = round(seq(-60,60,20),digits = 1),
                     breaks = round(seq(-60,60,20),digits = 1))+
  theme_bw()+
  labs(x=expression("AGB"["ALS_area_MCH"]~"(Mg/ha)"), 
       y=expression("Bias of AGB"["SAS"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))



#plot(reg)

#on-orbit GEDI L4A
mean(rsqd_orbit)
sd(rsqd_orbit)
mean(bias_orbit)
sd(bias_orbit)
mean(RMSE_orbit)
sd(RMSE_orbit)

p2 <- ggplot(testing_dataset, aes(x=Combo, y=AGBD))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=0,y=87,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_orbit),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(mean(RMSE_orbit),3)," Mg/ha",
             "\n" , " Bias: ", round(mean(bias_orbit),3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(a)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("AGB"["ALS_area_MCH"]~"(Mg/ha)"), 
       y=expression("AGB"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


breakbin = round(seq(0,60,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

df_bias_group_orbit<- do.call(rbind, bias_group_orbit)
df_bias_group_orbit <- df_bias_group_orbit[-1,]
colMeans(df_bias_group_orbit, na.rm = FALSE, dims = 1)
colSds(df_bias_group_orbit, na.rm = FALSE, dims = 1)
df_rbias_group_orbit <- do.call(rbind, rbias_group_orbit)
df_rbias_group_orbit <- df_rbias_group_orbit[-1,]
colMeans(df_rbias_group_orbit, na.rm = FALSE, dims = 1)
colSds(df_rbias_group_orbit, na.rm = FALSE, dims = 1)
testing_dataset$biases <- testing_dataset$AGBD - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)

b2<- ggplot(testing_dataset, aes(x=group, y=biases, group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  annotate("text",x=10,y=50,hjust = 0,size = 15,family= "A", label= "(a)") + 
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-60, 60))+
  scale_y_continuous(minor_breaks = round(seq(-60,60,20),digits = 1),
                     breaks = round(seq(-60,60,20),digits = 1))+
  theme_bw()+
  labs(x=expression("AGB"["ALS_area_MCH"]~"(Mg/ha)"), 
       y=expression("Bias of AGB"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


#DBT Africa model
mean(rsqd_dbt)
sd(rsqd_dbt)
mean(bias_dbt)
sd(bias_dbt)
mean(RMSE_dbt)
sd(RMSE_dbt)

p3 <- ggplot(testing_dataset, aes(x=Combo, y=AGBD_dbt))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=0,y=87,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_dbt),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(mean(RMSE_dbt),3)," Mg/ha",
             "\n" , " Bias: ", round(mean(bias_dbt),3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(b)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("AGB"["ALS_area_MCH"]~"(Mg/ha)"), 
       y=expression("AGB"["DBT_Africa"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



df_bias_group_dbt <- do.call(rbind, bias_group_dbt)
df_bias_group_dbt <- df_bias_group_dbt[-1,]
colMeans(df_bias_group_dbt, na.rm = FALSE, dims = 1)
colSds(df_bias_group_dbt, na.rm = FALSE, dims = 1)
df_rbias_group_dbt <- do.call(rbind, rbias_group_dbt)
df_rbias_group_dbt <- df_rbias_group_dbt[-1,]
colMeans(df_rbias_group_dbt, na.rm = FALSE, dims = 1)
colSds(df_rbias_group_dbt, na.rm = FALSE, dims = 1)
testing_dataset$biases <- testing_dataset$AGBD_dbt - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)

b3<- ggplot(testing_dataset, aes(x=group, y=biases, group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  annotate("text",x=10,y=50,hjust = 0,size = 15,family= "A", label= "(b)") + 
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-60, 60))+
  scale_y_continuous(minor_breaks = round(seq(-60,60,20),digits = 1),
                     breaks = round(seq(-60,60,20),digits = 1))+
  theme_bw()+
  labs(x=expression("AGB"["ALS_area_MCH"]~"(Mg/ha)"), 
       y=expression("Bias of AGB"["DBT_Africa"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))



#original
ggarrange(p2,p3,p1,ncol=3)
out = "E:\\Biomass\\CSIR\\figure\\three model comparison scatterplot.jpg"
ggsave(out,height=9, width=27, dpi=600)

ggarrange(b2,b3,b1,ncol=3)
out = "E:\\Biomass\\CSIR\\figure\\three model comparison biasplot.jpg"
ggsave(out,height=9, width=33, dpi=600)



# ------------------------------------------------------------------------------------------------ #
#Test PCA
library(FactoMineR)
library(ggcorrplot)
library(devtools)
library(factoextra)
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- Data[Data$sensitivity_a1 > 0.95,]
Data$Combo <-  4*Data$MEAN^1.19
Data <- Data[Data$Combo >0 & Data$Combo < 60,]

Data$RH1_98 <- sqrt(Data$RH1_98)
Data$RH1_50 <- sqrt(Data$RH1_50+5)
Data$Combo <- sqrt(Data$Combo)

Data_sel <- Data[,c("RH1_10","RH1_30","RH1_50","RH1_98","GEDI_Cover",
                    "GEDI_FHD","GEDI_PAI","GEDI_PAVD")]
Data_sel <- na.omit(Data_sel)

PCA_result <- princomp(Data_sel)
summary(PCA_result)
fviz_eig(PCA_result, addlabels = TRUE)
PCA_result$loadings[, 1:2]
fviz_pca_var(PCA_result, col.var = "black")
fviz_cos2(PCA_result, choice = "var", axes = 1:2)

Data_sel<- Data[,c("Combo")]
pcs <- as.data.frame(PCA_result$scores[,1:2])
Result_df <- cbind(Data_sel, pcs)
cor(Data$RH1_98,Data$GEDI_FHD)


random_sample <- createDataPartition(Result_df$Data_sel, p = 0.7, list = FALSE)
training_dataset  <- Result_df[random_sample, ]
testing_dataset <- Result_df[-random_sample, ]

repeat_cv <- trainControl(method='cv', number=5)
rf_model <- train(
  Data_sel~., 
  data=training_dataset, 
  method='rf', 
  ntree=100,
  nodesize = 2,
  trControl=repeat_cv)


rf_model$finalModel


preds <- predict(object=rf_model, 
                 newdata=testing_dataset[, ])

testing_dataset <- cbind(testing_dataset, preds)

testing_dataset$preds<- (testing_dataset$preds)^2
testing_dataset$Data_sel<- (testing_dataset$Data_sel)^2


r2 <- round(cor(testing_dataset$preds,testing_dataset$Data_sel)^2,3)
bias <- mean(testing_dataset$preds - testing_dataset$Data_sel)
RMSE <- sqrt(mean((testing_dataset$preds - testing_dataset$Data_sel)^2))
n<- nrow(testing_dataset)
ggplot(testing_dataset, aes(x=Data_sel, y=preds))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 60),ylim =  c(0, 60))+
  annotate("text",x=0,y=55,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(r2),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=48,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(mean(RMSE),3)," Mg/ha",
             "\n" , " Bias: ", round(mean(bias),3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Area-based (H only) AGB (Mg/ha)"), 
       y=expression("PCA+RF GEDI AGB (Mg/ha) (No correction)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# ------------------------------------------------------------------------------------------------ #
#test residual
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
#Data <- Data[Data$sensitivity_a1 > 0.95,]
Data$Combo <- 1.386052*Data$MEAN^1.996775
#Data$Combo <- Data$Combo*490.625/625
Data <- Data[Data$Combo >0,]
Data$RH1_50 <- Data$RH1_50 +100
#Data <- Data %>% group_by(gr) %>% slice_sample(n=50)
random_sample <- createDataPartition(Data$Combo, p = 0.7, list = FALSE)
training_dataset  <- Data[random_sample, ]
testing_dataset <- Data[-random_sample, ]



reg <- glm(formula = Combo ~ GEDI_PAVD+RH1_50+RH1_98,
           data = training_dataset,Gamma(link=identity), start=c(0,0.1,100,5))

#Weighted Least Squares Regression
wt <- 1 / lm(abs(reg$residuals) ~ reg$fitted.values)$fitted.values^2
reg <- lm(formula = Combo ~ GEDI_PAVD+RH1_50+RH1_98,
          data = training_dataset, weights=wt)


testing_dataset$preds <- predict(reg,newdata=testing_dataset)
testing_dataset$res <-  sqrt(abs(testing_dataset$preds - testing_dataset$Combo))
summary(reg)

p1<- ggplot(testing_dataset, aes(x=round(preds,3), y=round(res,3)))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 5),ylim = c(0, 10))+
  scale_color_viridis(direction = 1)+
  annotate("text",x=60,y=50,hjust = 0,size = 15,family= "A", label= "(c)") + 
  labs(x=expression("AGB"["SAS"]~"(Mg/ha)"), 
       y=expression("Sqrt of absolute residuals of AGB"["SAS"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))




testing_dataset$res_l4 <-  sqrt(abs(testing_dataset$AGBD - testing_dataset$Combo))

p2<- ggplot(testing_dataset, aes(x=round(AGBD,3), y=round(res_l4,3)))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 75),ylim = c(0, 10))+
  scale_color_viridis(direction = 1)+
  annotate("text",x=60,y=50,hjust = 0,size = 15,family= "A", label= "(a)") + 
  labs(x=expression("AGB"["L4A"]~"(Mg/ha)"), 
       y=expression("Sqrt of absolute residuals of AGB"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


testing_dataset$AGBD_dbt <- 1.092*(-118.408 + 1.957*sqrt(testing_dataset$RH1_50) 
                                   + 9.962 * sqrt(testing_dataset$RH1_98 + 100))^2

testing_dataset$res_dbt <-  sqrt(abs(testing_dataset$AGBD_dbt - testing_dataset$Combo))


p3<- ggplot(testing_dataset, aes(x=round(AGBD_dbt,3), y=round(res_dbt,3)))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 75),ylim = c(0, 10))+
  scale_color_viridis(direction = 1)+
  annotate("text",x=60,y=50,hjust = 0,size = 15,family= "A", label= "(b)") + 
  labs(x=expression("AGB"["DBT_Africa"]~"(Mg/ha)"), 
       y=expression("Sqrt of absolute residuals of AGB"["DBT_Africa"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


#original
ggarrange(p2,p3,p1,ncol=3)
out = "E:\\Biomass\\CSIR\\figure\\sqrt abs residual plot.jpg"
ggsave(out,height=9, width=27, dpi=600)



# ------------------------------------------------------------------------------------------------ #
#test heteroscedasticity -- robust standard error
library(lmtest)
library(sandwich)
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
#Data <- Data[Data$sensitivity_a1 > 0.95,]
Data$Combo <- 1.386052*Data$MEAN^1.996775
#Data$Combo <- Data$Combo*490.625/625
Data <- Data[Data$Combo >0,]
Data$RH1_50 <- Data$RH1_50 +100
#Data <- Data %>% group_by(gr) %>% slice_sample(n=50)
random_sample <- createDataPartition(Data$Combo, p = 0.7, list = FALSE)
training_dataset  <- Data[random_sample, ]
testing_dataset <- Data[-random_sample, ]

#local Updated
#reg <- lm(sqrt(Combo) ~ sqrt(RH1_98), data = training_dataset)
reg <- glm(formula = Combo ~ GEDI_PAVD+RH1_50+RH1_98,
           data = training_dataset)
testing_dataset$preds <- predict(reg,newdata=testing_dataset)
testing_dataset$res <-  testing_dataset$Combo - testing_dataset$preds
summary(reg)

reg <- glm(formula = Combo ~ GEDI_PAVD+RH1_50+RH1_98,
           data = training_dataset)
coeftest(reg, vcov = vcovHC(reg, type = "HC"))



# ------------------------------------------------------------------------------------------------ #
#three model comparison
#local updated model
dir <- "E:\\Biomass\\CSIR\\result\\All_02272023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- Data[Data$sensitivity_a1 > 0.95,]
Data$Combo <- 1.386052*Data$MEAN^1.996775
#Data$Combo <- Data$Combo*490.625/625
Data <- Data[Data$Combo >0,]

#hist(Data$Combo)

Data <- Data[,c("AGBD","Combo","MEAN","GEDI_PAVD", 
                                      "GEDI_PAI", "GEDI_FHD", "GEDI_Cover", "RH1_98", "RH1_50")]
Data$Combo <- round(Data$Combo,3)
Data$MEAN <- round(Data$MEAN,3)
write.table(Data,file='C:\\Users\\Shawn\\Desktop\\Output.csv',sep = ",",col.names=TRUE)





# ------------------------------------------------------------------------------------------------ #
#biomass 2018-2012
#boxplot
dir <- "E:\\Biomass\\CSIR\\ALS_25\\Justicia_25_2018_2012.csv"
Data = read.csv(dir)
Data <- Data[Data$Diff <80 & Data$Diff >-80,]

Data$group_AGBD <- cut(Data$Diff,
                       breaks = round(seq(-80,80,10),2),
                       dig.lab=3)

ggplot(data=Data) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_histogram(aes(x=group_AGBD,y = stat(count) / sum(count)*100),
                 color =  "white",fill = "blue",stat="count",alpha = 0.4) +
  coord_cartesian(ylim =  c(0, 70))+
  xlab("Difference between ALS derived 2018-2012 AGBD (Mg/ha)") + 
  ylab("Frequency (%)")+
  scale_y_continuous(minor_breaks = round(seq(0,70,10),digits = 1),
                     breaks = round(seq(0,70,10),digits = 1))+
  theme(text=element_text(family="A"))+
  theme(text=element_text(size=25),axis.text=element_text(size=25))




# ------------------------------------------------------------------------------------------------ #
# 1m CHM 2018 vs 2012
dir <- "E:\\Biomass\\CSIR\\ALS_1\\Welverdiendt.csv"
Data = read.csv(dir)
Data <- Data[Data$grid_code < 30 & Data$RASTERVALU <30,]
RMSE <- sqrt(mean((Data$grid_code - Data$RASTERVALU)^2))
MD <- mean(Data$grid_code - Data$RASTERVALU)
R2 <- round(cor(Data$grid_code, Data$RASTERVALU, method = "pearson")^2,3)



# ------------------------------------------------------------------------------------------------ #
# 25m CHM 2018 vs 2012
dir <- "E:\\Biomass\\CSIR\\Result_0606\\zonal_table\\Welverdiendt_index.csv"
Data_index = read.csv(dir)
Data_index = subset(Data_index, select=c(1))
colnames(Data_index)[1] <-  "ID"

dir <- "E:\\Biomass\\CSIR\\Result_0606\\zonal_table\\Welverdiendt_2012_H.csv"
Data = read.csv(dir)
Data_sel_2012_H = subset(Data, select=c(2,5))
colnames(Data_sel_2012_H)[1] <-  "ID"

dir <- "E:\\Biomass\\CSIR\\Result_0606\\zonal_table\\Welverdiendt_2018_H.csv"
Data = read.csv(dir)
Data_sel_2018_H = subset(Data, select=c(2,5))
colnames(Data_sel_2018_H)[1] <-  "ID"

dir <- "E:\\Biomass\\CSIR\\Result_0606\\zonal_table\\Welverdiendt_2012_AGBD.csv"
Data = read.csv(dir)
Data_sel_2012_AGBD = subset(Data, select=c(2,5))
colnames(Data_sel_2012_AGBD)[1] <-  "ID"

Data_index <- merge(Data_index, Data_sel_2012_H, by.x = "ID", by.y = "ID")
Data_index <- merge(Data_index, Data_sel_2018_H, by.x = "ID", by.y = "ID")
Data_index <- merge(Data_index, Data_sel_2012_AGBD, by.x = "ID", by.y = "ID")

new_column_names <- c("Index","Data_sel_2012_H", "Data_sel_2018_H", "Data_sel_2012_AGBD")
colnames(Data_index) <- new_column_names

Data_index$Data_sel_2018_AGBD <- round(1.386052*Data_index$Data_sel_2018_H^1.996775,3)

Data_index <- Data_index[Data_index$Data_sel_2018_H < 30 & Data_index$Data_sel_2012_H <30,]
Data_index <- Data_index[Data_index$Data_sel_2018_H > 0.5 & Data_index$Data_sel_2012_H <30,]

Data_index <- Data_index[Data_index$Data_sel_2018_AGBD < 100 & Data_index$Data_sel_2012_AGBD < 100,]

Data_Justicia <- Data_index
Data_Agincourt <- Data_index
Data_Ireagh <- Data_index
Data_Welverdiendt <- Data_index

Data_index <- rbind(Data_Agincourt,Data_Ireagh,Data_Justicia,Data_Welverdiendt)

RMSE <- sqrt(mean((Data_index$Data_sel_2018_H - Data_index$Data_sel_2012_H)^2))
MD <- mean(Data_index$Data_sel_2018_H - Data_index$Data_sel_2012_H)
R2 <- round(cor(Data_index$Data_sel_2018_H, Data_index$Data_sel_2012_H, method = "pearson")^2,3)

RMSE <- sqrt(mean((Data_index$Data_sel_2018_AGBD - Data_index$Data_sel_2012_AGBD)^2))
MD <- mean(Data_index$Data_sel_2018_AGBD - Data_index$Data_sel_2012_AGBD)
R2 <- round(cor(Data_index$Data_sel_2018_AGBD, Data_index$Data_sel_2012_AGBD, method = "pearson")^2,3)

Data_index <- Data_index %>% group_by(gr=cut(Data_sel_2012_AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1
Data_index <- Data_index %>% group_by(gr) %>% slice_sample(n=100)

ggplot(Data_index, aes(x=Data_sel_2012_H, y=Data_sel_2018_H))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_viridis(direction = 1)+
  labs(x="Mean of ALS CHM 2012 (m)", 
       y="Mean of ALS CHM 2018 (m)")+
  scale_x_continuous(limits = c(0,15),minor_breaks = seq(0,15,2))+
  scale_y_continuous(limits = c(0,15),minor_breaks = seq(0,15,2))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  #theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(Data_index, aes(x=Data_sel_2012_AGBD, y=Data_sel_2018_AGBD))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_viridis(direction = 1)+
  labs(x="Mean of ALS AGBD 2012 (m)", 
       y="Mean of ALS AGBD 2018 (m)")+
  scale_x_continuous(limits = c(0,100),minor_breaks = seq(0,100,10))+
  scale_y_continuous(limits = c(0,100),minor_breaks = seq(0,100,10))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  #theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))




# ------------------------------------------------------------------------------------------------ #
# 25m CHM 2018 vs 2012 difference
dir <- "E:\\Biomass\\CSIR\\Result_0606\\zonal_table\\Ireagh_index.csv"
Data_index = read.csv(dir)
Data_index = subset(Data_index, select=c(1))
colnames(Data_index)[1] <-  "ID"

dir <- "E:\\Biomass\\CSIR\\Result_0606\\zonal_table\\Ireagh_2012_H.csv"
Data = read.csv(dir)
Data_sel_2012_H = subset(Data, select=c(2,5))
colnames(Data_sel_2012_H)[1] <-  "ID"

dir <- "E:\\Biomass\\CSIR\\Result_0606\\zonal_table\\Ireagh_2018_H.csv"
Data = read.csv(dir)
Data_sel_2018_H = subset(Data, select=c(2,5))
colnames(Data_sel_2018_H)[1] <-  "ID"

Data_index <- merge(Data_index, Data_sel_2012_H, by.x = "ID", by.y = "ID")
Data_index <- merge(Data_index, Data_sel_2018_H, by.x = "ID", by.y = "ID")

new_column_names <- c("Index","Data_sel_2012_H", "Data_sel_2018_H")
colnames(Data_index) <- new_column_names

Data_index$diff <- Data_index$Data_sel_2018_H - Data_index$Data_sel_2012_H
Data_index = subset(Data_index, select=c(1,4))


Data_index$Index <- Data_index$Index  + 1

write.table(Data_index,
            file='E:\\Biomass\\CSIR\\Result_0606\\result\\Ireagh_diff.csv',
            sep = ",",col.names=TRUE,row.names=F)





# ------------------------------------------------------------------------------------------------ #
#SAS2012 vs. SAS 2018
dir <- "E:\\Biomass\\CSIR\\Result_0615\\Ireagh_2012_AGB.csv"
Data = read.csv(dir)
Data_sel_2012_AGBD = subset(Data, select=c(3))
colnames(Data_sel_2012_AGBD)[1] <-  "SAS_2012"

dir <- "E:\\Biomass\\CSIR\\Result_0615\\Ireagh_2018_AGB.csv"
Data = read.csv(dir)
Data_sel_2018_AGBD = subset(Data, select=c(3))
colnames(Data_sel_2018_AGBD)[1] <-  "SAS_2018"

Data_index <- cbind(Data_sel_2012_AGBD,Data_sel_2018_AGBD)
Data_index <- na.omit(Data_index)

Data_index <- Data_index[Data_index$SAS_2012 < 100 & Data_index$SAS_2018 < 100,]

RMSE <- sqrt(mean((Data_index$SAS_2018 - Data_index$SAS_2012)^2))
MD <- mean(Data_index$SAS_2018 - Data_index$SAS_2012)
R2 <- round(cor(Data_index$SAS_2018, Data_index$SAS_2012, method = "pearson")^2,3)


ggplot(Data_index, aes(x=SAS_2012, y=SAS_2018))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_viridis(direction = 1)+
  labs(x="SAS AGBD 2012 (m)", 
       y="SAS AGBD 2018 (m)")+
  annotate("text",x=10,y=90,hjust = 0,size = 8,family= "A",
           label= paste(" R2: ",round(R2,3),
             "\n" , "RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , "Bias: ", round(MD,3),"Mg/ha")) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  scale_x_continuous(limits = c(0,100),minor_breaks = seq(0,100,10))+
  scale_y_continuous(limits = c(0,100),minor_breaks = seq(0,100,10))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  #theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))

out = "E:\\Biomass\\CSIR\\Result_0615\\Ireagh.jpg"
ggsave(out,height=15, width=15, dpi=600)


# ------------------------------------------------------------------------------------------------ #
#SAS2012 vs. SAS 2018
dir <- "E:\\Biomass\\CSIR\\Result_0615"
dir_index <- file.path(dir, "Agincourt_2018_AGB.csv")
Data_index = read.csv(dir_index)
colnames(Data_index)[3] <- "AGB"

list_r2 <- list()
listt <- list("ScanSAR", "R2", "RMSE")
df_listt <- data.frame(listt)
colnames(df_listt) <- df_listt

files <- list.files(path = dir, full.names = TRUE, pattern = "Agincourt_SAR")

for (i in files){
  print(i)
  Data_sar = read.csv(i)
  Data_sar = subset(Data_sar, select=c(3))
  colnames(Data_sar)[1] <- "SAR"
  Data <- cbind(Data_index, Data_sar)
  Data <- Data[Data$SAR > -35,]
  Data <- Data[Data$AGB < 100,]

  reg <- lm(log(AGB) ~ SAR, data = Data)
  summary(reg)
  x<-Data$SAR
  y<-Data$AGB
  
  mult_nls <- nls(y ~ a*exp(b*x), start = list(a=1000, b=0.1))
  pred <- coef(mult_nls)[1]*exp(coef(mult_nls)[2]*x)
  RMSE<- sqrt(mean((pred - y)^2))
  
  list_r2[[1]] <- file_path_sans_ext(basename(i))
  list_r2[[2]] <- round(cor(log(Data$agbd), Data$SAR, method = "pearson")^2,3)
  list_r2[[3]] <- RMSE
  listt <- mapply(c, listt, list_r2, SIMPLIFY=FALSE)
  df1 <- data.frame(list_r2)
  colnames(df1) <- listt <- list("ScanSAR", "R2", "RMSE")
  df_listt <- rbind(df_listt,df1)
}


# ------------------------------------------------------------------------------------------------ #
#AGBD 2012 vs. SAR 2014 2015
dir <- "E:\\ScanSAR\\ScanSAR\\site_75m"
list_dir <- c("Agincourt","Ireagh","Justicia","Welverdiendt")
for (l in list_dir){
  print(l)
  dir_index <- file.path(dir,l,"index.csv")
  Data_index = read.csv(dir_index)
  colnames(Data_index)[1] <- "ID"
  i <- file.path(dir,l,paste0("CHM_2012_",l,"_75m.xlsx"))
  Data_sar = read_excel(i)
  Data_sar$AGBD <- 0.5956*Data_sar$MEAN^2.5357
  Data_sar = subset(Data_sar, select=c(2,6))
  colnames(Data_sar)[2] <-  str_sub(basename(i), end=-6)
  write.csv(Data_sar, file.path(dir,l,paste0("AGBD_2012_",l,"_75m.csv")), row.names = F)
  #Data_index <- merge(Data_index, Data_sar, by.x = "ID", by.y = "FID")
}


#stats
sar_date <- "_20140906_75m.csv"
sar_date <- "_20141129_75m.csv"
sar_date <- "_20150110_75m.csv"
sar_date <- "_20150221_75m.csv"
sar_date <- "_20150404_75m.csv"

#list_dir <- c("Agincourt","Ireagh","Justicia","Welverdiendt")
dir <- "E:\\ScanSAR\\ScanSAR\\site_75m"
l <- "Agincourt"
dir_index <- file.path(dir,l,"index.csv")
Data_index = read.csv(dir_index)
colnames(Data_index)[1] <- "ID"

dir_chm <- file.path(dir,l,paste0("AGBD_2012_",l,"_75m.csv"))
Data_chm = read.csv(dir_chm)
Data_index <- merge(Data_index, Data_chm, by.x = "ID", by.y = "FID")


dir_sar <- file.path(dir,l,paste0("ScanSAR_",l,sar_date))
Data_sar = read.csv(dir_sar)
Data_index <- merge(Data_index, Data_sar, by.x = "ID", by.y = "FID")

Data <- na.omit(Data_index)
colnames(Data) <- c("ID","AGB","SAR")

Data <- Data[is.finite(rowSums(Data)),]
Data <- Data[Data$AGB < 50 & Data$AGB > 0,]
Data <- Data[Data$SAR > -30 & Data$SAR < -10,]

reg <- lm(log(AGB) ~ SAR, data = Data)
summary(reg)



ggplot(Data, aes(x=AGB, y=SAR))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_viridis(direction = 1)+
  labs(x="AGB", 
       y="SAR")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  #theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))



x<-Data$SAR
y<-Data$AGB

mult_nls <- nls(y ~ a*exp(b*x), start = list(a=1000, b=0.1))
pred <- coef(mult_nls)[1]*exp(coef(mult_nls)[2]*x)
RMSE<- sqrt(mean((pred - y)^2))





# ------------------------------------------------------------------------------------------------ #
#test agb histogram per ALS site Colgan 2012 and SAS 2012
folder_path <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012"
files <- list.files(folder_path,pattern = "\\.tif$",  full.names = TRUE)
plot_list <- list()
for (i in 1:8){
  r <- raster(files[i])
  r[r > 100] <- NA
  year <- basename(files[i])
  df <- data.frame(rasterToPoints(r))
  new_column_names <- c("x","y", "agb")
  colnames(df) <- new_column_names
  df$group_agb <- cut(df$agb,breaks = seq(0,100, 10),dig.lab = 5)
  
  p <- ggplot(df, aes(x=group_agb)) +     
    theme_bw()+
    coord_cartesian(ylim =  c(0, 0.8))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = year,
         y = "Frequency")+
    geom_histogram(aes(y = stat(count) / sum(count)),
                   color =  "blue",fill = "red",stat="count",alpha = 0.8) +
    scale_x_discrete(labels=c("[0,10]","[10,20]","[20,30]","[30,40]",
                              "[40,50]","[50,60]","[60,70]","[70,80]",
                              "[80,90]","[90,100]"))+
    theme(legend.title = element_blank())+
    theme(text=element_text(size=20))
  
  plot_list[[i]] = p
}

ggarrange(plotlist=plot_list[7:8],nrow=1,ncol=2)
out <- "E:\\Biomass\\CSIR\\Result_0616\\Welverdiendt.jpg"
ggsave(out,height=10, width=20, dpi=600)




# ------------------------------------------------------------------------------------------------ #
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2012\\Agincourt_AGB_2012_Naidoo.tif"
r <- raster(dir)
r[r > 100] <- NA
year <- basename(dir)
df <- data.frame(rasterToPoints(r))
new_column_names <- c("x","y", "agb")

colnames(df) <- new_column_names
df$group_agb <- cut(df$agb,breaks = seq(0,100, 10),dig.lab = 5)

ggplot(df, aes(x=group_agb)) +     
  theme_bw()+
  coord_cartesian(ylim =  c(0, 0.8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = year,
       y = "Frequency")+
  geom_histogram(aes(y = stat(count) / sum(count)),
                 color =  "blue",fill = "red",stat="count",alpha = 0.8) +
  scale_x_discrete(labels=c("[0,10]","[10,20]","[20,30]","[30,40]",
                            "[40,50]","[50,60]","[60,70]","[70,80]",
                            "[80,90]","[90,100]"))+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=20))





# ------------------------------------------------------------------------------------------------ #
#all to ind csv
sites <- c('Agincourt', 'DNyala', 'Welverdiendt','Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh','Venetia')
dir_GEDI <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\All_02272023.csv"
for (ss in sites){
  Data = read.csv(dir_GEDI,header=T)
  Data <- na.omit(Data)
  Data <- Data[Data$status == "Leaf-on",]
  Data <- Data[Data$ql2a ==1,]
  Data <- Data[Data$ql4 ==1,]
  Data <- Data[Data$ql2a ==1,]
  Data <- Data[Data$ql2b ==1,]
  Data <- Data[Data$site ==ss,]
  colnames(Data)[1] <- "index"
  output <- file.path( "E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(ss,".csv"))
  write.table(Data,file=output,sep = ",",col.names=TRUE,row.names=F)
}


# ------------------------------------------------------------------------------------------------ #
#scatterplot 06192023
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
#Data <- Data[Data$RHS_98 <=2.35,]
#Data <- Data[Data$solar_elevation <0,]
#Data <- Data[Data$status == "Leaf-off",]

#not good
Data$SAS <- 18.786*Data$MEAN*Data$SUM/490 - 11.477
#relatively good
Data$SAS <- 10.35*Data$MEAN*Data$SUM/490 - 5.9236

Data <- Data[Data$SAS >0 & Data$SAS < 100,]
Data <- Data[Data$AGBD >0 & Data$AGBD < 100,]

r2 = round(cor(Data$SAS, Data$AGBD,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$SAS - Data$AGBD)^2))
MD <- mean(Data$SAS - Data$AGBD)



ggplot(Data, aes(x=round(SAS,3), y=round(AGBD,3)))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=87,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(r2),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(mean(RMSE),3)," Mg/ha",
             "\n" , " Bias: ", round(mean(MD),3)," Mg/ha")) + 
  #annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(c)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=expression("AGB"["SAS"]~"(Mg/ha)"), 
       y=expression("AGB"["GEDI"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))



Data$group_GEDI_AGB <- cut(Data$AGBD,breaks = seq(0,100, 10),dig.lab = 5)

p1<- ggplot(Data, aes(x=group_GEDI_AGB)) +     
  theme_bw()+
  coord_cartesian(ylim =  c(0, 0.8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "GEDI AGB (Mg/ha)",
       y = "Frequency")+
  geom_histogram(aes(y = stat(count) / sum(count)),
                 color =  "blue",fill = "red",stat="count",alpha = 0.8) +
  scale_x_discrete(labels=c("[0,10]","[10,20]","[20,30]","[30,40]",
                            "[40,50]","[50,60]","[60,70]","[70,80]",
                            "[80,90]","[90,100]"))+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=20))


Data$group_SAS <- cut(Data$SAS,breaks = seq(0,100, 10),dig.lab = 5)

p2<- ggplot(Data, aes(x=group_SAS)) +     
  theme_bw()+
  coord_cartesian(ylim =  c(0, 0.8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "SAS AGB (Mg/ha)",
       y = "Frequency")+
  geom_histogram(aes(y = stat(count) / sum(count)),
                 color =  "blue",fill = "red",stat="count",alpha = 0.8) +
  scale_x_discrete(labels=c("[0,10]","[10,20]","[20,30]","[30,40]",
                            "[40,50]","[50,60]","[60,70]","[70,80]",
                            "[80,90]","[90,100]"))+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=20))

ggarrange(p1,p2)
out <- 'E:\\Biomass\\CSIR\\Result_0618\\2018\\Histogram.jpg'
ggsave(out,height=12, width=24, dpi=600)



Data$bias <- Data$AGBD - Data$SAS
Data <- Data[Data$SAS >0 & Data$SAS < 60,]
Data$group_SAS <- cut(Data$SAS,breaks = seq(0,60, 5),dig.lab = 5)

ggplot(Data, aes(x=group_SAS, y=bias,group = group_SAS)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 80))+
  scale_y_continuous(minor_breaks = round(seq(-80,80,20),digits = 1),
                     breaks = round(seq(-80,80,20),digits = 1))+
  theme_bw()+
  labs(x = "25m GEDI L4A AGB (Mg/ha)",
       y="Bias of predicted GEDI AGB (H*CC) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


train_control <- trainControl(method = "cv", number = 10, savePredictions=TRUE)
model <- train(AGBD ~SAS, data = Data, method = "lm", trControl = train_control)
pred <- model$pred
pred <- pred[pred$obs < 60,]
pred$biases = pred$pred - pred$obs
pred$group <- cut(pred$obs,breaks = round(seq(0,60,5),2),dig.lab=5)
table(pred$group)
table(pred$group)/nrow(pred)
aggregate(pred$biases, list(pred$group), mean)
aggregate(pred$biases, list(pred$group), mean)[,2]/
  aggregate(pred$obs, list(pred$group), mean)[,2]*100






# ------------------------------------------------------------------------------------------------ #
#model comparison 06192023
#Local updated model
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]


suf <- "ALS_area_based_H_05"
Data$Combo <- 0.3502*Data$MEAN^2.5405 


suf <- "ALS_area_based_H_CC_linear"
Data$Combo <- 10.35*Data$MEAN*Data$SUM/490 - 5.9236



Data <- Data[Data$AGBD > 0 & Data$AGBD < 100,]
Data <- Data[Data$Combo > 0 & Data$Combo < 100,]
Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(Combo, breaks= seq(0, 100, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

rsqd_local <- 0
bias_local <- 0
RMSE_local <- 0
bias_group_local <- list(rep(0, 11))
rbias_group_local <- list(rep(0, 11))

rsqd_orbit <- 0
bias_orbit <- 0
RMSE_orbit <- 0
bias_group_orbit <- list(rep(0, 11))
rbias_group_orbit <- list(rep(0, 11))



#subsampling
for (i in 1:100){
  
  #Data <- Data %>% group_by(gr) %>% slice_sample(n=50)
  random_sample <- createDataPartition(Data$Combo, p = 0.7, list = FALSE)
  training_dataset  <- Data[random_sample, ]
  testing_dataset <- Data[-random_sample, ]
  
  #Local updated 
  reg <- glm(formula = Combo ~ GEDI_PAVD+RH1_50+RH1_98,data = training_dataset)
  testing_dataset$preds <- predict(reg,newdata=testing_dataset)

  rsqd_local[i] <- round(cor(testing_dataset$Combo,testing_dataset$preds)^2,3)
  bias_local[i] <- mean(testing_dataset$preds - testing_dataset$Combo)
  RMSE_local[i] <- sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))
  

  testing_dataset$biases = testing_dataset$preds - testing_dataset$Combo
  breakbin = round(seq(0,100,5))
  testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)
  bias_group_local <- append(bias_group_local,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]))
  rbias_group_local <- append(rbias_group_local,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]/
                                                       aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100))
  
  #On-orbit GEDI L4A
  rsqd_orbit[i] = round(cor(testing_dataset$AGBD, testing_dataset$Combo,method = "pearson")^2,3)
  RMSE_orbit[i] <- sqrt(mean((testing_dataset$AGBD - testing_dataset$Combo)^2))
  bias_orbit[i] <- mean(testing_dataset$AGBD - testing_dataset$Combo)
  
  testing_dataset$biases = testing_dataset$AGBD - testing_dataset$Combo
  bias_group_orbit <- append(bias_group_orbit,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]))
  rbias_group_orbit <- append(rbias_group_orbit,list(aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]/
                                                       aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100))
}


mean(rsqd_local)
mean(bias_local)
mean(RMSE_local)
n<- nrow(testing_dataset)


p1<- ggplot(testing_dataset, aes(x=Combo, y=preds))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=87,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_local),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(mean(RMSE_local),3)," Mg/ha",
             "\n" , " Bias: ", round(mean(bias_local),3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(b)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=paste0("AGBD_",suf, " (Mg/ha)"), 
       y=expression("AGBD"["SAS"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

table(testing_dataset$group)/nrow(testing_dataset)*100
df_bias_group_local <- do.call(rbind, bias_group_local)
df_bias_group_local <- df_bias_group_local[-1,]
colMeans(df_bias_group_local, na.rm = FALSE, dims = 1)
df_rbias_group_local <- do.call(rbind, rbias_group_local)
df_rbias_group_local <- df_rbias_group_local[-1,]
colMeans(df_rbias_group_local, na.rm = FALSE, dims = 1)
testing_dataset$biases <- testing_dataset$preds - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)
table(testing_dataset$group)/nrow(testing_dataset)

b1<- ggplot(testing_dataset, aes(x=group, y=biases, group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-60, 60),xlim = c(0, 12))+
  annotate("text",x=10,y=50,hjust = 0,size = 15,family= "A", label= "(b)") + 
  scale_y_continuous(minor_breaks = round(seq(-60,60,20),digits = 1),
                     breaks = round(seq(-60,60,20),digits = 1))+
  theme_bw()+
  labs(x=paste0("AGBD_",suf, " (Mg/ha)"), 
       y=expression("Bias of AGBD"["SAS"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


#on-orbit GEDI L4A
mean(rsqd_orbit)
mean(bias_orbit)
mean(RMSE_orbit)


p2 <- ggplot(testing_dataset, aes(x=Combo, y=AGBD))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=0,y=87,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_orbit),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(mean(RMSE_orbit),3)," Mg/ha",
             "\n" , " Bias: ", round(mean(bias_orbit),3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(a)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=paste0("AGBD_",suf, " (Mg/ha)"), 
       y=expression("AGBD"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

df_bias_group_orbit<- do.call(rbind, bias_group_orbit)
df_bias_group_orbit <- df_bias_group_orbit[-1,]
colMeans(df_bias_group_orbit, na.rm = FALSE, dims = 1)
df_rbias_group_orbit <- do.call(rbind, rbias_group_orbit)
df_rbias_group_orbit <- df_rbias_group_orbit[-1,]
colMeans(df_rbias_group_orbit, na.rm = FALSE, dims = 1)
testing_dataset$biases <- testing_dataset$AGBD - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)
table(testing_dataset$group)/nrow(testing_dataset)



b2<- ggplot(testing_dataset, aes(x=group, y=biases, group = group)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  annotate("text",x=10,y=50,hjust = 0,size = 15,family= "A", label= "(a)") + 
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-60, 60),xlim = c(0, 12))+
  scale_y_continuous(minor_breaks = round(seq(-60,60,20),digits = 1),
                     breaks = round(seq(-60,60,20),digits = 1))+
  theme_bw()+
  labs(x=paste0("AGBD_",suf, " (Mg/ha)"), 
       y=expression("Bias of AGBD"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


#original
ggarrange(p2,p1,ncol=2)
out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018\\AGB",paste0(suf,"_scatterplot.jpg"))
ggsave(out,height=9, width=18, dpi=600)

ggarrange(b2,b1,ncol=2)
out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018\\AGB",paste0(suf,"_boxplot.jpg"))
ggsave(out,height=9, width=27, dpi=600)





# ------------------------------------------------------------------------------------------------ #
#histogram of height and biomass
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 0.3502*Data$MEAN^2.5405 

Data <- Data[Data$AGBD > 0 & Data$AGBD < 100,]
Data <- Data[Data$Combo > 0 & Data$Combo <= 100,]



Data$group_Combo <- cut(Data$Combo,
                        breaks = c(0,5,10,15,20,25,30,35,40,45,50,55,60,200),
                        dig.lab=3)

summary(Data$Combo)

p1<- ggplot(data=Data) + 
  theme_bw()+
  annotate("text",x=10,y=50,hjust = 0,size = 15,family= "A", label= "(a)") + 
  theme(text=element_text(size=25),axis.text=element_text(size=25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_histogram(aes(x=group_Combo,y = stat(count) / sum(count)*100),
                 color =  "blue",fill = "red",stat="count",alpha = 0.8) +
  coord_cartesian(ylim =  c(0, 70))+
  labs(x=expression("AGB"["ALS_area_MCH_2018"]~"(Mg/ha)"))+
  labs(y="Frequency (%)")+
  scale_x_discrete(labels=c("0-5","5-10","10-15","15-20","20-25",
                            "25-30","30-35","35-40","40-45","45-50",
                            "50-55","55-60",">60"))+
  scale_y_continuous(minor_breaks = round(seq(0,70,10),digits = 1),
                     breaks = round(seq(0,70,10),digits = 1))+
  theme(text=element_text(family="A"))



dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

Data <- Data[Data$AGBD > 0 & Data$AGBD < 100,]
Data <- Data[Data$Combo > 0 & Data$Combo <= 60,]

Data$group_Combo <- cut(Data$Combo,
                        breaks = c(0,5,10,15,20,25,30,35,40,45,50,55,60,200),
                        dig.lab=3)

summary(Data$Combo)

p2<- ggplot(data=Data) + 
  theme_bw()+
  annotate("text",x=10,y=50,hjust = 0,size = 15,family= "A", label= "(b)") + 
  theme(text=element_text(size=25),axis.text=element_text(size=25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_histogram(aes(x=group_Combo,y = stat(count) / sum(count)*100),
                 color =  "blue",fill = "red",stat="count",alpha = 0.8) +
  coord_cartesian(ylim =  c(0, 40))+
  labs(x=expression("AGB"["ALS_area_MCH_CC_2018"]~"(Mg/ha)"))+
  labs(y="Frequency (%)")+
  scale_x_discrete(labels=c("0-5","5-10","10-15","15-20","20-25",
                            "25-30","30-35","35-40","40-45","45-50",
                            "50-55","55-60",">60"))+
  scale_y_continuous(minor_breaks = round(seq(0,70,10),digits = 1),
                     breaks = round(seq(0,70,10),digits = 1))+
  theme(text=element_text(family="A"))


dir <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012\\Agincourt_Colgan_2012.tif"
r <- raster(dir)
r[r > 100] <- NA
df1 <- data.frame(rasterToPoints(r))

dir <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012\\Ireagh_Colgan_2012.tif"
r <- raster(dir)
r[r > 100] <- NA
df2 <- data.frame(rasterToPoints(r))

dir <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012\\Justicia_Colgan_2012.tif"
r <- raster(dir)
r[r > 100] <- NA
df3 <- data.frame(rasterToPoints(r))

dir <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012\\Welverdiendt_Colgan_2012.tif"
r <- raster(dir)
r[r > 100] <- NA
df4 <- data.frame(rasterToPoints(r))

new_column_names <- c("x","y", "agb")
df <- rbind(df1,df2,df3,df4)
colnames(df) <- new_column_names

df <- df[df$agb > 0 & df$agb <= 60,]

df$group_agb <- cut(df$agb,
                        breaks = c(0,5,10,15,20,25,30,35,40,45,50,55,60,200),
                        dig.lab=3)
year <- basename(dir)

summary(df$agb)


p3<-ggplot(data=df) + 
  theme_bw()+
  annotate("text",x=10,y=50,hjust = 0,size = 15,family= "A", label= "(b)") + 
  theme(text=element_text(size=25),axis.text=element_text(size=25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_histogram(aes(x=group_agb,y = stat(count) / sum(count)*100),
                 color =  "blue",fill = "red",stat="count",alpha = 0.8) +
  coord_cartesian(ylim =  c(0, 40))+
  labs(x=expression("AGB"["ALS_area_MCH_CC_2012"]~"(Mg/ha)"))+
  labs(y="Frequency (%)")+
  scale_x_discrete(labels=c("0-5","5-10","10-15","15-20","20-25",
                            "25-30","30-35","35-40","40-45","45-50",
                            "50-55","55-60",">60"))+
  scale_y_continuous(minor_breaks = round(seq(0,70,10),digits = 1),
                     breaks = round(seq(0,70,10),digits = 1))+
  theme(text=element_text(family="A"))


ggarrange(p1,p2,p3,ncol=3)
out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\histogram_als_agb.jpg"
ggsave(out,height=15, width=40, dpi=600)


ggarrange(p2,p3,ncol=2)
out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\histogram_als_agb_2.jpg"
ggsave(out,height=20, width=40, dpi=600)



#------------------------------------------------------------------------------------------------#
#06302023
#SAS applied to all GEDI footprints -- glm
#Local updated model
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

suf <- "H_CC_linear_0_glm"
#Data$Combo <- 10.35*Data$MEAN*Data$SUM/490 - 5.9236
Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

#Data <- Data[Data$Combo > 0,]
Data <- Data[Data$AGBD > 0,]
Data <- Data[Data$SUM > 24,]

Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

#mean(Data$Combo)
#Data <- Data %>% group_by(gr) %>% slice_sample(n=100)
#mean(Data$Combo)

Data = Data[!duplicated(Data$shot_numbe_x),]
Data_50 <-Data[Data$AGBD < 40,]
Data_50_100 <- Data[Data$AGBD >= 40 & Data$AGBD < 100,]
random_sample <- createDataPartition(Data_50$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[-random_sample, ]
training_dataset <- rbind(training_dataset, Data_50_100)
testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset = testing_dataset[!duplicated(testing_dataset$shot_numbe_x),]
training_dataset = training_dataset[!duplicated(training_dataset$shot_numbe_x),]

#Local updated 
reg <- glm(formula = Combo ~ GEDI_PAVD+RH1_50+RH1_98, 
           data = training_dataset,
           family = Gamma(link = "identity"),
           start = c(0.5,0.5,102,5))
testing_dataset$preds <- predict(reg,newdata=testing_dataset)

rsqd_local <- round(cor(testing_dataset$Combo,testing_dataset$preds)^2,3)
bias_local <- mean(testing_dataset$preds - testing_dataset$Combo)
RMSE_local <- sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))
rRMSE_local <- 100*sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)


testing_dataset$biases = testing_dataset$preds - testing_dataset$Combo
breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)
bias_group_local <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_local <- aggregate(testing_dataset$biases, 
                               list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

#On-orbit GEDI L4A
rsqd_orbit = round(cor(testing_dataset$AGBD, testing_dataset$Combo,method = "pearson")^2,3)
RMSE_orbit <- sqrt(mean((testing_dataset$AGBD - testing_dataset$Combo)^2))
bias_orbit <- mean(testing_dataset$AGBD - testing_dataset$Combo)
rRMSE_orbit <- 100*sqrt(mean((testing_dataset$AGBD - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)

testing_dataset$biases = testing_dataset$AGBD - testing_dataset$Combo
bias_group_orbit <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_orbit <- aggregate(testing_dataset$biases, 
                               list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

bias_group_orbit
rbias_group_orbit

rsqd_local
bias_local
RMSE_local
n <- 1654
p1<- ggplot(testing_dataset, aes(x=Combo, y=preds))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_local),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE_local,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE_local,3)," %",
             "\n" , " Bias: ", round(bias_local,3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(b)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["SAS_GLM"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_local
rbias_group_local
testing_dataset$biases <- testing_dataset$preds - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_local<- round(bias_group_local[0:10],1)
c_local<- round(rbias_group_local[0:10],1)
c_local <- paste0(c_local, "%")
d_local <- a_local
d_local[d_local>0]=0
num <- table(testing_dataset$group)[0:10]

l <- cbind(x,a_local,num)
df <- data.frame(l)

b1<- ggplot(df, aes(x=x+2.5, y=a_local,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_local,y=d_local-2), size = 6)+
  geom_text(aes(label = c_local,y=d_local-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15,family= "A", label= "(b)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("Bias of AGBD"["SAS_GLM"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



#on-orbit GEDI L4A
rsqd_orbit
bias_orbit
RMSE_orbit
n <- 1654
p2 <- ggplot(testing_dataset, aes(x=Combo, y=AGBD))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_orbit),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE_orbit,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE_orbit,3)," %",
             "\n" , " Bias: ", round(bias_orbit,3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(a)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_orbit
rbias_group_orbit
testing_dataset$biases <- testing_dataset$AGBD - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_orbit<- round(bias_group_orbit[0:10],1)
c_orbit<- round(rbias_group_orbit[0:10],1)
c_orbit <- paste0(c_orbit, "%")
d_orbit <- a_orbit
d_orbit[d_orbit>0]=0
#num <- table(testing_dataset$group)[0:10]

l <- cbind(x,a_orbit,num)
df <- data.frame(l)


b2<- ggplot(df, aes(x=x+2.5, y=a_orbit,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_orbit,y=d_orbit-2), size = 6)+
  geom_text(aes(label = c_orbit,y=d_orbit-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15,family= "A", label= "(a)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("Bias of AGBD"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



# ------------------------------------------------------------------------------------------------ #
#rf
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

suf <- "H_CC_linear_0_RF"
#Data$Combo <- 10.35*Data$MEAN*Data$SUM/490 - 5.9236
Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

#Data <- Data[Data$Combo > 0,]
Data <- Data[Data$AGBD >= 0,]
#Data <- Data[Data$SUM > 24,]

Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

Data_sel <- Data[,c("RH1_50","RH1_75","RH1_98",
                    "GEDI_Cover","GEDI_FHD")]
Data_sel <- na.omit(Data_sel)

#PCA_result <- princomp(Data_sel)
#summary(PCA_result)
#fviz_eig(PCA_result, addlabels = TRUE)
#PCA_result$loadings[, 1:2]
#fviz_pca_var(PCA_result, col.var = "black")
#fviz_cos2(PCA_result, choice = "var", axes = 1:2)
#Data_sel<- Data[,c("Combo")]
#pcs <- as.data.frame(PCA_result$scores[,1:2])
#Data <- cbind(Data, pcs)

#Data <- Data %>% group_by(gr) %>% slice_sample(n=100)

Data = Data[!duplicated(Data$shot_numbe_x),]
Data_50 <-Data[Data$AGBD < 40,]
Data_50_100 <- Data[Data$AGBD >= 40 & Data$AGBD < 100,]
random_sample <- createDataPartition(Data_50$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[-random_sample, ]
training_dataset <- rbind(training_dataset, Data_50_100)
testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset = testing_dataset[!duplicated(testing_dataset$shot_numbe_x),]
training_dataset = training_dataset[!duplicated(training_dataset$shot_numbe_x),]

#Local updated 
repeat_cv <- trainControl(method='cv', number=5)
reg <- train(
  Combo~sensitivity_a1+
    GEDI_Cover + GEDI_FHD + 
    RH1_50 + RH1_75 + RH1_90 + RH1_98,
  data=training_dataset, 
  method='rf', 
  ntree=100,
  nodesize = 4,
  trControl=repeat_cv)

testing_dataset$preds <- predict(object=reg, 
                                 newdata=testing_dataset[, ])

rsqd_local2 <- round(cor(testing_dataset$Combo,testing_dataset$preds)^2,3)
bias_local2 <- mean(testing_dataset$preds - testing_dataset$Combo)
RMSE_local2 <- sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))
rRMSE_local2 <- 100*sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)
n <- 1654
mean(testing_dataset$Combo)
testing_dataset$biases = testing_dataset$preds - testing_dataset$Combo
breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)
bias_group_local2 <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_local2 <- aggregate(testing_dataset$biases, 
                               list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

p4<- ggplot(testing_dataset, aes(x=Combo, y=preds))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(rsqd_local2,3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE_local2,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE_local2,3),"%",
             "\n" , " Bias: ", round(bias_local2,3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(c)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["SAS_RF"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))

#suf = "RF_PCA_"
#out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_scatterplot.jpg"))
#ggsave(out,height=9, width=9, dpi=600)

breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_local2
rbias_group_local2
testing_dataset$biases <- testing_dataset$preds - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_local2<- round(bias_group_local2[0:10],1)
c_local2<- round(rbias_group_local2[0:10],1)
c_local2 <- paste0(c_local2, "%")
d_local2 <- a_local2
d_local2[d_local2>0]=0
#num <- table(testing_dataset$group)[0:10]

l <- cbind(x,a_local2,num)
df <- data.frame(l)

b4<- ggplot(df, aes(x=x+2.5, y=a_local2, group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_local2,y=d_local2-2), size = 6)+
  geom_text(aes(label = c_local2,y=d_local2-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15,family= "A", label= "(c)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("Bias of AGBD"["SAS_RF"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



#original
suf = "model4_"
ggarrange(p2,p1,p4,nrow = 1, ncol = 3)
out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_scatterplot.jpg"))
ggsave(out,height=10, width=30, dpi=600)

ggarrange(b2,b1,b4,nrow = 1, ncol = 3)
out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_boxplot.jpg"))
ggsave(out,height=10, width=30, dpi=600)





# ------------------------------------------------------------------------------------------------ #
#07282023

#0 intercept OLS model (c)
i <- "E:\\Biomass\\CSIR\\Biomass_2018\\CSIR_25m.xlsx"
Data = read_excel(i)

reg <- lm(agbd_ha~HCC + 0, data=Data)
summary(reg)

x<-Data$MEAN05
y<-Data$agbd_ha

mult_nls <- nls(y ~ a*exp(b*x), start = list(a=100, b=1))
coef(mult_nls)[1]
coef(mult_nls)[2]


Data$preds <- Data$HCC*coef(reg)[1]

rsqd_local <- round(cor(Data$agbd_ha,Data$preds)^2,3)
bias_local <- mean(Data$preds - Data$agbd_ha)
RMSE_local <- sqrt(mean((Data$preds - Data$agbd_ha)^2))



#nls power H only (b)
i <- "E:\\Biomass\\CSIR\\Biomass_2018\\CSIR_25m.xlsx"
Data = read_excel(i)
Data_50 <-Data[Data$CC15 < 1,]
Data_50 <-Data[Data$agbd_ha < 40,]
Data_50_100 <- Data[Data$agbd_ha >= 40 & Data$agbd_ha < 200,]
random_sample <- createDataPartition(Data_50$agbd_ha, p = 0.7, list = FALSE)
random_sample1 <- createDataPartition(Data$agbd_ha, p = 0.3, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[random_sample1, ]
training_dataset <- rbind(training_dataset, Data_50_100)
#testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset <- na.omit(testing_dataset)
training_dataset <- na.omit(training_dataset)

reg <- nls(agbd_ha~ a*MEAN05^b, start = list(a=1, b=1),data=training_dataset)
testing_dataset$preds <- coef(reg)[1]*testing_dataset$MEAN05^coef(reg)[2]

rsqd_local <- round(cor(testing_dataset$agbd_ha,testing_dataset$preds)^2,3)
bias_local <- mean(testing_dataset$preds - testing_dataset$agbd_ha)
RMSE_local <- sqrt(mean((testing_dataset$preds - testing_dataset$agbd_ha)^2))







# ------------------------------------------------------------------------------------------------ #
#scatterplot 25m csir vs. area based and ind biomass
filedir <- "E:\\Biomass\\CSIR\\result\\CSIR_25m\\biom_CSIR_25m_combo_1_h_cc.csv"
Data = read.csv(filedir)
Data <- na.omit(Data)
Data$ALS <-  Data$Ind
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np13"),]
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np15"),]
Data <- Data[Data$CSIR >0 ,]
#stats
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
MD <- mean(Data$ALS - Data$CSIR)

p1<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  annotate("text",x=0,y=118,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=110,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  annotate("text",x=100,y=120,hjust = 0,size = 15,family= "A", label= "(a)") + 
  coord_cartesian(xlim = c(0, 120),ylim =  c(0, 120))+
  labs(x=expression("AGB"["field_25m"]~"(Mg/ha)"), 
       y=expression("25 m AGB"["ALS_tree"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# ------------------------------------------------------------------------------------------------ #
#updated biomass model -- area-based H test
filedir <- "E:\\Biomass\\CSIR\\Biomass_2018\\CSIR_25m.xlsx"
Data = read_excel(filedir)
Data <- na.omit(Data)

Data<-Data[!(Data$plot=="Ven_S63_SP_5_np13"),]
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np15"),]
Data<-Data[!(Data$plot=="Ven_S68_SP_5_np10"),]
Data <- Data[Data$agbd_ha >0,]

i <- 1
rsqd_ <- 0
bias_ <- 0
RMSE_ <- 0
a_ <- 0 
b_ <- 0
for (i in 1:100){
  random_sample <- createDataPartition(Data$agbd_ha, p = 0.7, list = FALSE)
  training_dataset  <- Data[random_sample, ]
  testing_dataset <- Data[-random_sample, ]
  x <- training_dataset$MEAN05
  y <- training_dataset$agbd_ha
  reg <- nls(y~ a*x^b, start = list(a=1, b=1))
  testing_dataset$preds <- coef(reg)[1]*testing_dataset$MEAN05^coef(reg)[2]
  a_[i] <- coef(reg)[1] 
  b_[i] <- coef(reg)[2]
  rsqd_[i] <- round(cor(testing_dataset$agbd_ha,testing_dataset$preds)^2,3)
  bias_[i] <- mean(testing_dataset$preds - testing_dataset$agbd_ha)
  RMSE_[i] <- sqrt(mean((testing_dataset$preds - testing_dataset$agbd_ha)^2))
}

r2 <- mean(rsqd_)
MD <- mean(bias_)
RMSE <- mean(RMSE_)

data <- data.frame(cbind(rsqd_,a_,b_))
row_sel_rsqd <- data[data$rsqd_ == max(data$rsqd_),]
row_sel_rsqd


# ------------------------------------------------------------------------------------------------ #
#updated biomass model -- area-based H*CC test
filedir <- "E:\\Biomass\\CSIR\\Biomass_2018\\CSIR_25m.xlsx"
Data = read_excel(filedir)
Data <- na.omit(Data)

Data<-Data[!(Data$plot=="Ven_S63_SP_5_np13"),]
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np15"),]
Data<-Data[!(Data$plot=="Ven_S68_SP_5_np10"),]
Data <- Data[Data$agbd_ha >0,]

i <- 1
rsqd_ <- 0
bias_ <- 0
RMSE_ <- 0
a_ <- 0 
#b_ <- 0
for (i in 1:100){
  random_sample <- createDataPartition(Data$agbd_ha, p = 0.7, list = FALSE)
  training_dataset  <- Data[random_sample, ]
  testing_dataset <- Data[-random_sample, ]
  reg <- lm(agbd_ha~ HCC + 0,data=training_dataset)
  testing_dataset$preds <- coef(reg)[1]*testing_dataset$HCC
  a_[i] <- coef(reg)[1] 
  #b_[i] <- coef(reg)[2]
  rsqd_[i] <- round(cor(testing_dataset$agbd_ha,testing_dataset$preds)^2,3)
  bias_[i] <- mean(testing_dataset$preds - testing_dataset$agbd_ha)
  RMSE_[i] <- sqrt(mean((testing_dataset$preds - testing_dataset$agbd_ha)^2))
}

r2 <- mean(rsqd_)
MD <- mean(bias_)
RMSE <- mean(RMSE_)

data <- data.frame(cbind(rsqd_,a_))
row_sel_rsqd <- data[data$rsqd_ == max(data$rsqd_),]
row_sel_rsqd





# ------------------------------------------------------------------------------------------------ #
#histogram of GEDI L4A
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)

Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490


Data <- Data[Data$Combo > 0 & Data$Combo < 100,]

Data <- na.omit(Data)

Data$group_Combo <- cut(Data$Combo,
                        breaks = c(0,5,10,15,20,25,30,35,40,45,50,55,60),
                        dig.lab=3)




dir <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012\\Agincourt_Colgan_2012.tif"
r <- raster(dir)
r[r > 100] <- NA
df1 <- data.frame(rasterToPoints(r))

dir <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012\\Ireagh_Colgan_2012.tif"
r <- raster(dir)
r[r > 100] <- NA
df2 <- data.frame(rasterToPoints(r))

dir <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012\\Justicia_Colgan_2012.tif"
r <- raster(dir)
r[r > 100] <- NA
df3 <- data.frame(rasterToPoints(r))

dir <- "E:\\Biomass\\CSIR\\Result_0616\\SAS_Colgan_2012\\Welverdiendt_Colgan_2012.tif"
r <- raster(dir)
r[r > 100] <- NA
df4 <- data.frame(rasterToPoints(r))

new_column_names <- c("x","y", "agb")
df <- rbind(df1,df2,df3,df4)
colnames(df) <- new_column_names

df <- df[df$agb > 0 & df$agb <= 60,]

df$group_agb <- cut(df$agb,
                    breaks = c(0,5,10,15,20,25,30,35,40,45,50,55,60),
                    dig.lab=3)



df_merge <- cbind(table(Data$group_Combo)/nrow(Data)*100,table(df$group_agb)/nrow(df)*100)
df_merge<- data.frame(df_merge)
df_merge$bins <- rownames(df_merge)
colnames(df_merge) <- c("HCC_2018", "HCC_2012", "Bins")
df_merge <- melt(df_merge,id = "Bins")
colnames(df_merge) <- c("Bins", "AGBD", "freq")

ggplot(data=df_merge) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_bar(aes(x=reorder(Bins, -freq), y=freq, fill = AGBD), 
           stat = "identity", position = "dodge") +
  coord_cartesian(ylim =  c(0, 40))+
  xlab("ALS-based AGB (MCH*CC) (Mg/ha)") + 
  ylab("Frequency (%)")+
  theme(legend.position=c(.8,.3))+
  theme(text=element_text(family="A"))+
  theme(text=element_text(size=25),axis.text=element_text(size=25))



out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\GEDI_L4A_ALS_histogram_2.jpg"
ggsave(out,height=10, width=20, dpi=600)




# ------------------------------------------------------------------------------------------------ #
#histogram of height and biomass
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 0.3502*Data$MEAN^2.5405 

Data <- Data[Data$AGBD > 0 & Data$AGBD < 100,]
Data <- Data[Data$Combo > 0 & Data$Combo <= 100,]



Data$group_Combo <- cut(Data$Combo,
                        breaks = c(0,5,10,15,20,25,30,35,40,45,50),
                        dig.lab=3)



Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

Data <- Data[Data$AGBD > 0 & Data$AGBD < 100,]
Data <- Data[Data$Combo > 0 & Data$Combo <= 100,]

Data$group_Combo_HCC <- cut(Data$Combo,
                        breaks = c(0,5,10,15,20,25,30,35,40,45,50),
                        dig.lab=3)

df_merge <- cbind(table(Data$group_Combo)/nrow(Data)*100,table(Data$group_Combo_HCC)/nrow(Data)*100)
df_merge<- data.frame(df_merge)
df_merge$bins <- rownames(df_merge)
colnames(df_merge) <- c("ALS_based_H", "ALS_based_HCC", "Bins")
df_merge <- melt(df_merge,id = "Bins")
colnames(df_merge) <- c("Bins", "AGBD", "freq")

ggplot(data=df_merge) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_bar(aes(x=reorder(Bins, -freq), y=freq, fill = AGBD), 
           stat = "identity", position = "dodge") +
  coord_cartesian(ylim =  c(0, 70))+
  xlab("ALS-based AGBD (Mg/ha)") + 
  ylab("Frequency (%)")+
  theme(legend.position=c(.8,.3))+
  theme(text=element_text(family="A"))+
  theme(text=element_text(size=25),axis.text=element_text(size=25))



out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\histogram_als_agb_2.jpg"
ggsave(out,height=10, width=20, dpi=600)






dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo1 <- 0.3502*Data$MEAN^2.5405 
Data <- Data[Data$Combo1 > 0 & Data$Combo1 <= 100,]


Data$Combo2 <- 9.0665*Data$MEAN*Data$SUM/490
Data <- Data[Data$Combo2 > 0 & Data$Combo2 <= 100,]

ggplot(data=Data) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_density(aes(x=Combo1,fill="blue"), alpha=0.3,colour="red") +
  geom_density(aes(x=Combo2,fill="red"), alpha=0.3,colour="blue") +
  coord_cartesian(xlim =  c(0, 60),ylim =  c(0, 0.2))+
  xlab("ALS-based AGBD (Mg/ha)") + 
  ylab("Frequency (%)")+
  theme(legend.position=c(.7,.5))+
  theme(text=element_text(family="A"))+
  theme(text=element_text(size=30),axis.text=element_text(size=30))+
  scale_fill_discrete(name = "Legend",labels = c("ALS_based_H", "ALS_based_HCC"))


out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\density_als.jpg"
ggsave(out,height=10, width=10, dpi=600)



ggplot(data=Data) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_density(aes(x=AGBD,fill="blue"), alpha=0.3,colour="red") +
  geom_density(aes(x=Combo2,fill="red"), alpha=0.3,colour="blue") +
  coord_cartesian(xlim =  c(0, 60),ylim =  c(0, 0.2))+
  xlab("AGBD (Mg/ha)") + 
  ylab("Frequency (%)")+
  theme(legend.position=c(.7,.5))+
  theme(text=element_text(family="A"))+
  theme(text=element_text(size=30),axis.text=element_text(size=30))+
  scale_fill_discrete(name = "Legend",labels = c("GEDI", "ALS"))


out = "E:\\Biomass\\CSIR\\Result_0618\\2018\\density_als_gedi.jpg"
ggsave(out,height=10, width=10, dpi=600)




# ------------------------------------------------------------------------------------------------ #
#manually scatterplot and bias plot
#08222023
#SAS applied to all GEDI footprints -- glm
#Local updated model
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

suf <- "H_CC_linear_0_glm"
#Data$Combo <- 10.35*Data$MEAN*Data$SUM/490 - 5.9236
Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

#Data <- Data[Data$Combo > 0,]
Data <- Data[Data$AGBD > 0,]
Data <- Data[Data$SUM > 24,]

Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

#write.csv(Data, "C:\\Users\\Shawn\\Desktop\\Misha_0829.csv", row.names=FALSE)
#mean(Data$Combo)
#Data <- Data %>% group_by(gr) %>% slice_sample(n=100)
#mean(Data$Combo)

Data = Data[!duplicated(Data$shot_numbe_x),]
Data_50 <-Data[Data$AGBD < 40,]
Data_50_100 <- Data[Data$AGBD >= 40 & Data$AGBD < 100,]
random_sample <- createDataPartition(Data_50$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[-random_sample, ]
training_dataset <- rbind(training_dataset, Data_50_100)
testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset = testing_dataset[!duplicated(testing_dataset$shot_numbe_x),]
training_dataset = training_dataset[!duplicated(training_dataset$shot_numbe_x),]

#Local updated 
repeat_cv <- trainControl(method='cv', number=5)
training_dataset <-training_dataset[training_dataset$Combo > 0,]
reg <- train(
  data=training_dataset, 
  method='glm', 
  family = Gamma(link = "identity"),
  start=c(-600,-2,6,5),
  trControl=repeat_cv)

testing_dataset$preds <- predict(reg,newdata=testing_dataset)

rsqd_local <- round(cor(testing_dataset$Combo,testing_dataset$preds)^2,3)
bias_local <- mean(testing_dataset$preds - testing_dataset$Combo)
RMSE_local <- sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))
rRMSE_local <- 100*sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)


testing_dataset$biases = testing_dataset$preds - testing_dataset$Combo
breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)
bias_group_local <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_local <- aggregate(testing_dataset$biases, 
                               list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

#On-orbit GEDI L4A
rsqd_orbit = round(cor(testing_dataset$AGBD, testing_dataset$Combo,method = "pearson")^2,3)
RMSE_orbit <- sqrt(mean((testing_dataset$AGBD - testing_dataset$Combo)^2))
bias_orbit <- mean(testing_dataset$AGBD - testing_dataset$Combo)
rRMSE_orbit <- 100*sqrt(mean((testing_dataset$AGBD - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)

testing_dataset$biases = testing_dataset$AGBD - testing_dataset$Combo
bias_group_orbit <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_orbit <- aggregate(testing_dataset$biases, 
                               list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

#DBT Africa GEDI L4A
testing_dataset$AGBD_dbt <- 1.092*(-118.408 + 1.957*sqrt(testing_dataset$RH1_50) 
                                   + 9.962 * sqrt(testing_dataset$RH1_98 + 100))^2
rsqd_dbt = round(cor(testing_dataset$AGBD_dbt, testing_dataset$Combo,method = "pearson")^2,3)
RMSE_dbt <- sqrt(mean((testing_dataset$AGBD_dbt - testing_dataset$Combo)^2))
bias_dbt <- mean(testing_dataset$AGBD_dbt - testing_dataset$Combo)
rRMSE_dbt <- 100*sqrt(mean((testing_dataset$AGBD_dbt - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)


testing_dataset$biases = testing_dataset$AGBD_dbt - testing_dataset$Combo
bias_group_dbt <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_dbt <- aggregate(testing_dataset$biases, 
                             list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

n<- nrow(testing_dataset)

rsqd_local <- 0.536
bias_local
RMSE_local

p1<- ggplot(testing_dataset, aes(x=Combo, y=preds))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_local),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: 10.2 Mg/ha",
             "\n" , " rRMSE: 67.2 %",
             "\n" , " MSD: 1.39 Mg/ha",
             "\n" , " Sample size: 1654")) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(b)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["SAS_GLM"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))


breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_local
rbias_group_local
testing_dataset$biases <- testing_dataset$preds - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_local<- c(5.2,3.3,1,-0.8,-2.6,-6.6,-7.2,-9.9,-11.4,-16.2)
c_local<- c(197.5,44.4,8.1,-4.5,-11.6,-23.9,-22.3,-26.5,-27,-34.2)
c_local <- paste0(c_local, "%")
d_local <- a_local
d_local[d_local>0]=0
num <- c(409,334,271,180,130,101,73,54,30,23)

l <- cbind(x,a_local,num)
df <- data.frame(l)

b1<- ggplot(df, aes(x=x+2.5, y=a_local,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_local,y=d_local-2), size = 6)+
  geom_text(aes(label = c_local,y=d_local-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15, label= "(b)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("Bias of AGBD"["SAS_GLM"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



#on-orbit GEDI L4A
rsqd_orbit <-0.428
bias_orbit
RMSE_orbit

p2 <- ggplot(testing_dataset, aes(x=Combo, y=AGBD))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_orbit),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: 12.07 Mg/ha",
             "\n" , " rRMSE: 79.5 %",
             "\n" , " MSD: -3.87 Mg/ha",
             "\n" , " Sample size: 1654")) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(a)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_orbit
rbias_group_orbit
testing_dataset$biases <- testing_dataset$AGBD - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_orbit<- c(0.7,-0.9,-3.7,-5.9,-7.6,-10.4,-10.5,-9.9,-11.6,-16.5)
c_orbit<- c(27.9,-12.1,-29.9,-33.8,-34.2,-38.3,-32.2,-26.3,-27.6,-34.6)
c_orbit <- paste0(c_orbit, "%")
d_orbit <- a_orbit
d_orbit[d_orbit>0]=0
num <- c(409,334,271,180,130,101,73,54,30,23)

l <- cbind(x,a_orbit,num)
df <- data.frame(l)


b2<- ggplot(df, aes(x=x+2.5, y=a_orbit,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_orbit,y=d_orbit-2), size = 6)+
  geom_text(aes(label = c_orbit,y=d_orbit-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15, label= "(a)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("Bias of AGBD"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

suf <- "H_CC_linear_0_RF"
#Data$Combo <- 10.35*Data$MEAN*Data$SUM/490 - 5.9236
Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

#Data <- Data[Data$Combo > 0,]
Data <- Data[Data$AGBD >= 0,]
#Data <- Data[Data$SUM > 24,]

Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

Data_sel <- Data[,c("RH1_50","RH1_75","RH1_90","RH1_98",
                    "GEDI_Cover","GEDI_FHD","GEDI_PAI","GEDI_PAVD")]
Data_sel <- na.omit(Data_sel)

#PCA_result <- princomp(Data_sel)
#summary(PCA_result)
#fviz_eig(PCA_result, addlabels = TRUE)
#PCA_result$loadings[, 1:2]
#fviz_pca_var(PCA_result, col.var = "black")
#fviz_cos2(PCA_result, choice = "var", axes = 1:2)

#Data_sel<- Data[,c("Combo")]
#pcs <- as.data.frame(PCA_result$scores[,1:2])
#Data <- cbind(Data, pcs)


#Data <- Data %>% group_by(gr) %>% slice_sample(n=100)

Data = Data[!duplicated(Data$shot_numbe_x),]
Data_50 <-Data[Data$AGBD < 40,]
Data_50_100 <- Data[Data$AGBD >= 40 & Data$AGBD < 100,]
random_sample <- createDataPartition(Data_50$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[-random_sample, ]
training_dataset <- rbind(training_dataset, Data_50_100)
testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset = testing_dataset[!duplicated(testing_dataset$shot_numbe_x),]
training_dataset = training_dataset[!duplicated(training_dataset$shot_numbe_x),]

#Local updated 
repeat_cv <- trainControl(method='cv', number=5)
reg <- train(
  Combo~sensitivity_a1+
    GEDI_PAVD + GEDI_Cover + GEDI_PAI + GEDI_FHD + 
    RH1_50 + RH1_75 + RH1_90 + RH1_98,
  data=training_dataset, 
  method='rf', 
  ntree=100,
  nodesize = 4,
  trControl=repeat_cv)

rf <- randomForest(Combo~sensitivity_a1+
                     GEDI_Cover + GEDI_FHD + 
                     RH1_50 + RH1_75 + RH1_90 + RH1_98,
                     data=training_dataset,ntree = 100,
                   importance=TRUE, nodesize = 4)
importance(rf)
varImpPlot(rf)
varImp(rf)


testing_dataset$preds <- predict(object=reg, 
                                 newdata=testing_dataset[, ])

rsqd_local2 <- 0.712
bias_local2 <- mean(testing_dataset$preds - testing_dataset$Combo)
RMSE_local2 <- sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))
rRMSE_local2 <- 100*sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)
mean(testing_dataset$Combo)
testing_dataset$biases = testing_dataset$preds - testing_dataset$Combo
breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)
bias_group_local2 <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_local2 <- aggregate(testing_dataset$biases, 
                                list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100


p4<- ggplot(testing_dataset, aes(x=Combo, y=preds))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_local2),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: 7.3 Mg/ha",
             "\n" , " rRMSE: 53.3 %",
             "\n" , " MSD: 0.06 Mg/ha",
             "\n" , " Sample size: 1654")) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(c)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["SAS_RF"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))

#suf = "RF_PCA_"
#out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_scatterplot.jpg"))
#ggsave(out,height=9, width=9, dpi=600)

breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_local2
rbias_group_local2
testing_dataset$biases <- testing_dataset$preds - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_local2<- c(3.1,2.65,1.53,-0.9,-1.94,-4.19,-6.66,-8.77,-7.49,-11.7)
c_local2<- c(156.4,35.8,12.5,-5.2,-8.6,-15.1,-20.4,-23.5,-17.7,-24.8)
c_local2 <- paste0(c_local2, "%")
d_local2 <- a_local2
d_local2[d_local2>0]=0
num <- c(409,334,271,180,130,101,73,54,30,23)


l <- cbind(x,a_local2,num)
df <- data.frame(l)

b4<- ggplot(df, aes(x=x+2.5, y=a_local2, group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_local2,y=d_local2-2), size = 6)+
  geom_text(aes(label = c_local2,y=d_local2-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15, label= "(c)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("Bias of AGBD"["SAS_RF"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



#original
suf = "model4_"
ggarrange(p2,p1,p4,nrow = 1, ncol = 3)
out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_scatterplot.jpg"))
ggsave(out,height=10, width=30, dpi=600)

ggarrange(b2,b1,b4,nrow = 1, ncol = 3)
out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_boxplot.jpg"))
ggsave(out,height=10, width=30, dpi=600)





# ------------------------------------------------------------------------------------------------ #
#data setup 
i <- "E:\\Biomass\\CSIR\\Biomass_2018\\CSIR_25m.xlsx"
Data = read_excel(i)

height <- Data$MEAN05
cover <- Data$CC15


#Monte Carlo simulation for error propogation -- ALS - field modeling
biomass_estimate <- function(height,cover) {
  biomass <- 9.0665*height*cover
  return(biomass)
}

num_simulations <- 1000  # Number of simulation runs

# Initialize an empty data frame to store results
results <- data.frame()

for (i in 1:num_simulations) {
  # Randomly sample values from the distributions
  height_sample <- sample(height, 1)
  cover_sample <- sample(cover, 1)
  
  # Calculate biomass estimate for this simulation run
  biomass_result <- biomass_estimate(height_sample, cover_sample)
  
  # Store the results in the data frame
  results <- bind_rows(results, data.frame(Simulation = i, Biomass = biomass_result))
}

# Summary statistics
summary_stats <- results %>%
  summarise(Mean_Biomass = mean(Biomass),
            SD_Biomass = sd(Biomass),
            Lower_CI = quantile(Biomass, 0.025),
            Upper_CI = quantile(Biomass, 0.975))

# Histogram of biomass estimates
ggplot(results, aes(x = Biomass)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  labs(title = "Monte Carlo Biomass Estimates",
       x = "Biomass",
       y = "Frequency")

# Print summary statistics
print(summary_stats)





dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data <- Data[c("sensitivity_a1", "GEDI_PAVD", "GEDI_Cover", "GEDI_PAI",
               "GEDI_FHD", "RH1_50", "RH1_75", "RH1_90", "RH1_98")]
stats_cov <- cov(Data)
write.table(stats_cov,file='C:\\Users\\Shawn\\Desktop\\stats_cov.csv',sep = ",",col.names=TRUE)
write.table(Data,file='C:\\Users\\Shawn\\Desktop\\dataframe.csv',sep = ",",col.names=TRUE)






#------------------------------------------------------------------------------------------------#
#09252023
#SAS applied to all GEDI footprints -- glm
#Local updated model
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

suf <- "H_CC_linear_0_glm"
Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

Data <- Data[Data$AGBD > 0,]
Data <- Data[Data$SUM > 24,]

Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

Data = Data[!duplicated(Data$shot_numbe_x),]
Data_50 <-Data[Data$AGBD < 40,]
Data_50_100 <- Data[Data$AGBD >= 40 & Data$AGBD < 100,]
random_sample <- createDataPartition(Data_50$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[-random_sample, ]
training_dataset <- rbind(training_dataset, Data_50_100)
testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset = testing_dataset[!duplicated(testing_dataset$shot_numbe_x),]
training_dataset = training_dataset[!duplicated(training_dataset$shot_numbe_x),]

repeat_cv <- trainControl(method='cv', number=5)
reg <- train(
  Combo ~ GEDI_FHD + RH1_50 + RH1_98, 
  data=training_dataset, 
  method='glm', 
  family = Gamma(link = "identity"),
  start=c(0,0.1,100,5),
  trControl=repeat_cv)
testing_dataset$preds <- predict(object=reg, 
                                 newdata=testing_dataset[, ])

rsqd_local <- round(cor(testing_dataset$Combo,testing_dataset$preds)^2,3)
bias_local <- mean(testing_dataset$preds - testing_dataset$Combo)
RMSE_local <- sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))
rRMSE_local <- 100*sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)


testing_dataset$biases = testing_dataset$preds - testing_dataset$Combo
breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)
bias_group_local <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_local <- aggregate(testing_dataset$biases, 
                               list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

#On-orbit GEDI L4A
rsqd_orbit = round(cor(testing_dataset$AGBD, testing_dataset$Combo,method = "pearson")^2,3)
RMSE_orbit <- sqrt(mean((testing_dataset$AGBD - testing_dataset$Combo)^2))
bias_orbit <- mean(testing_dataset$AGBD - testing_dataset$Combo)
rRMSE_orbit <- 100*sqrt(mean((testing_dataset$AGBD - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)

testing_dataset$biases = testing_dataset$AGBD - testing_dataset$Combo
bias_group_orbit <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_orbit <- aggregate(testing_dataset$biases, 
                               list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

bias_group_orbit
rbias_group_orbit

rsqd_local
bias_local
RMSE_local
n <- 1654
p1<- ggplot(testing_dataset, aes(x=Combo, y=preds))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_local),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE_local,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE_local,3)," %",
             "\n" , " Bias: ", round(bias_local,3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(b)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["SAS_GLM"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_local
rbias_group_local
testing_dataset$biases <- testing_dataset$preds - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_local<- round(bias_group_local[0:10],1)
c_local<- round(rbias_group_local[0:10],1)
c_local <- paste0(c_local, "%")
d_local <- a_local
d_local[d_local>0]=0
num <- table(testing_dataset$group)[0:10]

l <- cbind(x,a_local,num)
df <- data.frame(l)

b1<- ggplot(df, aes(x=x+2.5, y=a_local,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_local,y=d_local-2), size = 6)+
  geom_text(aes(label = c_local,y=d_local-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15,family= "A", label= "(b)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("MSD of AGBD"["SAS_GLM"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))

#on-orbit GEDI L4A
rsqd_orbit
bias_orbit
RMSE_orbit
n <- 1654
p2 <- ggplot(testing_dataset, aes(x=Combo, y=AGBD))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(mean(rsqd_orbit),3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE_orbit,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE_orbit,3)," %",
             "\n" , " Bias: ", round(bias_orbit,3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(a)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_orbit
rbias_group_orbit
testing_dataset$biases <- testing_dataset$AGBD - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_orbit<- round(bias_group_orbit[0:10],1)
c_orbit<- round(rbias_group_orbit[0:10],1)
c_orbit <- paste0(c_orbit, "%")
d_orbit <- a_orbit
d_orbit[d_orbit>0]=0

l <- cbind(x,a_orbit,num)
df <- data.frame(l)


b2<- ggplot(df, aes(x=x+2.5, y=a_orbit,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_orbit,y=d_orbit-2), size = 6)+
  geom_text(aes(label = c_orbit,y=d_orbit-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15,family= "A", label= "(a)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("MSD of AGBD"["L4A"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



#rf
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

suf <- "H_CC_linear_0_RF"
Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

Data <- Data[Data$AGBD >= 0,]

Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

Data_sel <- Data[,c("RH1_50","RH1_75","RH1_98",
                    "GEDI_Cover","GEDI_FHD")]
Data_sel <- na.omit(Data_sel)

Data = Data[!duplicated(Data$shot_numbe_x),]
Data_50 <-Data[Data$AGBD < 40,]
Data_50_100 <- Data[Data$AGBD >= 40 & Data$AGBD < 100,]
random_sample <- createDataPartition(Data_50$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[-random_sample, ]
training_dataset <- rbind(training_dataset, Data_50_100)
testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset = testing_dataset[!duplicated(testing_dataset$shot_numbe_x),]
training_dataset = training_dataset[!duplicated(training_dataset$shot_numbe_x),]

#Local updated 
repeat_cv <- trainControl(method='cv', number=5)
reg <- train(
  Combo~sensitivity_a1+
    GEDI_Cover + GEDI_FHD + 
    RH1_50 + RH1_75 + RH1_90 + RH1_98,
  data=training_dataset, 
  method='rf', 
  ntree=100,
  nodesize = 4,
  trControl=repeat_cv)

rf <- randomForest(Combo ~ 
                     sensitivity_a1+
                     GEDI_FHD + GEDI_Cover +
                     RH1_50 + RH1_75 + RH1_90 + RH1_98, 
                   data = training_dataset,
                   ntree = 100,
                   importance=TRUE,
                   nodesize = 4)

testing_dataset$preds <- predict(object=reg, 
                                 newdata=testing_dataset[, ])
round(cor(testing_dataset$Combo,testing_dataset$preds)^2,3)

rsqd_local2 <- round(cor(testing_dataset$Combo,testing_dataset$preds)^2,3)
bias_local2 <- mean(testing_dataset$preds - testing_dataset$Combo)
RMSE_local2 <- sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))
rRMSE_local2 <- 100*sqrt(mean((testing_dataset$preds - testing_dataset$Combo)^2))/mean(testing_dataset$Combo)
n <- 1654
mean(testing_dataset$Combo)
testing_dataset$biases = testing_dataset$preds - testing_dataset$Combo
breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)
bias_group_local2 <- aggregate(testing_dataset$biases, list(testing_dataset$group), mean)[,2]
rbias_group_local2 <- aggregate(testing_dataset$biases, 
                                list(testing_dataset$group), mean)[,2]/
  aggregate(testing_dataset$Combo, list(testing_dataset$group), mean)[,2]*100

p4<- ggplot(testing_dataset, aes(x=Combo, y=preds))+ 
  geom_pointdensity()+
  theme_bw()+
  annotate("text",x=0,y=91,hjust = 0,size = 10,family= "A",
           label= paste(expression(" "~R^2),": ",round(rsqd_local2,3) 
           ),parse=TRUE) + 
  annotate("text",x=0,y=77,hjust = 0,size = 10,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE_local2,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE_local2,3),"%",
             "\n" , " Bias: ", round(bias_local2,3)," Mg/ha",
             "\n" , " Sample size: ", n)) + 
  annotate("text",x=80,y=100,hjust = 0,size = 15,family= "A", label= "(c)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("AGBD"["SAS_RF"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))

#suf = "RF_PCA_"
#out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_scatterplot.jpg"))
#ggsave(out,height=9, width=9, dpi=600)

breakbin = round(seq(0,100,5))
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=5)

bias_group_local2
rbias_group_local2
testing_dataset$biases <- testing_dataset$preds - testing_dataset$Combo
testing_dataset <- na.omit(testing_dataset)
table(testing_dataset$group)

x <- round(seq(0,45,5))
a_local2<- round(bias_group_local2[0:10],1)
c_local2<- round(rbias_group_local2[0:10],1)
c_local2 <- paste0(c_local2, "%")
d_local2 <- a_local2
d_local2[d_local2>0]=0
#num <- table(testing_dataset$group)[0:10]

l <- cbind(x,a_local2,num)
df <- data.frame(l)

b4<- ggplot(df, aes(x=x+2.5, y=a_local2, group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=9), size = 6)+
  geom_text(aes(label = a_local2,y=d_local2-2), size = 6)+
  geom_text(aes(label = c_local2,y=d_local2-4), size = 6)+
  theme_bw()+
  coord_cartesian(ylim = c(-25, 10))+
  annotate("text",x=45,y=5,hjust = 0,size = 15,family= "A", label= "(c)") + 
  scale_y_continuous(minor_breaks = seq(-30, 10, 3),breaks =  seq(-30, 10, 3))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y=expression("MSD of AGBD"["SAS_RF"]~"(Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



#original
suf = "model4_"
ggarrange(p2,p1,p4,nrow = 1, ncol = 3)
out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_scatterplot_glm_updated.jpg"))
ggsave(out,height=10, width=30, dpi=600)

ggarrange(b2,b1,b4,nrow = 1, ncol = 3)
out = file.path("E:\\Biomass\\CSIR\\Result_0618\\2018",paste0(suf,"_boxplot_glm_updated.jpg"))
ggsave(out,height=10, width=30, dpi=600)





#violin plots -- ALS AGBD, L4A, SAS GLM, SAS RF
#------------------------------------------------------------------------------------------------#
#09252023
#SAS applied to all GEDI footprints -- glm
#Local updated model
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

suf <- "H_CC_linear_0_glm"
Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

Data <- Data[Data$AGBD > 0,]
Data <- Data[Data$SUM > 24,]

Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

Data = Data[!duplicated(Data$shot_numbe_x),]
Data_50 <-Data[Data$AGBD < 40,]
Data_50_100 <- Data[Data$AGBD >= 40 & Data$AGBD < 100,]
random_sample <- createDataPartition(Data_50$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[-random_sample, ]
training_dataset <- rbind(training_dataset, Data_50_100)
testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset = testing_dataset[!duplicated(testing_dataset$shot_numbe_x),]
training_dataset = training_dataset[!duplicated(training_dataset$shot_numbe_x),]

repeat_cv <- trainControl(method='cv', number=5)
reg <- train(
  Combo ~ GEDI_FHD + RH1_50 + RH1_98, 
  data=training_dataset, 
  method='glm', 
  family = Gamma(link = "identity"),
  start=c(0,0.1,100,5),
  trControl=repeat_cv)
testing_dataset$preds <- predict(object=reg, 
                                 newdata=testing_dataset[, ])

breakbin = seq(0,60,10)
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=6)
testing_dataset <- na.omit(testing_dataset)

#ALS AGBD 
p0<- ggplot(testing_dataset, aes(x="", y=Combo))+ 
  geom_violin(trim=FALSE, fill="gray")+
  theme_bw()+
  coord_cartesian(ylim = c(0, 80))+
  scale_y_continuous(minor_breaks = seq(0, 80, 10),breaks =  seq(0, 80, 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression("AGBD"["ALS_area_MCH_CC"]~"(Mg/ha)"), 
       y="AGBD (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))

#glm
p1<- ggplot(testing_dataset, aes(x="", y=preds))+ 
  geom_violin(trim=FALSE, fill="gray")+
  theme_bw()+
  coord_cartesian(ylim = c(0, 80))+
  scale_y_continuous(minor_breaks = seq(0, 80, 10),breaks =  seq(0, 80, 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression("AGBD"["SAS_GLM"]), 
       y="AGBD (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


#on-orbit GEDI L4A
p2<- ggplot(testing_dataset, aes(x="", y=AGBD))+ 
  geom_violin(trim=FALSE, fill="gray")+
  theme_bw()+
  coord_cartesian(ylim = c(0, 80))+
  scale_y_continuous(minor_breaks = seq(0, 80, 10),breaks =  seq(0, 80, 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression("AGBD"["L4A"]), 
       y="AGBD (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


#rf
dir <- "E:\\Biomass\\CSIR\\Result_0618\\2018\\csv\\All_06192023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

suf <- "H_CC_linear_0_RF"
Data$Combo <- 9.0665*Data$MEAN*Data$SUM/490

Data <- Data[Data$AGBD >= 0,]

Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1

Data_sel <- Data[,c("RH1_50","RH1_75","RH1_98","GEDI_Cover","GEDI_FHD")]
Data_sel <- na.omit(Data_sel)

Data = Data[!duplicated(Data$shot_numbe_x),]
Data_50 <-Data[Data$AGBD < 40,]
Data_50_100 <- Data[Data$AGBD >= 40 & Data$AGBD < 100,]
random_sample <- createDataPartition(Data_50$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data_50[random_sample, ]
testing_dataset <- Data_50[-random_sample, ]
training_dataset <- rbind(training_dataset, Data_50_100)
testing_dataset <- rbind(testing_dataset, Data_50_100)
testing_dataset = testing_dataset[!duplicated(testing_dataset$shot_numbe_x),]
training_dataset = training_dataset[!duplicated(training_dataset$shot_numbe_x),]

#rf
repeat_cv <- trainControl(method='cv', number=5)
reg <- train(
  Combo~sensitivity_a1+
    GEDI_Cover + GEDI_FHD + 
    RH1_50 + RH1_75 + RH1_90 + RH1_98,
  data=training_dataset, 
  method='rf', 
  ntree=100,
  nodesize = 4,
  trControl=repeat_cv)

testing_dataset$preds <- predict(object=reg, 
                                 newdata=testing_dataset[, ])

breakbin = seq(0,60,10)
testing_dataset$group <- cut(testing_dataset$Combo,breaks = breakbin,dig.lab=6)
testing_dataset <- na.omit(testing_dataset)
p3<- ggplot(testing_dataset, aes(x="", y=preds))+ 
  geom_violin(trim=FALSE, fill="gray")+
  theme_bw()+
  coord_cartesian(ylim = c(0, 80))+
  scale_y_continuous(minor_breaks = seq(0, 80, 10),breaks =  seq(0, 80, 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x=expression("AGBD"["SAS_RF"]), 
       y="AGBD (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(p0,p2,p1,p3,ncol=4)
out = file.path("C:\\Users\\Shawn\\Desktop\\figures\\violin.jpg")
ggsave(out,height=10, width=30, dpi=600)





# ------------------------------------------------------------------------------------------------ #
#Field-ALS-SAR AGBD vs. ALS AGBD
Data_merge <- 0

dir <- "E:\\Biomass\\CSIR\\Result_01292024"
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 'DNyala',
              'Venetia',"Limpopo1", "Limpopo2", "Limpopo3")


for (l in list_dir){
  print(l)
  dir_chm <- file.path(dir,paste0(l,"_MCH.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(1,5))
  
  dir_cover <- file.path(dir,paste0(l,"_CC.csv"))
  Data = read.csv(dir_cover)
  Data_cover = subset(Data, select=c(1,5))
  
  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_chm, Data_cover, by.x = "OID_", by.y = "OID_")
  colnames(Data_index) <- c("FID","MCH","CC")
  Data_index <- na.omit(Data_index)
  Data_index$CC <- Data_index$CC/625
  Data_index <- Data_index[Data_index$MCH > 0 & Data_index$MCH < 15,]
  Data_index <- Data_index[Data_index$CC > 0 & Data_index$CC < 1,]
  Data_merge <- rbind(Data_merge, Data_index)
}



file_path <- file.path("E:\\Biomass\\CSIR\\Result_01292024\\merge.csv")
write.csv(Data_merge, file = file_path)





# ------------------------------------------------------------------------------------------------ #
#I2 output - 20m
Data_merge = 0
dir <- "C:\\Users\\Shawn\\Desktop\\GEDI_0130\\I2"
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 'DNyala',
              'Venetia',"Limpopo1", "Limpopo2", "Limpopo3")


for (l in list_dir){
  print(l)
  dir_indix <- file.path(dir,paste0("I2_20m.xlsx"))
  Data_I2 = read_excel(dir_indix)
  
  dir_chm <- file.path(dir,paste0(l,"_I2_20m_MCH.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  
  dir_cover <- file.path(dir,paste0(l,"_I2_20m_CC.csv"))
  Data = read.csv(dir_cover)
  Data <- Data[Data$COUNT > 255,]
  Data_cover = subset(Data, select=c(2,5))

  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_chm, Data_cover, by.x = "FID", by.y = "FID")
  colnames(Data_index) <- c("FID","MCH","CC")
  Data_index <- na.omit(Data_index)
  Data_index$CC <- Data_index$CC/260
  Data_index <- Data_index[Data_index$MCH > 0 & Data_index$MCH < 15,]
  Data_index <- Data_index[Data_index$CC > 0 & Data_index$CC <= 1,]
  Data_I2 <- merge(Data_I2, Data_index, by.x = "FID", by.y = "FID")
  Data_merge <- rbind(Data_merge, Data_I2)
}

file_path <- file.path("C:\\Users\\Shawn\\Desktop\\GEDI_0130\\I2_20m.csv")
write.csv(Data_merge, file = file_path)



#I2 output - 100m
Data_merge = 0
dir <- "C:\\Users\\Shawn\\Desktop\\GEDI_0130\\I2"
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 'DNyala',
              'Venetia',"Limpopo1", "Limpopo2", "Limpopo3")


for (l in list_dir){
  print(l)
  dir_indix <- file.path(dir,paste0("I2_100m.xlsx"))
  Data_I2 = read_excel(dir_indix)
  
  dir_chm <- file.path(dir,paste0(l,"_I2_100m_MCH.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  
  dir_cover <- file.path(dir,paste0(l,"_I2_100m_CC.csv"))
  Data = read.csv(dir_cover)
  Data <- Data[Data$COUNT > 1200,]
  Data_cover = subset(Data, select=c(2,5))
  
  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_chm, Data_cover, by.x = "FID", by.y = "FID")
  colnames(Data_index) <- c("FID","MCH","CC")
  Data_index <- na.omit(Data_index)
  Data_index$CC <- Data_index$CC/1300
  Data_index <- Data_index[Data_index$MCH > 0 & Data_index$MCH < 15,]
  Data_index <- Data_index[Data_index$CC > 0 & Data_index$CC <= 1,]
  Data_I2 <- merge(Data_I2, Data_index, by.x = "FID", by.y = "FID")
  Data_merge <- rbind(Data_merge, Data_I2)
}

file_path <- file.path("C:\\Users\\Shawn\\Desktop\\GEDI_0130\\I2_100m.csv")
write.csv(Data_merge, file = file_path)


