#Author: Xiaoxuan Li
#08302022
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
windowsFonts(A = windowsFont("Times New Roman"))



# ------------------------------------------------------------------------------------------------ #
#recalculate CSIR biomass 
dir <- "E:\\Biomass\\Biomass.xlsx"
data <- read_excel(dir, sheet = "Biomass")

#3-5 cm
Select <- data[(data$SP == 1.10) | (data$SP == 2.10)  | (data$SP == 3.10) | (data$SP == 4.10), ]

Select2 <- Select[(Select$Diameter >= 3) &  (Select$Diameter < 5), ]
Select2 <- Select2[!grepl(paste(c('Plot','plot', 'Just', 'border', 'outside','Pushed',
                                  'Next','chopped','Dead','dead','Bent','elephants',
                                  'boundary','Edge','GPS','previous','Subploy'), 
                                collapse="|"), Select2$Notes),]
aggregate(Select2[, "Biomass (Kg)"], list(Select2$Site), sum)

#5-10 cm
Select <- data[(data$SP == 1.25) | (data$SP == 2.25)  | 
                 (data$SP == 3.25) | (data$SP == 4.25), ]

Select2 <- Select[(Select$Diameter >= 5) &  (Select$Diameter < 10), ]
Select2 <- Select2[!grepl(paste(c('Plot','plot', 'Just', 'border', 'outside','Pushed',
                                  'Next','chopped','Dead','dead','Bent','elephants',
                                  'boundary','Edge','GPS','previous','Subploy'), 
                              collapse="|"), Select2$Notes),]
#Select3 <- Select2[(Select2$Site == "S70") , ]

aggregate(Select2[, "Biomass (Kg)"], list(Select2$Site), sum)


#>10 cm
#Select <- data[(data$SP == 5), ]
Select <- data[(data$Diameter >= 10), ]
Select <- Select[!grepl(paste(c('Plot', 'plot', 'Just', 'border', 'outside','Pushed',
                                'Next','chopped','Dead','dead','Bent','elephants',
                                'boundary','Edge','GPS','previous','Subploy'), 
                              collapse="|"), Select$Notes),]
#Select3 <- Select[(Select$Site == "S12") , ]

aggregate(Select[, "Biomass (Kg)"], list(Select$Site), sum)
  



#individual tree segementation -- Colgan 2013---improved
chm <- raster("E:\\Biomass\\ras\\DN_S8_SP_5_DNyala_h.tif")
intens <- raster("E:\\Biomass\\ras\\DN_S8_SP_5_DNyala_i.tif")
plot(chm)
#height cv
fsd_h <- focal(chm, w=matrix(1,3,3), fun=sd, na.rm=TRUE, pad=TRUE)
fcv_h <- fsd_h/chm

#intensity cv
fsd_i <- focal(intens, w=matrix(1,3,3), fun=sd,na.rm=TRUE, pad=TRUE)
fcv_i <- fsd_i/intens

fcv_h <- crop(fcv_h,fcv_i)
chm <- crop(chm,fcv_i)
fcv_i <- crop(fcv_i,chm)

#test
cv1 <- 0.3
#ground
chm[] <- ifelse(chm[]<=0.5, 0, chm[])
#tree
chm[] <- ifelse(chm[]>=1.5 & fcv_h[]<cv1 , -1, chm[])
#shurb
chm[] <- ifelse(chm[]<=1.5 & chm[]>=0.5 & fcv_h[]<cv1 , -2, chm[])
#edge
chm[] <- ifelse(chm[] != 0 & chm[] != -1 & chm[] != -2, -3, chm[])

myColor <- c("black", "red", "green", "white")
plot(chm,col=myColor)



chm[] <- ifelse(chm[]>=1.5 & chm[]<=5 &fcv_h[]>=0.3, 0, chm[])
plot(chm)
crowns = dalponte2016(chm, ttops, 
                      max_cr = 11, 
                      th_seed = 0.5,
                      th_cr = 0.5,
                      th_tree = 3)()

myColor <- randomcoloR::distinctColorPalette(k = 100)
plot(chm)
plot(crowns,col=myColor,add=T)
plot(ttops,add=T)

#second round
chm[!is.na(crowns[])] <- NA
ttops2 <- find_trees(chm, lmf(13, 1.5,shape = "circular"))

#algorithms
crowns2 = dalponte2016(chm, ttops2, 
                       max_cr = 11, 
                       th_seed = 0.55,
                       th_cr = 0.55,
                       th_tree = 1.5)()

plot(chm)
plot(crowns2,col=myColor,add=T)
plot(ttops2,add=T)

writeRaster(crowns, 'E:\\Biomass\\lidR_crown\\DN_S8_SP_5_DNyala_h.tif',overwrite=TRUE)





# ------------------------------------------------------------------------------------------------ #

#new method ---buld process
chm <- raster("E:\\Biomass\\ras\\DN_S8_SP_5_DNyala_h.tif")
#find tree
smooth <- focal(chm, w=matrix(1,3,3), fun=mean)
#ttops <- find_trees(smooth, lmf(10, 5,shape = "circular"))
ttops <- find_trees(smooth, lmf(30, 8,shape = "circular"))
cv_test <- 0.3

plot(chm)
plot(ttops,add=T)
df_tree <- data.frame(ttops)[,1:2]
#grow region algorithm
#parameter: coefficient of variation of height 
fsd_h <- focal(chm, w=matrix(1,3,3), fun=sd, na.rm=TRUE, pad=TRUE)
fcv_h <- fsd_h/chm
#extract the single tree pixel from original chm
chm_s <- rasterize(ttops, chm, field = "treeID")
plot(chm_s)
#identify cell number of identified tree pixel 
cn <- cellFromXY(chm_s, ttops@coords)
#identify cell numbers of all neighborhood pixels around the region/tree pixel
ad <- adjacent(chm, cell = c(cn), directions=8, sorted=TRUE) 

for (k in 1: 12){
  #extract chm/cv values corresponding to identified cell numbers
  chm_ext <- chm[c(ad[,2])]
  cv_ext <- fcv_h[c(ad[,2])]
  Tree <- data.frame(cbind(ad[,2],chm_ext,cv_ext))
  print(k)
  combo$index <- 0
  combo <- unique(combo)
  for (i in 1: nrow(combo)){
    combo[,4][i] <- which.min(abs(abs(colFromCell(chm_s, c(combo[,1][i])) - 
                                        colFromCell(chm_s, c(cn)))  + 
                                    abs(rowFromCell(chm_s, c(combo[,1][i])) - 
                                          rowFromCell(chm_s, c(cn))) ))
  }
  chm_s[c(combo[,1])] <- c(combo[,4])
  #identify cell numbers of all neighborhood pixels around the region/tree pixel
  combo <- merge(x = combo, y = df_tree, by.x = "index", by.y = "treeID")
  #query chm value, and cv value
  combo <- combo[combo$chm_ext > 1.5 & 
                   combo$chm_ext >= combo$Z/2 & 
                   combo$cv_ext<=cv_test,]
  #combo <- combo[-1,]
  plot(chm_s)
  ad <- adjacent(chm, cell = c(combo[,2]), directions=8, sorted=TRUE) 
  
}






# ------------------------------------------------------------------------------------------------ #
#new method ---individual tree process
filename <- "Test_0929.tif"
suf<- sapply(strsplit(filename, "_"), "[", 2)
out_dir <- file.path("E:\\Biomass\\CSIR\\lidR_crown", file_path_sans_ext(filename))
dir.create(out_dir, showWarnings = FALSE)
chm <- raster(file.path("E:\\Biomass\\CSIR\\ras",filename))
chm[is.na(chm[])] <- 0
plot(chm)
d_list <- c(40, 40, 30, 30, 20, 20)
h_dist <- c(10, 10, 8,  8,  6,  6)
k_list <- c(15, 15, 10, 10, 10, 10)
cv_list<- c(0.2,0.2,0.2,0.2,0.2)
hs_list <- c(0.45,0.45,0.45,0.45,0.45)
for (x in 1:6){
  print(x)
  try({
  #find tree
  if (x <=2){
      chm <- focal(chm, w=matrix(1,3,3), fun=mean)
      fsd_h <- focal(chm, w=matrix(1,3,3), fun=sd, na.rm=TRUE, pad=TRUE)
      fcv_h <- fsd_h/chm
    }
  ttops <- find_trees(chm, lmf(d_list[x], h_dist[x], shape = "circular"))
  chm_append <- rasterize(ttops, chm, field = "treeID")
  chm_append[is.na(chm_append[])] <- 0
  chm_append[!is.na(chm_append[])] <- 0
  df_tree <- data.frame(ttops)[,1:2]

  for (t in 1: nrow(ttops)){
      sel <- ttops[ttops$treeID == t,]
      #extract the single tree pixel from original chm
      chm_s <- rasterize(sel, chm, field = "treeID")
      chm_s[is.na(chm_s)] <- 0
      #identify cell number of identified tree pixel 
      cn <- cellFromXY(chm_s, sel@coords)
      chm_s[cn] <- ttops[[1]][t]
      #identify cell numbers of all neighborhood pixels around the region/tree pixel
      ad <- adjacent(chm, cell = c(cn), directions=8, sorted=TRUE) 
      for (k in 1: k_list[x]){
        combo<-0
        #extract chm/cv values corresponding to identified cell numbers
        chm_ext <- chm[c(ad[,2])]
        cv_ext <- fcv_h[c(ad[,2])]
        combo <- data.frame(cbind(ad[,2],chm_ext,cv_ext))
        #combo <- unique(combo)
        combo$index <- ttops[[1]][t]
        combo <- combo[combo$chm_ext > 0.5 & 
                         combo$chm_ext >= ttops[[2]][t]*hs_list[x] & 
                         combo$cv_ext<=cv_list[x],]
        chm_s[c(combo[,1])] <- ttops[[1]][t]
        ad <- adjacent(chm, cell = c(combo[,1]), directions=8, sorted=TRUE)
        ad <- unique(ad)

      }
      combo <- data.frame(cbind(ad[,2],chm_ext))
      combo <- combo[combo$chm_ext > 1.5,]
      chm_s[c(combo[,1])] <- ttops[[1]][t]
      chm_append <- chm_append + chm_s
  }

  chm_append[chm_append==0] <- NA
  writeRaster(chm_append+100*x, 
              file=file.path(out_dir,
                             paste0(x,"_res.tif")),
                             overwrite=TRUE)
  chm[chm_append > 0] <- 0
  })
  plot(chm_append)
  plot(chm)
  plot(ttops,add=T)
}

#dalponte2016 algorithms find short trees
ttops <- find_trees(chm, lmf(7, 3, shape = "circular"))
myColor <- randomcoloR::distinctColorPalette(k = 100)
crowns = dalponte2016(chm, ttops, 
                       max_cr = 5, 
                       th_seed = 0.75,
                       th_cr = 0.75,
                       th_tree = 1.5)()
plot(chm)
plot(ttops,add=T)
plot(crowns,col=myColor,add=T)
writeRaster(crowns, file.path(out_dir,paste0("10_res.tif")),overwrite=TRUE)

a <- overlay(chm, fcv_h, fun = function(x, y) {x[x >=0.5  & x<=1.5 & y<=0.2] <- 1; x})
a[a != 1] <- NA
plot(a)
plot(chm)
biomass_shurb <- (freq(a)[1,][2]*3.04-111)/1000
biomass_shurb
writeRaster(a, file.path(out_dir,paste0("shurb.tif")),overwrite=TRUE)


lst_res <- list.files(path = out_dir, full.names = TRUE, pattern = "\\.tif$")
rlist <- lapply(lst_res, raster)
rlist$fun <- max
rlist$na.rm <- TRUE
x <- do.call(mosaic, rlist) 
plot(x)
writeRaster(x, file.path(out_dir,paste0("result.tif")),overwrite=TRUE)




#for 10 field plot
filename <- "DN_S10_SP_5_DNyala_h.tif"
chm <- raster(file.path("E:\\Biomass\\ras",filename))
chm[is.na(chm[])] <- 0
out_dir <- file.path("E:\\Biomass\\lidR_crown", file_path_sans_ext(filename))
dir.create(out_dir, showWarnings = FALSE)
chm <- focal(chm, w=matrix(1,3,3), fun=mean)

ttops <- find_trees(chm, lmf(7, 3, shape = "circular"))
myColor <- randomcoloR::distinctColorPalette(k = 100)
crowns = dalponte2016(chm, ttops, 
                      max_cr = 15, 
                      th_seed = 0.75,
                      th_cr = 0.75,
                      th_tree = 1.5)()
plot(chm)
plot(ttops,add=T)
plot(crowns,col=myColor,add=T)
writeRaster(crowns, file.path(out_dir,paste0("result.tif")),overwrite=TRUE)


#25m GEDI footprint based biom
dir = "E:\\Biomass\\CSIR\\ras\\CSIR"
f <- list.files(path=dir, pattern='.tif$', full.names=TRUE)
for (filename in f){
  print(filename)
  out_dir <- file.path("E:\\Biomass\\CSIR\\lidR_crown", file_path_sans_ext(basename(filename)))
  dir.create(out_dir, showWarnings = FALSE)
  chm <- raster(filename)
  chm[is.na(chm[])] <- 0
  plot(chm)
  d_list <- c(10, 10, 10, 10, 10, 10, 10)
  h_dist <- c(10, 9, 8, 7, 6, 5, 4)
  k_list <- c(10, 9, 8, 7, 6, 5, 4)
  cv_list<- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
  hs_list <- c(0.65, 0.65, 0.65, 0.65, 0.65, 0.65)
  for (x in 1:6){
    print(x)
    try({
      #find tree
      if (x <=2){
        chm <- focal(chm, w=matrix(1,3,3), fun=mean)
        fsd_h <- focal(chm, w=matrix(1,3,3), fun=sd, na.rm=TRUE, pad=TRUE)
        fcv_h <- fsd_h/chm
      }
      ttops <- find_trees(chm, lmf(d_list[x], h_dist[x], shape = "circular"))
      chm_append <- rasterize(ttops, chm, field = "treeID")
      chm_append[is.na(chm_append[])] <- 0
      chm_append[!is.na(chm_append[])] <- 0
      df_tree <- data.frame(ttops)[,1:2]
      for (t in 1: nrow(ttops)){
        sel <- ttops[ttops$treeID == t,]
        #extract the single tree pixel from original chm
        chm_s <- rasterize(sel, chm, field = "treeID")
        chm_s[is.na(chm_s)] <- 0
        #identify cell number of identified tree pixel 
        cn <- cellFromXY(chm_s, sel@coords)
        chm_s[cn] <- ttops[[1]][t]
        #identify cell numbers of all neighborhood pixels around the region/tree pixel
        ad <- adjacent(chm, cell = c(cn), directions=8, sorted=TRUE) 
        for (k in 1: k_list[x]){
          combo<-0
          #extract chm/cv values corresponding to identified cell numbers
          chm_ext <- chm[c(ad[,2])]
          cv_ext <- fcv_h[c(ad[,2])]
          combo <- data.frame(cbind(ad[,2],chm_ext,cv_ext))
          #combo <- unique(combo)
          combo$index <- ttops[[1]][t]
          combo <- combo[combo$chm_ext > 0.5 & 
                           combo$chm_ext >= ttops[[2]][t]*hs_list[x] & 
                           combo$cv_ext<=cv_list[x],]
          chm_s[c(combo[,1])] <- ttops[[1]][t]
          ad <- adjacent(chm, cell = c(combo[,1]), directions=8, sorted=TRUE)
          ad <- unique(ad)
          
        }
        combo <- data.frame(cbind(ad[,2],chm_ext))
        combo <- combo[combo$chm_ext > 1.5,]
        chm_s[c(combo[,1])] <- ttops[[1]][t]
        chm_append <- chm_append + chm_s
        }
        chm_append[chm_append==0] <- NA
        writeRaster(chm_append+100*x,  
                    file=file.path(out_dir,
                                   paste0(x,"_res.tif")),
                    overwrite=TRUE)
        chm[chm_append > 0] <- 0

        plot(chm_append)
        plot(chm)
        plot(ttops,add=T)
    })
    }
  try({
    ttops <- find_trees(chm, lmf(10, 3, shape = "circular"))
    myColor <- randomcoloR::distinctColorPalette(k = 100)
    crowns = dalponte2016(chm, ttops, 
                          max_cr = 7, 
                          th_seed = 0.65,
                          th_cr = 0.65,
                          th_tree = 1.5)()
    plot(chm)
    plot(ttops,add=T)
    plot(crowns,col=myColor,add=T)
    writeRaster(crowns, file.path(out_dir,paste0("10_res.tif")),overwrite=TRUE)
    
    a <- overlay(chm, fcv_h, fun = function(x, y) {x[x >=0.5  & x<=1.5 & y<=0.2] <- 1000; x})
    a[a != 1000] <- NA
    plot(a)
    plot(chm)
    biomass_shurb <- (freq(a)[1,][2]*3.04-111)/1000
    biomass_shurb
    #if (biomass_shurb >0){writeRaster(a, file.path(out_dir,paste0("shurb.tif")),overwrite=TRUE)}
    
    lst_res <- list.files(path = out_dir, full.names = TRUE, pattern = "\\.tif$")
    rlist <- lapply(lst_res, raster)
    rlist$fun <- max
    rlist$na.rm <- TRUE
    x <- do.call(mosaic, rlist) 
    plot(x)
    writeRaster(x, file.path(out_dir,paste0("result.tif")),overwrite=TRUE)
  })
}



# ------------------------------------------------------------------------------------------------ #
#new method ---individual tree process --- test single ras 10/03/2022
dir = "E:\\Biomass\\CSIR\\ras"
f <- list.files(path=dir, pattern='.tif$', full.names=TRUE)
list_shurb = 0
for (filename in f){
  print(filename)
  suf <- sapply(strsplit(filename, "_"), "[", 2)
  out_dir <- file.path("E:\\Biomass\\CSIR\\lidR_crown", file_path_sans_ext(basename(filename)))
  dir.create(out_dir, showWarnings = FALSE)
  chm <- raster(filename)
  chm[is.na(chm[])] <- 0
  chm <- focal(chm, w=matrix(1,3,3), fun=mean)
  fsd_h <- focal(chm, w=matrix(1,3,3), fun=sd, na.rm=TRUE, pad=TRUE)
  fcv_h <- fsd_h/chm
  ttops <- find_trees(chm, lmf(5, 3, shape = "circular"))
  chm_append <- rasterize(ttops, chm, field = "treeID")
  chm_append[is.na(chm_append[])] <- 0
  chm_append[!is.na(chm_append[])] <- 0
  df_tree <- data.frame(ttops)[,1:2]
  plot(chm)
  plot(ttops,add=T)
  for (t in 1: nrow(ttops)){
    try({
      sel <- ttops[ttops$treeID == t,]
      chm_s <- rasterize(sel, chm, field = "treeID")
      chm_s[is.na(chm_s)] <- 0
      cn <- cellFromXY(chm_s, sel@coords)
      chm_s[cn] <- ttops[[1]][t]
      ad <- adjacent(chm, cell = c(cn), directions=8, sorted=TRUE)
      for (k in 1: 15){
        combo <- 0
        chm_ext <- chm[c(ad[,2])]
        cv_ext <- fcv_h[c(ad[,2])]
        combo <- data.frame(cbind(ad[,2],chm_ext,cv_ext))
        combo$index <- ttops[[1]][t]
        combo <- combo[combo$chm_ext > 1.5 & 
                         combo$chm_ext >= ttops[[2]][t]*0.6 & 
                         combo$cv_ext <= 0.1,]
        chm_s[c(combo[,1])] <- ttops[[1]][t]
        ad <- adjacent(chm, cell = c(combo[,1]), directions=8, sorted=TRUE)
        ad <- unique(ad)
        }
        combo <- data.frame(cbind(ad[,2],chm_ext))
        chm_s[c(combo[,1])] <- ttops[[1]][t]
        chm_append <- chm_append + chm_s
        chm[chm_append > 0] <- 0
        ad <- adjacent(chm, cell = c(combo[,1]), directions=8, sorted=TRUE)
        s <- unique(ad[,1])
        m <- unique(ad[,2])
        m<- m[!(m %in% s)]
        chm[c(m)] <- 0
    })
    }
    chm_append[chm_append == 0] <- NA
    plot(chm_append,add=T)
    writeRaster(chm_append+100,
                file=file.path(out_dir,paste0("tree_1.tif")),
                overwrite=TRUE)
    ttops <- find_trees(chm, lmf(5, 3, shape = "circular"))
    plot(chm)
    plot(ttops,add=T)
    chm_append <- rasterize(ttops, chm, field = "treeID")
    chm_append[is.na(chm_append[])] <- 0
    chm_append[!is.na(chm_append[])] <- 0
    df_tree <- data.frame(ttops)[,1:2]
    myColor <- randomcoloR::distinctColorPalette(k = 100)
    crowns = dalponte2016(chm, ttops,
                          max_cr = 15,
                          th_seed = 0.7,
                          th_cr = 0.7,
                          th_tree = 1.5)()
    plot(chm)
    plot(ttops,add=T)
    plot(crowns,col=myColor,add=T)
    writeRaster(crowns, file.path(out_dir,paste0("tree_2.tif")),overwrite=TRUE)
    #result
    lst_res <- list.files(path = out_dir, full.names = TRUE, pattern = "\\.tif$")
    rlist <- lapply(lst_res, raster)
    rlist$fun <- max
    rlist$na.rm <- TRUE
    x <- do.call(mosaic, rlist)
    plot(x)
    writeRaster(x, file.path(out_dir,paste0("result.tif")),overwrite=TRUE)
    #shurb
    a <- overlay(chm, fun = function(x, y) {x[x >=0.5 & x<=1.5] <- 1; x})
    a[a != 1] <- NA
    biomass_shurb <- (freq(a)[1,][2]*3.04-111)/1000
    list_shurb = c(list_shurb, biomass_shurb, suf)
}
list_shurb
write.csv(list_shurb, file = "C:\\Users\\Shawn\\Desktop\\biom_shurb.csv")










# ------------------------------------------------------------------------------------------------ #
#new method ---individual tree process --- test all site ras 10/04/2022
dir = "E:\\Biomass\\CSIR\\ras"
sites <- list.files(path=dir, full.names=FALSE)
#sites[sites != "Welverdiendt"]
for (site in sites){
  print(site)
  new_dir <- file.path(dir,site)
  out_dir <- file.path("E:\\Biomass\\CSIR\\lidR_crown",site)
  list_shurb = 0
  list_filename = 0
  f <- list.files(path=new_dir, pattern='.tif$', full.names=TRUE)
  for (filename in f){
    print(filename)
    chm <- raster(filename)
    chm[is.na(chm[])] <- 0
    chm <- focal(chm, w=matrix(1,3,3), fun=mean)
    ttops <- find_trees(chm, lmf(5, 3, shape = "circular"))
    crowns = dalponte2016(chm, ttops,
                          max_cr = 15,
                          th_seed = 0.7,
                          th_cr = 0.7,
                          th_tree = 1.5)()
    writeRaster(crowns, file.path(out_dir, paste0(file_path_sans_ext(basename(filename)), ".tif")), overwrite=TRUE)
    #shurb
    chm[chm >=0.5 & chm<=1.5] <- 1
    chm[chm != 1] <- NA
    biomass_shurb <- (freq(chm)[1,][2]*3.04-111)/1000
    if (biomass_shurb<0){
      biomass_shurb=0
    }
    list_shurb = c(list_shurb, biomass_shurb)
    list_filename = c(list_filename, file_path_sans_ext(basename(filename)))
  }
  
  list = cbind(list_shurb,list_filename)
  write.csv(list, file = file.path("E:\\Biomass\\CSIR\\result",site,"biom_shurb.csv"))
}




# ------------------------------------------------------------------------------------------------ #
#new method ---individual tree process --- test single ras 10/07/2022, CSIR 1ha
dir = "E:\\Biomass\\CSIR\\ras\\CSIR"
out_dir <- file.path("E:\\Biomass\\CSIR\\lidR_crown\\CSIR")
list_shurb = 0
list_filename = 0
f <- list.files(path=dir, pattern='_h.tif$', full.names=TRUE)
for (filename in f){
  print(filename)
  chm <- raster(filename)
  chm[is.na(chm[])] <- 0
  chm <- focal(chm, w=matrix(1,3,3), fun=mean)
  ttops <- find_trees(chm, lmf(5, 3, shape = "circular"))
  crowns = dalponte2016(chm, ttops,
                        max_cr = 15,
                        th_seed = 0.7,
                        th_cr = 0.7,
                        th_tree = 1.5)()
  writeRaster(crowns, file.path(out_dir, paste0(file_path_sans_ext(basename(filename)), ".tif")), overwrite=TRUE)
  #shurb
  chm[chm >=0 & chm<=1.5] <- 1
  chm[chm != 1] <- NA
  biomass_shurb <- (freq(chm)[1,][2]*3.04-111)/1000
  if (biomass_shurb<0){
    biomass_shurb=0
  }
  list_shurb = c(list_shurb, biomass_shurb)
  list_filename = c(list_filename, file_path_sans_ext(basename(filename)))
}

list = cbind(list_shurb,list_filename)
write.csv(list, file = "E:\\Biomass\\CSIR\\result\\CSIR\\biom_shurb_CSIR.csv")




# ------------------------------------------------------------------------------------------------ #
#new method ---individual tree process --- test single ras 1/28/2023, CSIR 25m
dir = "E:\\Biomass\\CSIR\\ras\\CSIR_25m"
out_dir <- file.path("E:\\Biomass\\CSIR\\lidR_crown\\CSIR_25m")
list_shurb = 0
list_filename = 0
f <- list.files(path=dir, pattern=glob2rx('*_h.tif$*'), full.names=TRUE)

f <- grep('_Ven_', f, value = TRUE)

for (filename in f){
  print(filename)
  chm <- raster(filename)
  chm[is.na(chm[])] <- 0
  chm <- focal(chm, w=matrix(1,3,3), fun=mean)
  ttops <- find_trees(chm, lmf(5, 3, shape = "circular"))
  crowns = dalponte2016(chm, ttops,
                        max_cr = 15,
                        th_seed = 0.7,
                        th_cr = 0.7,
                        th_tree = 1.5)()
  writeRaster(crowns, file.path(out_dir, paste0(file_path_sans_ext(basename(filename)), ".tif")), overwrite=TRUE)
  #shurb
  chm[chm >=0.5 & chm<=1.5] <- 1
  chm[chm != 1] <- NA
  biomass_shurb <- (freq(chm)[1,][2]*3.04-111)/1000
  if (biomass_shurb<0){
    biomass_shurb=0
  }
  list_shurb = c(list_shurb, biomass_shurb)
  list_filename = c(list_filename, file_path_sans_ext(basename(filename)))
}

list = cbind(list_shurb,list_filename)
write.csv(list, file = "E:\\Biomass\\CSIR\\result\\CSIR_25m\\biom_shurb_CSIR_25m_BBR.csv")



# ------------------------------------------------------------------------------------------------ #
#new method ---individual tree process --- test single ras 10/11/2022, individual tree GPS
filename = "E:\\Biomass\\CSIR\\ras_CSIR\\Ven_S60_SP_5_Venetia_h.tif"
out_dir <- "E:\\Biomass\\CSIR\\lidR_crown"
shp = "E:\\Biomass\\CSIR\\individual_tree\\Venetia_Trees\\Ven_S60_Trees.shp"
xx <- readOGR(shp)
shp_utm <- spTransform(xx, crs(proj4string(chm)))
print(filename)
chm <- raster(filename)
chm[is.na(chm[])] <- 0
chm <- focal(chm, w=matrix(1,3,3), fun=mean)
plot(chm)
plot(shp_utm ,add=T)
ttops <- find_trees(chm, lmf(5, 3, shape = "circular"))
crowns = dalponte2016(chm, shp_utm,
                      max_cr = 5,
                      th_seed = 0.7,
                      th_cr = 0.7,
                      th_tree = 1.5)()
plot(crowns)
writeRaster(crowns, file.path(out_dir, paste0(file_path_sans_ext(basename(filename)), "_2.tif")), overwrite=TRUE)
#shurb
chm[chm >=0.5 & chm<=1.5] <- 1
chm[chm != 1] <- NA
biomass_shurb <- (freq(chm)[1,][2]*3.04-111)/1000
if (biomass_shurb<0){
  biomass_shurb=0
}
list_shurb = c(list_shurb, biomass_shurb)
list_filename = c(list_filename, file_path_sans_ext(basename(filename)))


list = cbind(list_shurb,list_filename)
write.csv(list, file = "E:\\Biomass\\CSIR\\result\\biom_shurb_CSIR.csv")




#------------------------------------------------------------------------------------------------ #
#multi-linear regression 
dir <- "E:\\Biomass\\CSIR\\result\\All_biom_l2b.csv"
Data = read.csv(dir,header=T)
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- Data[Data$Combo !=0,]
Data <- Data[Data$status == "Leaf-on",]


#GEDI AGBD vs. GEDI metrics
reg <- lm(formula = AGBD_1 ~ RH1_98, data = Data)
summary(reg)
reg <- lm(formula = AGBD_1 ~ RH1_50 + RH1_50*RH1_75 + RH1_98, data = Data)
summary(reg)

#multi-variate regression -- Laura method
#GEDI metrics derived AGBD vs. CHM reference AGBD
reg <- lm(formula = Combo ~ RH1_50 + RH1_50*RH1_75 + RH1_98, data = Data)
summary(reg)
#multi-variate regression -- more parameters added
#GEDI metrics derived AGBD vs. CHM reference AGBD
reg <- lm(formula = Combo ~ site+orbit+Beam+sensitivity_a1+solar_elevation+
            leaf_on_doy +leaf_on_cycle + leaf_off_doy+leaf_off_flag+
            pgap_theta+pgap_theta_z+
            landsat_treecover+modis_treecover+modis_nonvegetated+
            GEDI_PAI_1 + GEDI_PAI_2 + GEDI_PAI_3 + 
            GEDI_PAVD_1 + GEDI_PAVD_2 + GEDI_PAVD_3 + 
            GEDI_cc_1 + GEDI_cc_2 + GEDI_cc_3 + GEDI_FHD +
            RH1_50 + RH1_50*RH1_75 + RH1_95 + RH1_98 + RH1_100, data = Data)
summary(reg)

Data$Combo <- Data$Combo/0.049
Data$pred <- predict(reg, Data)/0.049
Data <- na.omit(Data)
#stats
r2_1 = round(cor(Data$Combo, Data$pred,method = "pearson")^2,3)
RMSE_1 <- sqrt(mean((Data$Combo - Data$pred)^2))
rRMSE_1 <- 100*sqrt(mean((Data$pred - Data$Combo)^2))/mean(Data$Combo)
r2_2 = round(cor(Data$Combo, Data$AGBD_1,method = "pearson")^2,3)
RMSE_2 <- sqrt(mean((Data$Combo - Data$AGBD_1)^2))
rRMSE_2 <- 100*sqrt(mean((Data$AGBD_1 - Data$Combo)^2))/mean(Data$Combo)


#random forest
dir <- "E:\\Biomass\\CSIR\\result\\All_biom_l2b.csv"
Data = read.csv(dir,header=T)
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data <- Data[Data$Combo !=0,]
Data <- Data[Data$status == "Leaf-on",]
Data$Combo <- sqrt(Data$Combo)
Data$RH1_50 <- sqrt(Data$RH1_50)
Data$rh1_75 <- sqrt(Data$RH1_75)
Data$rh1_98 <- sqrt(Data$RH1_98)
spit = c(train = .8, test = .2)
gs = sample(cut(seq(nrow(Data)), nrow(Data)*cumsum(c(0,spit)),labels = names(spit)))
res = split(Data, gs)
sapply(res, nrow)/nrow(Data)
rf <- randomForest(Combo ~ 
                     site+Beam+sensitivity_a1+elev_lowestmode_a1+
                     landsat_treecover+modis_treecover+modis_nonvegetated+
                     GEDI_PAI_3 + GEDI_PAVD_1 + GEDI_PAVD_2 +
                     GEDI_cc_1 + GEDI_cc_2 + GEDI_FHD +
                     RH1_50 + RH1_50*RH1_75 + RH1_98,
                   data = res$train,ntree = 100,importance=TRUE, nodesize = 20)
pred <- predict(rf, res$test)
round(cor(pred, res$test$Combo,method = "pearson")^2,3)
importance(rf)
varImpPlot(rf)
varImp(rf)



# ------------------------------------------------------------------------------------------------ #
#scatterplot --- biomass pred vs chm 
dir <- "E:\\Biomass\\CSIR\\result\\All_biom_z_L2B.csv"

Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
#Data <- Data[Data$Combo !=0,]
#Data$Combo <- Data$Tree
Data <- Data[Data$status == "Leaf-on",]

r2 = round(cor(Data$RHS_98, Data$GEDI_sim_rh_98,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$RHS_98 - Data$GEDI_sim_rh_98)^2))
rRMSE <- 100*sqrt(mean((Data$RHS_98 - Data$GEDI_sim_rh_98)^2))/mean(Data$GEDI_sim_rh_98)
MD <- mean(Data$RHS_98 - Data$GEDI_sim_rh_98)
RB <- 100*MD/mean(Data$GEDI_sim_rh_98)

p3<- ggplot(Data, aes(x=GEDI_sim_rh_98, y=RHS_98))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 15),ylim =  c(0, 15))+
  annotate("text",x=1,y=14,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=1,y=12,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"m",             
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " Bias: ", round(MD,3),"m",
             "\n" , " %Bias: ", round(RB,3),"%",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("Simulated GEDI RH98 (Leaf-on) (m)"), 
       y=expression("On-orbit GEDI RH98 (m)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#AGBD
Data <- Data[Data$ql4 ==1,]
Data$Combo <- Data$Combo/0.049
#Data$Combo <- 9.8441*(Data$MEAN*Data$CC)
#Data$Combo <- 6.85*Data$MEAN^0.95
Data <- Data[Data$Combo >0,]
reg <- lm(formula = sqrt(Combo) ~ 
            site+Beam+sensitivity_a1+elev_lowestmode_a1+
            landsat_treecover+modis_treecover+modis_nonvegetated+
            GEDI_PAI_3 + GEDI_PAVD_1 + GEDI_PAVD_2 +
            GEDI_cc_1 + GEDI_cc_2 + GEDI_FHD +
            sqrt(RH1_50) + sqrt(RH1_50)*sqrt(RH1_75) + sqrt(RH1_98), 
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



dir <- "E:\\Biomass\\CSIR\\result\\All_biom_z_L2B.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
#Data <- Data[Data$Combo !=0,]
#Data$Combo <- Data$Tree
Data <- Data[Data$status == "Leaf-off",]

r2 = round(cor(Data$RHS_98, Data$GEDI_sim_rh_98,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$RHS_98 - Data$GEDI_sim_rh_98)^2))
rRMSE <- 100*sqrt(mean((Data$RHS_98 - Data$GEDI_sim_rh_98)^2))/mean(Data$GEDI_sim_rh_98)
MD <- mean(Data$RHS_98 - Data$GEDI_sim_rh_98)
RB <- 100*MD/mean(Data$GEDI_sim_rh_98)

p4<- ggplot(Data, aes(x=GEDI_sim_rh_98, y=RHS_98))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 15),ylim =  c(0, 15))+
  annotate("text",x=1,y=14,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=1,y=12,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"m",             
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " Bias: ", round(MD,3),"m",
             "\n" , " %Bias: ", round(RB,3),"%",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("Simulated GEDI RH98 (Leaf-off) (m)"), 
       y=expression("On-orbit GEDI RH98 (m)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#AGBD
Data <- Data[Data$ql4 ==1,]
Data$Combo <- Data$Combo/0.049
#Data$Combo <- 9.8441*(Data$MEAN*Data$CC)
#Data$Combo <- 6.85*Data$MEAN^0.95
Data <- Data[Data$Combo >0,]
reg <- lm(formula = sqrt(Combo) ~ site+Beam+sensitivity_a1+elev_lowestmode_a1+
            landsat_treecover+modis_treecover+modis_nonvegetated+
            GEDI_PAI_3 + GEDI_PAVD_1 + GEDI_PAVD_2 +
            GEDI_cc_1 + GEDI_cc_2 + GEDI_FHD +
            sqrt(RH1_50) + sqrt(RH1_50)*sqrt(RH1_75) + sqrt(RH1_98), data = Data)

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
#out = "E:\\Biomass\\CSIR\\figure\\agbd_chm_mean.jpg"
#out = "E:\\Biomass\\CSIR\\figure\\agbd_H_CC.jpg"
out = "E:\\Biomass\\CSIR\\figure\\agbd_individual.jpg"
ggsave(out,height=12, width=24, dpi=600)


ggarrange(p3,p4)
out = "E:\\Biomass\\CSIR\\figure\\rh.jpg"
ggsave(out,height=12, width=24, dpi=600)



# ------------------------------------------------------------------------------------------------ #
#GEDI L4A vs. ALS CHM
dir <- "E:\\Biomass\\CSIR\\result\\All_01292023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
#Data$Combo <- Data$Tree
Data <- Data[Data$solar_elevation <0,]
Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

#Data$Combo <- 6.85*Data$MEAN^0.95
#Data$Combo <- Data$Combo/0.049
#Data$Combo <- 9.8441*(Data$MEAN*Data$CC15)
#Data$Combo <- (9.8441*(Data$MEAN*Data$CC15)/10000)*16
#Data$Combo <- 25.61*(Data$MEAN*Data$CC05)+3.95
#Data$Combo <- 23.2*(Data$MEAN*Data$CC05)-7.15
Data$Combo <-  25.8*(Data$MEAN*Data$CC15)-11.5
#Data$Combo <- 14.222*(Data$MEAN*Data$CC15)+2.0921

Data <- Data[Data$Combo >0,]
#Data$Combo <- Data$Combo*490.625/706.5

r2 = round(cor(Data$AGBD, Data$Combo,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$AGBD - Data$Combo)^2))
MD <- mean(Data$AGBD - Data$Combo)

p1<- ggplot(Data, aes(x=Combo, y=AGBD))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=15,y=80,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2
           ),parse=TRUE) + 
  annotate("text",x=15,y=70,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("Footprint level CHM derived AGB (Leaf-on) (Mg/ha)"), 
       y=expression("On-orbit GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dir <- "E:\\Biomass\\CSIR\\result\\All_biom_z_L2B.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
#Data <- Data[Data$Combo !=0,]
#Data$Combo <- Data$Tree
Data <- Data[Data$status == "Leaf-off",]
Data <- Data[Data$ql4 ==1,]
#Data$Combo <- Data$Combo/0.049
Data$Combo <- 9.8441*(Data$MEAN*Data$CC)
#Data$Combo <- 6.85*Data$MEAN^0.95
Data <- Data[Data$Combo >0,]
r2 = round(cor(Data$AGBD_1, Data$Combo,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$AGBD_1 - Data$Combo)^2))
rRMSE <- 100*sqrt(mean((Data$AGBD_1 - Data$Combo)^2))/mean(Data$Combo)
p2<- ggplot(Data, aes(x=Combo, y=AGBD_1))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 100),ylim =  c(0, 100))+
  annotate("text",x=15,y=80,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2
           ),parse=TRUE) + 
  annotate("text",x=15,y=70,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " %RMSE: ", round(rRMSE,3),"%",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("Footprint level CHM derived AGB (Leaf-off) (Mg/ha)"), 
       y=expression("On-orbit GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggarrange(p1,p2)


#out = "E:\\Biomass\\CSIR\\figure\\agbd_L4A_chm_mean.jpg"
out = "E:\\Biomass\\CSIR\\figure\\agbd_L4A_H_CC.jpg"
#out = "E:\\Biomass\\CSIR\\figure\\agbd_L4A_individual.jpg"

ggsave(out,height=12, width=24, dpi=600)



# ------------------------------------------------------------------------------------------------ #
#boxplot
#histogram of GEDI L4A
dir <- "E:\\Biomass\\CSIR\\result\\All_biom_l2b.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data = Data[Data$ql4 == 1,]
Data <- Data[Data$Combo !=0,]
Data <- Data[Data$status == "Leaf-on",]
Data$Combo <- Data$Combo/0.049
reg <- lm(formula = sqrt(Combo) ~ site+Beam+sensitivity_a1+elev_lowestmode_a1+
            landsat_treecover+modis_treecover+modis_nonvegetated+
            GEDI_PAI_3 + GEDI_PAVD_1 + GEDI_PAVD_2 +
            GEDI_cc_1 + GEDI_cc_2 + GEDI_FHD +
            sqrt(RH1_50) + sqrt(RH1_50)*sqrt(RH1_75) + sqrt(RH1_98), data = Data)
summary(reg)
Data$pred <- predict(reg, Data)^2
Data$diff <- Data$pred - Data$Combo
Data <- Data[Data$Combo <= 90,]
Data$group_RH <- cut(Data$Combo,breaks = seq(0,90, 10),
                     dig.lab = 5)
a2 <- tapply(Data$diff, cut(Data$Combo,breaks = seq(0,90, 10),
                            dig.lab = 5), mean)
b2 <- tapply(Data$Combo, cut(Data$Combo,breaks = seq(0,90, 10),
                             dig.lab=1), mean)
RB2 <- a2/b2*100
#write.table(a2,file="C:\\Users\\Shawn\\Desktop\\a2_on.csv",sep=" ")
table(Data$group_RH)
p1<- ggplot(Data, aes(x=group_RH, y=diff, group=group_RH)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 40))+
  scale_y_continuous(minor_breaks = seq(-100,40, 10),breaks = seq(-100,40, 10))+
  theme_bw()+
  labs(x = "Footprint level CHM derived AGB (Leaf-on) (Mg/ha)",
       y="Bias of multi-variate predicted GEDI AGB and CHM derived AGB (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))

dir <- "E:\\Biomass\\CSIR\\result\\All_biom_l2b.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data = Data[Data$ql4 == 1,]
Data <- Data[Data$Combo !=0,]
Data <- Data[Data$status == "Leaf-off",]
Data$Combo <- Data$Combo/0.049
reg <- lm(formula = sqrt(Combo) ~ site+Beam+sensitivity_a1+elev_lowestmode_a1+
            landsat_treecover+modis_treecover+modis_nonvegetated+
            GEDI_PAI_3 + GEDI_PAVD_1 + GEDI_PAVD_2 +
            GEDI_cc_1 + GEDI_cc_2 + GEDI_FHD +
            sqrt(RH1_50) + sqrt(RH1_50)*sqrt(RH1_75) + sqrt(RH1_98), data = Data)
summary(reg)
Data$pred <- predict(reg, Data)^2
Data$diff <- Data$pred - Data$Combo
Data <- Data[Data$Combo <= 90,]
Data$group_RH <- cut(Data$Combo,breaks = seq(0,90, 10),
                     dig.lab = 5)
a2 <- tapply(Data$diff, cut(Data$Combo,breaks = seq(0,90, 10),
                            dig.lab = 5), mean)
b2 <- tapply(Data$Combo, cut(Data$Combo,breaks = seq(0,90, 10),
                             dig.lab=1), mean)
RB2 <- a2/b2*100
#write.table(a2,file="C:\\Users\\Shawn\\Desktop\\a2_on.csv",sep=" ")
table(Data$group_RH)
p2<- ggplot(Data, aes(x=group_RH, y=diff, group=group_RH)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 40))+
  scale_y_continuous(minor_breaks = seq(-100,40, 10),breaks = seq(-100,40, 10))+
  theme_bw()+
  labs(x = "Footprint level CHM derived AGB (Leaf-off) (Mg/ha)",
       y="Bias of multi-variate predicted GEDI AGB and CHM derived AGB (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(p1,p2)
out = "E:\\Biomass\\CSIR\\figure\\Biom_bias.jpg"
ggsave(out,height=12, width=30, dpi=600)



# ------------------------------------------------------------------------------------------------ #
#Boxplot -- GEDI L4A vs. ALS CHM
dir <- "E:\\Biomass\\CSIR\\result\\All_biom_l2b.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data = Data[Data$ql4 == 1,]
Data <- Data[Data$Combo !=0,]
Data <- Data[Data$status == "Leaf-on",]
Data$Combo <- Data$Combo/0.049
Data$diff <- Data$AGBD_1 - Data$Combo
Data <- Data[Data$Combo <= 90,]
Data$group_RH <- cut(Data$Combo,breaks = seq(0,90, 10),
                     dig.lab = 5)
a2 <- tapply(Data$diff, cut(Data$Combo,breaks = seq(0,90, 10),
                            dig.lab = 5), mean)
b2 <- tapply(Data$Combo, cut(Data$Combo,breaks = seq(0,90, 10),
                             dig.lab=1), mean)
RB2 <- a2/b2*100
#write.table(a2,file="C:\\Users\\Shawn\\Desktop\\a2_on.csv",sep=" ")
table(Data$group_RH)
p1<- ggplot(Data, aes(x=group_RH, y=diff, group=group_RH)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 40))+
  scale_y_continuous(minor_breaks = seq(-100,40, 10),breaks = seq(-100,40, 10))+
  theme_bw()+
  labs(x = "Footprint level CHM derived AGB (Leaf-on) (Mg/ha)",
       y="Bias of GEDI L4A AGB and CHM derived AGB (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


dir <- "E:\\Biomass\\CSIR\\result\\All_biom_l2b.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data = Data[Data$ql4 == 1,]
Data <- Data[Data$Combo !=0,]
Data <- Data[Data$status == "Leaf-off",]
Data$Combo <- Data$Combo/0.049
Data$diff <- Data$AGBD_1 - Data$Combo
Data <- Data[Data$Combo <= 90,]
Data$group_RH <- cut(Data$Combo,breaks = seq(0,90, 10),
                     dig.lab = 5)
a2 <- tapply(Data$diff, cut(Data$Combo,breaks = seq(0,90, 10),
                            dig.lab = 5), mean)
b2 <- tapply(Data$Combo, cut(Data$Combo,breaks = seq(0,90, 10),
                             dig.lab=1), mean)
RB2 <- a2/b2*100
#write.table(a2,file="C:\\Users\\Shawn\\Desktop\\a2_on.csv",sep=" ")
table(Data$group_RH)
p2<- ggplot(Data, aes(x=group_RH, y=diff, group=group_RH)) + 
  stat_boxplot(geom ='errorbar', width = 0.3,position = position_dodge(width = 0.75)) +
  geom_boxplot(width = 0.5)+
  coord_cartesian(ylim = c(-80, 40))+
  scale_y_continuous(minor_breaks = seq(-100,40, 10),breaks = seq(-100,40, 10))+
  theme_bw()+
  labs(x = "Footprint level CHM derived AGB (Leaf-off) (Mg/ha)",
       y="Bias of GEDI L4A AGB and CHM derived AGB (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


ggarrange(p1,p2)
out = "E:\\Biomass\\CSIR\\figure\\Biom_bias_L4A.jpg"
ggsave(out,height=12, width=30, dpi=600)



# ------------------------------------------------------------------------------------------------ #
#height distribution
dir <- "E:\\Biomass\\CSIR\\result\\All_biom_z_L2B.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]
Data = Data[Data$ql4 == 1,]
#Data <- Data[Data$Combo !=0,]
Data$Combo <- Data$Combo/0.049
#Data$Combo <- 6.85*Data$MEAN^0.95
#Data$Combo <- 9.8441*(Data$MEAN*Data$CC)

Data <- Data[Data$Combo >0,]
Data <- Data[Data$status == "Leaf-on",]
Data$group_RH <- cut(Data$Combo,breaks = seq(0,90, 10),
                     dig.lab = 5)
table(Data$group_RH)
ggplot(Data, aes(x=group_RH)) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Footprint level CHM derived AGB (Mg/ha)",
       y = "Frequency (%)")+
  geom_histogram(aes(y = stat(count) / sum(count)),
                 color =  "white",fill = "grey",stat="count",alpha = 0.8) +
  scale_x_discrete(labels=c("[0,10]","[10,20]","[20,30]","[30,40]",
                            "[40,50]","[50,60]","[60,70]","[70,80]",
                            "[80,90]","[90,100]"))+
  scale_y_continuous(labels = scales::percent)+
  theme(text=element_text(family="A"))+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25),axis.text=element_text(size=25))

out = "E:\\Biomass\\CSIR\\figure\\Biom_hist.jpg"
out = "E:\\Biomass\\CSIR\\figure\\Biom_hist_chm_mean.jpg"
out = "E:\\Biomass\\CSIR\\figure\\Biom_hist_H_CC.jpg"

ggsave(out,height=12, width=30, dpi=600)





#1/26/2023
# ------------------------------------------------------------------------------------------------ #
#CSIR field vs. ALS CHM AGB
filedir <- "E:\\Biomass\\CSIR\\result\\CSIR\\biom_CSIR_combo.csv"
Data = read.csv(filedir)
Data<-Data[!(Data$Name1=="Ven_S68_SP_5_Venetia_h.tif"),]
Data$ALS = Data$Tree
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
Mean <- mean(Data$CSIR)       
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
MD <- mean(Data$ALS - Data$CSIR)
p1<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  annotate("text",x=0,y=148,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=128,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " Bias: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  coord_cartesian(xlim = c(0, 150),ylim =  c(0, 150))+
  labs(x=expression("1 ha CSIR AGB (Mg/ha)"), 
       y=expression("Individual tree based AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


filedir <- "E:\\Biomass\\CSIR\\result\\CSIR\\biom_CSIR_combo.csv"
Data = read.csv(filedir)
Data<-Data[!(Data$Name1=="Ven_S68_SP_5_Venetia_h.tif"),]
Data$ALS = Data$Area
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
Mean <- mean(Data$CSIR)       
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
MD <- mean(Data$ALS - Data$CSIR)

p2<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  annotate("text",x=0,y=148,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=128,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " Bias: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  coord_cartesian(xlim = c(0, 150),ylim =  c(0, 150))+
  labs(x=expression("1 ha CSIR AGB (Mg/ha)"), 
       y=expression("Area based AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#scatterplot csir 25m vs. ALS MCH biomass
filedir <- "E:\\Biomass\\CSIR\\result\\CSIR_25m\\biom_CSIR_25m_combo_1_h_cc.csv"
Data = read.csv(filedir)
Data <- na.omit(Data)
Data$ALS <- Data$Ind
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
MD <- mean(Data$ALS - Data$CSIR)

p3<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  annotate("text",x=0,y=148,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=128,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " Bias: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  coord_cartesian(xlim = c(0, 150),ylim =  c(0, 150))+
  labs(x=expression("25 m CSIR AGB (Mg/ha)"), 
       y=expression("Individual tree based AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())





filedir <- "E:\\Biomass\\CSIR\\result\\CSIR_25m\\biom_CSIR_25m_combo_1_h_cc.csv"
Data = read.csv(filedir)
Data <- na.omit(Data)
#Data$ALS <- Data$MCH*625/706.5
#Data$ALS <-  6.85*Data$MEAN^0.95
Data$ALS <-  Data$Ind
Data$ALS <-  5.5*Data$MEAN*Data$CC + 7.434
Data$ALS <-  5.2*Data$MEAN^1.06
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np13"),]
Data<-Data[!(Data$plot=="Ven_S63_SP_5_np15"),]
Data <- Data[Data$CSIR >0 & Data$CSIR < 60,]
#stats
r2 = round(cor(Data$ALS, Data$CSIR,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$ALS - Data$CSIR)^2))
MD <- mean(Data$ALS - Data$CSIR)

p4<- ggplot(Data, aes(x=CSIR, y=ALS))+ 
  geom_abline(intercept = 0, slope = 1,color="gray", linetype="dashed", size=1.5)+
  geom_pointdensity(color="blue",shape=19, size=3, alpha=0.5)+
  theme_bw()+
  annotate("text",x=0,y=70,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=0,y=63,hjust = 0,size = 8,family= "A",
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " Bias: ", round(MD,3)," Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  coord_cartesian(xlim = c(0, 70),ylim =  c(0, 70))+
  labs(x=expression("25 m CSIR AGB (Mg/ha)"), 
       y=expression("Area based AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



ggarrange(p1,p3,p2,p4,nrow = 2,ncol = 2)
out = "E:\\Biomass\\CSIR\\figure\\CSIR_field_CHM.jpg"
ggsave(out,height=12, width=12, dpi=600)



#1/29/2023
# ------------------------------------------------------------------------------------------------ #
#GEDI L4A vs. ALS CHM -- individual tree method
dir <- "E:\\Biomass\\CSIR\\result\\All_01292023.csv"
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
r2 = round(cor(Data$AGBD, Data$Combo,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$AGBD - Data$Combo)^2))
MD <- mean(Data$AGBD - Data$Combo)

p1<- ggplot(Data, aes(x=Combo, y=AGBD))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 50),ylim =  c(0, 50))+
  annotate("text",x=15,y=50,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=15,y=45,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Individual-based CHM-derived AGB (Mg/ha)"), 
       y=expression("GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



#GEDI L4A vs. ALS CHM -- area based method 
dir <- "E:\\Biomass\\CSIR\\result\\All_01292023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)

#Data <- Data[Data$solar_elevation <0,]
#Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

#Data$Combo <- 2.02*Data$MEAN^1.8
Data$Combo <- 6.85*Data$MEAN^0.95

#Data$Combo <- 25.8*Data$MEAN*Data$CC15-11.5
#Data$Combo <- 11.72*Data$MEAN*Data$CC15-8.165
#Data <- Data[Data$CC15 > 0.1,]

Data <- Data[Data$Combo >0,]
Data$Combo <- Data$Combo*490.625/625

r2 = round(cor(Data$AGBD, Data$Combo,method = "pearson")^2,3)
RMSE <- sqrt(mean((Data$AGBD - Data$Combo)^2))
MD <- mean(Data$AGBD - Data$Combo)

p2 <- ggplot(Data, aes(x=Combo, y=AGBD))+ 
  geom_pointdensity()+
  scale_color_viridis()+
  theme_bw()+
  coord_cartesian(xlim = c(0, 50),ylim =  c(0, 50))+
  annotate("text",x=15,y=50,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",r2 
           ),parse=TRUE) + 
  annotate("text",x=15,y=45,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Area-based CHM-derived AGB (Mg/ha)"), 
       y=expression("GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggarrange(p1,p2)

out = "E:\\Biomass\\CSIR\\figure\\L4A_ALS_original.jpg"
ggsave(out,height=12, width=24, dpi=600)





# ------------------------------------------------------------------------------------------------ #
#improved GEDI L4A vs. ALS CHM -- individual tree method 
dir <- "E:\\Biomass\\CSIR\\result\\All_01292023.csv"
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
model <- train(sqrt(Combo) ~
                 fhd_normal+pgap_theta_z+elev_lowestmode_a1+
                 GEDI_cc_1 + GEDI_cc_2 + GEDI_cc_3 + 
                 sqrt(RH1_98), 
               data = Data, 
               method = "lm", 
               trControl = train_control)

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
  coord_cartesian(xlim = c(0, 50),ylim =  c(0, 50))+
  annotate("text",x=15,y=50,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(r2,3) 
           ),parse=TRUE) + 
  annotate("text",x=15,y=45,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Area-based CHM-derived AGB (Mg/ha)"), 
       y=expression("GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




#test multicollinearity
reg <- lm(formula = sqrt(Combo) ~ 
            fhd_normal+pgap_theta_z+elev_lowestmode_a1+
            GEDI_cc_1 + GEDI_cc_2 + GEDI_cc_3 + 
            sqrt(RH1_98), 
          data = Data)
summary(reg)
vif(reg)
p = Data[, c("fhd_normal", "pgap_theta_z", "elev_lowestmode_a1", 
             "GEDI_cc_1", "GEDI_cc_2", "GEDI_cc_3", "RH1_98")]
cor(p)


#test importance
spit = c(train = .8, test = .2)
gs = sample(cut(seq(nrow(Data)), nrow(Data)*cumsum(c(0,spit)),labels = names(spit)))
res = split(Data, gs)
sapply(res, nrow)/nrow(Data)
rf <- randomForest(Combo ~ 
                     fhd_normal+pgap_theta_z+elev_lowestmode_a1+
                     GEDI_cc_1 + GEDI_cc_2 + GEDI_cc_3 + 
                     RHS_98, 
                   data = res$train,ntree = 100,importance=TRUE, nodesize = 20)
importance(rf)
varImpPlot(rf)
varImp(rf)




#GEDI L4A vs. ALS CHM -- area based method 
dir <- "E:\\Biomass\\CSIR\\result\\All_01292023.csv"
Data = read.csv(dir,header=T)
Data <- na.omit(Data)

#Data <- Data[Data$solar_elevation <0,]
#Data <- Data[Data$status == "Leaf-on",]
Data <- Data[Data$ql4 ==1,]
Data <- Data[Data$ql2a ==1,]
Data <- Data[Data$ql2b ==1,]

Data$Combo <- 6.85*Data$MEAN^0.95

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

model <- train(sqrt(Combo) ~
                 fhd_normal+pgap_theta_z+elev_lowestmode_a1+
                 GEDI_cc_1 + GEDI_cc_2 + GEDI_cc_3 + 
                 sqrt(RH1_98), 
               data = Data, 
               method = "lm", 
               trControl = train_control)

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
  coord_cartesian(xlim = c(0, 50),ylim =  c(0, 50))+
  annotate("text",x=15,y=50,hjust = 0,size = 8,family= "A",
           label= paste(expression(" "~R^2),": ",round(r2,3) 
           ),parse=TRUE) + 
  annotate("text",x=15,y=45,hjust = 0,size = 8,family= "A",
           label= paste(
             "  RMSE: ", round(RMSE,3),"Mg/ha",             
             "\n" , " Bias: ", round(MD,3),"Mg/ha",
             "\n" , " Sample size: ", nrow(Data))) + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(x=expression("25 m Area-based CHM-derived AGB (Mg/ha)"), 
       y=expression("GEDI L4A AGB (Mg/ha)"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



ggarrange(p1,p2)

out = "E:\\Biomass\\CSIR\\figure\\L4A_ALS_improved.jpg"
ggsave(out,height=12, width=24, dpi=600)

