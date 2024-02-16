library(lidR)
library(raster)

# ------------------------------------------------------------------------------------------------ #
#new method ---individual tree process --- test single ras 1/28/2023, CSIR 25m
dir = ""
out_dir <- file.path("")
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
  #writeRaster(crowns, file.path(out_dir, paste0(file_path_sans_ext(basename(filename)), ".tif")), overwrite=TRUE)
}
