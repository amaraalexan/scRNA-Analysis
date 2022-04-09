###Reading in Data -----

options(stringsAsFactors = F) 
#Keep strings intact
setwd("ZFBrain_Analysis/Input/ZFdata")
#personalize

folders <- list.files("ZFBrain_Analysis/Input/ZFdata")
library(Seurat)
library(dplyr)
zfbrainList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
zfbrain.combined <- merge(zfbrainList[[1]], 
                 y = c(zfbrainList[[2]],zfbrainList[[3]],zfbrainList[[4]],zfbrainList[[5]],
                       zfbrainList[[6]],zfbrainList[[7]],zfbrainList[[8]],zfbrainList[[9]],
                       zfbrainList[[10]],zfbrainList[[11]],zfbrainList[[12]],zfbrainList[[13]],
                       zfbrainList[[14]],zfbrainList[[15]]), 
                 add.cell.ids = folders, 
                 project = "zebrafish")
