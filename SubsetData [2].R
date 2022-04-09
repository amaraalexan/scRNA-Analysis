 library(Seurat)
 library(dplyr)
#Subset
pdyn <- subset(x = zfbrain.combined, pdyn >0)

#Subset Further
crhb_only <- subset(pdyn, crhb >0 & avp ==0)
crhb_avp <- subset(pdyn, crhb >0 & avp > 0 )
avp_only <- subset(pdyn, crhb == 0 & avp >0)
pdyn_only <- subset(pdyn, avp ==0 & crhb == 0)