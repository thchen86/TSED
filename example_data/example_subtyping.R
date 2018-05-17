source("D:/Tools/TSED/expression_data_tumor_subtyping.R")
source("D:/Tools/TSED/centroid_clustering.R")
setwd("D:/Tools/TSED/example_data")


indication <- "OV"
inFile<-paste0("D:/Tools/TSED/example_data/OV_log2TPM_for_classification.txt")
datalib_dir<- "D:/Tools/TSED"

tumor_subtyping(indication=indication,inFile=inFile,prefix=NA,
                  datalib_dir=datalib_dir)



