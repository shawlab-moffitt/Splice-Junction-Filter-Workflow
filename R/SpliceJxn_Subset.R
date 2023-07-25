
####---- User Input ----####

Junction_File <- "TCGA_CHOL_SpliceJxn_noGTEx_noBM_GRTR0.txt"

Meta_File <- "TCGA_CHOL_SpliceJxn_Meta.txt"

Column_Name <- "tcga.cgc_sample_sample_type"





####---- Run Script ----####

library(readr)

jxn <- as.data.frame(read_delim(Junction_File, delim = '\t', col_names = T))
met <- as.data.frame(read_delim(Meta_File, delim = '\t', col_names = T))
  
jxnFilePre <- gsub(".txt","",Junction_File)
MetFilePre <- gsub(".txt","",Meta_File)
  
for (j in unique(met[,Column_Name])) {
  type <- gsub(" ","",j)
  metSub <- met[which(met[,Column_Name] == j),]
  jxnSub <- jxn[,c("Junction",metSub[,"external_id"]), drop = F]
  write_delim(jxnSub,paste0(jxnFilePre,"_",type,".txt"), delim = '\t')
  write_delim(metSub,paste0(MetFilePre,"_",type,".txt"), delim = '\t')
}


