
####---- User Input ----####

JxnFile <- "Example_Data/TCGA_CHOL_SpliceJxn_noGTEx_noBM_GRTR0.txt"

AnnoFile <- "Example_Data/TCGA_CHOL_SpliceJxn_noGTEx_noBM_GRTR0_Anno.txt"

OutFile <- "Example_Data/TCGA_CHOL_SpliceJxn_noGTEx_noBM_GRTR0"




####---- Run Script ----####

library(readr)
library(tidyr)
library(stringr)


jxn <- as.data.frame(read_delim(JxnFile,delim = '\t', col_names = T))
nJxnCol <- ncol(jxn)

if (exists("AnnoFile")) {
  if (file.exists(AnnoFile)) {
    anno <- as.data.frame(read_delim(AnnoFile,delim = '\t', col_names = T))
  } else {anno <- NULL}
} else {anno <- NULL}

bed <- jxn[,1,drop = F]
bed2 <- str_split_fixed(bed$Junction,":",3)
bed3 <- str_split_fixed(bed2[,2],"-",2)
bed_split <- data.frame(Junction = bed[,1],
                        chr = bed2[,1],
                        start = bed3[,1],
                        stop = bed3[,2],
                        strand = bed2[,3])
bed <- bed_split

if (!is.null(anno)) {
  bed <- merge(bed,anno[,c(1,3)], sort = F)
  bed$annotated <- ifelse(1,"known","unknown")
} else {bed$annotated <- "known"}

ExprFilter <- function(vec,start,stop) {
  vec <- vec[c(start:stop)]
  vec <- as.numeric(vec)
  meet <- sum(vec>0)
  return(meet)
}

jxn$AverageExpr <- rowMeans(jxn[,-1], na.rm = T)
jxn$SumExpr <- apply(jxn,1,function(x) ExprFilter(x,2,nJxnCol))

bed_AvgExpr <- merge(bed,jxn[,c("Junction","AverageExpr")], sort = F)
bed_SumExpr <- merge(bed,jxn[,c("Junction","SumExpr")], sort = F)

bed_AvgExpr_out <- bed_AvgExpr[,-1]
bed_SumExpr_out <- bed_SumExpr[,-1]

if (exists("OutFile")) {
  if (grepl(".bed$",OutFile)) {
    gsub(".bed$","",OutFile)
  }
} else {OutFile <- gsub(".txt$","",JxnFile)}

write_delim(bed_AvgExpr_out,paste0(OutFile,"_AvgExpr.bed"), col_names = F, delim = '\t')
write_delim(bed_SumExpr_out,paste0(OutFile,"_TotalSampExpr.bed"), col_names = F, delim = '\t')








