

####---- User Input ----####

Project_Directory <- "Example_Data/"

Project_Name_List <- "CHOL,ACC"

FilterOut_Below <- 1

SampleN_FilterMin <- 10






####---- Run Script ----####

####---- R Libraries ----####
library(recount3)
library(readr)
library(dplyr)
library(data.table)
library(slam)
library(stringr)


options(timeout = max(1000, getOption("timeout")))

if (length(Project_Name_List) == 1) {
Project_Name_List <- strsplit(Project_Name_List,",")[[1]]
}

options("recount3_organism_human_project_homes_URL_http://duffel.rail.bio/recount3" = c("data_sources/sra", "data_sources/gtex", "data_sources/tcga"))

load(url("http://shawlab.science/shiny/GTEx_AMLBM_JunctionList/GTEx_AMLBM_JxnPresent_Vector.RData"))

recount3_proj <- available_projects()

last_char <- str_sub(Project_Directory,-1,-1)
if (last_char != "/") {
  Project_Directory <- paste(Project_Directory,"/",sep = "")
}

Head_Proj_Dir <- Project_Directory
if (!dir.exists(Head_Proj_Dir)) {
  dir.create(Head_Proj_Dir,recursive = T)
}

Cache_dir <- paste0(Head_Proj_Dir,"/Recount3_Cache/")
if (!dir.exists(Cache_dir)) {
  dir.create(Cache_dir,recursive = T)
}

if (exists("FilterOut_Below")) {
  if (any(FilterOut_Below == 0, is.na(FilterOut_Below), is.null(FilterOut_Below))) {
    FilterOut_Below <- NULL
  }
} else {FilterOut_Below <- NULL}

if (exists("SampleN_FilterMin")) {
  if (any(is.na(SampleN_FilterMin), is.null(SampleN_FilterMin))) {
    FilterOut_Below <- NULL
  }
} else {FilterOut_Below <- NULL}



dgTMatrix_norm <- function(m) {
  b <- rle(m@j)
  dt <- data.frame(number = b$values, lengths = b$lengths)
  dt$end <- cumsum(dt$lengths)
  dt$start <- dt$end - dt$lengths + 1
  xnew <- c()
  for (i in unique(dt$number)) {
    idxStart <- dt[which(dt$number==i),"start"]
    idxEnd <- dt[which(dt$number==i),"end"]
    xSub <- m@x[idxStart:idxEnd]
    xSub_norm <- xSub/sum(xSub)
    xSub_normx <- xSub_norm*1000000
    xnew <- c(xnew,xSub_normx)
  }
  m@x <- xnew
  return(m)
}

dgTMatrix_filter <- function(m,criteria,number) {
  JxnExpr <- function(vec,criteria,number) {
    meet <- sum(vec>criteria)
    if (meet > number) {
      return(TRUE)
    } else { return(FALSE) }
  }
  if (!is.null(criteria) & !is.null(number)) {
    rowsFiltered <- rowapply_simple_triplet_matrix(as.simple_triplet_matrix(m), function(x) JxnExpr(x,criteria,number))
    m_filter <- m[names(rowsFiltered)[which(rowsFiltered == TRUE)],]
  }
  return(m_filter)
}

Jxn2Bed <- function(JxnFile,AnnoFile,OutFile) {
  jxn <- JxnFile
  nJxnCol <- ncol(jxn)
  
  anno <- AnnoFile
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
  
  write_delim(bed_AvgExpr_out,paste0(OutFile,"_AvgExpr.bed"), col_names = F, delim = '\t')
  write_delim(bed_SumExpr_out,paste0(OutFile,"_TotalSampExpr.bed"), col_names = F, delim = '\t')
}

for (i in Project_Name_List) {
  
  Proj_home <- recount3_proj[which(recount3_proj$project == i),"project_home"]
  Proj_org <- recount3_proj[which(recount3_proj$project == i),"organism"]
  Proj_source <- toupper(recount3_proj[which(recount3_proj$project == i),"file_source"])
  File_prefix <- paste0(Proj_source,"_",i)
  
  ## Prep output folder
  Sub_Proj_dir <- paste0(Head_Proj_Dir,File_prefix,"/")
  if (!dir.exists(Sub_Proj_dir)) {
    dir.create(Sub_Proj_dir,recursive = T)
  }
  write(paste(i,"Junction Filtering Log File"),
        file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
  
  ## Load in Recount3 Data
  print(paste0("Loading Project ",File_prefix))
  temp_rse_jxn <- recount3::create_rse_manual(
    project = i,
    project_home = Proj_home,
    organism = Proj_org,
    type = "jxn",
    #annotation = "gencode_v29",
    verbose = FALSE,
    bfc = recount3_cache(Cache_dir)
  )
  
  ## Junction Data
  mat <- assay(temp_rse_jxn)
  temp_meta <- as.data.frame(t(do.call(rbind,temp_rse_jxn@colData@listData)))
  print(paste0(Sys.time(),": Sarting Junction Number: ",nrow(mat)))
  write(paste0(Sys.time(),": Sarting junction number: ",nrow(mat)),
        file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
  
  if (SampleN_FilterMin == 10) {
    if (nrow(temp_meta) <= SampleN_FilterMin) {
      SampleN_FilterMin <- 3
    }
  }
  
  JunctionData <- as.data.frame(t(do.call(rbind, temp_rse_jxn@rowRanges@elementMetadata@listData)))
  rownames(JunctionData) <- rownames(mat)
  
  ## Meta Data
  write(paste0(Sys.time(),": Writing out Clinical Data: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_Meta.txt")),
        file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
  write_delim(temp_meta,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_Meta.txt"), delim = '\t')
  
  ## Normalizing
  print(paste0(Sys.time(), ": Normalizing"))
  mat_norm <- dgTMatrix_norm(mat)
  
  ## Filtering based on GTEx
  print(paste0(Sys.time(), ": Filtering based on presense in GTEx"))
  KeepRows <- row.names(mat_norm)[which(!row.names(mat_norm) %in% Normal_Jxn_list)]
  
  if (length(KeepRows) > 0) {
    
    mat_norm_noGTEx <- mat_norm[KeepRows,]
    print(paste0(Sys.time(),": Number of junctions not present in GTEx: ",nrow(mat_norm_noGTEx)))
    write(paste0(Sys.time(),": Number of junctions not present in GTEx: ",nrow(mat_norm_noGTEx)),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    mat_norm_noGTEx_df <- as.data.frame(as.matrix(mat_norm_noGTEx))
    mat_norm_noGTEx_df <- tibble::rownames_to_column(mat_norm_noGTEx_df,var = "Junction")
    write(paste0(Sys.time(),": Writing out Junctions filtered by GTEx: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM.txt")),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    write_delim(mat_norm_noGTEx_df,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM.txt"), delim = '\t')
    
    JunctionData_filtered <- JunctionData[rownames(mat_norm_noGTEx),]
    JunctionData_filtered <- tibble::rownames_to_column(JunctionData_filtered,var = "Junction")
    write(paste0(Sys.time(),": Writing out filtered junction annotation: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_Anno.txt")),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    write_delim(JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_Anno.txt"), delim = '\t')
    
    #Jxn2Bed(mat_norm_noGTEx_df,JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM"))
    
    ## Remove Junctions of all 0 counts
    print(paste0(Sys.time(), ": Filtering based on row sum > 0"))
    mat_norm_noGTEx_nz <- mat_norm_noGTEx[which(rowSums(mat_norm_noGTEx) > 0),]
    mat_norm_noGTEx_nz_df <- as.data.frame(as.matrix(mat_norm_noGTEx_nz))
    mat_norm_noGTEx_nz_df <- tibble::rownames_to_column(mat_norm_noGTEx_nz_df,var = "Junction")
    print(paste0(Sys.time(),": Number of junctions with a row sum > 0: ",nrow(mat_norm_noGTEx_nz_df)))
    write(paste0(Sys.time(),": Number of junctions with a row sum > 0: ",nrow(mat_norm_noGTEx_nz_df)),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    write(paste0(Sys.time(),": Writing out Junctions with a row sum > 0: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR0.txt")),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    write_delim(mat_norm_noGTEx_nz_df,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR0.txt"), delim = '\t')
    
    JunctionData_filtered <- JunctionData[rownames(mat_norm_noGTEx_nz),]
    JunctionData_filtered <- tibble::rownames_to_column(JunctionData_filtered,var = "Junction")
    write(paste0(Sys.time(),": Writing out filtered junction annotation: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR0_Anno.txt")),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    write_delim(JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR0_Anno.txt"), delim = '\t')
    
    #Jxn2Bed(mat_norm_noGTEx_nz_df,JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR0"))
    
    ## Filtering based on criteria
    if (!is.null(FilterOut_Below) & !is.null(SampleN_FilterMin)) {
      print(paste0(Sys.time(), ": Filtering based on criteria"))
      mat_norm_noGTEx_nz_filter <- dgTMatrix_filter(mat_norm_noGTEx_nz,FilterOut_Below,SampleN_FilterMin)
      if (nrow(mat_norm_noGTEx_nz_filter > 0)) {
        mat_norm_noGTEx_nz_filter_df <- as.data.frame(as.matrix(mat_norm_noGTEx_nz_filter))
        mat_norm_noGTEx_nz_filter_df <- tibble::rownames_to_column(mat_norm_noGTEx_nz_filter_df,var = "Junction")
        print(paste0(Sys.time(),": Number of junctions with a CPM > ",FilterOut_Below," in more than ",SampleN_FilterMin," samples: ",nrow(mat_norm_noGTEx_nz_filter_df)))
        write(paste0(Sys.time(),": Number of junctions with a CPM > ",FilterOut_Below," in more than ",SampleN_FilterMin," samples: ",nrow(mat_norm_noGTEx_nz_filter_df)),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
        write(paste0(Sys.time(),": Writing out Junctions filtered by criteria: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR",FilterOut_Below,"in",SampleN_FilterMin,".txt")),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
        write_delim(mat_norm_noGTEx_nz_filter_df,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR",FilterOut_Below,"in",SampleN_FilterMin,".txt"), delim = '\t')
        
        JunctionData_filtered <- JunctionData[rownames(mat_norm_noGTEx_nz_filter),]
        JunctionData_filtered <- tibble::rownames_to_column(JunctionData_filtered,var = "Junction")
        write(paste0(Sys.time(),": Writing out filtered junction annotation: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR",FilterOut_Below,"in",SampleN_FilterMin,"_anno.txt")),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
        write_delim(JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR",FilterOut_Below,"in",SampleN_FilterMin,"_anno.txt"), delim = '\t')
        
        #Jxn2Bed(mat_norm_noGTEx_nz_filter_df,JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_noGTEx_noBM_GRTR",FilterOut_Below,"in",SampleN_FilterMin))
      } else {
        print(paste0(Sys.time(),": Number of junctions with a CPM > ",FilterOut_Below," in more than ",SampleN_FilterMin," samples: ",nrow(mat_norm_noGTEx_nz_filter_df)))
        write(paste0(Sys.time(),": Number of junctions with a CPM > ",FilterOut_Below," in more than ",SampleN_FilterMin," samples: ",nrow(mat_norm_noGTEx_nz_filter_df)),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
        print(paste0(Sys.time(),": No Junctions written to file"))
        write(paste0(Sys.time(),": No Junctions written to file"),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
      }
      
    }
  } else {
    
    ## Let user know all junctions found in GTEx
    print(paste0(Sys.time(),": Number of junctions not present in GTEx: 0"))
    write(paste0(Sys.time(),": Number of junctions not present in GTEx: 0"),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    
    ## Remove Junctions of all 0 counts
    print(paste0(Sys.time(), ": Filtering based on row sum > 0"))
    mat_norm_nz <- mat_norm[which(rowSums(mat_norm) > 0),]
    mat_norm_nz_df <- as.data.frame(as.matrix(mat_norm_nz))
    mat_norm_nz_df <- tibble::rownames_to_column(mat_norm_nz_df,var = "Junction")
    print(paste0(Sys.time(),": Number of junctions with a row sum > 0: ",nrow(mat_norm_nz_df)))
    write(paste0(Sys.time(),": Number of junctions with a row sum > 0: ",nrow(mat_norm_nz_df)),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    write(paste0(Sys.time(),": Writing out Junctions with a row sum > 0: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR0.txt")),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    write_delim(mat_norm_nz_df,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR0.txt"), delim = '\t')
    
    JunctionData_filtered <- JunctionData[rownames(mat_norm_nz),]
    JunctionData_filtered <- tibble::rownames_to_column(JunctionData_filtered,var = "Junction")
    write(paste0(Sys.time(),": Writing out filtered junction annotation: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR0_Anno.txt")),
          file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
    write_delim(JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR0_Anno.txt"), delim = '\t')
    
    #Jxn2Bed(mat_norm_nz_df,JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR0"))
    
    ## Filtering based on criteria
    if (!is.null(FilterOut_Below) & !is.null(SampleN_FilterMin)) {
      print(paste0(Sys.time(), ": Filtering based on criteria"))
      mat_norm_nz_filter <- dgTMatrix_filter(mat_norm_nz,FilterOut_Below,SampleN_FilterMin)
      if (nrow(mat_norm_nz_filter > 0)) {
        mat_norm_nz_filter_df <- as.data.frame(as.matrix(mat_norm_nz_filter))
        mat_norm_nz_filter_df <- tibble::rownames_to_column(mat_norm_nz_filter_df,var = "Junction")
        print(paste0(Sys.time(),": Number of junctions with a CPM > ",FilterOut_Below," in more than ",SampleN_FilterMin," samples: ",nrow(mat_norm_nz_filter_df)))
        write(paste0(Sys.time(),": Number of junctions with a CPM > ",FilterOut_Below," in more than ",SampleN_FilterMin," samples: ",nrow(mat_norm_nz_filter_df)),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
        write(paste0(Sys.time(),": Writing out Junctions filtered by criteria: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR",FilterOut_Below,"in",SampleN_FilterMin,".txt")),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
        write_delim(mat_norm_nz_filter_df,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR",FilterOut_Below,"in",SampleN_FilterMin,".txt"), delim = '\t')
        
        JunctionData_filtered <- JunctionData[rownames(mat_norm_nz_filter),]
        JunctionData_filtered <- tibble::rownames_to_column(JunctionData_filtered,var = "Junction")
        write(paste0(Sys.time(),": Writing out filtered junction annotation: ",paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR",FilterOut_Below,"in",SampleN_FilterMin,"_Anno.txt")),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
        write_delim(JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR",FilterOut_Below,"in",SampleN_FilterMin,"_Anno.txt"), delim = '\t')
        
        Jxn2Bed(mat_norm_nz_filter_df,JunctionData_filtered,paste0(Sub_Proj_dir,File_prefix,"_SpliceJxn_GRTR",FilterOut_Below,"in",SampleN_FilterMin))
      } else {
        print(paste0(Sys.time(),": Number of junctions with a CPM > ",FilterOut_Below," in more than ",SampleN_FilterMin," samples: ",nrow(mat_norm_nz_filter_df)))
        write(paste0(Sys.time(),": Number of junctions with a CPM > ",FilterOut_Below," in more than ",SampleN_FilterMin," samples: ",nrow(mat_norm_nz_filter_df)),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
        print(paste0(Sys.time(),": No Junctions written to file"))
        write(paste0(Sys.time(),": No Junctions written to file"),
              file = paste0(Sub_Proj_dir,File_prefix,"_JunctionFilter_Log.txt"),append = T, ncolumns = 1)
      }
      
    }
  }
  
  
}


