###############################################
#### calculation: [min_dist]-based scores  ####
#### - using SCimilarity's output          ####
#### - score name: [ cluster_label_score ] ####
###############################################
# version: 0.1.3 (2024-12-13)
# by hyojinsong
# ---

## [SCimilarity] - by Genentech ##
## - source: [https://github.com/Genentech/scimilarity] #GitHub
## - source: [https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html] #tutorial

#### directory ####
## set working directory ####
## - note: need to customize the directory below (i.e., path to save "final" output files)
mydir <- "/Users/your_username/" #need-to-customize
setwd(mydir)

#### project ####
## - note: [prj] variable will be included in the final output files
prj <- "TEST" #need-to-customize

#### input file ####
## - note: either "absolute" or "relative" file path can be used
fn_obs <- paste0(mydir,"obs.csv") #[relative] file path
fn_obs <- "~/obs.csv" #[absolute] file path

# ===
########################################
#### function: [calc_score_minDist] ####
########################################
calc_score_minDist <- function(prj, fn_obs, save_txt, save_rds) {
  print("################################")
  print(paste0("#### [ calc_score_minDist ] ####"))
  print("################################")
  
  print(paste0("=== [BEGIN] run 'calc_score_minDist': ", Sys.time()," ==="))
  
  #### packages ####
  mypackages <- c(
    "plyr", "dplyr"
  )

  # install.packages("plyr", dependencies = T) #if needed
  # install.packages("dplyr", dependencies = T) #if needed
  
  lapply(mypackages, library, character.only = TRUE)
  # ---
  
  
  # ---
  ## pre-step | extract <obs.csv> file from the <.h5ad> file SCimilarity "constrained" prediction ####
  # adata.write_csvs(dirname=OUTDIR_CSV, skip_data=True, sep=',') #python
  # ---
  
  
  #---
  ## step 0 | load the <obs.csv> file as a data frame ####
  ## - note: required to use the <obs.csv> file generated from the "unconstrained" prediction - as this step will calculate the QC metrics including the [min_dist] column.
  fn_obs_c <- fn_obs
  obs_c <- read.csv(fn_obs_c, header = T, row.names = 1)
  print(head(obs_c,3))
  print(paste0("=== INPUT [obs_c]: nRow: ",nrow(obs_c)," | nColumn: ",ncol(obs_c)," ===")) #checkpoint
  print(":")
  # ---
  
  # ---
  ## step 1 | summarize the [min_dist] values -  by: [seurat_clusters],[sample],[expr_GCAR],[celltype_hint] ####
  ## - required columns: [min_dist], [seurat_clusters], [celltype_hint], [cellProp_byClst] ##
  print(paste0("[BEGIN] run 'step 1' summarize [min_dist] values: ", Sys.time()))
  sum_mindist_byClst <- ddply(obs_c, 
                              c("seurat_clusters","celltype_hint"), 
                              summarise,
                              # sum_minDist=sum(min_dist), #optional
                              max_minDist=max(min_dist),
                              # min_minDist=min(min_dist), 
                              median_minDist=median(min_dist)
                              # mean_minDist=mean(min_dist) #optional
  )
  print(paste0("--- [obs_c]: nRow: ",nrow(obs_c)," | nColumn: ",ncol(obs_c)," ---")) #checkpoint
  print(paste("[END] run 'step 1' summarize [min_dist] values: ", Sys.time(), sep = ""))
  print(":")
  # ---
  
  # ---
  ## step 2 | calculate key variables "required" to calculate the [minDistScore] ####
  ## - the key variables include: [nCell],[totCell_byClst],[cellProp_byClst]
  print(paste("[BEGIN] run 'step 2' calculate key variables [min_dist] values: ", Sys.time(), sep = ""))
  
  ## step 2a: make a summary cell prediction data frame based on the <obs_c> data frame ####
  sum_cellpred_byClst <- as.data.frame(table(obs_c$celltype_hint,obs_c$seurat_clusters))
  print(head(sum_cellpred_byClst),3)
  print(dim(sum_cellpred_byClst))
  ## exclude rows w/ [Freq]=zero
  sum_cellpred_byClst <- subset(sum_cellpred_byClst, Freq>0)
  print(dim(sum_cellpred_byClst))
  ## rename the columns
  colnames(sum_cellpred_byClst) <- c("celltype_hint","seurat_clusters","nCell")
  print(head(sum_cellpred_byClst),3)
  
  ## step 2b: calculate the "total number of cells" per cell cluster: [totCell_byClst] ####
  sum_cellpred_byClst <- sum_cellpred_byClst %>%
    group_by(seurat_clusters) %>%
    mutate(
      totCell_byClst=sum(nCell)
    ) %>%
    as.data.frame()
  head(sum_cellpred_byClst,3)
  tail(sum_cellpred_byClst,3)
  
  ## step 2c: calculate the cell proportions of each "predicted" cell types per cell cluster: [cellProp_byClst] ####
  sum_cellpred_byClst <- sum_cellpred_byClst %>%
    group_by(seurat_clusters) %>%
    mutate(
      cellProp_byClst=(nCell/totCell_byClst)*100
    ) %>%
    as.data.frame()
  head(sum_cellpred_byClst,3)
  tail(sum_cellpred_byClst,3)
  # subset(sum_cellpred_byClst, seurat_clusters==18) #checkpoint
  
  ## step 2d: add an identifier column: [id__cluster_celltype] ####
  sum_cellpred_byClst$id__cluster_celltype <- paste(sum_cellpred_byClst$seurat_clusters,sum_cellpred_byClst$celltype_hint,sep="_")
  
  print(paste("[DONE] run 'step 2' calculate key variables [min_dist] values: ", Sys.time(), sep = ""))
  print(":")
  # ---
  
  # ---
  ## step 3 | calculate the [score1_minDist_median] & [score1a_minDist_median_cellProp] ####
  print(paste("[BEGIN] run 'step 3' calculate [min_dist]-based scores: ", Sys.time(), sep = ""))
  
  ## step 3a: set the "maximum" threshold of the [min_dist] - generated via SCimilarity's constrained prediction approach
  max_minDist <- 0.1
  
  ## step 3b: add a calculation-based column: [score1_minDist_median] ####
  sum_mindist_byClst <- sum_mindist_byClst %>%
    group_by(seurat_clusters) %>%
    mutate(
      score1_minDist_median=(max_minDist-median_minDist)
    ) %>%
    as.data.frame()
  print(head(sum_mindist_byClst,3))
  
  ## step 3c: add an identifier column: [id__cluster_celltype] ####
  sum_mindist_byClst$id__cluster_celltype <- paste(sum_mindist_byClst$seurat_clusters,sum_mindist_byClst$celltype_hint,sep="_")
  print(head(sum_mindist_byClst,3))
  
  ## step 3d: add a column from the preliminary summary table (generated from the [step 2c]): [cellProp_byClst] ####
  sum_mindist_byClst$cellProp_byClst <- sum_cellpred_byClst$cellProp_byClst[match(sum_mindist_byClst$id__cluster_celltype, sum_cellpred_byClst$id__cluster_celltype)]
  
  ## step 3e: add a calculation-based column: [score1a_minDist_median_cellProp] ####
  sum_mindist_byClst <- sum_mindist_byClst %>%
    group_by(seurat_clusters) %>%
    mutate(
      score1a_minDist_median_cellProp=((max_minDist-median_minDist)*cellProp_byClst)
    ) %>%
    arrange(-score1a_minDist_median_cellProp) %>%
    as.data.frame()
  
  ## checkpoints:
  print(":")
  print("--- [CHECKPOINT] first 10 rows:")
  print(head(sum_mindist_byClst,10))
  print(":")
  print("--- [CHECKPOINT] last 10 rows:")
  print(tail(sum_mindist_byClst,10))
  print(":")
  
  print(paste0("=== OUTPUT [sum_mindist_byClst]: nRow: ",nrow(sum_mindist_byClst)," | nColumn: ",ncol(sum_mindist_byClst)," ===")) #checkpoint
  print(paste("[END] run 'step 3' calculate [min_dist]-based scores: ", Sys.time(), sep = ""))
  
  print(":")
  print(paste0("=== [DONE] run 'calc_score_minDist': ", Sys.time()," ==="))
  # ---
  
  # ---
  ## save as TXT: [sum_mindist_byClst] ####
  if (isTRUE(save_txt)) {
    fn_txt <- paste0(mydir,"sum_mindist_byClst","_wMinDistScore_",prj,".txt")
    if (!file.exists(fn_txt)) {
      write.table(x = sum_mindist_byClst, file = fn_txt, sep = "\t", col.names = T, row.names = F, quote = F)
      print(paste0("[DONE] save as TXT: ",fn_txt))
    } else {
      print(paste0("[SKIP] TXT file exists: ",fn_txt))
    }
  } #closer: if (isTRUE(save_txt)) {
  # ---
  
  # ---
  ## save as RDS: [sum_mindist_byClst] ####
  if (isTRUE(save_rds)) {
    fn_rds <- paste0(mydir,"sum_mindist_byClst","_wMinDistScore_",prj,".rds")
    if (!file.exists(fn_rds)) {
      saveRDS(object = sum_mindist_byClst, file = fn_rds)
      print(paste0("[DONE] save as RDS: ",fn_rds))
    } else {
      print(paste0("[SKIP] RDS file exists: ",fn_rds))
    }
  } #closer: if (isTRUE(save_rds)) {
  # ---
}
# ===


## run the function: [calc_score_minDist] ####
calc_score_minDist(prj = "TEST", fn_obs = fn_obs, save_txt = T, save_rds = T)
# ---
